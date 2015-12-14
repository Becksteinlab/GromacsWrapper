
#
# JD: adopted from 'pytopol' from https://github.com/resal81/PyTopol
#

import logging
from . import blocks
from collections import OrderedDict as odict

module_logger = logging.getLogger('gromacs.fileformats')


class TOP(blocks.System):
    def __init__(self, fname):

        super(TOP, self).__init__()

        self.lgr = logging.getLogger('gromacs.fileformats.TOP')
        self.fname = fname

        self.defaults = {
            'nbfunc': None, 'comb-rule':None, 'gen-pairs':None, 'fudgeLJ':None, 'fudgeQQ':None,
        }

        self.dict_molname_mol = odict()   # contains molname:mol
        self.found_sections   = []
        self.forcefield       = 'gromacs'

        self.molecules = []
        self._parse(fname)
        self.molecules = tuple(self.molecules)

    def write(self, filename):
        SystemToGroTop(self, filename)

    def __repr__(self):
        moltypenames = list(self.dict_molname_mol.keys())
        moltypenames.sort()

        data = []
        data.append('\n')

        main_items = set(['atomtypes', 'pairtypes', 'bondtypes', 'angletypes', 'dihedraltypes'])
        other_items = ['%s (%d)' % (m, len(self.information[m])) for m in list(self.information.keys()) if m not in main_items]
        other_items = ' '.join(other_items)
        nattype = len(self.atomtypes)
        nprtype = len(self.pairtypes)
        nbndtype= len(self.bondtypes)
        nangtype= len(self.angletypes)
        ndihtype= len(self.dihedraltypes)
        nimptype= len(self.impropertypes)
        data.append('{:>20s}  {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s}'.format('Param types:', 'atom', 'pair', 'bond', 'ang', 'dih', 'imp'))
        msg = '{:20s}  {:7d} {:7d} {:7d} {:7d} {:7d} {:7d}    {:s}'.format('', nattype, nprtype, nbndtype, nangtype, ndihtype, nimptype, other_items)
        data.append('=' * 69)
        data.append(msg)
        data.append('\n')


        main_items = set(['atoms', 'pairs', 'bonds', 'angles', 'dihedrals'])
        data.append('{:>20s}  {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s}'.format('Params:', 'atom', 'pair', 'bond', 'ang', 'dih', 'imp'))
        data.append('=' * 69)
        for mname in moltypenames:
            mol = self.dict_molname_mol[mname]
            other_items = ['%s (%d)' % (m, len(mol.information[m])) for m in list(mol.information.keys()) if m not in main_items]
            other_items = ' '.join(other_items)

            natoms = len(mol.atoms)
            npairs = len(mol.pairs)
            nbonds = len(mol.bonds)
            nangles= len(mol.angles)
            ndih   = len(mol.dihedrals)
            nimp   = len(mol.impropers)
            msg = '{:20s}  {:7d} {:7d} {:7d} {:7d} {:7d} {:7d}    {:s}'.format(mol.name, natoms, npairs, nbonds, nangles, ndih, nimp, other_items)
            data.append(msg)




        return '\n'.join(data)




    def _parse(self, fname):

        _find_section = lambda line: line.strip('[').strip(']').strip()

        def _add_info(sys_or_mol, section, container):
            # like (mol, 'atomtypes', mol.atomtypes)
            if sys_or_mol.information.get(section, False) is False:
                sys_or_mol.information[section] = container

        mol        = None   # to hold the current mol
        curr_sec   = None
        cmap_lines = []

        with open(fname) as f:
            for i_line, line in enumerate(f):

                # trimming
                if ';' in line:
                    line = line[0:line.index(';')]
                line = line.strip()

                if line == '':
                    continue

                if line[0] == '*':
                    continue

                # the topology must be stand-alone (i.e. no includes)
                if line.startswith('#include'):
                    msg = 'The topology file has "#include" statements.'
                    msg+= ' You must provide a processed topology file that grompp creates.'
                    raise ValueError(msg)

                # find sections
                if line[0] == '[':
                    curr_sec = _find_section(line)
                    self.found_sections.append(curr_sec)
                    continue

                fields = line.split()

                if curr_sec == 'defaults':
                    '''
                    # ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
                    #1               2               yes             0.5     0.8333
                    '''

                    assert len(fields) == 5

                    self.defaults['nbfunc']    = int(fields[0])
                    self.defaults['comb-rule'] = int(fields[1])
                    self.defaults['gen-pairs'] = fields[2]
                    self.defaults['fudgeLJ']   = float(fields[3])
                    self.defaults['fudgeQQ']   = float(fields[4])

                elif curr_sec == 'atomtypes':
                    '''
                    # ;name               at.num    mass         charge    ptype  sigma   epsilon
                    # ;name   bond_type   at.num    mass         charge    ptype  sigma   epsilon
                    '''

                    if len(fields) not in (7,8):
                        print('skipping atomtype line with neither 7 or 8 fields: \n %s' % line)
                        continue

                    shift = 0 if len(fields) == 7 else 1
                    at = blocks.AtomType('gromacs')
                    at.atype = fields[0]
                    if shift == 1: at.bond_type = fields[1]

                    at.mass  = float(fields[2+shift])
                    at.charge= float(fields[3+shift])

                    particletype = fields[4+shift]
                    assert particletype in ('A', 'S', 'V', 'D')
                    if particletype not in ('A',):
                        print('warning: non-atom particletype: "%s"' % line)

                    sig = float(fields[5+shift])
                    eps = float(fields[6+shift])

                    at.gromacs= {'param': {'lje':eps, 'ljl':sig, 'lje14':None, 'ljl14':None} }

                    self.atomtypes.append(at)

                    _add_info(self, curr_sec, self.atomtypes)


                # extend system.molecules
                elif curr_sec == 'moleculetype':
                    assert len(fields) == 2

                    mol = blocks.Molecule()

                    mol.name = fields[0]
                    mol.exclusion_numb = int(fields[1])

                    self.dict_molname_mol[mol.name] = mol


                elif curr_sec == 'atoms':
                    '''
                    #id    at_type     res_nr  residu_name at_name  cg_nr  charge   mass  typeB    chargeB      massB
                    # 1       OC          1       OH          O1       1      -1.32

                    OR

                    [ atoms ]
                    ; id   at type  res nr  residu name at name     cg nr   charge
                    1       OT      1       SOL              OW             1       -0.834

                    '''

                    aserial = int(fields[0])
                    atype   = fields[1]
                    resnumb = int(fields[2])
                    resname = fields[3]
                    aname   = fields[4]
                    cgnr    = int(fields[5])
                    charge  = float(fields[6])
                    rest = fields[7:]

                    atom         = blocks.Atom()
                    atom.name    = aname
                    atom.atomtype= atype
                    atom.number  = aserial
                    atom.resname = resname
                    atom.resnumb = resnumb
                    atom.charge  = charge

                    if len(rest) >= 1:
                        mass = float(rest[0])
                        atom.mass = mass

                    mol.atoms.append(atom)

                    _add_info(mol, curr_sec, mol.atoms)

                elif curr_sec in ('pairtypes', 'pairs', 'pairs_nb'):
                    '''
                    section     #at     fu      #param
                    ---------------------------------
                    pairs       2       1       V,W
                    pairs       2       2       fudgeQQ, qi, qj, V, W
                    pairs_nb    2       1       qi, qj, V, W

                    '''

                    ai, aj = fields[:2]
                    fu     = int(fields[2])
                    assert fu in (1,2)

                    pair = blocks.InteractionType('gromacs')
                    if fu == 1:
                        if curr_sec=='pairtypes':
                            pair.atype1 = ai
                            pair.atype2 = aj
                            v, w = list(map(float, fields[3:5]))
                            pair.gromacs = {'param': {'lje':None, 'ljl':None, 'lje14':w, 'ljl14':v}, 'func':fu }

                            self.pairtypes.append(pair)
                            _add_info(self, curr_sec, self.pairtypes)

                        elif curr_sec == 'pairs':
                            ai, aj = list( map(int, [ai,aj]) )
                            pair.atom1 = mol.atoms[ai-1]
                            pair.atom2 = mol.atoms[aj-1]
                            pair.gromacs['func'] = fu

                            mol.pairs.append(pair)
                            _add_info(mol, curr_sec, mol.pairs)

                        else:
                            raise ValueError

                    else:
                        raise NotImplementedError('%s with functiontype %d is not supported' % (curr_sec,fu))

                elif curr_sec == 'nonbond_params':
                    '''
                    ; typei typej  f.type sigma   epsilon
                    ; f.type=1 means LJ (not buckingham)
                    ; sigma&eps since mixing-rule = 2          
                    '''

                    assert len(fields) == 5
                    ai, aj = fields[:2]
                    fu     = int(fields[2])

                    assert fu == 1
                    sig    = float(fields[3])
                    eps    = float(fields[4])

                    nonbond_param = blocks.NonbondedParamType('gromacs')
                    nonbond_param.atype1 = ai
                    nonbond_param.atype2 = aj
                    nonbond_param.gromacs['func'] = fu
                    nonbond_param.gromacs['param'] = {'eps': eps, 'sig': sig} 

                    self.nonbond_params.append(nonbond_param)
                    _add_info(self, curr_sec, self.nonbond_params)

                elif curr_sec in ('bondtypes', 'bonds'):
                    '''
                    section     #at     fu      #param
                    ----------------------------------
                    bonds       2       1       2
                    bonds       2       2       2
                    bonds       2       3       3
                    bonds       2       4       2
                    bonds       2       5       ??
                    bonds       2       6       2
                    bonds       2       7       2
                    bonds       2       8       ??
                    bonds       2       9       ??
                    bonds       2       10      4
                    '''

                    ai, aj = fields[:2]
                    fu     = int(fields[2])
                    assert fu in (1,2,3,4,5,6,7,8,9,10)

                    if fu != 1:
                        raise NotImplementedError('function %d is not yet supported' % fu)

                    bond = blocks.BondType('gromacs')

                    if fu == 1:
                        if curr_sec == 'bondtypes':
                            bond.atype1 = ai
                            bond.atype2 = aj
                            b0, kb = list(map(float, fields[3:5]))

                            bond.gromacs = {'param':{'kb':kb, 'b0':b0}, 'func':fu}

                            self.bondtypes.append(bond)
                            _add_info(self, curr_sec, self.bondtypes)

                        elif curr_sec == 'bonds':
                            ai, aj = list(map(int, [ai, aj]))
                            bond.atom1 = mol.atoms[ai-1]
                            bond.atom2 = mol.atoms[aj-1]
                            bond.gromacs['func'] = fu

                            mol.bonds.append(bond)
                            _add_info(mol, curr_sec, mol.bonds)

                    else:
                        raise NotImplementedError

                elif curr_sec in ('angletypes', 'angles'):
                    '''
                    section     #at     fu      #param
                    ----------------------------------
                    angles      3       1       2
                    angles      3       2       2
                    angles      3       3       3
                    angles      3       4       4
                    angles      3       5       4
                    angles      3       6       6
                    angles      3       8       ??
                    '''

                    ai, aj , ak = fields[:3]
                    fu          = int(fields[3])
                    assert fu in (1,2,3,4,5,6,8)  # no 7

                    if fu not in (1,5):
                        raise NotImplementedError('function %d is not yet supported' % fu)

                    ang = blocks.AngleType('gromacs')
                    if fu == 1:
                        if curr_sec == 'angletypes':
                            ang.atype1 = ai
                            ang.atype2 = aj
                            ang.atype3 = ak

                            tetha0, ktetha = list(map(float, fields[4:6]))
                            ang.gromacs = {'param':{'ktetha':ktetha, 'tetha0':tetha0, 'kub':None, 's0':None}, 'func':fu}

                            self.angletypes.append(ang)
                            _add_info(self, curr_sec, self.angletypes)

                        elif curr_sec == 'angles':
                            ai, aj, ak = list(map(int, [ai, aj, ak]))
                            ang.atom1 = mol.atoms[ai-1]
                            ang.atom2 = mol.atoms[aj-1]
                            ang.atom3 = mol.atoms[ak-1]
                            ang.gromacs['func'] = fu

                            mol.angles.append(ang)
                            _add_info(mol, curr_sec, mol.angles)

                        else:
                            raise ValueError

                    elif fu == 5:
                        if curr_sec == 'angletypes':
                            ang.atype1 = ai
                            ang.atype2 = aj
                            ang.atype3 = ak
                            tetha0, ktetha, s0, kub = list(map(float, fields[4:8]))

                            ang.gromacs = {'param':{'ktetha':ktetha, 'tetha0':tetha0, 'kub':kub, 's0':s0}, 'func':fu}

                            self.angletypes.append(ang)
                            _add_info(self, curr_sec, self.angletypes)

                        elif curr_sec == 'angles':
                            ai, aj, ak = list(map(int, [ai, aj, ak]))
                            ang.atom1 = mol.atoms[ai-1]
                            ang.atom2 = mol.atoms[aj-1]
                            ang.atom3 = mol.atoms[ak-1]
                            ang.gromacs['func'] = fu

                            mol.angles.append(ang)
                            _add_info(mol, curr_sec, mol.angles)

                        else:
                            raise ValueError

                    else:
                        raise NotImplementedError


                elif curr_sec in  ('dihedraltypes', 'dihedrals'):
                    '''
                    section     #at     fu      #param
                    ----------------------------------
                    dihedrals   4       1       3
                    dihedrals   4       2       2
                    dihedrals   4       3       6
                    dihedrals   4       4       3
                    dihedrals   4       5       4
                    dihedrals   4       8       ??
                    dihedrals   4       9       3
                    '''

                    if curr_sec == 'dihedraltypes' and len(fields) == 6:
                        # in oplsaa - quartz parameters
                        fields.insert(2, 'X')
                        fields.insert(0, 'X')

                    ai, aj, ak, am = fields[:4]
                    fu = int(fields[4])
                    assert fu in (1,2,3,4,5,8,9)

                    if fu not in (1,2,3,4,9):
                        raise NotImplementedError('function %d is not yet supported' % fu)

                    dih = blocks.DihedralType('gromacs')
                    imp = blocks.ImproperType('gromacs')

                    if fu in (1,3,4,9):
                        if curr_sec == 'dihedraltypes':
                            dih.atype1 = ai
                            dih.atype2 = aj
                            dih.atype3 = ak
                            dih.atype4 = am

                            dih.line = i_line + 1  

                            if fu == 1:
                                delta, kchi, n = list(map(float, fields[5:8]))
                                dih.gromacs['param'].append({'kchi':kchi, 'n':n, 'delta':delta})
                            elif fu == 3:
                                c0, c1, c2, c3, c4, c5 = list(map(float, fields[5:11]))
                                m = dict(c0=c0, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5)
                                dih.gromacs['param'].append(m)
                            elif fu == 4:
                                delta, kchi, n = list(map(float, fields[5:8]))
                                dih.gromacs['param'].append({'kchi':kchi, 'n':int(n), 'delta':delta})                                                                                               
                            elif fu == 9:
                                delta, kchi, n = list(map(float, fields[5:8]))
                                dih.gromacs['param'].append({'kchi':kchi, 'n':int(n), 'delta':delta})
                            else:
                                raise ValueError

                            dih.gromacs['func'] = fu
                            self.dihedraltypes.append(dih)
                            _add_info(self, curr_sec, self.dihedraltypes)

                        elif curr_sec == 'dihedrals':
                            ai, aj, ak, am = list(map(int, fields[:4]))
                            dih.atom1 = mol.atoms[ai-1]
                            dih.atom2 = mol.atoms[aj-1]
                            dih.atom3 = mol.atoms[ak-1]
                            dih.atom4 = mol.atoms[am-1]
                            dih.gromacs['func'] = fu

                            dih.line = i_line + 1

                            if fu == 1:
                                pass
                            elif fu == 3:
                                pass
                            elif fu == 4:
                                pass                                
                            elif fu == 9:
                                if len(fields[5:8]) == 3:
                                    delta, kchi, n = list(map(float, fields[5:8]))
                                    dih.gromacs['param'].append({'kchi':kchi, 'n':int(n), 'delta':delta})
                            else:
                                raise ValueError

                            mol.dihedrals.append(dih)
                            _add_info(mol, curr_sec, mol.dihedrals)

                        else:
                            raise ValueError

                    elif fu in (2,4):
                        if curr_sec == 'dihedraltypes':
                            imp.atype1 = ai
                            imp.atype2 = aj
                            imp.atype3 = ak
                            imp.atype4 = am

                            imp.line = i_line + 1

                            if fu == 2:
                                psi0 , kpsi = list(map(float, fields[5:7]))
                                imp.gromacs['param'].append({'kpsi':kpsi, 'psi0': psi0})
                            elif fu == 4:
                                psi0 , kpsi, n = list(map(float, fields[5:8]))
                                imp.gromacs['param'].append({'kpsi':kpsi, 'psi0': psi0, 'n':n})
                            else:
                                raise ValueError

                            imp.gromacs['func'] = fu
                            self.impropertypes.append(imp)
                            _add_info(self, curr_sec, self.impropertypes)

                        elif curr_sec == 'dihedrals':
                            ai, aj, ak, am = list(map(int, fields[:4]))
                            imp.atom1 = mol.atoms[ai-1]
                            imp.atom2 = mol.atoms[aj-1]
                            imp.atom3 = mol.atoms[ak-1]
                            imp.atom4 = mol.atoms[am-1]
                            imp.gromacs['func'] = fu

                            imp.line = i_line + 1

                            if fu == 2:
                                pass
                            elif fu == 4:
                                pass
                            else:
                                raise ValueError

                            mol.impropers.append(imp)
                            _add_info(mol, curr_sec, mol.impropers)

                        else:
                            raise ValueError

                    else:
                        raise NotImplementedError


                elif curr_sec in ('cmaptypes', 'cmap'):

                    cmap = blocks.CMapType('gromacs')
                    if curr_sec == 'cmaptypes':
                        cmap_lines.append(line)
                        _add_info(self, curr_sec, self.cmaptypes)
                    else:
                        ai, aj, ak, am, an = list(map(int, fields[:5]))
                        fu = int(fields[5])
                        assert fu == 1
                        cmap.atom1 = mol.atoms[ai-1]
                        cmap.atom2 = mol.atoms[aj-1]
                        cmap.atom3 = mol.atoms[ak-1]
                        cmap.atom4 = mol.atoms[am-1]
                        cmap.atom8 = mol.atoms[an-1]
                        cmap.gromacs['func'] = fu

                        mol.cmaps.append(cmap)
                        _add_info(mol, curr_sec, mol.cmaps)


                elif curr_sec == 'settles':
                    '''
                    section     #at     fu      #param
                    ----------------------------------
                    '''

                    assert len(fields) == 4
                    ai = int(fields[0])
                    fu = int(fields[1])
                    assert fu == 1

                    settle = blocks.SettleType('gromacs')
                    settle.atom = mol.atoms[ai-1]
                    settle.dOH = float(fields[2])
                    settle.dHH = float(fields[3])

                    mol.settles.append(settle)
                    _add_info(mol, curr_sec, mol.settles)

                elif curr_sec == "virtual_sites3":
                    '''
                        ; Dummy from            funct   a       b
                        4   1   2   3   1   0.131937768 0.131937768                    
                    '''
                    assert len(fields) == 7
                    ai = int(fields[0])
                    aj = int(fields[1])
                    ak = int(fields[2])
                    al = int(fields[3])
                    fu = int(fields[4])
                    assert fu == 1
                    a = float(fields[5])
                    b = float(fields[6])

                    vs3 = blocks.VirtualSites3Type('gromacs')
                    vs3.atom1 = ai
                    vs3.atom2 = aj
                    vs3.atom3 = ak
                    vs3.atom4 = al
                    vs3.gromacs['func'] = fu
                    vs3.gromacs['param'] = { 'a': a, 'b':b }
                    mol.virtual_sites3.append(vs3)
                    _add_info(mol, curr_sec, mol.virtual_sites3)


                elif curr_sec in ('exclusions',):
                    ai = int(fields[0])
                    other = list(map(int, fields[1:]))

                    exc = blocks.Exclusion()
                    exc.main_atom  = mol.atoms[ai-1]
                    exc.other_atoms= [mol.atoms[k-1] for k in other]

                    mol.exclusions.append(exc)
                    _add_info(mol, curr_sec, mol.exclusions)


                elif curr_sec in ('constrainttypes', 'constraints'):
                    '''
                    section     #at     fu      #param
                    ----------------------------------
                    constraints 2       1       1
                    constraints 2       2       1
                    '''

                    ai, aj = fields[:2]
                    fu = int(fields[2])
                    assert fu in (1,2)

                    cons = blocks.ConstraintType('gromacs')

                    # TODO: what's different between 1 and 2
                    if fu == 1 or fu == 2:
                        if curr_sec == 'constrainttypes':
                            cons.atype1 = ai
                            cons.atype2 = aj
                            b0 = float(fields[3])
                            cons.gromacs = {'param':{'b0':b0}, 'func': fu}

                            self.constrainttypes.append(cons)
                            _add_info(self, curr_sec, self.constrainttypes)

                        elif curr_sec == 'constraints':
                            ai, aj = list(map(int, fields[:2]))
                            cons.atom1 = mol.atoms[ai-1]
                            cons.atom2 = mol.atoms[aj-1]
                            cons.gromacs['func'] = fu

                            mol.constraints.append(cons)
                            _add_info(mol, curr_sec, mol.constraints)

                        else:
                            raise ValueError
                    else:
                        raise ValueError

                elif curr_sec in ('position_restraints',
                                  'distance_restraints',
                                  'dihedral_restraints',
                                  'orientation_restraints',
                                  'angle_restraints',
                                  'angle_restraints_z'):
                    pass


                elif curr_sec in ('implicit_genborn_params',):
                    '''
                    attype   sar     st      pi      gbr      hct
                    '''
                    pass

                elif curr_sec == 'system':
                    #assert len(fields) == 1
                    self.name = fields[0]


                elif curr_sec == 'molecules':
                    assert len(fields) == 2
                    mname, nmol = fields[0], int(fields[1])

                    # if the number of a molecule is more than 1, add copies to system.molecules
                    for i in range(nmol):
                        self.molecules.append(self.dict_molname_mol[mname])

                else:
                    print('Uknown section in topology: %s' % curr_sec)
        
        # process cmap_lines
        curr_cons = None
        for line in cmap_lines:

            # cmaptype opening line
            if len(line.split()) == 8:
                cons = blocks.CMapType('gromacs')
                
                atype1, atype2, atype3, atype4, atype8, func, sizeX, sizeY = line.replace("\\","").split()
                func, sizeX, sizeY = int(func), int(sizeX), int(sizeY)
                cons.atype1 = atype1
                cons.atype2 = atype2
                cons.atype3 = atype3
                cons.atype4 = atype4
                cons.atype8 = atype8
                cons.gromacs = {'param':[], 'func': func}

                curr_cons = cons

            # cmap body
            elif len(line.split()) == 10:
                cmap_param = map(float, line.replace("\\","").split())
                cons.gromacs['param'] += cmap_param

            # cmaptype cloning line
            elif len(line.split()) == 6:
                cmap_param = map(float, line.replace("\\","").split())
                cons.gromacs['param'] += cmap_param
                self.cmaptypes.append(curr_cons)
            else:
                raise ValueError



class SystemToGroTop(object):

    formats = {
        'atomtypes'      : '{:<7s} {:3s} {:3d} {:>7.6f}   {:4.3f}   {:3s}     {:14.12f}     {:10.9f}  \n',
        'atoms'          : '{:6d} {:>10s} {:6d} {:6s} {:6s} {:6d} {:f} {:11.4f} \n',
        'atoms_nomass'   : '{:6d} {:>10s} {:6d} {:6s} {:6s} {:6d} {:f}\n',
        'nonbond_params' : '{:20s}  {:20s}  {:1d}  {:10.5f}  {:10.5f}\n',
        'bondtypes'      : '{:5s}  {:5s}  {:1d}  {:6.4f}  {:6.1f}\n',
        'bonds'          : '{:3d}  {:3d}   {:1d}\n',
        'settles'        : '{:3d}  {:3d}  {:11.5f} {:11.5f}\n',
        'virtual_sites3' : '{:3d}  {:3d}  {:3d}  {:3d}   {:1d}  {:12.10f}  {:12.10f}\n',
        'exclusions'     : '{:3d}  {:3d}  {:3d}\n',
        'pairtypes'      : '{:6s} {:6s}   {:d}    {:14.12f}     {:14.12f}    \n',
        'pairs'          : '{:3d} {:3d}   {:1d}\n',
        'angletypes_1'   : '{:>8s} {:>8s} {:>8s} {:1d}    {:8.4f}    {:10.5f}\n',
        'angletypes_5'   : '{:>8s} {:>8s} {:>8s} {:1d}    {:8.4f}    {:10.5f}    {:9.5f}    {:11.5f}\n',
        'constrainttypes': '{:6s} {:6s} {:1d}    {:8.6f}\n',
        'angles'         : '{:3d} {:3d} {:3d}   {:1d}\n',
        'dihedraltypes'  : '{:6s} {:6s} {:6s} {:6s}   {:1d}    {:6.2f}    {:f}    {:1d}\n',
        'dihedrals'      : '{:3d} {:3d} {:3d} {:3d}   {:1d}\n',
        'dihedrals_ext'  : '{:3d} {:3d} {:3d} {:3d}   {:1d}    {:6.2f}    {:f}    {:1d}\n',
        'impropertypes'  : '{:6s} {:6s} {:6s} {:6s}   {:1d} {:6.2f} {:8.4f} \n',
        'impropers'      : '{:3d} {:3d} {:3d} {:3d}   {:1d}\n',
        'impropers_ext'  : '{:3d} {:3d} {:3d} {:3d}   {:1d} {:6.2f} {:8.4f} \n',
    }


    toptemplate = ""
    toptemplate += "[ defaults ]      \n*DEFAULTS*    \n"
    toptemplate += "[ atomtypes ]      \n*ATOMTYPES*    \n"
    toptemplate += "[ nonbond_params ] \n*NONBOND_PARAM* \n"
    toptemplate += "[ pairtypes ]    \n*PAIRTYPES*    \n"
    toptemplate += "[ bondtypes ]    \n*BONDTYPES*    \n"
    toptemplate += "[ angletypes ]   \n*ANGLETYPES*   \n"
    toptemplate += "[ constrainttypes ]   \n*CONSTRAINTTYPES*   \n"
    toptemplate += "[ dihedraltypes ]\n*DIHEDRALTYPES*\n"
    toptemplate += "[ dihedraltypes ]\n*IMPROPERTYPES*\n"
    toptemplate += "[ cmaptypes ]    \n*CMAPTYPES*\n"

    itptemplate = ""
    itptemplate += "[ moleculetype ] \n*MOLECULETYPE* \n"
    itptemplate += "[ atoms ]        \n*ATOMS*        \n"
    itptemplate += "[ bonds ]        \n*BONDS*        \n"
    itptemplate += "[ pairs ]        \n*PAIRS*        \n"
    itptemplate += "[ settles ]        \n*SETTLES*        \n"
    itptemplate += "[ virtual_sites3 ]        \n*VIRTUAL_SITES3*        \n"
    itptemplate += "[ exclusions ]        \n*EXCLUSIONS*        \n"
    itptemplate += "[ angles ]       \n*ANGLES*       \n"
    itptemplate += "[ dihedrals ]    \n*DIHEDRALS*    \n"
    itptemplate += "[ dihedrals ]    \n*IMPROPERS*    \n"
    itptemplate += "[ cmap ]        \n*CMAPS*    \n"




    def __init__(self, psfsystem, outfile="top.top", multiple_output=False):
        self.lgr = logging.getLogger('gromacs.fileformats.SystemToGroTop')
        self.lgr.debug(">> entering SystemToGroTop")

        self.system   = psfsystem
        self.outfile = outfile
        self.multiple_output = multiple_output
        self.assemble_topology()

        self.lgr.debug("<< leaving SystemToGroTop")


    @staticmethod
    def _redefine_atomtypes(mol):

        i = 1

        for atom in mol.atoms:
            atom.atomtype = 'at%03d' % i
            i += 1


    def assemble_topology(self, redefine_atom_types = False):

        self.lgr.debug("starting to assemble topology...")

        top = '[ defaults ] ; \n'
        top += ';nbfunc    comb-rule    gen-pairs    fudgeLJ    fudgeQQ \n'

        if self.system.forcefield == 'charmm':
            top += '1          2           yes          1.0       1.0 \n'

        self.lgr.debug("making atom/pair/bond/angle/dihedral/improper types")
        top += self.toptemplate
        top = top.replace('*DEFAULTS*',       ''.join( self._make_defaults(self.system)) )
        top = top.replace('*ATOMTYPES*',      ''.join( self._make_atomtypes(self.system)) )
        top = top.replace('*NONBOND_PARAM*',  ''.join( self._make_nonbond_param(self.system)) )
        top = top.replace('*PAIRTYPES*',      ''.join( self._make_pairtypes(self.system)) )
        top = top.replace('*BONDTYPES*',      ''.join( self._make_bondtypes(self.system)) )
        top = top.replace('*CONSTRAINTTYPES*',''.join( self._make_constrainttypes(self.system)))
        top = top.replace('*ANGLETYPES*',     ''.join( self._make_angletypes(self.system)))
        top = top.replace('*DIHEDRALTYPES*',  ''.join( self._make_dihedraltypes(self.system)) )
        top = top.replace('*IMPROPERTYPES*',  ''.join( self._make_impropertypes(self.system)) )
        top = top.replace('*CMAPTYPES*',      ''.join( self._make_cmaptypes(self.system)) )


        for i,(molname,m) in enumerate(self.system.dict_molname_mol.items()):

            itp = self.itptemplate
            itp = itp.replace('*MOLECULETYPE*',  ''.join( self._make_moleculetype(m, molname))  )
            itp = itp.replace('*ATOMS*',         ''.join( self._make_atoms(m))  )
            itp = itp.replace('*BONDS*',         ''.join( self._make_bonds(m))  )
            itp = itp.replace('*PAIRS*',         ''.join( self._make_pairs(m))  )
            itp = itp.replace('*SETTLES*',       ''.join( self._make_settles(m))  )
            itp = itp.replace('*VIRTUAL_SITES3*',''.join( self._make_virtual_sites3(m))  )
            itp = itp.replace('*EXCLUSIONS*',    ''.join( self._make_exclusions(m))  )
            itp = itp.replace('*ANGLES*',        ''.join( self._make_angles(m)) )
            itp = itp.replace('*DIHEDRALS*',     ''.join( self._make_dihedrals(m)) )
            itp = itp.replace('*IMPROPERS*',     ''.join( self._make_impropers(m)) )
            itp = itp.replace('*CMAPS*',         ''.join( self._make_cmaps(m)) )
            if not self.multiple_output:
                top += itp
            else:
                outfile = "mol_{}.itp".format(molname)
                top += '#include "mol_%s.itp" \n' % molname                
                with open(outfile, "w") as f:
                    f.writelines([itp])

        top += '\n[system]  \nConvertedSystem\n\n'
        top += '[molecules] \n'
        molecules = [("", 0)]

        for m in self.system.molecules:
            if (molecules[-1][0] != m.name): 
                molecules.append([m.name, 0])
            if molecules[-1][0] == m.name:
                molecules[-1][1] += 1

        for molname, n in molecules[1:]:
            top += '%s     %s\n' % (molname, n)
        top += '\n'

        with open(self.outfile, 'w') as f:
            f.writelines([top])

        return
        for i, m in enumerate(self.system.molecules):
            molname = 'mol_%02d' % (i+1)
            top += '#include "itp_%s.itp" \n' % molname

        
        



        self.lgr.debug('writing top finished')


        self.lgr.debug("generating atom/pair/bond/angle/dihedral/improper for the itp files")

        return
        for i,m in enumerate(self.system.molecules):
            molname = 'mol_%02d' % (i+1)
            itp = self.itptemplate
            itp = itp.replace('*MOLECULETYPE*',  ''.join( self._make_moleculetype(m, molname))  )
            itp = itp.replace('*ATOMS*',         ''.join( self._make_atoms(m))  )
            itp = itp.replace('*BONDS*',         ''.join( self._make_bonds(m))  )
            itp = itp.replace('*PAIRS*',         ''.join( self._make_pairs(m))  )
            itp = itp.replace('*ANGLES*',        ''.join( self._make_angles(m)) )
            itp = itp.replace('*DIHEDRALS*',     ''.join( self._make_dihedrals(m)) )
            itp = itp.replace('*IMPROPERS*',     ''.join( self._make_impropers(m)) )
            itp = itp.replace('*CMAPS*',         ''.join( self._make_cmaps(m)) )

            #with open('itp_%s.itp' % molname, 'w') as f:
            #    f.writelines([itp])
            with open('top.top', 'a') as f:
                f.writelines([top])

        self.lgr.debug('writing %d itp files finished' % (i+1))

    def _make_defaults(self,m):
        return ['{:d}          {:d}           {}          {:.1f}       {:.1f} \n'.format(m.defaults['nbfunc'], m.defaults['comb-rule'], m.defaults['gen-pairs'] , m.defaults['fudgeLJ'], m.defaults['fudgeQQ'])]


    def _make_atomtypes(self,m):
        def get_prot(at):
            # TODO improve this
            _protons = {'C':6, 'H':1, 'N':7, 'O':8, 'S':16, 'P':15}
            if at[0] in list(_protons.keys()):
                return _protons[at[0]]
            else:
                return 0

        result = []
        for at in m.atomtypes:
            at.convert('gromacs')
            prot = get_prot(at.atype)
            ljl  = at.gromacs['param']['ljl']
            lje  = at.gromacs['param']['lje']
            line = self.formats['atomtypes'].format(at.atype, at.bond_type if at.bond_type else "", prot, at.mass, at.charge, 'A', ljl, lje)
            if at.comment : line += at.comment
            result.append(line)

        return result

    def _make_nonbond_param(self, m):
        result = []
        for pr in m.nonbond_params:
            at1 = pr.atype1
            at2 = pr.atype2

            #pr.convert('gromacs')
            eps = pr.gromacs['param']['eps']
            sig = pr.gromacs['param']['sig']

            fu = 1  # TODO
            line = self.formats['nonbond_params'].format(at1, at2, fu, sig, eps)
            result.append(line)

        return result

    def _make_pairtypes(self,m):

        result = []
        for pt in m.pairtypes:
            at1, at2 = pt.atype1, pt.atype2
            fu, l14, e14 = pt.gromacs['func'], pt.gromacs['param']['ljl14'], pt.gromacs['param']['lje14']
            line = self.formats['pairtypes'].format(at1, at2, fu, l14, e14)
            if pt.comment : line = line[:-1] + pt.comment
            result.append(line)

        return result



    def _make_bondtypes(self,m):
        result = []
        for bond in m.bondtypes:
            at1 = bond.atype1
            at2 = bond.atype2
            bond.convert('gromacs')

            kb = bond.gromacs['param']['kb']
            b0 = bond.gromacs['param']['b0']
            fu = bond.gromacs['func']

            line = self.formats['bondtypes'].format(at1, at2, fu, b0, kb)
            result.append(line)

        return result


    def _make_constrainttypes(self,m):
        result = []

        for con in m.constrainttypes:
            at1 = con.atype1
            at2 = con.atype2

            fu  = con.gromacs['func']
            b0  = con.gromacs['param']['b0']

            line = self.formats['constrainttypes'].format(at1, at2, fu, b0)
            result.append(line)

        return result


    def _make_angletypes(self,m):
        result = []
        for ang in m.angletypes:
            at1 = ang.atype1
            at2 = ang.atype2
            at3 = ang.atype3
            ang.convert('gromacs')

            ktetha = ang.gromacs['param']['ktetha']
            tetha0 = ang.gromacs['param']['tetha0']
            kub    = ang.gromacs['param']['kub']
            s0     = ang.gromacs['param']['s0']

            fu = ang.gromacs['func']

            angletypes = 'angletypes_{:d}'.format(fu)
            line = self.formats[angletypes].format(at1, at2, at3, fu, tetha0, ktetha, s0, kub)
            result.append(line)

        return result


    def _make_dihedraltypes(self,m):
        result = []
        for dih in m.dihedraltypes:
            at1 = dih.atype1
            at2 = dih.atype2
            at3 = dih.atype3
            at4 = dih.atype4

            dih.convert('gromacs')
            fu = dih.gromacs['func']

            for dpar in dih.gromacs['param']:
                kchi = dpar['kchi']
                n    = dpar['n']
                delta= dpar['delta']
                
                if not dih.disabled:
                    line = self.formats['dihedraltypes'].format(at1, at2, at3, at4, fu, delta, kchi, n)
                else:
                    line = self.formats['dihedraltypes'].format(at1, at2, at3, at4, fu, delta, kchi, n)
                    line = dih.comment + line
                result.append(line)

        return result

    def _make_impropertypes(self,m):
        result = []
        for imp in m.impropertypes:
            at1 = imp.atype1
            at2 = imp.atype2
            at3 = imp.atype3
            at4 = imp.atype4

            imp.convert('gromacs')
            fu = imp.gromacs['func']

            for ipar in imp.gromacs['param']:

                kpsi = ipar['kpsi']
                psi0 = ipar['psi0']

                if not imp.disabled:
                    line = self.formats['impropertypes'].format(at1, at2, at3, at4, fu, psi0, kpsi)
                else: 
                    line = self.formats['impropertypes'].format(at1, at2, at3, at4, fu, psi0, kpsi)
                    line = imp.comment + line
                result.append(line)

        return result

    def _make_cmaptypes(self, m):
        result = []
        for cmap in m.cmaptypes:
            at1 = cmap.atype1
            at2 = cmap.atype2
            at3 = cmap.atype3
            at4 = cmap.atype4
            #at5 = cmap.atype5
            #at6 = cmap.atype6
            #at7 = cmap.atype7
            at8 = cmap.atype8

            cmap.convert('gromacs')

            fu = cmap.gromacs['func']
            line = '%s %s %s %s %s %d 24 24' % (at1, at2, at3, at4, at8, fu)
            for i,c in enumerate(cmap.gromacs['param']):
                if i%10 == 0:
                    line += '\\\n'
                else:
                    line += ' '
                line += '%12.8f' % c

            line += '\n\n'
            result.append(line)

        return result

    def _make_moleculetype(self,m, molname):
        return ['; Name \t\t  nrexcl \n %s    3 \n' % molname]

    def _make_atoms(self,m):
        result = []
        #i = 1
        for atom in m.atoms:
            numb = cgnr = atom.number
            atype = atom.get_atomtype()
            
            assert atype!= False and hasattr(atom, 'charge') #and hasattr(atom, 'mass')

            if hasattr(atom, 'mass'):
                line = self.formats['atoms'].format(
                    numb, atype, atom.resnumb, atom.resname, atom.name, cgnr, atom.charge, atom.mass)
            else:
                line = self.formats['atoms_nomass'].format(
                    numb, atype, atom.resnumb, atom.resname, atom.name, cgnr, atom.charge)
            result.append(line)

        result.insert(0,'; %5d atoms\n' % len(result))
        return result

    def _make_pairs(self,m):

        result = []
        for pr in m.pairs:
            fu = 1
            p1 = pr.atom1.number
            p4 = pr.atom2.number

            line = self.formats['pairs'].format(p1, p4, fu)
            result.append(line)

        result.insert(0,'; %5d pairs\n' % len(result))
        return result


    def _make_bonds(self,m):
        result = []
        for bond in m.bonds:
            fu = 1
            line = self.formats['bonds'].format(bond.atom1.number, bond.atom2.number, fu)
            result.append(line)

        result.insert(0,'; %5d bonds\n' % len(result))
        return result

    def _make_angles(self,m):
        result = []
        for ang in m.angles:
            fu = 5
            line = self.formats['angles'].format(ang.atom1.number, ang.atom2.number, ang.atom3.number, fu)
            result.append(line)

        result.insert(0,'; %5d angles\n' % len(result))
        return result

    def _make_settles(self,m):
        result = []
        for st in m.settles:
            line = self.formats['settles'].format(st.atom.number, 1, st.dOH, st.dHH)
            result.append(line)

        result.insert(0,'; %5d settles\n' % len(result))
        return result


    def _make_virtual_sites3(self,m):
        result = []
        for vs in m.virtual_sites3:
            fu = 1
            line = self.formats['virtual_sites3'].format(vs.atom1, vs.atom2, vs.atom3, vs.atom4, fu, vs.gromacs['param']['a'], vs.gromacs['param']['b'])
            result.append(line)

        result.insert(0,'; %5d virtual_sites3\n' % len(result))
        return result


    def _make_exclusions(self,m):
        result = []
        for excl in m.exclusions:
            
            line = self.formats['exclusions'].format(excl.main_atom.number, excl.other_atoms[0].number, excl.other_atoms[1].number)
            result.append(line)

        result.insert(0,'; %5d exclusions\n' % len(result))
        return result

    def _make_dihedrals(self,m):
        result = []
        for dih in m.dihedrals:
            fu = 9
            
            if not len(dih.gromacs['param']):
                line = self.formats['dihedrals'].format(
                    dih.atom1.number, dih.atom2.number, dih.atom3.number, dih.atom4.number, fu)
                result.append(line)

            for dpar in dih.gromacs['param']:
                kchi = dpar['kchi']
                n    = dpar['n']
                delta= dpar['delta']

                line = self.formats['dihedrals_ext'].format(dih.atom1.number, dih.atom2.number, dih.atom3.number, dih.atom4.number, fu, delta, kchi, n)
                if dih.comment: line = dih.comment + line
                result.append(line)

        result.insert(0,'; %5d dihedrals\n' % len(result))
        return result

    def _make_impropers(self,m):
        result = []
        for imp in m.impropers:
            fu = 2

            if not len(imp.gromacs['param']):
                line = self.formats['impropers'].format(
                    imp.atom1.number, imp.atom2.number, imp.atom3.number, imp.atom4.number, fu)
                result.append(line)

            for ipar in imp.gromacs['param']:
                kpsi = ipar['kpsi']
                psi0 = ipar['psi0']

                line = self.formats['impropers_ext'].format(imp.atom1.number, imp.atom2.number, imp.atom3.number, imp.atom4.number, fu, psi0, kpsi)
                if imp.comment: line = imp.comment + line
                result.append(line)

        result.insert(0,'; %5d impropers\n' % len(result))
        return result

    def _make_cmaps(self, m):
        result = []

        for cmap in m.cmaps:
            fu = 1
            line = '%5d %5d %5d %5d %5d   %d\n' % (
                cmap.atom1.number, cmap.atom2.number, cmap.atom3.number, cmap.atom4.number,
                cmap.atom8.number, fu)
            result.append(line)

        result.insert(0,'; %5d cmaps\n' % len(result))
        return result



if __name__ == '__main__':
    import sys
    grotop = GroTop(sys.argv[1])
    print(grotop)



