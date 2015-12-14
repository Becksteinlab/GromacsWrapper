

#
# JD: adopted from 'pytopol' from https://github.com/resal81/PyTopol
#

class System(object):
    def __init__(self):
        self.molecules = tuple([])

        self.atomtypes        = []
        self.bondtypes        = []
        self.nonbond_params   = []
        self.angletypes       = []
        self.dihedraltypes    = []
        self.impropertypes    = []
        self.cmaptypes        = []
        self.interactiontypes = []
        self.pairtypes        = []
        self.constrainttypes  = []
        self.forcefield=  None

        self.information = {} # like 'atomtypes': self.atomtypes



class Molecule(object):

    def __init__(self):
        self.chains    = []
        self.atoms     = []
        self.residues  = []

        self.bonds     = []
        self.angles    = []
        self.dihedrals = []
        self.impropers = []

        self.cmaps     = []
        self.pairs     = []
        self.exclusion_numb = None  # 0, 1, 2, ..
        self.virtual_sites3 = []
        self.exclusions = []
        self.settles    = []
        self.constraints= []

        self.information = {}  # like 'atoms': self.atoms

        self.name = None

        self._anumb_to_atom = {}


    def anumb_to_atom(self, anumb):
        '''Returns the atom object corresponding to an atom number'''

        assert isinstance(anumb, int), "anumb must be integer"

        if len(self._anumb_to_atom) == 0:   # empty dictionary

            if len(self.atoms) != 0:
                for atom in self.atoms:
                    self._anumb_to_atom[atom.number] = atom
                return self._anumb_to_atom[anumb]
            else:
                print("no atoms in the molecule")
                return False

        else:
            if anumb in self._anumb_to_atom:
                return self._anumb_to_atom[anumb]
            else:
                print("no such atom number (%d) in the molecule" % (anumb))
                return False


    def renumber_atoms(self):

        if len(self.atoms) != 0:

            # reset the mapping
            self._anumb_to_atom = {}

            for i,atom in enumerate(self.atoms):
                atom.number = i+1   # starting from 1

        else:
            print("the number of atoms is zero - no renumbering")


class Chain(object):
    """
        name    = str,
        residues= list,
        molecule= Molecule
    """

    def __init__(self):
        self.residues = []


class Residue(object):
    """
        name    = str,
        number  = int,
        chain   = Chain,
        chain_name = str,
        atoms   = list,
    """

    def __init__(self):
        self.atoms  = []


class Atom(object):
    """
        name    = str,
        number  = int,
        flag    = str,        # HETATM
        coords  = list,
        residue = Residue,
        occup   = float,
        bfactor = float,
        altlocs = list,
        atomtype= str,
        radius  = float,
        charge  = float,
        mass    = float,
        chain   = str,
        resname = str,
        resnumb = int,
        altloc  = str,         # per atoms

    """

    def __init__(self):

        self.coords = []        # a list of coordinates (x,y,z) of models
        self.altlocs= []        # a list of (altloc_name, (x,y,z), occup, bfactor)




    def get_atomtype(self):
        if hasattr(self, 'atomtype'):
            return self.atomtype
        else:
            print("atom %s doesn't have atomtype" % self)
            return False


class Param(object):
    def __init__(self):
        self.comment = None
        self.line  = None
        self.disabled = False

    def convert(self, reqformat):
        assert reqformat in ('charmm', 'gromacs')

        if reqformat == self.format:
            if reqformat == 'charmm':
                return self.charmm
            elif reqformat == 'gromacs':
                return self.gromacs
            else:
                raise NotImplementedError



        if isinstance(self, AtomType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                self.gromacs['param']['lje'] = abs(self.charmm['param']['lje']) * 4.184
                self.gromacs['param']['ljl'] = self.charmm['param']['ljl'] * 2 * 0.1 / (2**(1.0/6.0))

                if self.charmm['param']['lje14'] is not None:
                    self.gromacs['param']['lje14'] = abs(self.charmm['param']['lje14']) * 4.184
                    self.gromacs['param']['ljl14'] = self.charmm['param']['ljl14'] * 2 * 0.1 / (2**(1.0/6.0))
                else:
                    self.gromacs['param']['lje14'] = None
                    self.gromacs['param']['ljl14'] = None
            else:
                raise NotImplementedError



        elif isinstance(self, BondType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                self.gromacs['param']['kb'] = self.charmm['param']['kb'] * 2 * 4.184 * (1.0 / 0.01)   # nm^2
                self.gromacs['param']['b0'] = self.charmm['param']['b0'] * 0.1
                self.gromacs['func'] = 1
            else:
                raise NotImplementedError



        elif isinstance(self, AngleType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                self.gromacs['param']['ktetha'] = self.charmm['param']['ktetha'] * 2 * 4.184
                self.gromacs['param']['tetha0'] = self.charmm['param']['tetha0']
                self.gromacs['param']['kub'] = self.charmm['param']['kub'] * 2 * 4.184 * 10 * 10
                self.gromacs['param']['s0'] = self.charmm['param']['s0'] * 0.1
                self.gromacs['func'] = 5
            else:
                raise NotImplementedError



        elif isinstance(self, DihedralType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                for dih in self.charmm['param']:
                    convdih = {}
                    convdih['kchi']  = dih['kchi'] * 4.184
                    convdih['n']     = dih['n']
                    convdih['delta'] = dih['delta']
                    self.gromacs['param'].append(convdih)
                self.gromacs['func'] = 9
            else:
                raise NotImplementedError



        elif isinstance(self, ImproperType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                for imp in self.charmm['param']:
                    convimp = {}
                    convimp['kpsi']  = imp['kpsi'] * 2 * 4.184
                    convimp['psi0']  = imp['psi0']
                    if imp.get('n', False):
                        convimp['n']     = imp['n']
                    self.gromacs['param'].append(convimp)
                self.gromacs['func'] = 2

                # self.gromacs['param']['kpsi'] = self.charmm['param']['kpsi'] * 2 * 4.184
                # self.gromacs['param']['psi0'] = self.charmm['param']['psi0']
                # self.gromacs['func'] = 2
            else:
                raise NotImplementedError



        elif isinstance(self, CMapType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                self.gromacs['param']= [n*4.184 for n in self.charmm['param']]
                self.gromacs['func'] = 1
            else:
                raise NotImplementedError



        elif isinstance(self, InteractionType):
            if reqformat == 'gromacs' and self.format == 'charmm':
                if self.charmm['param']['lje'] is not None:
                    self.gromacs['param']['lje'] = abs(self.charmm['param']['lje']) * 4.184
                    self.gromacs['param']['ljl'] = self.charmm['param']['ljl'] * 0.1 /  (2**(1.0/6.0))   # no *2
                else:
                    self.gromacs['param']['lje'] = None
                    self.gromacs['param']['ljl'] = None

                if self.charmm['param']['lje14'] is not None:
                    self.gromacs['param']['lje14'] = abs(self.charmm['param']['lje14']) * 4.184
                    self.gromacs['param']['ljl14'] = self.charmm['param']['ljl14'] * 0.1 / (2**(1.0/6.0))
                else:
                    self.gromacs['param']['lje14'] = None
                    self.gromacs['param']['ljl14'] = None
            else:
                raise NotImplementedError

        else:
            raise NotImplementedError


class AtomType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format
        self.atype  = None
        self.atnum  = None
        self.mass   = None
        self.charge = None
        self.bond_type = None

        self.charmm = {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }
        self.gromacs= {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }

        super(AtomType,self).__init__()

class BondType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format


        self.atom1 = None
        self.atom2 = None

        self.atype1 = None
        self.atype2 = None

        self.charmm = {'param': {'kb':None, 'b0':None} }
        self.gromacs= {'param': {'kb':None, 'b0':None}, 'func':None}

        super(BondType,self).__init__()


class AngleType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format


        self.atom1 = None
        self.atom2 = None
        self.atom3 = None

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None

        self.charmm = {'param':{'ktetha':None, 'tetha0':None, 'kub':None, 's0':None} }
        self.gromacs= {'param':{'ktetha':None, 'tetha0':None, 'kub':None, 's0':None}, 'func':None}

        super(AngleType,self).__init__()


class DihedralType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format


        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None
        self.atype4 = None

        self.charmm = {'param':[]}  # {kchi, n, delta}
        self.gromacs= {'param':[]}

        super(DihedralType,self).__init__()

class ImproperType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None
        self.atype4 = None

        self.charmm = {'param':[]}
        self.gromacs= {'param':[], 'func': None}  # {'kpsi': None, 'psi0':None}

        super(ImproperType,self).__init__()


class CMapType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format

        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None
        self.atom5 = None
        self.atom6 = None
        self.atom7 = None
        self.atom8 = None

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None
        self.atype4 = None
        self.atype5 = None
        self.atype6 = None
        self.atype7 = None
        self.atype8 = None

        self.charmm = {'param': []}
        self.gromacs= {'param': []}

        super(CMapType,self).__init__()


class InteractionType(Param):
    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format

        self.atom1  = None
        self.atom2  = None

        self.atype1 = None
        self.atype2 = None

        self.charmm = {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }
        self.gromacs= {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None}, 'func':None }

        super(InteractionType,self).__init__()


class SettleType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)
        self.atom = None
        self.dOH  = None
        self.dHH  = None

        super(SettleType,self).__init__()

class ConstraintType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)

        self.atom1 = None
        self.atom2 = None

        self.atype1 = None
        self.atype2 = None

        self.gromacs= {'param': {'b0':None}, 'func':None}

        super(ConstraintType,self).__init__()


class NonbondedParamType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)

        self.atype1 = None
        self.atype2 = None

        self.gromacs= {'param': {'esp':None, 'sig':None}, 'func':None}

        super(NonbondedParamType,self).__init__()


class VirtualSites3Type(Param):
    def __init__(self, format):
        assert format in ('gromacs',)

        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None

        self.gromacs= {'param': {'a':None, 'b': None}, 'func':None}

        super(VirtualSites3Type,self).__init__()

class Exclusion:
    def __init__(self):
        self.main_atom  = None
        self.other_atoms = []
