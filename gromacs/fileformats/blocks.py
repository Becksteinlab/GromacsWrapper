# GromacsWrapper: top.py
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Gromacs TOP - BLOCKS boiler-plate code
======================================


Classes
-------

.. autoclass:: System
    :members:

.. autoclass:: Molecule
    :members:

.. autoclass:: Atom
    :members:

.. autoclass:: Param
    :members:

.. autoclass:: AtomType
    :members:
.. autoclass:: BondType
    :members:
.. autoclass:: AngleType
    :members:
.. autoclass:: DihedralType
    :members:
.. autoclass:: ImproperType
    :members:
.. autoclass:: CMapType
    :members:
.. autoclass:: InteractionType
    :members:
.. autoclass:: SettleType
    :members:
.. autoclass:: ConstraintType
    :members:
.. autoclass:: NonbondedParamType
    :members:
.. autoclass:: VirtualSites3Type
    :members:

.. autoclass:: Exclusion
    :members:

History
-------

Sources adapted from code by Reza Salari https://github.com/resal81/PyTopol

"""

import logging

class System(object):
    """Top-level class containing molecule topology.

       Contains all the parameter types (AtomTypes, BondTypes, ... )
       and molecules.

    """
    logger = logging.getLogger('gromacs.formats.BLOCKS')

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
    """Class that represents a Molecule

    Contains all the molecule attributes: atoms, bonds, angles dihedrals.
    Also contains settle, cmap and exclusion sections, if present.

    .. Example input::

        ; Include water topology
        [ moleculetype ]
        ; molname      nrexcl
        SOL             2

        [ atoms ]
        ; id   at type  res nr  residu name     at name         cg nr   charge
        1       OWT3    1       SOL              OW             1       -0.834
        2       HWT3    1       SOL             HW1             1        0.417
        3       HWT3    1       SOL             HW2             1        0.417


        [ settles ]
        ; i       j     funct   length
        1         1     0.09572 0.15139

        [ exclusions ]
        1 2 3
        2 1 3
        3 1 2



    """
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

        if not self._anumb_to_atom:   # empty dictionary

            if self.atoms:
                for atom in self.atoms:
                    self._anumb_to_atom[atom.number] = atom
                return self._anumb_to_atom[anumb]
            else:
                self.logger("no atoms in the molecule")
                return False

        else:
            if anumb in self._anumb_to_atom:
                return self._anumb_to_atom[anumb]
            else:
                self.logger("no such atom number ({0:d}) in the molecule".format(anumb))
                return False


    def renumber_atoms(self):
        """Reset the molecule's atoms :attr:`number` to be 1-indexed"""
        if self.atoms:

            # reset the mapping
            self._anumb_to_atom = {}

            for i,atom in enumerate(self.atoms):
                atom.number = i+1   # starting from 1

        else:
            self.logger("the number of atoms is zero - no renumbering")

class Atom(object):
    """Class that represents an Atom

    Contains only the simplest atom attributes, that are contained like in
    section example below.

    :class:`Molecule` cantains an :attr:`atoms` that's a list-container for
    :class:`Atom` instances.

    .. Example input::

        [ atoms ]
        ;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
        ; residue  43 SER rtp SER  q  0.0
            22        NH1     43    SER      N     22      -0.47     14.007   ; qtot 0.53
            23          H     43    SER     HN     23       0.31      1.008   ; qtot 0.84
            24        CT1     43    SER     CA     24       0.07     12.011   ; qtot 0.91
            25         HB     43    SER     HA     25       0.09      1.008   ; qtot 1
            26        CT2     43    SER     CB     26       0.05     12.011   ; qtot 1.05
            27         HA     43    SER    HB1     27       0.09      1.008   ; qtot 1.14
            28         HA     43    SER    HB2     28       0.09      1.008   ; qtot 1.23
            29        OH1     43    SER     OG     29      -0.66     15.999   ; qtot 0.57
            30          H     43    SER    HG1     30       0.43      1.008   ; qtot 1
            31          C     43    SER      C     31       0.51     12.011   ; qtot 1.51
            32          O     43    SER      O     32      -0.51     15.999   ; qtot 1


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

        self.name     = None
        self.atomtype = None
        self.number   = None
        self.resname  = None
        self.resnumb  = None
        self.charge   = None

    def get_atomtype(self):
        if hasattr(self, 'atomtype'):
            return self.atomtype
        else:
            self.logger("atom {0} doesn't have atomtype".format(self))
            return False


class Param(object):
    """Class that represents an abstract Parameter.

    This class is the parent to AtomType, BondType and all the other parameter types.

    The class understands a parameter line and that a :attr:`comment` that may follow.
    CMapType is an exception (it's a multi-line parameter).

    :meth:`convert` provides a rudimentary support for parameter unit conversion between
    GROMACS and CHARMM notation: change kJ/mol into kcal/mol and nm into Angstrom.

    :attr:`disabled` for supressing output when writing-out to a file.
    """

    def __init__(self, format):
        assert format in ('charmm', 'gromacs')
        self.format = format

        self.comment = None
        self.line  = None
        self.disabled = False
        self.charmm = None
        self.gromacs = None

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

        super(AtomType,self).__init__(format)

        self.atype  = None
        self.atnum  = None
        self.mass   = None
        self.charge = None
        self.bond_type = None

        self.charmm = {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }
        self.gromacs= {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }

    def __eq__(self, other):
        return \
            self.atype == other.atype and \
            self.atnum == other.atnum and \
            self.mass == other.mass and \
            self.charge == other.charge and \
            self.bond_type == other.bond_type and \
            self.charmm == other.charmm

    def __repr__(self):
        return '<{0!s} {1!s} m={2:g} q={3:g} (gromacs:{4!s})>'.format(
            self.__class__.__name__, self.atype, self.mass, self.charge, self.gromacs)


class BondType(Param):
    def __init__(self, format):

        super(BondType,self).__init__(format)

        self.atom1 = None
        self.atom2 = None

        self.atype1 = None
        self.atype2 = None

        self.charmm = {'param': {'kb':None, 'b0':None} }
        self.gromacs= {'param': {'kb':None, 'b0':None}, 'func':None}

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm


class AngleType(Param):
    def __init__(self, format):

        super(AngleType,self).__init__(format)

        self.atom1 = None
        self.atom2 = None
        self.atom3 = None

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None

        self.charmm = {'param':{'ktetha':None, 'tetha0':None, 'kub':None, 's0':None} }
        self.gromacs= {'param':{'ktetha':None, 'tetha0':None, 'kub':None, 's0':None}, 'func':None}

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.atype3 == other.atype3 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm

class DihedralType(Param):
    def __init__(self, format):

        super(DihedralType,self).__init__(format)

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

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.atype3 == other.atype3 and \
            self.atype4 == other.atype4 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm

class ImproperType(Param):
    def __init__(self, format):

        super(ImproperType,self).__init__(format)

        self.atype1 = None
        self.atype2 = None
        self.atype3 = None
        self.atype4 = None

        self.charmm = {'param':[]}
        self.gromacs= {'param':[], 'func': None}  # {'kpsi': None, 'psi0':None}

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.atype3 == other.atype3 and \
            self.atype4 == other.atype4 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm

class CMapType(Param):
    def __init__(self, format):

        super(CMapType,self).__init__(format)

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

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.atype3 == other.atype3 and \
            self.atype4 == other.atype4 and \
            self.atype5 == other.atype5 and \
            self.atype6 == other.atype6 and \
            self.atype7 == other.atype7 and \
            self.atype8 == other.atype8 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm

class InteractionType(Param):
    def __init__(self, format):

        super(InteractionType,self).__init__(format)

        self.atom1  = None
        self.atom2  = None

        self.atype1 = None
        self.atype2 = None

        self.charmm = {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None} }
        self.gromacs= {'param': {'lje':None, 'ljl':None, 'lje14':None, 'ljl14':None}, 'func':None }

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm

    def __repr__(self):
        return '<{0!s} {1!s} {2!s} (gromacs:{3!s})>'.format(
            self.__class__.__name__, self.atype1, self.atype2, self.gromacs)


class SettleType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)
        super(SettleType,self).__init__(format)

        self.atom = None
        self.dOH  = None
        self.dHH  = None


class ConstraintType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)
        super(ConstraintType,self).__init__(format)

        self.atom1 = None
        self.atom2 = None

        self.atype1 = None
        self.atype2 = None

        self.gromacs= {'param': {'b0':None}, 'func':None}

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm


class NonbondedParamType(Param):
    def __init__(self, format):
        assert format in ('gromacs',)
        super(NonbondedParamType,self).__init__(format)

        self.atype1 = None
        self.atype2 = None

        self.gromacs= {'param': {'eps':None, 'sig':None}, 'func':None}

    def __eq__(self, other):
        return \
            self.atype1 == other.atype1 and \
            self.atype2 == other.atype2 and \
            self.gromacs == other.gromacs and \
            self.charmm == other.charmm


class VirtualSites3Type(Param):
    def __init__(self, format):
        assert format in ('gromacs',)
        super(VirtualSites3Type,self).__init__(format)

        self.atom1 = None
        self.atom2 = None
        self.atom3 = None
        self.atom4 = None

        self.gromacs= {'param': {'a':None, 'b': None}, 'func':None}


class Exclusion(object):
    """Class to define non-interacting pairs of atoms, or "exclusions".


    .. Note::

       Does not inherit from :class:`Param` unlike other classes in :mod:`blocks`
    """
    def __init__(self):
        self.main_atom  = None
        self.other_atoms = []
