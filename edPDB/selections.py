# edPDB.selections 
"""
:mod:`edPDB.selections`
=================

Extensions to Bio.PDB, some useful selections.

Partly published on http://biopython.org/wiki/Reading_large_PDB_files

License: like Biopython

Module content
--------------

.. automodule:: edPDB.selections
"""

import Bio.PDB
from Bio.PDB.PDBIO import Select
from Bio.PDB.Residue import Residue

from  utilities import asiterable


#: How to recognize a protein.
PROTEIN_RESNAMES = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
                    'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 
                    'HIS':'H', 'HSD':'H', 'HSE':'H', 'HSP':'H', 
                    'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 
                    'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
                    'TRP':'W', 'TYR':'Y', 'VAL':'V', 
                    'ALAD':'AA', 'CHO':'?', 'EAM':'?'}

def canonical(resname):
    """Return canonical representation of resname.
    
    space stripped and  upper case
    """
    return resname.strip().upper()

class ResnameSelect(Select):
    """Select all atoms that match *resnames*."""
    def __init__(self, resnames, complement=False):
        """Supply a *resname*, e.g. 'SOL' or 'PHE' or a list."""
        self.resnames = dict([(canonical(r),True) for r in asiterable(resnames)])
        if not complement:
            self.accept_residue = self._accept_residue
        else:
            self.accept_residue = self._accept_not_residue
    def _accept_residue(self,residue):
        # use a dict --- 'in' is probably faster on dict keys than on
        # lists ... TODO = check ;-) --- this seems to be a bottle neck
        return canonical(residue.resname) in self.resnames
    def _accept_not_residue(self,residue):
        return not canonical(residue.resname) in self.resnames


class ResidueSelect(Select):
    """Select all atoms that are in the *residues* list."""
    def __init__(self, residues, complement=False):
        """Supply a list of Bio.PDB residues for the search."""        
        self.residues = residues
        if not complement:
            self.accept_residue = self._accept_residue
        else:
            self.accept_residue = self._accept_not_residue
    def _accept_residue(self,residue):
        return residue in self.residues
    def _accept_not_residue(self,residue):
        return not residue in self.residues

class NotResidueSelect(ResidueSelect):
    """Select all atoms that are *not* in the *residues* list.

    (Same as :class:`ResidueSelect(residues, complement=True)`.)
    """
    def __init__(self, residues, complement=False):
        """Supply a list of Bio.PDB residues for the search."""        
        ResidueSelect.__init__(self, residues,complement=(not complement))

class ProteinSelect(Select):
    """Select all amino acid residues."""
    def __init__(self, complement=False):
        if not complement:
            self.accept_residue = self._accept_residue
        else:
            self.accept_residue = self._accept_not_residue
    def _accept_residue(self,residue):
        return canonical(residue.resname) in PROTEIN_RESNAMES
    def _accept_not_residue(self,residue):
        return not canonical(residue.resname) in PROTEIN_RESNAMES



def residues_by_resname(structure, resnames):
    """Return a list of residue instances that match *resnames*.

    *resnames* can be a single string or a list of strings.
    """
    #return [r for r in Bio.PDB.Selection.unfold_entities(structure, 'R')
    #        if r.resname.strip() == resname]
    return residues_by_selection(structure, ResnameSelect(resnames))

def residues_by_selection(structure, selection):
    """General residue selection: supply a Bio.PDB.PDBIO.Select instance."""
    return [r for r in Bio.PDB.Selection.unfold_entities(structure, 'R')
            if selection.accept_residue(r)]
    

def find_water(structure, ligand, radius=3.0, water='SOL'):
    """Find all water (SOL) molecules within radius of ligand.

    :Arguments:
        structure
            Bio.PDB structure of Mhp1 system with water
         ligand : list
            Bio.PDB list of atoms of the ligand (Bio.PDB.Atom.Atom
            instances)
         radius : float
            Find waters for which the ligand-atom - OW  distance is <
            radius.
         water : string
            resname of a water molecule

    :Returns: list of residue instances
    """
    
    # get all SOL (water)
    solvent = residues_by_resname(structure, water)
    # NOT working in script (but in ipython) ?!?!
    #solvent_OW = [a for a in r.get_list() for r in solvent if a.name == 'OW']
    solvent_OW = []
    for r in solvent:
        for a in r.get_list():
            if a.name == "OW":
                solvent_OW.append(a)
    # sanity check:
    assert len(solvent) == len(solvent_OW)
    
    # set up KDtree neighbour search (use the biggest group for the
    # tree, i.e. solvent not the ligand)
    ns = Bio.PDB.NeighborSearch(solvent_OW)

    names_centers = [(a.name, a.get_coord())
                     for a in Bio.PDB.Selection.unfold_entities(ligand, 'A')]
    water_shell = AtomGroup()
    for name,center in names_centers:
        _shell = AtomGroup(ns.search(center, radius))
        logger.debug("around %6r: %3d OW = " % (name, len(_shell)) + str(_shell))
        water_shell.update(_shell)  # keep unique residues only
    return sorted([a.parent for a in water_shell])

