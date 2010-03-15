# $Id: xpdb.py 467 2009-08-19 15:34:56Z oliver $
"""
:mod:`edPDB.xpdb`
=================

Extensions to Bio.PDB, such as handling of large pdb files and some useful selections.

Partly published on http://biopython.org/wiki/Reading_large_PDB_files

License: like Biopython

Module content
--------------

.. automodule:: edPDB.xpdb
"""

import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.PDBIO import Select
from Bio.PDB.Residue import Residue

from  utilities import asiterable

import logging
logger = logging.getLogger('edPDB.xpdb')

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
    """Select all atoms that match *resname*."""
    def __init__(self,resnames):
        """Supply a *resname*, e.g. 'SOL' or 'PHE' or a list."""
        self.resnames = dict([(canonical(r),True) for r in asiterable(resnames)])
    def accept_residue(self,residue):
        # use a dict --- 'in' is probably faster on dict keys than on
        # lists ... TODO = check ;-) --- this seems to be a bottle neck
        return canonical(residue.resname) in self.resnames

class ResidueSelect(Select):
    """Select all atoms that are in the residue list."""
    def __init__(self,residues):
        """Supply a list of Bio.PDB residues for the search."""        
        self.residues = residues
    def accept_residue(self,residue):
        return residue in self.residues

class NotResidueSelect(Select):
    """Select all atoms that are *not* in the residue list."""
    def __init__(self,residues):
        """Supply a list of Bio.PDB residues for the search."""        
        self.residues = residues
    def accept_residue(self,residue):
        return not residue in self.residues

class ProteinSelect(Select):
    """Select all amino acid residues."""
    def accept_residue(self,residue):
        return canonical(residue.resname) in PROTEIN_RESNAMES

class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Cope with resSeq < 10,000 limitation by just incrementing internally.

    # Q: What's wrong here??
    #   Some atoms or residues will be missing in the data structure.
    #   WARNING: Residue (' ', 8954, ' ') redefined at line 74803.
    #   PDBConstructionException: Blank altlocs in duplicate residue SOL (' ', 8954, ' ') at line 74803.
    #
    # A: resSeq only goes to 9999 --> goes back to 0 (PDB format is not really good here)
    """

    # NOTE/TODO:
    # - H and W records are probably not handled yet (don't have examples to test)

    def __init__(self,verbose=False):
        Bio.PDB.StructureBuilder.StructureBuilder.__init__(self)
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        """
        Initiate a new Residue object.
        
        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for 
            hetero residues, otherwise blanc.
        o resseq - int, sequence identifier
        o icode - string, insertion code
        """
        if field!=" ":
            if field=="H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field="H_"+resname 
        res_id=(field, resseq, icode)

        if resseq > self.max_resseq:
            self.max_resseq  = resseq
        
        if field==" ":
            fudged_resseq = False
            while (self.chain.has_id(res_id) or resseq == 0):
                # There already is a residue with the id (field, resseq, icode).
                # resseq == 0 catches already wrapped residue numbers which do not
                # trigger the has_id() test.
                # 
                # Be sloppy and just increment...
                # (This code will not leave gaps in resids... I think)
                #
                # XXX: shouldn't we also do this for hetero atoms and water??
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)    # use max_resseq!
                fudged_resseq = True
                
            if fudged_resseq and self.verbose:
                logger.debug("Residues are wrapping (Residue ('%s', %i, '%s') at line %i)." 
                             % (field, resseq, icode, self.line_counter) +
                             ".... assigning new resid %d.\n" % self.max_resseq)
        residue=Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue=residue

class SloppyPDBIO(Bio.PDB.PDBIO):
    """PDBIO class that can deal with large pdb files as used in MD simulations.
    
    - resSeq simply wrap and are printed modulo 10,000.
    - atom numbers wrap at 99,999 and are printed modulo 100,000    
    """
    # directly copied from PDBIO.py
    # (has to be copied because of the package layout it is not externally accessible)
    _ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname, 
        resseq, icode, chain_id, element="  ", charge="  "):
        """
        Returns an ATOM PDB string that is guaranteed to fit into the ATOM format.

        - Resid (resseq) is wrapped (modulo 10,000) to fit into %4i (4I) format
        - Atom number (atom_number) is wrapped (modulo 100,000) to fit into %4i (4I) format
        """
        if hetfield!=" ":
            record_type="HETATM"
        else:
            record_type="ATOM  "
        name=atom.get_fullname()
        altloc=atom.get_altloc()
        x, y, z=atom.get_coord()
        bfactor=atom.get_bfactor()
        occupancy=atom.get_occupancy()
        args=(record_type, atom_number % 100000, name, altloc, resname, chain_id,
            resseq % 10000, icode, x, y, z, occupancy, bfactor, segid,
            element, charge)
        return self._ATOM_FORMAT_STRING % args

    


sloppyparser = Bio.PDB.PDBParser(PERMISSIVE=True,structure_builder=SloppyStructureBuilder())

def get_structure(pdbfile,pdbid='system'):
    return sloppyparser.get_structure(pdbid,pdbfile)



class AtomGroup(set):
    def __init__(self, atoms=None):
        atoms = atoms or []
        super(AtomGroup, self).__init__(atoms)

    @property
    def identifiers(self):
        return [(a.parent.resname, a.parent.id[1], a.name) for a in self]

    def __repr__(self):
        return "[" + ", ".join(["%s%d:%s" % x for x in self.identifiers]) + "]"

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

writelogger = logging.getLogger('edPDB.write_pdb')
def write_pdb(structure, filename, exclusions=None, inclusions=None, **kwargs):
    """Write Bio.PDB molecule *structure* to *filename*.
    
    :Arguments:
       *structure*
         Bio.PDB structure instance
       *filename*
         pdb file
       *exclusions*
         list of **residue** instances that will *not* be included
       *inclusions*
         list of **residue** instances that will be included
       *chain*
         set the chain identifier for **all** atoms written; this can be
         useful to  simply to erase all chain ids by setting it to ' '

    Typical use is to supply a list of water molecules that should not
    be written or a ligand that should be include.

    .. Note:: Currently only either *exclusions* or *inclusions* is 
              supported.
    """
    if exclusions and inclusions:
        raise NotImplementedError("Sorry, only either exclusions OR "
                                  "inclusions are supported")
    if exclusions:
        selection = NotResidueSelect(exclusions)
    else:
        inclusions = inclusions or []
        selection = ResidueSelect(inclusions)

    chain = kwargs.pop('chain',None)
    if not chain is None:
        writelogger.info("Setting the chain id for ALL atoms to %(chain)r",  vars())
        for c in structure.get_chains():
            c.id = chain

    writelogger.debug("Starting PDBIO...")
    io = SloppyPDBIO()     # deals with resSeq > 9999
    io.set_structure(structure)
    io.save(filename, select=selection)

    writelogger.info("Wrote pdb %(filename)r." % vars())


