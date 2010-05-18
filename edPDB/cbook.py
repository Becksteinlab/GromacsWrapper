# edPDB.cbook
"""
:mod:`edPDB.cbook` -- Recipes for editing PDB files
===================================================

The cook book contains short python functions that demonstrate how to
implement basic PDB editing functionality. They do not do exhaustive
error checking and might have to be altered for your purpose.

.. autoclass:: PDB
   :members:
.. autofunction:: align_ligand
.. autofunction:: remove_overlap_water
.. autofunction:: extract_residue

"""

import os.path
from warnings import warn

import Bio.PDB

import xpdb
import selections

import logging
logger = logging.getLogger('edPDB.cbook')

def align_ligand(protein_struct, ligand_struct, ligand_resname, output='ligand_aligned.pdb'):
    """Align a ligand to the same ligand in a protein, based on the heavy atoms.

    This is useful when a new ligand was generated with hydrogens but
    the position in space changed.

    :Arguments:
      *protein_struct*
         protein + ligand pdb
      *ligand_struct*
         ligand pdb
      *ligand_resname*
         residue name of the ligand in the file *protein_struct*
      *ligand_aligned*
         output

    :Returns: RMSD of the fit in Angstroem.

    .. Warning:: Assumes only heavy atoms in PDB (I think... check source!) 
    """
    logger.debug("align_ligand(%(protein_struct)r, %(ligand_struct)r, "
                "%(ligand_resname)r, output=%(output)r)" % vars())

    prot = xpdb.get_structure('prot', protein_struct)
    lig = xpdb.get_structure('lig', ligand_struct)
    plig = selections.residues_by_resname(prot, ligand_resname)

    heavyMoving = [a for a in lig.get_atoms() if a.name[0] != 'H']
    heavyFixed = Bio.PDB.Selection.unfold_entities(plig, 'A')  # only heavies in pdb

    S = Bio.PDB.Superimposer()
    S.set_atoms(heavyFixed, heavyMoving)    # get rotation matrix
    S.apply(lig.get_atoms())                # rotate ligand atoms

    io = Bio.PDB.PDBIO()
    io.set_structure(lig)
    io.save(output)

    logger.info("Wrote aligned ligand  %r (RMSD = %g A)" % (ligand_aligned, S.rms))
    return S.rms


def remove_overlap_water(pdbname, output, ligand_resname, distance=3.0, water="SOL", **kwargs):
    """Remove water (SOL) molecules overlapping with ligand.

    :Arguments:
      pdbname
        pdb file that contains the ligand and the water
      output
        pdb output filename
      ligand_resname
        name of the ligand residue(s) in the pdb
      distance
        overlap is defined as a centre-centre distance of any solvent
        OW atom with any ligand atom of less than *distance*

    .. note:: The residue and atom numbering will be fairly
              meaningless in the final PDB because it wraps at 100,000
              or 10,000.

              Also make sure that there are either consistent chain
              identifiers or none (blank) because otherwise the
              residue blocks migh become reordered. (This is due to the
              way the Bio.PDB.PDBIO writes files.)
    """
    logger.debug("remove_overlap_water(%(pdbname)r, %(output)r, %(ligand_resname)r, "
                "distance=%(distance)r)" % vars())

    structure = xpdb.get_structure(pdbname)
    ligand = selections.residues_by_resname(structure, ligand_resname)

    w = xpdb.find_water(structure, ligand, radius=distance, water=water)
    logger.debug("waters found: %r" % w)
    logger.info("removed %d %r molecules overlapping %r", 
                len(w), water, ligand_resname)
    xpdb.write_pdb(structure, output, exclusions=w, **kwargs)

# the following are deprecated

def extract_resnames(pdbname, output, resnames, **kwargs):
    """Write a pdb file with *resname* extracted.

    """
    warn("extract_resname() will be removed", category=DeprecationWarning)
    logger.debug("extract_residue(%(pdbname)r, %(output)r, %(resnames)r)" % vars())
    structure = xpdb.get_structure(pdbname)
    residues = selections.residues_by_resname(structure, resnames)
    xpdb.write_pdb(structure, output, inclusions=residues, **kwargs)

def extract_protein(pdbname, output, **kwargs):
    """Write a pdb file with the protein (i.e. all amino acids) extracted.
    """
    warn("extract_protein() will be removed", category=DeprecationWarning)
    logger.debug("extract_protein(%(pdbname)r, %(output)r)" % vars())
    structure = xpdb.get_structure(pdbname)
    residues = selections.residues_by_selection(structure, xpdb.ProteinSelect())
    xpdb.write_pdb(structure, output, inclusions=residues, **kwargs)

def extract_lipids(pdbname, output, lipid_resnames='POPC|POPG|POPE|DMPC|DPPE|DOPE', **kwargs):
    """Write a pdb file with the lipids extracted.

    Note that resnames are also tried truncated to the first three
    characters, which means that POPE and POPG are identical and
    cannot be distinguished.
    """
    warn("extract_lipids() will be removed", category=DeprecationWarning)
    logger.debug("extract_lipids(%(pdbname)r, %(output)r, lipid_resnames=%(lipid_resnames)r)" % vars())
    resnames = lipid_resnames.split('|')
    resnames.extend([r[:3] for r in resnames])

    structure = xpdb.get_structure(pdbname)
    residues = selections.residues_by_selection(structure, xpdb.ResnameSelect(resnames))
    xpdb.write_pdb(structure, output, inclusions=residues, **kwargs)


class PDB(object):
    """Class that represents a PDB file and allows extractions of interesting parts.
    
    The structure itself is never changed. In order to extract
    sub-parts of a structure one selects these parts and writes them
    as new pdb file.

    The advantage over a simple :program:`grep` is that you will be able to
    read any odd pdb file and you will also able to do things like
    :meth:`~edPDB.cbook.PDB.extract_protein` or
    :meth:`~edPDB.cbook.PDB.extract_lipids`.
    """

    def __init__(self, pdbname):
        """Load structure from file *pdbname*."""
        self.pdbname = pdbname
        self.structure = xpdb.get_structure(pdbname)
        self.logger = logging.getLogger('edPDB.PDB')
        self.logger.info("Loaded pdb file %(pdbname)r." % vars())

    def write(self, filename, **kwargs):
        """Write pdbfile which includes or excludes residues.

        :Arguments:
            *filename*
                output pdb filename
            *inclusions*
                list of residues to include
            *exclusions*
                list of residues to exclude
            *chain*
                relabel the selection with a new chain identifier


         Residues must be BioPDB residues as returned by, for
         instance, :meth:`~edPDB.cbook.PDB.residues_by_resname`.

         .. Note:: Currently only either *inclusions* or *exclusions* can be
                   supplied, not both.
        """
        self.logger.debug("write(): file %r, args %r", filename, kwargs.keys())
        xpdb.write_pdb(self.structure, filename, **kwargs)
        
    def residues_by_resname(self, resnames, **kwargs):
        """Return a list of BioPDB residues that match *resnames*.

        *resnames* can be a string or a list.
        """
        return selections.residues_by_resname(self.structure, resnames)

    def residues_by_selection(self, selection):
        """Return  a list of BioPDB residues that are selected by *selection*.

        *selection* must be BioPDB.PDB.PDBIO.Select instance (see for
        example :class:`edPDB.xpdb.ProteinSelect`).
        """
        return selections.residues_by_selection(self.structure, selection)

    def extract_resnames(self, filename, resnames, **kwargs):
        """Write a pdb file with *resnames* extracted."""
        self.logger.info("extract_resnames(%(filename)r, %(resnames)r)" % vars())
        kwargs['selection'] = selections.ResnameSelect(resnames)
        self.write(filename, **kwargs)

    def extract_protein(self, filename, **kwargs):
        """Write a pdb file with the protein (i.e. all amino acids) extracted."""
        self.logger.info("extract_protein(%(filename)r)" % vars())
        kwargs['selection'] = selections.ProteinSelect()
        self.write(filename, **kwargs)

    def extract_notprotein(self, filename, **kwargs):
        """Write a pdb file without any amino acids extracted."""
        self.logger.info("extract_notprotein(%(filename)r)" % vars())
        kwargs['selection'] = selections.NotProteinSelect()
        self.write(filename, **kwargs)

    def extract_lipids(self, filename, lipid_resnames='POPC|POPG|POPE|DMPC|DPPE|DOPE', **kwargs):
        """Write a pdb file with the lipids extracted.

        Note that resnames are also tried truncated to the first three
        characters, which means that POPE and POPG are identical and
        cannot be distinguished.
        """
        self.logger.info("extract_lipids(%(filename)r, lipid_resnames=%(lipid_resnames)r)" % vars())
        resnames = lipid_resnames.split('|')
        resnames.extend([r[:3] for r in resnames])
        kwargs['selection'] = selections.ResnameSelect(resnames)
        self.write(filename, **kwargs)
        

    
