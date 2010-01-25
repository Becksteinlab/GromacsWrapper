# build module
"""
:mod:`edPDB.cbook` -- Recipes for editing PDB files
===================================================

The cook book contains short python functions that deomstrate how to
implemen basic PDB editing functionality. They do not do exhaustive
error checking and might have to be altered for your purpose.


.. autofunction:: align_ligand
.. autofunction:: remove_overlap_water

"""

import os.path

import Bio.PDB

import xpdb
import config


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

    # XXX: untested since I changed the code to use xpdb
    prot = xpdb.get_structure('prot', protein_struct)
    lig = xpdb.get_structure('lig', ligand_struct)
    plig = xpdb.residues_by_resname(prot, ligand_resname)

    heavyMoving = [a for a in lig.get_atoms() if a.name[0] != 'H']
    heavyFixed = Bio.PDB.Selection.unfold_entities(plig, 'A')  # only heavies in pdb

    S = Bio.PDB.Superimposer()
    S.set_atoms(heavyFixed, heavyMoving)    # get rotation matrix
    S.apply(lig.get_atoms())                # rotate ligand atoms

    io = Bio.PDB.PDBIO()
    io.set_structure(lig)
    io.save(output)

    print "Wrote aligned ligand  %r (RMSD = %g A)" % (ligand_aligned, S.rms)
    return S.rms


def remove_overlap_water(pdbname, output, ligand_resname, distance=3.0):
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

    structure = xpdb.get_structure(pdbname)
    ligand = xpdb.residues_by_resname(structure, ligand_resname)

    w = xpdb.find_water(structure, ligand, radius=distance)
    print "w = %r" % w

    xpdb.write_pdb(structure, output, exclusions=w)


    
    
