#!/usr/bin/env python

desc = """
Merge n topology files. Useful for cases when one wishes to create a bond
between two molecules (protein and lipid anchor).

The script will merge in the topologies, in the order provided, by changing 
atom ids (and corresponding bond, angle and dihedral sections) to create 
a conigous, 1-based topology file.  

Example use
-----------
- building up a lipid from chemical fragments - two tails and a head-group,
- adding a lipid anchor to a protein topology 

No terms will be added to the topology - it is up to the user to modify the 
file and add any bonded terms.

Example case
------------
Two molecules, each with two atoms and one bond:
  # moleculeA.itp
    [ molecule ]
    mol_A
    [ atoms ]
    1 <typeA>
    2 <typeB>
    [ bonds ]
    1 2 <bond term>

  #  moleculeB.itp
    [ molecule ]
    mol_B
    [ atoms ]
    1 <typeC>
    2 <typeD>
    [ bonds ]
    1 2 <bond term>

Upon merging (order matters) as -p moleculeA.itp moleculeB.itp -o topolout.itp
become:
  # topolout.itp
    [ molecule ]
    toplogy
    [ atoms ]
    1 <typeA>
    2 <typeB>
    3 <typeC>
    4 <typeD>
    [ bonds ]
    1 2 <bond term>
    3 4 <bond term>

Angle, dihedral and pair terms will be corrected the same way as the [ bonds ]
section. 

NOTE: it is up to the user to create bonded terms between the two merged
molecules. 

Author: Jan Domanski (jandom@gmail.com)
Year: 2012
License: WTFPL http://en.wikipedia.org/wiki/WTFPL
"""
from argparse import RawTextHelpFormatter
import gromacs.fileformats
import argparse
import numpy as np


def get_no_atoms(molecule):
  atoms = molecule.sections['header'].sections['moleculetype'].sections['atoms'].data.atomnr
  # The first atom has to be '1' 
  assert np.min(atoms) == 1
  # Check if atomnumbering is contigous
  assert not (range(np.min(atoms), np.max(atoms)+1) - atoms).all()
  return len(atoms)
  

def add_comment(molecule, blocks, comment):
  for block in blocks:
     molecule.sections['header'].sections['moleculetype'].sections[block].data.comment[0] += comment

def merge_topologies(topol_list, topolout, blocks, name="MOL", dirty=True):
  molecules = []
  for topol in topol_list:
    molecule = gromacs.fileformats.itp.ITP(topol)
    molecules.append(molecule)
  
  no_atoms = 0
  #
  # * Shift all the atom numbering in all sections
  #
  for i, molecule in enumerate(molecules):
    if dirty:
      add_comment(molecule, blocks, " [from %s]" % topol_list[i])

    # Adjust everything to the first molecule (i == 0), shift the atoms of the next molecules by this amount
    if not i: 
      no_atoms =  get_no_atoms(molecule)
      continue
    no_atoms0 = get_no_atoms(molecule)
    # JD The setting could be generalized using :meth:`getattr()` but here, I've decided
    # to keep some locality/redundancy, in case of a major blunder on my side. 
    molecule = molecule.sections['header'].sections['moleculetype']
    
    molecule.sections['atoms'].data.atomnr += no_atoms
    
    molecule.sections['bonds'].data.ai += no_atoms
    molecule.sections['bonds'].data.aj += no_atoms

    molecule.sections['angles'].data.ai += no_atoms
    molecule.sections['angles'].data.aj += no_atoms
    molecule.sections['angles'].data.ak += no_atoms
    
    molecule.sections['dihedrals'].data.ai += no_atoms
    molecule.sections['dihedrals'].data.aj += no_atoms
    molecule.sections['dihedrals'].data.ak += no_atoms
    molecule.sections['dihedrals'].data.al += no_atoms

    molecule.sections['pairs'].data.ai += no_atoms
    molecule.sections['pairs'].data.aj += no_atoms
   
    no_atoms +=  no_atoms0
  #
  # * Merege all the topology sections.data 
  #
  atoms = np.concatenate([m.sections['header'].sections['moleculetype'].sections['atoms'].data for m in molecules])
  bonds = np.concatenate([m.sections['header'].sections['moleculetype'].sections['bonds'].data for m in molecules])
  angle = np.concatenate([m.sections['header'].sections['moleculetype'].sections['angles'].data for m in molecules])
  dihed = np.concatenate([m.sections['header'].sections['moleculetype'].sections['dihedrals'].data for m in molecules])
  pairs = np.concatenate([m.sections['header'].sections['moleculetype'].sections['pairs'].data for m in molecules])
  
  mol = molecules[0]
  molecule = mol.sections['header'].sections['moleculetype']
  
  molecule.sections['atoms'].set_data(atoms)
  molecule.sections['bonds'].set_data(bonds)
  molecule.sections['angles'].set_data(angle)
  molecule.sections['dihedrals'].set_data(dihed)
  molecule.sections['pairs'].set_data(pairs)
  
  molecule.data['name'] = name
  
  mol.write(topolout)

def main():

  parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
  parser.add_argument('-n', dest='name', metavar="MOL", default="MOL", help="Molecule name")
  parser.add_argument('-d', dest='dirty',metavar="dirty", type=bool, default=True, help="Comment first atoms/terms with the input file they come from")  
  parser.add_argument('-p', dest='topol', metavar="topol.itp",default=["topol.itp"], nargs='*')
  parser.add_argument('-b', dest='blocks', metavar="block", default=['atoms', 'bonds', 'angles', 'dihedrals', 'pairs'], nargs='*')
  parser.add_argument('-po', dest='topolout', metavar="topolout.itp", default="topolout.itp")
  args = parser.parse_args()
  
  merge_topologies(args.topol, args.topolout, args.blocks, args.name, args.dirty)


if __name__ == "__main__":
    main()
