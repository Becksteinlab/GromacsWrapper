About the tests
---------------

These are somewhat non-standard tests for the gromacs.fileformats.TOP module

The philosophy behing the test is to read-in a processed.top file, parse it, 
then write it out as output.top, and read-in and parse output.top. This cycle 
of parse-save-parse (deserialize-serialize-deserialize, if you wish) is a 
minimal design to check if the read/write procedure introduces errors. 

Tests are of two types:

- GROMACS-based, check the energies with `mdrun -rerun`; however, they may not 
  reveal discrepancies to unused force-field terms/parameters (since only a 
  subset of the force-field parameters is used in a protein/water system)

- pure-python tests that have no concept of energies but are more careful in 
  checking if all terms are accounted for

Running tests
-------------

Ther are 3 directories 'amber03star', 'amber03w' and 'charmm22st', each one 
contains the same chemical system - a protein in a water box. 

In each directory running 

  $ bash ../test.sh 

Should return no errors, if gromacs energies match (see test_energies.py, 
for how the comparision is done). 

To run python-only test, execute

  $ bash test-python.sh


Future work
-----------

1. Tests for modifications, a modification should introduce only the desired changes. 
2. Re-write both the GROMACS test and the python tests as unittest test (i.e. standard 
   way in which python handles testing)


