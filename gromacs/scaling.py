# GromacsWrapper: scaling.py
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""

:mod:`gromacs.scaling` -- Partial tempering
===========================================

:Author: Jan Domanski, @jandom

.. versionadded:: 0.5.0

Helper functions for scaling gromacs topologies; useful for setting up
simulations with Hamiltonian replicate exchange and partial tempering
(REST2).

.. autofunction:: scale_dihedrals
.. autofunction:: scale_impropers
.. autofunction:: partial_tempering

"""
from __future__ import absolute_import

import math
import copy

import numpy as np

from .fileformats import TOP
from .fileformats import blocks


def scale_dihedrals(mol, dihedrals, scale, banned_lines=None):
        """Scale dihedral angles"""

        if banned_lines is None:
                banned_lines = []
        new_dihedrals = []
        for dh in mol.dihedrals:
                atypes = dh.atom1.get_atomtype(), dh.atom2.get_atomtype(), dh.atom3.get_atomtype(), dh.atom4.get_atomtype()
                atypes = [a.replace("_", "").replace("=","") for a in atypes]
                for iswitch in range(32):
                        if  (iswitch%2==0 ):
                                a1=atypes[0]; a2=atypes[1]; a3=atypes[2]; a4=atypes[3]
                        else:
                                a1=atypes[3]; a2=atypes[2]; a3=atypes[1]; a4=atypes[0]

                        if((iswitch//2)%2==1): a1="X";
                        if((iswitch//4)%2==1): a2="X";
                        if((iswitch//8)%2==1): a3="X";
                        if((iswitch//16)%2==1): a4="X";
                        key = "{0}-{1}-{2}-{3}-{4}".format(a1, a2, a3, a4, dh.gromacs['func'])
                        if (key in dihedrals):
                                for i, dt in enumerate(dihedrals[key]):
                                        dhA = copy.deepcopy(dh)
                                        param = copy.deepcopy(dt.gromacs['param'])
                                        # Only check the first dihedral in a list
                                        if not dihedrals[key][0].line in banned_lines:
                                                for p in param: p['kchi'] *= scale
                                        dhA.gromacs['param'] = param
                                        #if key == "CT3-C-NH1-CT1-9": print i, dt, key
                                        if i == 0:
                                                dhA.comment = "; banned lines {0} found={1}\n".format(" ".join(map(str, banned_lines)), 1 if dt.line in banned_lines else 0)
                                                dhA.comment += "; parameters for types {}-{}-{}-{}-9 at LINE({})\n".format(dhA.atom1.atomtype, dhA.atom2.atomtype, dhA.atom3.atomtype, dhA.atom4.atomtype, dt.line).replace("_","")
                                        name = "{}-{}-{}-{}-9".format(dhA.atom1.atomtype, dhA.atom2.atomtype, dhA.atom3.atomtype, dhA.atom4.atomtype).replace("_","")
                                        #if name == "CL-CTL2-CTL2-HAL2-9": print dihedrals[key], key
                                        new_dihedrals.append(dhA)
                                break


        mol.dihedrals = new_dihedrals
        #assert(len(mol.dihedrals) == new_dihedrals)
        return mol

def scale_impropers(mol, impropers, scale, banned_lines=None):
        """Scale improper dihedrals"""
        if banned_lines is None:
                banned_lines = []
        new_impropers = []
        for im in mol.impropers:
                atypes = im.atom1.get_atomtype(), im.atom2.get_atomtype(), im.atom3.get_atomtype(), im.atom4.get_atomtype()
                atypes = [a.replace("_", "").replace("=","") for a in atypes]
                for iswitch in range(32):
                        if  (iswitch%2==0 ):
                                a1=atypes[0]; a2=atypes[1]; a3=atypes[2]; a4=atypes[3];
                        else:
                                a1=atypes[3]; a2=atypes[2]; a3=atypes[1]; a4=atypes[0];
                        if((iswitch/2)%2==1): a1="X";
                        if((iswitch/4)%2==1): a2="X";
                        if((iswitch/8)%2==1): a3="X";
                        if((iswitch/16)%2==1): a4="X";
                        key = "{0}-{1}-{2}-{3}-{4}".format(a1, a2, a3, a4, im.gromacs['func'])
                        if (key in impropers):
                                for i, imt in enumerate(impropers[key]):
                                        imA = copy.deepcopy(im)
                                        param = copy.deepcopy(imt.gromacs['param'])
                                        # Only check the first dihedral in a list
                                        if not impropers[key][0].line in banned_lines:
                                                for p in param: p['kpsi'] *= scale
                                        imA.gromacs['param'] = param
                                        if i == 0: imA.comment = "; banned lines {0} found={1}\n ; parameters for types {2}-{3}-{4}-{5}-9 at LINE({6})\n".format(" ".join(map(str, banned_lines)), 1 if imt.line in banned_lines else 0 , imt.atype1, imt.atype2, imt.atype3, imt.atype4, imt.line)
                                        new_impropers.append(imA)
                                break
        #assert(len(mol.impropers) == new_impropers)
        mol.impropers = new_impropers
        return mol


def partial_tempering(args):
        """Set up topology for partial tempering (REST2) replica exchange"""

        banned_lines = map(int, args.banned_lines.split())
        top = TOP(args.input)
        groups = [("_", args.scale_protein), ("=", args.scale_lipids)]

        #
        # CMAPTYPES
        #
        cmaptypes = []
        for ct in top.cmaptypes:
                cmaptypes.append(ct)
                for gr, scale in groups:
                        ctA = copy.deepcopy(ct)
                        ctA.atype1 += gr; ctA.atype2 += gr; ctA.atype3 += gr; ctA.atype4 += gr;  ctA.atype8 += gr;
                        ctA.gromacs['param'] = [ v*scale for v in ct.gromacs['param'] ]
                        cmaptypes.append(ctA)
        print("cmaptypes was {0}, is {1}".format(len(top.cmaptypes), len(cmaptypes)))
        top.cmaptypes = cmaptypes


        #
        # ATOMTYPES
        #
        atomtypes = []
        for at in top.atomtypes:
                atomtypes.append(at)
                for gr, scale in groups:
                        atA = copy.deepcopy(at)
                        atA.atnum = atA.atype
                        atA.atype += gr
                        atA.gromacs['param']['lje'] *= scale
                        atomtypes.append(atA)
        top.atomtypes = atomtypes

        #
        # PAIRTYPES
        #
        pairtypes = []
        for pt in top.pairtypes:
                pairtypes.append(pt)
                for gr, scale in groups:
                        ptA = copy.deepcopy(pt)
                        ptA.atype1 += gr; ptA.atype2 += gr
                        ptA.gromacs['param']['lje14'] *= scale

                        pairtypes.append(ptA)
        top.pairtypes = pairtypes

        #
        # BONDTYPES
        #
        bondtypes = []
        for bt in top.bondtypes:
                #break
                bondtypes.append(bt)
                for gr, scale in groups:
                        btA = copy.deepcopy(bt)
                        btA.atype1 += gr; btA.atype2 += gr
                        bondtypes.append(btA)
        top.bondtypes = bondtypes


        #
        # ANGLETYPES
        #
        angletypes = []
        for at in top.angletypes:
                #break
                angletypes.append(at)
                for gr, scale in groups:
                        atA = copy.deepcopy(at)
                        atA.atype1 += gr; atA.atype2 += gr;  atA.atype3 += gr;
                        angletypes.append(atA)
        top.angletypes = angletypes

        #
        # Build dihedral dictionary
        #
        dihedraltypes = {}
        for dt in top.dihedraltypes:
                dt.disabled = True
                dt.comment = "; type=%s-%s-%s-%s-9\n; LINE(%d) " % (dt.atype1, dt.atype2, dt.atype3, dt.atype4, dt.line)
                dt.comment = dt.comment.replace("_","")

                #if "X-CTL2-CTL2-X-9" in dt.comment: print dt
                name = "{0}-{1}-{2}-{3}-{4}".format(dt.atype1, dt.atype2, dt.atype3, dt.atype4, dt.gromacs['func'])
                if not name in dihedraltypes: dihedraltypes[name] = []
                dihedraltypes[name].append(dt)
        print("Build dihedraltypes dictionary with {0} entries".format(len(dihedraltypes)))

        #
        # Build improper dictionary
        #
        impropertypes = {}
        for it in top.impropertypes:
                it.disabled = True
                it.comment = "; LINE(%d) " % it.line
                name = "{0}-{1}-{2}-{3}-{4}".format(it.atype1, it.atype2, it.atype3, it.atype4, it.gromacs['func'])
                if not name in impropertypes: impropertypes[name] = []
                impropertypes[name].append(it)
        print("Build impropertypes dictionary with {0} entries".format(len(impropertypes)))

        if 'Protein' in top.dict_molname_mol:
                mol = top.dict_molname_mol['Protein']
                for at in mol.atoms: at.charge *= math.sqrt(args.scale_protein)
                mol = scale_dihedrals(mol, dihedraltypes, args.scale_protein, banned_lines)
                mol = scale_impropers(mol, impropertypes, 1.0, banned_lines)

        # Remove non-default moleculetypes
        for k in top.dict_molname_mol.keys():
                if k in ["Protein", "SOL", "Ion" ]: continue
                del top.dict_molname_mol[k]


        top.write(args.output)
