import gromacs as gmx
from unittest.mock import patch, mock_open
import enum


class AnglesFuncT(enum.Enum):
    Angle = 1
    G96_angle = 2
    cross_bond_bond = 3
    cross_bond_angle = 4
    urey_bradley = 5
    quartic_angle = 6
    # there is no funct = 7
    tabulated_angle = 8
    linear_angle = 9
    restricted_bending_potential = 10


angles_types = {
    AnglesFuncT.Angle: {
        "theta0",
    }
}


def create_topology_data() -> str:
    """
    https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#tab-topfile2
    """
    return """
[ moleculetype ]
; Name            nrexcl
Example     3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    cha
; residue   1 ASN rtp ASN  q +1.0
     1        NH3      1    ASN      N      1       -0.3     14.007   ; qtot -0.3
     2         HC      1    ASN     H1      2       0.33      1.008   ; qtot 0.03
     3         HC      1    ASN     H2      3       0.33      1.008   ; qtot 0.36

[ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3
    2     1     3     5 

[ molecules ]
; Compound        #mols
Example	1
"""


@patch("builtins.open", mock_open(read_data=create_topology_data()))
def test_angles():
    topol = gmx.fileformats.top.TOP("topol.top")
    topol.molecules[0].angles[0].gromacs
