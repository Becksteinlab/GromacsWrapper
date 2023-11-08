import gromacs as gmx
from unittest.mock import patch, mock_open
import random
import pytest
from gromacs.fileformats import blocks


def generate_topology_line(func_type: blocks.AngleFunctionType):
    # Generates random parameters for the different function types
    def random_params(num):
        return [random.uniform(0.0, 100.0) for _ in range(num)]

    # Define the pattern for each function type
    patterns = {
        blocks.AngleFunctionType.HARMONIC: "{:5d}{:5d}{:5d}    1{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.G96_ANGLE: "{:5d}{:5d}{:5d}    2{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.CROSS_BOND_BOND: "{:5d}{:5d}{:5d}    3{:10.4f}{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.CROSS_BOND_ANGLE: "{:5d}{:5d}{:5d}    4{:10.4f}{:10.4f}{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.UREY_BRADLEY: "{:5d}{:5d}{:5d}    5{:10.4f}{:10.4f}{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.QUARTIC_ANGLE: "{:5d}{:5d}{:5d}    6{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.TABULATED_ANGLE: "{:5d}{:5d}{:5d}    8{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.LINEAR_ANGLE: "{:5d}{:5d}{:5d}    9{:10.4f}{:10.4f}",
        blocks.AngleFunctionType.RESTRICTED_BENDING: "{:5d}{:5d}{:5d}   10{:10.4f}{:10.4f}",
    }

    atoms = [1, 2, 3]
    params = random_params(func_type.num_params)
    line = patterns[func_type].format(*atoms, *params)
    return line


def create_topology_data(func_type: blocks.AngleFunctionType):
    """
    https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#tab-topfile2
    """
    line = generate_topology_line(func_type)
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
{line}

[ molecules ]
; Compound        #mols
Example	1
""".format(
        line=line
    )


@pytest.mark.parametrize("func_type", blocks.AngleFunctionType)
def test_angles(func_type: blocks.AngleFunctionType):
    with patch(
        "builtins.open", mock_open(read_data=create_topology_data(func_type))
    ) as mock:
        topol = gmx.fileformats.top.TOP("topol.top")
        [molecule] = topol.molecules
        [angle] = molecule.angles
        assert len(angle.gromacs["param"].keys()) == func_type.num_params
