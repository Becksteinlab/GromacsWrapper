
from gromacs.scaling import partial_tempering
import numpy as np
import math
import copy, argparse

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("--scale_protein", type=float)
	parser.add_argument("--scale_lipids", type=float)
	parser.add_argument("--scale_protein_lipids", type=float)
	parser.add_argument("input")
	parser.add_argument("output")
	parser.add_argument("--banned_lines", default="")
	return parser.parse_args()

args = parse_args()
partial_tempering(args)