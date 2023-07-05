from bmapqml.chemxpl.random_walk import (
    randomized_change,
    randomized_cross_coupling,
    TrajectoryPoint,
    ChemGraph,
    add_heavy_atom_chain,
    replace_heavy_atom,
    remove_heavy_atom,
    change_bond_order,
)
from bmapqml.chemxpl.rdkit_utils import SMILES_to_egc
from bmapqml.chemxpl.minimized_functions.morfeus_quantity_estimates import (
    morfeus_FF_xTB_code_quants,
)
from bmapqml.utils import loadpkl, mkdir
import random, sys
import numpy as np


def solv_en(tp):
    if isinstance(tp, ChemGraph):
        true_tp = TrajectoryPoint(cg=tp)
    else:
        true_tp = tp
    quant = "solvation_energy"
    true_tp.calculated_data = {}
    output = morfeus_FF_xTB_code_quants(
        true_tp,
        num_conformers=16,
        num_attempts=16,
        ff_type="MMFF94",
        quantities=[quant],
        solvent="water",
        remaining_rho=0.9,
    )["mean"][quant]
    return output


SMILES_txt = sys.argv[1]

SMILES_input = open(SMILES_txt, "r")

SMILES_lines = SMILES_input.readlines()

SMILES_input.close()

for s in SMILES_lines:
    try:
        print(int(s))
    except ValueError:
        egc = SMILES_to_egc(s)
        tp = TrajectoryPoint(egc=egc)
        print(s[:-1], solv_en(tp))
