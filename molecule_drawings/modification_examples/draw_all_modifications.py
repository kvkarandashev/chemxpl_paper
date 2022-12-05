from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_modification_possibilities,
    full_change_list,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.utils import mkdir
import numpy as np
import os

chemgraph_strings = [
    "6#3@1:8@2:6#3",
    "6#3@1:6@2:7@3:6#1@4:6@1@5:9",
    "6#3@1:15#1@2:6#3",
    "6#3@1:15@2@3:6#3:9",
    "6#3@1:6#1@2@3:6#3:6#3",
]

kwargs = {
    "size": (300, 200),
    "highlightAtomRadius": 0.4,
    "color_change": LIGHTRED,
    "color_change_neighbors": LIGHTBLUE,
    "randomized_change_params": {
        "change_prob_dict": full_change_list,
        "possible_elements": ["B", "C", "O"],
        "added_bond_orders": [1, 2],
        "chain_addition_tuple_possibilities": False,
        "bond_order_changes": [-1, 1],
        "bond_order_valence_changes": [-2, 2],
        "max_fragment_num": 1,
        "added_bond_orders_val_change": [1, 2],
    },
}

for str_id, cg_str in enumerate(chemgraph_strings):
    dump_dir = "mod_possibilities_" + str(str_id)
    mkdir("mod_possibilities_" + str(str_id))
    os.chdir(dump_dir)

    cg = str2ChemGraph(cg_str)

    file_prefix = "test_mod_cg_" + str(str_id) + "_"

    draw_all_modification_possibilities(cg, file_prefix, **kwargs)

    os.chdir("..")
