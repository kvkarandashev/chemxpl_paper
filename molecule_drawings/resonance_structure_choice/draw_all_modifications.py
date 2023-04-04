from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_modification_possibilities,
    draw_all_possible_resonance_structures,
    full_change_list,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.utils import mkdir
import os

chemgraph_string = "16#1@1:6@2:15@3@4:6#3:6#3"

highlight_radius = 0.4

kwargs = {
    "highlightAtomRadius": 0.4,
    "color_change_main": LIGHTRED,
    "color_change_minor": LIGHTBLUE,
    "color_change_special": LIGHTGREEN,
    "randomized_change_params": {
        "change_prob_dict": full_change_list,
        "possible_elements": ["S", "P", "C"],
        "added_bond_orders": [1, 2],
        "chain_addition_tuple_possibilities": False,
        "bond_order_changes": [-1, 1],
        "bond_order_valence_changes": [-2, 2],
        "max_fragment_num": 1,
        "added_bond_orders_val_change": [1, 2],
    },
}


rotation = {"rotated": 90, "standard": None}

for rotate_label in rotation.keys():
    mkdir(rotate_label)
    os.chdir(rotate_label)
    dump_dir = "mod_possibilities"
    mkdir(dump_dir)
    os.chdir(dump_dir)

    cg = str2ChemGraph(chemgraph_string)

    file_prefix = "test_mod_cg_"

    centreMoleculesBeforeDrawing = True
    padding = 0.1
    size = (300, 200)

    if rotate_label == "rotated":
        size = size[::-1]

    draw_all_modification_possibilities(
        cg,
        file_prefix,
        size=size,
        rotate=rotation[rotate_label],
        centreMoleculesBeforeDrawing=centreMoleculesBeforeDrawing,
        padding=padding,
        **kwargs
    )

    os.chdir("..")
    res_struct_dir = "resonance_structures"
    mkdir(res_struct_dir)
    os.chdir(res_struct_dir)
    draw_all_possible_resonance_structures(
        cg,
        "res_struct_",
        size=size,
        rotate=rotation[rotate_label],
        centreMoleculesBeforeDrawing=centreMoleculesBeforeDrawing,
        padding=padding,
        highlightAtomColors={3: LIGHTRED, 4: LIGHTRED, 2: LIGHTGREEN},
        highlightAtomRadius=highlight_radius,
        highlightBondTupleColors={(2, 3): LIGHTRED, (2, 4): LIGHTRED},
    )
    os.chdir("..")
    os.chdir("..")
