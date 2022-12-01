from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_chemgraph_to_file,
    LIGHTRED,
    LIGHTBLUE,
)

chemgraph_str = "7#1@1:15#1@2:7#1"

where_add = [0, 2]

post_added_bonds = {"pre_move": None, "post_move": [(*where_add, 1)]}
highlightBondTuples = {"pre_move": [], "post_move": [tuple(where_add)]}

cg = str2ChemGraph(chemgraph_str)
for image in ["pre_move", "post_move"]:
    draw_chemgraph_to_file(
        cg,
        "illustration_" + image + ".png",
        size=(300, 200),
        highlightAtoms=where_add,
        highlightAtomColor=LIGHTBLUE,
        highlightBondTuples=highlightBondTuples[image],
        highlightBondTupleColor=LIGHTRED,
        highlightAtomRadius=0.4,
        post_added_bonds=post_added_bonds[image],
    )

other_cg = str2ChemGraph("7@1@2:15#1@2:7")
draw_chemgraph_to_file(
    other_cg,
    "illustration_correct.png",
    size=(300, 200),
    highlightAtoms=where_add,
    highlightAtomColor=LIGHTBLUE,
    highlightBondTuples=highlightBondTuples["post_move"],
    highlightBondTupleColor=LIGHTRED,
    highlightAtomRadius=0.4,
)
