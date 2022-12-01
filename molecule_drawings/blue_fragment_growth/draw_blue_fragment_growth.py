from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_possible_resonance_structures,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit
from bmapqml.chemxpl.modify import FragmentPair
import numpy as np
import random

random.seed(1)
np.random.seed(1)

chemgraph_strings = [
    "6#3@1:6@2:7@3:6#1@4:6@1@5:9",
    "9@1:6@2@3:6#1@4:6#1@5:6#1@5:7@6:6#2@7:6#3",
]

kwargs = {
    "size": (300, 200),
    "highlightAtomRadius": 0.4,
    "bw_palette": True,
    "abbrevs": None,
}

added_member_atoms_list = [[[0], [1], [2]], [[4], [5, 2]]]

for i_cg, chemgraph_string in enumerate(chemgraph_strings):
    cg = str2ChemGraph(chemgraph_string)

    membership_vector = np.zeros((cg.nhatoms(),), dtype=int)

    for i, added_member_atoms in enumerate(added_member_atoms_list[i_cg]):

        for ama in added_member_atoms:
            membership_vector[ama] = 1

        file_prefix_cg = "resonance_" + str(i_cg) + "_" + str(i) + "_"

        fp = FragmentPair(cg, membership_vector)

        draw_all_possible_resonance_structures(fp, file_prefix_cg, **kwargs)
