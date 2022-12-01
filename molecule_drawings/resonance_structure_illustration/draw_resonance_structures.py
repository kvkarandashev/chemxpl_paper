from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_possible_resonance_structures,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit
import numpy as np
import random

random.seed(1)
np.random.seed(1)

chemgraph_strings = [
    "6#3@1:6@2:7@3:6#1@4:6@1@5:9",
    "6@1@5@6@7:16@2:6@3:16@4:6@8@9@10:9:9:9:9:9:9",
]

kwargs = {
    "size": (300, 200),
    "highlightAtomRadius": 0.4,
    "bw_palette": True,
    "abbrevs": True,
}

for str_id, cg_str in enumerate(chemgraph_strings):
    cg = str2ChemGraph(cg_str)

    SMILES = chemgraph_to_canonical_rdkit(cg, SMILES_only=True)

    print(SMILES)

    file_prefix_cg = "resonance_" + str(str_id) + "_"

    draw_all_possible_resonance_structures(cg, file_prefix_cg, **kwargs)
