from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_possible_fragment_pairs,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.utils import mkdir
import os
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

for i, chemgraph_str in enumerate(chemgraph_strings):
    cg=str2ChemGraph(chemgraph_str)
    output_dir="frags_"+str(i)
    mkdir(output_dir)
    os.chdir(output_dir)
    draw_all_possible_fragment_pairs(cg)
    os.chdir("..")