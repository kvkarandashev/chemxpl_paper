from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_successful_random_coupling,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit
import numpy as np
import random, os
from bmapqml.utils import mkdir

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
}

cgs = [str2ChemGraph(cg_str) for cg_str in chemgraph_strings]

for i in range(8):
    folder_name = "cross_couplings/" + str(i)
    mkdir(folder_name)
    os.chdir(folder_name)

    file_prefix_cg = "resonance_" + str(i) + "_"

    draw_successful_random_coupling(cgs, (3, 3), **kwargs)
    os.chdir("../..")
