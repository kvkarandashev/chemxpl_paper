from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_cross_couplings,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
import os
from bmapqml.utils import mkdir

chemgraph_strings = [
    "6#3@1:6@2:7@3:6#1@4:6@1@5:9",
    "9@1:6@2@3:6#1@4:6#1@5:6#1@5:7@6:6#2@7:6#3",
]

kwargs = {
    "nhatoms_range" : [1, 9],
    "size": (300, 200),
    "highlightAtomRadius": 0.4,
    "bw_palette": True,
}

cgs = [str2ChemGraph(cg_str) for cg_str in chemgraph_strings]

folder_name = "cross_couplings"
mkdir(folder_name)
os.chdir(folder_name)

draw_all_cross_couplings(cgs, fragment_ratio_range=[.3, .7], **kwargs)
os.chdir("..")
