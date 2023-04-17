from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_utils import SMILES_to_egc
from bmapqml.chemxpl.cross_coupling import FragmentPair
from bmapqml.chemxpl.rdkit_draw_utils import (
    draw_all_possible_fragment_pairs,
    draw_all_possible_resonance_structures,
    LIGHTBLUE,
    LIGHTRED,
    LIGHTGREEN,
)
import os
from bmapqml.utils import mkdir

# chemgraph_strings = [
#    "6#3@1:6@2:7@3:6#1@4:6@1@5:9",
#    "9@1:6@2@3:6#1@4:7@5:6#1@6:6#1@6:6@7:8#1",
#    "9@1:6@2:6@3:8#1",
#    ]

# origin_sizes=[(5, 1), (0, 6)]

SMILES = ["N1=C(C)C(F)=C1", "C1=C(O)C=CC(F)=N1", "FC#CO", "C1=C(F)N=CC=NC(C)=C1"]

kwargs = {
    "nhatoms_range": [1, 20],
    "highlightAtomRadius": 0.4,
    "bw_palette": True,
}

size = (300, 200)
rotate = None

for i, s in enumerate(SMILES):
    dump_dir = "all_frags_" + str(i)
    cg = SMILES_to_egc(s).chemgraph  # str2ChemGraph(s)
    mkdir(dump_dir)
    os.chdir(dump_dir)
    draw_all_possible_fragment_pairs(
        cg, size=size, rotate=rotate, filename_prefix="fragment_pair_" + str(i) + "_"
    )
    os.chdir("..")
