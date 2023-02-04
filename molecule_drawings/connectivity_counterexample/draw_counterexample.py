from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import draw_chemgraph_to_file

chemgraph_strs = {
    "egp_connectivity_counterexample": "6#3@1:15@2@3@4@5:6#2@6:6#2@6:6#2@6:6#2@6:15@7:6#3"
}

for regime, chemgraph_str in chemgraph_strs.items():
    cg = str2ChemGraph(chemgraph_str)

    draw_chemgraph_to_file(cg, regime + ".png", size=(300, 200), bw_palette=True)
