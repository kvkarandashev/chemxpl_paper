from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.rdkit_draw_utils import draw_chemgraph_to_file

chemgraph_strs = {
    "qm9_weak_dipole": "8@7:7#2@8:7@5@6@8:7@6@7:7@7@8:6#3:6#1:6:6",
    "qm9_weak_solvation_energy": "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
    "qm9_strong_dipole": "8@7:8@8:7#1@7@8:7#1@5@7:7#1@6@8:6#2@6:6#2:6:6",
    "qm9_strong_solvation_energy": "8@7:8@8:7#1@7@8:7#1@5@7:7#1@6@8:6#2@6:6#2:6:6",
    "egp_weak_dipole": "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
    "egp_weak_solvation_energy": "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
    "egp_strong_dipole": "16@2@4@5@8:16@3@6@7@9:9:9:8:8:8:8:6#2@9:6#2",
    "egp_strong_solvation_energy": "16@2@3@6@11:16@4@5@7@12:8:8:8:8:7#2:7#2:6#2@10@11:6#2@10@12:6#2:6#2:6#2",
}

for regime, chemgraph_str in chemgraph_strs.items():
    cg = str2ChemGraph(chemgraph_str)

    draw_chemgraph_to_file(cg, regime + ".png", size=(300, 200), bw_palette=True)
