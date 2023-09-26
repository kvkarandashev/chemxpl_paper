import os, glob, subprocess


def all_xyz_vals(xyz_file):
    i = open(xyz_file, "r")
    lines = i.readlines()
    if len(lines) == 1:
        return None
    l = lines[1].split()
    output = {}
    for q in l:
        q_spl = q.split("=")
        quant_name = q_spl[0]
        output[quant_name] = q.split(quant_name + "=")[1]
    return output


folder = "/store/konst/chemxpl_related/minimization_runs_xTB_dipole_solvation_cheap_post_sub/qm9_randomized_init_none_strong_solvation_energy"

candidate_SMILES = {}

best_cand_files = glob.glob(folder + "/data_*/best_candidates/best_candidate_0.xyz")

for best_cand_file in best_cand_files:
    xyz_vals = all_xyz_vals(best_cand_file)
    cur_SMILES = xyz_vals["SMILES"]
    if cur_SMILES not in candidate_SMILES:
        candidate_SMILES[cur_SMILES] = 0
    candidate_SMILES[cur_SMILES] += 1

print("Best candidate SMILES:", candidate_SMILES)
