import os, sys, glob, subprocess
from bmapqml.utils import mkdir
from bmapqml.chemxpl.utils import SMILES_to_egc
import numpy as np
from misc_procedures import all_xyz_vals, extract_hist_size
from tabulate import tabulate


def average_overconv_cheap_deviation_norm(data_dir):
    xyzs = glob.glob(data_dir + "/best_candidates/best_candidate_*.xyz")
    prop_coeff = None
    cheap_val_label = None
    diff_values = []
    final_quant_values = []
    cheap_quant_values = []
    for xyz in xyzs:
        best_cand_data = all_xyz_vals(xyz)
        if best_cand_data is None:
            continue
        if cheap_val_label is None:
            for k in best_cand_data.keys():
                if "min_" in k:
                    cheap_val_label = k
                    break
        cheap_func_val = float(best_cand_data[cheap_val_label])
        overconv_quant_val = float(best_cand_data["overconv_quant_mean"])
        if prop_coeff is None:
            prop_coeff = overconv_quant_val / float(best_cand_data["final_minfunc_val"])
        cheap_quant_val = cheap_func_val * prop_coeff
        diff_values.append(abs(cheap_quant_val - overconv_quant_val))
        final_quant_values.append(overconv_quant_val)
        cheap_quant_values.append(cheap_quant_val)
    return np.mean(diff_values) / np.std(final_quant_values), np.mean(
        diff_values
    ) / np.std(cheap_quant_values)


parent_folder = sys.argv[1]

summary_dir = os.getcwd() + "/summary_" + os.path.basename(parent_folder)

mkdir(summary_dir)

os.chdir(summary_dir)

biases = ["none", "weak", "stronger"]

gap_constraints = ["weak", "strong"]

quantities = ["dipole", "solvation_energy"]

STD_RMSE_coeff = (
    0.25  # Since for the overconverged quantity we take the average over 16 attempts.
)

output = open(summary_dir + "/all_output.txt", "w")

for quantity in quantities:
    for gap_constraint in gap_constraints:
        print(
            "Quantity ",
            quantity,
            "constraint ",
            gap_constraint,
            " :",
            file=output,
        )
        #        print_ref_data_params(quantity=quantity, gap_constraint=gap_constraint, , file=output)
        headers = [
            "BIAS STRENGTH",
            "TOTAL BEST VALUE",
            "AGREEMENT",
            "TOTAL BEST SMILES",
            "MAX BEST RMSE",
            "AV BEST VALUE",
            "STDDEV ...",
            "AV REQ CONSIDERED TPS",
            "STDDEV ...",
            "AV TOT CONSIDERED TPS",
            "STDDEV ...",
            "NORM FINAL MIN NOISE",
            "NORM CHEAP ...",
        ]
        table_values = []
        for bias in biases:
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            folder = parent_folder + bulk_name
            best_cand_files = glob.glob(
                folder + "/data_*/best_candidates/best_candidate_0.xyz"
            )
            min_data_dir = None
            min_val = None
            min_SMILES = None
            all_SMILES = []
            max_RMSE = 0.0
            vals = []
            needed_considered_tps = []
            tot_considered_tps = []
            final_norm_noise_quants = []
            cheap_norm_noise_quants = []
            for data_dir in glob.glob(folder + "/data_*"):
                best_cand_file = data_dir + "/best_candidates/best_candidate_0.xyz"
                xyz_vals = all_xyz_vals(best_cand_file)

                cur_val = float(xyz_vals["overconv_quant_mean"])

                vals.append(cur_val)

                max_RMSE = max(
                    float(xyz_vals["overconv_quant_std"]) * STD_RMSE_coeff, max_RMSE
                )

                tot_considered_tps.append(extract_hist_size(data_dir))

                needed_considered_tps.append(float(xyz_vals["req_num_tps"]))

                (
                    final_noise_quant,
                    cheap_noise_quant,
                ) = average_overconv_cheap_deviation_norm(data_dir)
                final_norm_noise_quants.append(final_noise_quant)
                cheap_norm_noise_quants.append(cheap_noise_quant)

                SMILES = xyz_vals["SMILES"]
                all_SMILES.append(SMILES)
                if (min_val is None) or (abs(min_val) < abs(cur_val)):
                    min_val = cur_val
                    min_data_dir = data_dir
                    min_SMILES = SMILES
            if min_data_dir is None:
                print("Incomplete data for:", folder)
                quit()
            for ending in ["xyz", "png"]:
                subprocess.run(
                    [
                        "cp",
                        min_data_dir + "/best_candidates/best_candidate_0." + ending,
                        summary_dir + "/best" + bulk_name + "." + ending,
                    ]
                )
            egcs = [SMILES_to_egc(s) for s in all_SMILES]
            min_egc = SMILES_to_egc(min_SMILES)
            table_values.append(
                [
                    bias,
                    min_val,
                    sum(egc == min_egc for egc in egcs),
                    min_SMILES,
                    max_RMSE,
                    np.mean(vals),
                    np.std(vals),
                    np.mean(needed_considered_tps),
                    np.std(needed_considered_tps),
                    np.mean(tot_considered_tps),
                    np.std(tot_considered_tps),
                    np.mean(final_norm_noise_quants),
                    np.mean(cheap_norm_noise_quants),
                ]
            )
        print(tabulate(table_values, headers=["# " + s for s in headers]), file=output)

output.close()
