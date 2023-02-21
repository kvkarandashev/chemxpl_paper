import os, sys, glob, subprocess
from bmapqml.utils import mkdir
from bmapqml.chemxpl.utils import SMILES_to_egc
from bmapqml.chemxpl.rdkit_draw_utils import draw_all_possible_resonance_structures
import numpy as np
from misc_procedures import (
    all_xyz_vals,
    extract_hist_size,
    table_to_csv,
    table_to_latex,
)
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

latex_quantity_name = {"dipole": "\\dipole", "solvation_energy": "\\dEsolv"}

best = {"dipole": "max.", "solvation_energy": "min."}
# Since for the overconverged quantity we take the average over 16 attempts.
STD_RMSE_coeff = 0.25

output = open(summary_dir + "/all_output.txt", "w")

candidate_SMILES = {}

best_SMILES = []

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
            "TOTAL BEST SMILES",
            "AGREEMENT",
            "MAX BEST RMSE",
            "AV BEST VALUE",
            "STDDEV ...",
            "AV REQ CONSIDERED TPS",
            "STDDEV ...",
            "AV TOT CONSIDERED TPS",
            "STDDEV ...",
            "AV STEP BEST FOUND",
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
            global_step_first_encounters = []
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

                global_step_first_encounters.append(
                    float(xyz_vals["first_global_MC_step_encounter"])
                )

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
            # Also replot SMILES in a format fitting for the paper.
            png_dir = "best" + bulk_name + "_pngs"
            mkdir(png_dir)
            os.chdir(png_dir)
            draw_all_possible_resonance_structures(
                min_egc.chemgraph,
                "best" + bulk_name + "_rot_",
                size=(200, 300),
                rotate=90,
                bw_palette=True,
            )
            draw_all_possible_resonance_structures(
                min_egc.chemgraph,
                "best" + bulk_name + "_",
                size=(300, 200),
                rotate=None,
                bw_palette=True,
            )
            os.chdir("..")
            best_SMILES.append(min_SMILES)
            for SMILES in all_SMILES:
                if SMILES not in candidate_SMILES:
                    candidate_SMILES[SMILES] = {}
                if bulk_name not in candidate_SMILES[SMILES]:
                    candidate_SMILES[SMILES][bulk_name] = 0
                candidate_SMILES[SMILES][bulk_name] += 1

            table_values.append(
                [
                    bias,
                    min_val,
                    min_SMILES,
                    sum(egc == min_egc for egc in egcs),
                    max_RMSE,
                    np.mean(vals),
                    np.std(vals),
                    np.mean(needed_considered_tps),
                    np.std(needed_considered_tps),
                    np.mean(tot_considered_tps),
                    np.std(tot_considered_tps),
                    np.mean(global_step_first_encounters),
                    np.std(global_step_first_encounters),
                    np.mean(final_norm_noise_quants),
                    np.mean(cheap_norm_noise_quants),
                ]
            )
        print(tabulate(table_values, headers=["# " + s for s in headers]), file=output)
        summary_file_prefix = (
            summary_dir
            + "/summary_gap_constraint_"
            + gap_constraint
            + "_quant_"
            + quantity
        )
        table_to_csv(
            table_values,
            headers,
            summary_file_prefix + ".csv",
        )
        fin_res_label = latex_quantity_name[quantity] + "^{\mathrm{best}}"
        latex_headers = [
            "bias strength",
            best[quantity] + " $" + fin_res_label + "$, a.u.",
            None,
            "agreement",
            "max. RMSE ($" + fin_res_label + "$), a.u.",
            "$\overline{" + fin_res_label + "}$, a.u.",
            "$\sigma(" + fin_res_label + ")$, a.u.",
            "$\overline{\tpreq}$",
            "$\sigma(\tpreq)$",
            "$\overline{\tottp}$",
            "$\sigma(\tottp)$",
            "$\overline{\\beststepfound}$",
            "$\sigma(\\beststepfound)$",
            "$\cheapquantnoise$",
            None,
        ]
        for transpose in [True, False]:
            if transpose:
                ending = "_tr.tex"
            else:
                ending = ".tex"
            table_to_latex(
                table_values,
                latex_headers,
                summary_file_prefix + ending,
                transpose=transpose,
            )


output.close()

runner_ups = "runner_ups"

mkdir(runner_ups)
os.chdir(runner_ups)
runner_up_output = open("runner_ups.txt", "w")
better_runnerup_ordering = [1, 2, 4, 0, 3, 5]
better_candidate_ordering = [2, 3, 0, 1]

unordered_best = []
unordered_runnerup = []

candS = list(candidate_SMILES.keys())
for S in candS:
    if S in best_SMILES:
        unordered_best.append(S)
    else:
        unordered_runnerup.append(S)

ordered_best = [unordered_best[i] for i in better_candidate_ordering]
ordered_runnerup = [unordered_runnerup[i] for i in better_runnerup_ordering]

image_label_prefix = {True: "C", False: "R"}

LaTeX_label_prefix = {True: "\\bestcand", False: "\runnerupmol"}

all_printed_SMILES = []

for is_best, l in zip([True, False], [ordered_best, ordered_runnerup]):
    all_printed_SMILES += [(i, S, is_best) for i, S in enumerate(l)]


def image_label(i, is_best):
    return image_label_prefix[is_best] + str(i + 1)


def LaTeX_label(i, is_best):
    return LaTeX_label_prefix[is_best] + "{" + str(i + 1) + "}"


runner_up_table = []

reordered_SMILES = []

for i, S, is_best in all_printed_SMILES:
    image_name = image_label(i, is_best)
    enc_data = candidate_SMILES[S]
    print(image_name, S, ":", file=runner_up_output)
    for prob, enc in enc_data.items():
        print(prob, enc, file=runner_up_output)
    draw_all_possible_resonance_structures(
        SMILES_to_egc(S).chemgraph,
        image_name + "_",
        size=(200, 300),
        rotate=90,
        bw_palette=True,
    )
    checked_SMILES = ""
    for c in S:
        if c == "#":
            checked_SMILES += "\\" + c
        else:
            checked_SMILES += c
    runner_up_table.append([LaTeX_label(i, is_best), checked_SMILES])
    reordered_SMILES.append(S)

runner_up_output.close()

phantom = "\phantom{\_}"

# def brute_multcol(s, ncols):
#    return "\multicolumn{"+str(ncols)+"}"+"{l}{"+s+"}"

# angle_header=("optimized quantity", "gap constraint", "biasing")

# runner_up_headers=[angle_header, tuple(phantom for _ in range(len(angle_header)))]


def multrow(s, nrows):
    return tuple(
        ["\multirow{" + str(nrows) + "}{*}{" + s + "}"]
        + [phantom for _ in range(nrows - 1)]
    )


runner_up_headers = [multrow(s, 3) for s in ["molecule", "SMILES"]]

bulk_names = []


def rotatebox(s):
    #    return "\rotatebox[origin=l]{90}{"+bias+"}")
    return "\STAB{\rotatebox[origin=l]{90}{" + s + "}}"


for quantity in quantities[::-1]:
    for gap_constraint in gap_constraints:
        for bias in biases:
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            bulk_names.append(bulk_name)
            runner_up_headers.append(
                (
                    "opt. $" + latex_quantity_name[quantity] + "$",
                    gap_constraint,
                    rotatebox(bias),
                )
            )
#            runner_up_headers.append(("opt. $"+latex_quantity_name[quantity]+"$", gap_constraint, bias[0]))

for row_id, s in enumerate(reordered_SMILES):
    for bulk_name in bulk_names:
        num_dict = candidate_SMILES[s]
        if bulk_name in num_dict:
            num = num_dict[bulk_name]
        else:
            num = 0
        runner_up_table[row_id].append(num)


table_to_latex(
    runner_up_table,
    runner_up_headers,
    "runner_up_table.tex",
    transpose=False,
    multicolumn=True,
    multirow=True,
)
