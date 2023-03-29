import os, sys, glob, subprocess
from bmapqml.utils import mkdir
from bmapqml.chemxpl.utils import SMILES_to_egc
from bmapqml.chemxpl.rdkit_draw_utils import draw_all_possible_resonance_structures
import numpy as np
import pandas as pd
from copy import deepcopy
from misc_procedures import (
    all_xyz_vals,
    extract_hist_size,
    table_to_csv,
    table_to_latex,
)
from tabulate import tabulate

biases = ["none", "weak", "stronger"]

bias_strength = {"none": "0.0", "weak": "0.2", "stronger": "0.4"}

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

best = {"dipole": "max.", "solvation_energy": "min."}

latex_quantity_name = {"dipole": "\\dipole", "solvation_energy": "\\dEsolv"}

phantom = "\phantom{\_}"


class OptProbLabel:
    def __init__(self, quant, gap_constr):
        self.quant = quant
        self.gap_constr = gap_constr

    def __gt__(self, other):
        if self.quant != other.quant:
            return quantities.index(self.quant) > quantities.index(other.quant)
        if self.gap_constr != other.gap_constr:
            return gap_constraints.index(self.gap_constr) > gap_constraints.index(
                other.gap_constr
            )
        return False

    def __eq__(self, other):
        if other is None:
            return False
        return (self.quant == other.quant) and (self.gap_constr == other.gap_constr)

    def __lt__(self, other):
        if self == other:
            return False
        return not (self > other)

    def __str__(self):
        return self.quant + "_" + self.gap_constr

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    def tex_label(self):
        return (
            best[self.quant]
            + " $"
            + latex_quantity_name[self.quant]
            + "$/"
            + self.gap_constr
            + " $\gap$ constr."
        )

    def opt_name_row(self):
        return ["\\midrule\multicolumn{6}{c}{" + self.tex_label() + "}"] + [
            phantom for _ in range(5)
        ]


def better(val1, val2, quantity):
    if quantity == "solvation_energy":
        return val1 < val2
    if quantity == "dipole":
        return val1 > val2
    raise Exception


# For LaTeX conversion.
def preexp_power(n):
    s = "{:0.3e}".format(n)
    parts = s.split("e")
    return parts[0], int(parts[1])


def adjusted_error(n, power):
    npe, np = preexp_power(n)
    power_diff = power - np
    if power_diff < 0:
        raise Exception("something is wrong")
    if power_diff == 0:
        return npe
    numbers_split = npe.split(".")
    numbers = numbers_split[0] + numbers_split[1]
    if power_diff > 1:
        for _ in range(power_diff - 1):
            numbers = "0" + numbers
    numbers = "0." + numbers

    return numbers[:5]


class EncDict:
    def __init__(self):
        self.dict = {}
        for bias in biases:
            self.dict[bias] = 0
        self.opt_prob_label = None
        self.best_quant_est = None
        self.best_quant_est_RMSE = None

    def add(self, new_opt_prob, bias, quant_est, quant_est_RMSE):
        if self.opt_prob_label is None:
            self.opt_prob_label = new_opt_prob
        else:
            if self.opt_prob_label != new_opt_prob:
                raise Exception()
        self.dict[bias] += 1
        if self.best_quant_est is None:
            update_best_quant_est = True
        else:
            update_best_quant_est = better(
                quant_est, self.best_quant_est, self.opt_prob_label.quant
            )
        if update_best_quant_est:
            self.best_quant_est = quant_est
            self.best_quant_est_RMSE = quant_est_RMSE

    def __gt__(self, other):
        if self.opt_prob_label != other.opt_prob_label:
            return self.opt_prob_label > other.opt_prob_label
        return better(
            other.best_quant_est, self.best_quant_est, self.opt_prob_label.quant
        )

    def __repr__(self):
        return str(self)

    def __str__(self):
        return (
            str(self.opt_prob_label)
            + "_"
            + str(self.dict)
            + ":"
            + str(self.best_quant_est)
            + " pm "
            + str(self.best_quant_est_RMSE)
        )

    def latex_val(self, phantom_minus_alignment=False):
        est_mean_pe, est_mean_power = preexp_power(self.best_quant_est)
        est_err_num = adjusted_error(self.best_quant_est_RMSE, est_mean_power)
        output = est_mean_pe + " \pm " + est_err_num
        if est_mean_power == 0:
            output = "\phantom{(}" + output
        else:
            output = "(" + output + ")\cdot 10^{" + str(est_mean_power) + "}"
        if phantom_minus_alignment:
            if self.best_quant_est > 0.0:
                output = "\phantom{-}" + output
        return "$" + output + "$"


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
        opt_label = OptProbLabel(quantity, gap_constraint)
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

                cur_RMSE = float(xyz_vals["overconv_quant_std"]) * STD_RMSE_coeff

                vals.append(cur_val)

                max_RMSE = max(cur_RMSE, max_RMSE)

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
                if SMILES not in candidate_SMILES:
                    candidate_SMILES[SMILES] = EncDict()
                candidate_SMILES[SMILES].add(
                    OptProbLabel(quantity, gap_constraint), bias, cur_val, cur_RMSE
                )
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

            table_values.append(
                [
                    bias_strength[bias],
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
            "$\\biasprop$",
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

display_ordering = list(range(len(best_SMILES)))

unordered = []

enc_ids = []

for S, enc_data in candidate_SMILES.items():
    enc_ids.append(enc_data)

enc_ids = [(enc, i) for i, enc in enumerate(enc_ids)]

enc_ids.sort()
better_ordering = [t[1] for t in enc_ids]

ordered = [list(candidate_SMILES.keys())[i] for i in better_ordering]

image_label_prefix = {True: "C", False: "R"}

LaTeX_label_prefix = {True: "\\bestcand", False: "\runnerupmol"}


def image_label(i, is_best):
    return image_label_prefix[is_best] + str(i)


def LaTeX_label(i, is_best):
    return LaTeX_label_prefix[is_best] + "{" + str(i) + "}"


def multrow(s, nrows, add_cline=None):
    phantoms = [phantom for _ in range(nrows - 1)]
    if (add_cline is not None) and (nrows > 1):
        phantoms[-1] = "\cline{" + add_cline + "}" + phantoms[-1]
    return tuple(["\multirow{" + str(nrows) + "}{*}{" + s + "}"] + phantoms)


runner_up_headers = []
for i, h in enumerate(["molecule", "SMILES", "$\phantom{-(}$opt. quant."]):
    if i == 0:
        add_cline = "4-6"
    else:
        add_cline = None
    runner_up_headers.append(multrow(h, 2, add_cline=add_cline))


def rotatebox(s):
    #    return "\rotatebox[origin=l]{90}{"+bias+"}")
    return "\STAB{\rotatebox[origin=l]{90}{" + s + "}}"


# diff_bias_headers = [
#    (
#        "\\tworowcell{num. enc.}{with bias:}",
#        rotatebox(bias),
#    )
#    for bias in biases
# ]

diff_bias_headers = [("enc. with $\\biasprop$", bias_strength[bias]) for bias in biases]

all_runner_up_headers = runner_up_headers + diff_bias_headers

runner_up_table = []

individual_runner_up_tables = {}

best_id = 0
runner_up_id = 0

cur_opt_prob = None

max_colwidth = 0

for SMILES in ordered:
    is_best = SMILES in best_SMILES
    if is_best:
        best_id += 1
        cur_id = best_id
    else:
        runner_up_id += 1
        cur_id = runner_up_id
    image_name = image_label(cur_id, is_best)
    enc_data = candidate_SMILES[SMILES]
    cg = SMILES_to_egc(SMILES).chemgraph
    draw_all_possible_resonance_structures(
        cg,
        image_name + "_",
        size=(200, 300),
        rotate=90,
        bw_palette=True,
    )
    draw_all_possible_resonance_structures(
        cg,
        image_name + "_rot_",
        size=(300, 200),
        bw_palette=True,
    )

    if enc_data.opt_prob_label != cur_opt_prob:
        cur_opt_prob = enc_data.opt_prob_label
        extra_row = cur_opt_prob.opt_name_row()
        max_colwidth = max(max_colwidth, len(extra_row[0]))
        runner_up_table.append(cur_opt_prob.opt_name_row())
        individual_runner_up_tables[cur_opt_prob] = []

    checked_SMILES = ""
    for c in SMILES:
        if c == "#":
            checked_SMILES += "\\" + c
        else:
            checked_SMILES += c

    cur_enc_data = candidate_SMILES[SMILES]

    new_row = [
        LaTeX_label(cur_id, is_best),
        checked_SMILES,
        cur_enc_data.latex_val(phantom_minus_alignment=True),
    ]
    new_row1 = [
        LaTeX_label(cur_id, is_best),
        checked_SMILES,
        cur_enc_data.latex_val(phantom_minus_alignment=True),
    ]

    for bias in biases:
        new_row.append(cur_enc_data.dict[bias])
        new_row1.append(cur_enc_data.dict[bias])

    runner_up_table.append(new_row)
    individual_runner_up_tables[cur_opt_prob].append(new_row1)

with pd.option_context("max_colwidth", max_colwidth):
    table_to_latex(
        runner_up_table,
        all_runner_up_headers,
        "runner_up_table.tex",
        transpose=False,
        multicolumn=True,
        multirow=True,
    )
    for opt_prob, tab in individual_runner_up_tables.items():
        cur_runner_up_table = deepcopy(all_runner_up_headers)
        latex_opt = latex_quantity_name[opt_prob.quant] + "^{\\mathrm{conv.}}"
        #        if opt_prob.quant == "solvation_energy":
        #            latex_opt = "\phantom{-}" + latex_opt
        cur_runner_up_table[2] = multrow("$\phantom{-(}" + latex_opt + "$", 2)
        table_to_latex(
            tab,
            cur_runner_up_table,
            "runner_up_table_" + str(opt_prob) + ".tex",
            transpose=False,
            multicolumn=True,
            multirow=True,
        )
        ruStxt = "runner_up_SMILES_" + str(opt_prob) + ".txt"
        with open(ruStxt, "w") as SMILES_txt:
            for row in tab[::-1]:
                print(row[1], file=SMILES_txt)
        subprocess.run(["sed", "-i", "s/\\\\#/#/g", ruStxt])

subprocess.run(["../../final_summary_postproc.sh", "runner_up_table.tex"])
