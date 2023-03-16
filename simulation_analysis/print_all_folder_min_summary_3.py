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

biases = ["none", "weak", "stronger"]

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

best = {"dipole": "max.", "solvation_energy": "min."}

latex_quantity_name = {"dipole": "\\dipole", "solvation_energy": "\\dEsolv"}


class EncLabel:
    def __init__(self, quant, gap_constr, bias):
        self.quant = quant
        self.gap_constr = gap_constr
        self.bias = bias

    def __gt__(self, other):
        if self.quant != other.quant:
            return quantities.index(self.quant) > quantities.index(other.quant)
        if self.gap_constr != other.gap_constr:
            return gap_constraints.index(self.gap_constr) > gap_constraints.index(
                other.gap_constr
            )
        if self.bias != other.bias:
            return biases.index(self.bias) > biases.index(other.bias)
        return False

    def __eq__(self, other):
        return (
            (self.quant == other.quant)
            and (self.gap_constr == other.gap_constr)
            and (self.bias == other.bias)
        )

    def __lt__(self, other):
        if self == other:
            return False
        return not (self > other)

    def __str__(self):
        return self.quant + "_" + self.gap_constr + "_" + self.bias

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)


def EncLaTeXLabel(el):
    return (
        best[el.quant]
        + " $"
        + latex_quantity_name[el.quant]
        + "$/"
        + el.gap_constr
        + " $\gap$ constr."
    )


class EncDict:
    def __init__(self):
        self.dict = {}

    def add(self, enc_label):
        if enc_label not in self.dict:
            self.dict[enc_label] = 0
        self.dict[enc_label] += 1

    def comp_list(self):
        return sorted([(h, -v) for h, v in self.dict.items()])

    def __gt__(self, other):
        cl1 = self.comp_list()
        cl2 = other.comp_list()
        for t1, t2 in zip(cl1, cl2):
            if t1 != t2:
                return t1 > t2
        return len(cl1) > len(cl2)

    def __str__(self):
        return str(self.comp_list())

    def __repr__(self):
        return str(self)


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
        for bias in biases:
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            bulk_label = EncLabel(quantity, gap_constraint, bias)
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
                    candidate_SMILES[SMILES] = EncDict()
                candidate_SMILES[SMILES].add(bulk_label)

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

better_candidate_ordering = list(range(len(best_SMILES)))

unordered_best = []
unordered_runnerup = []

runnerup_enc_ids = []

for S, enc_data in candidate_SMILES.items():
    if S in best_SMILES:
        unordered_best.append(S)
    else:
        unordered_runnerup.append(S)
        runnerup_enc_ids.append(enc_data)

runnerup_enc_ids = [(enc, i) for i, enc in enumerate(runnerup_enc_ids)]

runnerup_enc_ids.sort()
better_runnerup_ordering = [t[1] for t in runnerup_enc_ids]

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
    for prob, enc in enc_data.dict.items():
        print(prob, enc, file=runner_up_output)
    draw_all_possible_resonance_structures(
        SMILES_to_egc(S).chemgraph,
        image_name + "_",
        size=(200, 300),
        rotate=90,
        bw_palette=True,
    )
    draw_all_possible_resonance_structures(
        SMILES_to_egc(S).chemgraph,
        image_name + "_rot_",
        size=(300, 200),
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


runner_up_headers = [multrow(s, 2) for s in ["molecule", "SMILES", "opt. problem"]]


def rotatebox(s):
    #    return "\rotatebox[origin=l]{90}{"+bias+"}")
    return "\STAB{\rotatebox[origin=l]{90}{" + s + "}}"


diff_bias_headers = [
    (
        "\\tworowcell{num. enc.}{with bias:}",
        rotatebox(bias),
    )
    for bias in biases
]

all_runner_up_headers = runner_up_headers + diff_bias_headers

#            runner_up_headers.append(("opt. $"+latex_quantity_name[quantity]+"$", gap_constraint, bias[0]))

for row_id, s in enumerate(reordered_SMILES):
    cur_tuple = [0 for _ in biases]
    s_ll = None
    for bulk_lable, num in candidate_SMILES[s].dict.items():
        ll = EncLaTeXLabel(bulk_lable)
        if s_ll is None:
            s_ll = ll
        else:
            if s_ll != ll:
                raise Exception()
        cur_tuple[biases.index(bulk_lable.bias)] = num
    runner_up_table[row_id].append(ll)
    runner_up_table[row_id] += cur_tuple

table_to_latex(
    runner_up_table,
    all_runner_up_headers,
    "runner_up_table.tex",
    transpose=False,
    multicolumn=True,
    multirow=True,
)


table_to_latex(
    runner_up_table,
    runner_up_headers,
    "runner_up_table_no_last_columns.tex",
    transpose=False,
    multicolumn=True,
    multirow=True,
)
