import os, sys, glob
from bmapqml.utils import mkdir
import numpy as np
import pandas as pd
from misc_procedures import (
    all_xyz_vals,
    quant_STD,
    brute_table_write,
    best_ref_vals,
)

biases = ["none", "weak", "stronger"]

bias_strength = {"none": "0.0", "weak": "0.2", "stronger": "0.4"}

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

best = {"dipole": "max.", "solvation_energy": "min."}

latex_quantity_name = {"dipole": "\\dipole", "solvation_energy": "\\dEsolv"}

quant_conversion_coefficients = {"solvation_energy": 2625.5002, "dipole": 2.541746}

quant_measure_name = {"solvation_energy": "kJ/mol", "dipole": "debye"}

phantom = "\phantom{\_}"

parent_folder = sys.argv[1]

dataset_name = sys.argv[2]

if dataset_name == "QM9":
    improvement_numbers_to_the_left = 1
    improvement_numbers_to_the_right = 3
    improvement_RMSE_numbers_to_the_left = 1
    improvement_RMSE_numbers_to_the_right = 3
else:
    improvement_numbers_to_the_left = 3
    improvement_numbers_to_the_right = 2
    improvement_RMSE_numbers_to_the_left = 1
    improvement_RMSE_numbers_to_the_right = 2

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


def latex_single_number_form(val, phantom_minus_alignment=False):
    pe, power = preexp_power(val)
    output = pe
    if power != 0:
        output += "\cdot 10^{" + str(power) + "}"
    if phantom_minus_alignment:
        if val > 0.0:
            output = "\phantom{-}" + output
        output = "\phantom{(}" + output
    return output


def phantom_zero_pad(n):
    s = ""
    for _ in range(abs(n)):
        s += "\phantom{0}"
    return s


def phantom_zero_padded_value_str(
    val_str, numbers_to_the_left=1, numbers_to_the_right=3
):
    val_str_split = val_str.split(".")
    cur_to_the_left = len(val_str_split[0])
    cur_to_the_right = len(val_str_split[1])
    return (
        phantom_zero_pad(numbers_to_the_left - cur_to_the_left)
        + val_str
        + phantom_zero_pad(numbers_to_the_right - cur_to_the_right)
    )


def special_latex_form(
    val,
    RMSE,
    numbers_to_the_left=1,
    numbers_to_the_right=3,
    RMSE_numbers_to_the_left=1,
    RMSE_numbers_to_the_right=3,
):
    _, power = preexp_power(val)
    form_str = (
        "{:."
        + str(numbers_to_the_right + numbers_to_the_left - 2 - max(0, power))
        + "f}"
    )

    val_str = form_str.format(val)
    val_str = phantom_zero_padded_value_str(
        val_str,
        numbers_to_the_left=numbers_to_the_left,
        numbers_to_the_right=numbers_to_the_right,
    )

    RMSE_str = form_str.format(RMSE)
    RMSE_str = phantom_zero_padded_value_str(
        RMSE_str,
        numbers_to_the_left=RMSE_numbers_to_the_left,
        numbers_to_the_right=RMSE_numbers_to_the_right,
    )

    return "$" + val_str + " \pm " + RMSE_str + "$"


def latex_single_number_winline(val, phantom_minus_alignment=False):
    return (
        "$"
        + latex_single_number_form(val, phantom_minus_alignment=phantom_minus_alignment)
        + "$"
    )


def latex_form(mean, RMSE, phantom_minus_alignment=False, quant=None):
    if quant is not None:
        conversion_coefficient = quant_conversion_coefficients[quant]
        used_mean = mean * conversion_coefficient
        used_RMSE = RMSE * conversion_coefficient
    else:
        used_mean = mean
        used_RMSE = RMSE
    est_mean_pe, est_mean_power = preexp_power(used_mean)
    est_err_num = adjusted_error(used_RMSE, est_mean_power)
    output = est_mean_pe + " \pm " + est_err_num
    if est_mean_power == 0:
        output = "\phantom{(}" + output
    else:
        output = "(" + output + ")\cdot 10^{" + str(est_mean_power) + "}"
    if phantom_minus_alignment:
        if used_mean > 0.0:
            output = "\phantom{-}" + output
    return "$" + output + "$"


def opt_prob_quant_str(quant):
    return best[quant] + " $" + latex_quantity_name[quant] + "$"


def opt_prob_str(quant, gap_constr):
    return opt_prob_quant_str(quant) + " (" + gap_constr + " $\gap$ constr.)"


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
        return opt_prob_str(self.quant, self.gap_constr)

    def best_val_tuple(self):
        return best_ref_vals[dataset_name][self.gap_constr][self.quant]

    def STD(self):
        return quant_STD[dataset_name][self.gap_constr][self.quant]

    def opt_name_row(self, ncols):
        if self.quant == "solvation_energy" and self.gap_constr == "weak":
            midrule_str = ""
        else:
            midrule_str = "\\midrule"
        return [
            midrule_str
            + "\multicolumn{"
            + str(ncols)
            + "}{c}{"
            + self.tex_label()
            + "}"
        ] + [phantom for _ in range(ncols - 1)]


def better(val1, val2, quantity):
    if quantity == "solvation_energy":
        return val1 < val2
    if quantity == "dipole":
        return val1 > val2
    raise Exception


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
        return latex_form(
            self.best_quant_est,
            self.best_quant_est_RMSE,
            phantom_minus_alignment=phantom_minus_alignment,
            quant=self.opt_prob_label.quant,
        )

    def latex_improve_val(self):
        dataset_best_tuple = self.opt_prob_label.best_val_tuple()
        dataset_quant_STD = self.opt_prob_label.STD()
        dataset_best_RMSE = dataset_best_tuple[1] * 0.25
        rel_improvement_RMSE = (
            np.sqrt(self.best_quant_est_RMSE**2 + dataset_best_RMSE**2)
            / dataset_quant_STD
        )
        rel_improvement = (
            abs(self.best_quant_est - dataset_best_tuple[0]) / dataset_quant_STD
        )
        return special_latex_form(
            rel_improvement,
            rel_improvement_RMSE,
            numbers_to_the_left=improvement_numbers_to_the_left,
            numbers_to_the_right=improvement_numbers_to_the_right,
            RMSE_numbers_to_the_left=improvement_RMSE_numbers_to_the_left,
            RMSE_numbers_to_the_right=improvement_RMSE_numbers_to_the_right,
        )


summary_dir = os.getcwd() + "/summary_post_sub_" + dataset_name

mkdir(summary_dir)

os.chdir(summary_dir)

# Since for the overconverged quantity we take the average over 16 attempts.
STD_RMSE_coeff = 0.25


def step_comp_index(gap_constraint, bias):
    return (gap_constraint, bias_strength[bias])


def bias_strength_adj(bias):
    return "\phantom{(}$\phantom{-}$" + bias_strength[bias]


del_str = "delete"


def centered_cell(s):
    return r"\multicolumn{1}{c}{" + s + "}"


opt_quant_val_label = "opt. quant. val."
rel_improv_label = "rel. improv."
best_label = centered_cell("best")
worst_label = centered_cell("worst")


def multrow2(string):
    return (r"\multirow{2}{*}{" + string + "}", r"\phantom{\_}")


# opt_quant_label=list(multrow2("opt. quant."))
# opt_quant_label[1]="\cline{3-4}\cline{5-6}"+opt_quant_label[1]
# opt_quant_label=tuple(opt_quant_label)
# gap_constr_label=multrow2(r"$\Delta\epsilon$ constr.")

opt_quant_label = (
    centered_cell("opt."),
    "\cline{3-4}\cline{6-7}" + centered_cell("quant."),
)
gap_constr_label = (centered_cell(r"$\Delta\epsilon$"), centered_cell("constr."))

minmax_candidate_table = {
    opt_quant_label: [],
    gap_constr_label: [],
    (opt_quant_val_label, best_label): [],
    (opt_quant_val_label, worst_label): [],
    ("\phantom{.}", ""): ["" for _ in range(4)],
    (rel_improv_label, best_label): [],
    (rel_improv_label, worst_label): [],
}


for quantity in quantities:
    minmax_candidate_table[opt_quant_label] += list(
        multrow2("$" + latex_quantity_name[quantity] + "$")
    )

    if quantity == "dipole":
        minmax_candidate_table[opt_quant_label][-2] = (
            "\cline{1-2}" + minmax_candidate_table[opt_quant_label][-2]
        )
    for gap_constraint in gap_constraints:
        minmax_candidate_table[gap_constr_label].append(gap_constraint)

        candidate_SMILES = {}
        opt_label = OptProbLabel(quantity, gap_constraint)
        # Find best and worst molecule over all bias values.
        for bias in biases:
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            folder = parent_folder + bulk_name
            best_cand_files = glob.glob(
                folder + "/data_*/best_candidates/best_candidate_0.xyz"
            )
            max_RMSE = 0.0
            for data_dir in glob.glob(folder + "/data_*"):
                best_cand_file = data_dir + "/best_candidates/best_candidate_0.xyz"
                xyz_vals = all_xyz_vals(best_cand_file)

                cur_val = float(xyz_vals["overconv_quant_mean"])

                cur_RMSE = float(xyz_vals["overconv_quant_std"]) * STD_RMSE_coeff

                max_RMSE = max(cur_RMSE, max_RMSE)

                SMILES = xyz_vals["SMILES"]
                if SMILES not in candidate_SMILES:
                    candidate_SMILES[SMILES] = EncDict()
                candidate_SMILES[SMILES].add(
                    OptProbLabel(quantity, gap_constraint), bias, cur_val, cur_RMSE
                )

        enc_dicts = list(candidate_SMILES.values())
        enc_dicts.sort()

        best_enc = enc_dicts[0]
        best_val = best_enc.latex_val(phantom_minus_alignment=True)
        best_improv_val = best_enc.latex_improve_val()

        minmax_candidate_table[(opt_quant_val_label, best_label)].append(best_val)
        minmax_candidate_table[(rel_improv_label, best_label)].append(best_improv_val)

        if len(enc_dicts) == 1:
            # The best candidate is the only one.
            for type_label in [opt_quant_val_label, rel_improv_label]:
                minmax_candidate_table[(type_label, worst_label)].append(
                    centered_cell(r"\_")
                )
            continue
        # Get the worst candidate.
        worst_enc = enc_dicts[-1]
        worst_val = worst_enc.latex_val(phantom_minus_alignment=True)
        worst_improv_val = worst_enc.latex_improve_val()

        minmax_candidate_table[(opt_quant_val_label, worst_label)].append(worst_val)
        minmax_candidate_table[(rel_improv_label, worst_label)].append(worst_improv_val)


brute_table_write(
    minmax_candidate_table,
    summary_dir + "/best_worst_candidates_" + dataset_name + ".tex",
    del_str=None,
)
