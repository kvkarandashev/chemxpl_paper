from print_dataset_size_table_1 import (
    reference_datasets,
    ref_dataset_egcs,
    compliant_egc_list,
)
from bmapqml.chemxpl.utils import chemgraph_to_canonical_rdkit
import pandas as pd
import numpy as np


quant_conversion_coefficients = {
    "solvation_energy": 2625.5002,
    "dipole": 2.541746,
    "HOMO_LUMO_gap": 27.211324570273,
}

phantom = "\phantom{\_}"

# For LaTeX conversion.
def preexp_power(n):
    s = "{:0.3e}".format(n)
    parts = s.split("e")
    return parts[0], int(parts[1])


def adjusted_error(n, power):
    format_str = "{:0.3f}"
    return format_str.format(n * 10 ** (-power))


def latex_scientific(number_in):
    def_sci = "{:0.3e}".format(number_in)
    parts = def_sci.split("e")
    output = parts[0]
    if len(parts) > 1:
        extra = parts[1]
        if extra != "+00":
            output += "\cdot 10^{" + str(int(extra)) + "}"
    return output


class latex_table_format:
    def __init__(self, present_negative=False):
        self.present_negative = present_negative

    def __call__(self, number_in, wmath=True):
        if isinstance(number_in, str):
            return number_in
        if isinstance(number_in, int):
            output = str(number_in)
        else:
            #            output="{:0.3e}".format(number_in)
            output = latex_scientific(number_in)
        if self.present_negative and number_in >= 0:
            output = "\phantom{-}" + output
        if wmath:
            output = "$" + output + "$"
        return output


def multrow(s, nrows):
    return ["\multirow{" + str(nrows) + "}{*}{" + str(s) + "}"] + [
        phantom for _ in range(nrows - 1)
    ]


def multcol(s, ncols):
    return "\multicolumn{" + str(ncols) + "}{c}{" + str(s) + "}"


std_RMSE_coeff = 0.25  # because averaging is done over 16 instances.


def better(val1, val2, find_max):
    if val2 is None:
        return True
    if find_max:
        return val1 > val2
    else:
        return val1 < val2


def egc2checked_SMILES(egc):
    SMILES = chemgraph_to_canonical_rdkit(egc.chemgraph, SMILES_only=True)
    res = ""
    for c in SMILES:
        if c == "#":
            res += "\\" + c
        else:
            res += c
    return res


def find_extrema(egcs, mean_or_std, quant_name, find_max):
    conversion_coefficient = quant_conversion_coefficients[quant_name]

    best_val = None
    best_RMSE = None
    best_SMILES = None
    for egc in egcs:
        ad = egc.additional_data
        val = ad[mean_or_std][quant_name]
        if val is None:
            continue
        if mean_or_std == "std":
            val *= std_RMSE_coeff
        if better(val, best_val, find_max):
            best_SMILES = egc2checked_SMILES(egc)
            best_val = val
            if mean_or_std == "mean":
                best_RMSE = ad["std"][quant_name] * std_RMSE_coeff
    # Now write all this in proper format.
    if best_val is not None:
        best_val *= conversion_coefficient
    if best_RMSE is not None:
        best_RMSE *= conversion_coefficient
    if mean_or_std == "std":
        val_str = "\phantom{-(}" + latex_scientific(best_val)
    else:
        mean_preexp, mean_power = preexp_power(best_val)
        err_str = adjusted_error(best_RMSE, mean_power)
        val_str = mean_preexp + " \pm " + err_str
        if mean_power == 0:
            val_str = "\phantom{(}" + val_str
        else:
            val_str = "(" + val_str + ")\cdot 10^{" + str(mean_power) + "}"
        if best_val > 0:
            val_str = "\phantom{-}" + val_str
    return "$" + val_str + "$", "$\phantom{-(}$" + best_SMILES


def find_global_std(egcs, quant_name):
    #    format_func = latex_table_format(present_negative=True)
    val_arr = []
    for egc in egcs:
        ad = egc.additional_data
        val = ad["mean"][quant_name]
        if val is None:
            continue
        val_arr.append(val)
    output_number = np.std(val_arr)
    output_number *= quant_conversion_coefficients[quant_name]
    output_str = "$\phantom{-(}" + latex_scientific(output_number) + "$"
    print("STD FOR", quant_name, output_number)
    return output_str


hnum_mol = "num. mol."
hmin_dEsolv = "min. $\dEsolv$, kJ/mol"
hsmin_dEsolv = "min. $\dEsolv$ SMILES"
hmax_RMSE_dEsolv = "max. RMSE($\dEsolv$), kJ/mol"
hsmax_RMSE_dEsolv = "max. RMSE($\dEsolv$) SMILES"
hSTD_dEsolv = "STD $\dEsolv$, kJ/mol"
hmax_dipole = "max. $\dipole$, debye"
hsmax_dipole = "max. $\dipole$ SMILES"
hmax_RMSE_dipole = "max. RMSE($\dipole$), debye"
hsmax_RMSE_dipole = "max. RMSE($\dipole$) SMILES"
hSTD_dipole = "STD $\dipole$, debye"
hmax_RMSE_gap = "max. RMSE($\gap$), eV"
hsmax_RMSE_gap = "max. RMSE($\gap$) SMILES"


def egc_extrema_dictionnary(relevant_egcs):
    #    lf = latex_table_format(present_negative=True)
    #    output = {hnum_mol: lf(len(relevant_egcs))}
    output = {hnum_mol: "$\phantom{-(}$" + str(len(relevant_egcs))}
    output[hmin_dEsolv], output[hsmin_dEsolv] = find_extrema(
        relevant_egcs, "mean", "solvation_energy", False
    )
    output[hmax_RMSE_dEsolv], output[hsmax_RMSE_dEsolv] = find_extrema(
        relevant_egcs, "std", "solvation_energy", True
    )
    output[hmax_dipole], output[hsmax_dipole] = find_extrema(
        relevant_egcs, "mean", "dipole", True
    )
    output[hmax_RMSE_dipole], output[hsmax_RMSE_dipole] = find_extrema(
        relevant_egcs, "std", "dipole", True
    )
    output[hmax_RMSE_gap], output[hsmax_RMSE_gap] = find_extrema(
        relevant_egcs, "std", "HOMO_LUMO_gap", True
    )

    output[hSTD_dEsolv] = find_global_std(relevant_egcs, "solvation_energy")
    output[hSTD_dipole] = find_global_std(relevant_egcs, "dipole")

    return output


# Evaluated in molopt/examples/chemxpl/MC_minimization/long_run_scripts/xTB_dipole_solvation/post_submission
num_constrained_egcs = {"QM9": 89015, "EGP": 9437}


def dataset_subheader(name):
    if name == "QM9":
        midrule_str = ""
    else:
        midrule_str = "\midrule"
    dataset_line = (
        midrule_str
        + "\multicolumn{3}{c}{"
        + name
        + " and "
        + name
        + "* intersection}\\\\\midrule"
    )

    #    tot_num_line=" tot. num. mol. & \multicolumn{2}{c}{"+str(num_constrained_egcs[name])+"}\\\\"
    tot_num_line = (
        " tot. num. mol. & $\phantom{-(}"
        + str(num_constrained_egcs[name])
        + "$ & $ \phantom{-(}\_$\\\\"
    )
    return [dataset_line, tot_num_line]


def main():
    gap_constraint_vals = ["weak", "strong"]
    #    latex_quantity = {
    #        "HOMO_LUMO_gap": "$\gap$",
    #        "solvation_energy": "$\dEsolv$",
    #        "dipole": "$\dipole$",
    #    }

    extrema_headers = [
        hnum_mol,
        hmin_dEsolv,
        hsmin_dEsolv,
        hmax_RMSE_dEsolv,
        hsmax_RMSE_dEsolv,
        hSTD_dEsolv,
        hmax_dipole,
        hsmax_dipole,
        hmax_RMSE_dipole,
        hsmax_RMSE_dipole,
        hSTD_dipole,
        hmax_RMSE_gap,
        hsmax_RMSE_gap,
    ]

    for ref_dataset in reference_datasets:
        table = {"$\gap$ constr.": extrema_headers}
        all_egcs = ref_dataset_egcs(ref_dataset)

        for gap_constr in gap_constraint_vals:
            relevant_egcs = compliant_egc_list(
                all_egcs, gap_constraint=gap_constr, chemspace_constraint=ref_dataset
            )
            #            cur_h = (multcol(ref_dataset, 2), gap_constr)
            cur_h = "$\phantom{-(}$" + gap_constr
            print("For", ref_dataset, gap_constr, ":")
            extrema_dict = egc_extrema_dictionnary(relevant_egcs)
            table[cur_h] = [extrema_dict[eh] for eh in extrema_headers]

        #            for writeout_h in [hSTD_dEsolv, hSTD_dipole]:
        #                print(writeout_h, extrema_dict[writeout_h])

        df = pd.DataFrame(table)
        df.to_latex(
            "dataset_reference_table_" + ref_dataset + ".tex",
            index=False,
            escape=False,
            multicolumn=True,
        )

    # Combining two tables into one.
    QM9_lines = open("dataset_reference_table_QM9.tex", "r").readlines()
    EGP_lines = open("dataset_reference_table_EGP.tex", "r").readlines()
    bottom_line_id = 4
    upper_line_id = len(extrema_headers) + bottom_line_id
    top = QM9_lines[:bottom_line_id]
    QM9_subheader = dataset_subheader("QM9")
    QM9_body = QM9_lines[bottom_line_id:upper_line_id]
    EGP_subheader = dataset_subheader("EGP")
    EGP_body = EGP_lines[bottom_line_id:upper_line_id]
    bottom = EGP_lines[upper_line_id:]
    big_table = open("dataset_reference_table.tex", "w")
    for l in top + QM9_subheader + QM9_body + EGP_subheader + EGP_body + bottom:
        big_table.write(l)
    big_table.close()


if __name__ == "__main__":
    main()
