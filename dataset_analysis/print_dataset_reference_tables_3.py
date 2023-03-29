from print_dataset_size_table_1 import (
    reference_datasets,
    ref_dataset_egcs,
    compliant_egc_list,
)
from bmapqml.chemxpl.utils import chemgraph_to_canonical_rdkit
import pandas as pd
import numpy as np

phantom = "\phantom{\_}"

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
    output_str = "$\phantom{-(}" + latex_scientific(output_number) + "$"
    return output_str


hnum_mol = "num. mol."
hmin_dEsolv = "min. $\dEsolv^{\mathrm{conv.}}$, a.u."
hsmin_dEsolv = "min. $\dEsolv^{\mathrm{conv.}}$ SMILES"
hmax_RMSE_dEsolv = "max. RMSE($\dEsolv^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_dEsolv = "max. RMSE($\dEsolv^{\mathrm{conv.}}$) SMILES"
hSTD_dEsolv = "STD $\dEsolv^{\mathrm{conv.}}$, a.u."
hmax_dipole = "max. $\dipole^{\mathrm{conv.}}$, a.u."
hsmax_dipole = "max. $\dipole^{\mathrm{conv.}}$ SMILES"
hmax_RMSE_dipole = "max. RMSE($\dipole^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_dipole = "max. RMSE($\dipole^{\mathrm{conv.}}$) SMILES"
hSTD_dipole = "STD $\dipole^{\mathrm{conv.}}$, a.u."
hmax_RMSE_gap = "max. RMSE($\gap^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_gap = "max. RMSE($\gap^{\mathrm{conv.}}$) SMILES"


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
    table = {"$\gap$ constr.": extrema_headers}

    for ref_dataset in reference_datasets:
        table = {"$\gap$ constr.": extrema_headers}
        all_egcs = ref_dataset_egcs(ref_dataset)

        for gap_constr in gap_constraint_vals:
            relevant_egcs = compliant_egc_list(
                all_egcs, gap_constraint=gap_constr, chemspace_constraint=ref_dataset
            )
            #            cur_h = (multcol(ref_dataset, 2), gap_constr)
            cur_h = "$\phantom{-(}$" + gap_constr
            extrema_dict = egc_extrema_dictionnary(relevant_egcs)
            table[cur_h] = [extrema_dict[eh] for eh in extrema_headers]

        df = pd.DataFrame(table)
        df.to_latex(
            "dataset_reference_table_" + ref_dataset + ".tex",
            index=False,
            escape=False,
            multicolumn=True,
        )


if __name__ == "__main__":
    main()
