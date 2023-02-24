from print_dataset_size_table_1 import (
    reference_datasets,
    ref_dataset_egcs,
    compliant_egc_list,
)
from bmapqml.chemxpl.utils import chemgraph_to_canonical_rdkit
from bmapqml.utils import mkdir
import pandas as pd

phantom = "\phantom{\_}"


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

    def __call__(self, number_in):
        if isinstance(number_in, str):
            return number_in
        if isinstance(number_in, int):
            output = str(number_in)
        else:
            #            output="{:0.3e}".format(number_in)
            output = latex_scientific(number_in)
        if self.present_negative and number_in >= 0:
            output = "\phantom{-}" + output
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
    format_func = latex_table_format(present_negative=True)
    best_val = None
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
    return format_func(best_val), "$\phantom{-}$" + best_SMILES


hnum_mol = "num. mol."
hmin_dEsolv = "min. $\dEsolv^{\mathrm{conv.}}$, a.u."
hsmin_dEsolv = "min. $\dEsolv^{\mathrm{conv.}}$ SMILES"
hmax_RMSE_dEsolv = "max. RMSE($\dEsolv^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_dEsolv = "max. RMSE($\dEsolv^{\mathrm{conv.}}$) SMILES"
hmax_dipole = "max. $\dipole^{\mathrm{conv.}}$, a.u."
hsmax_dipole = "max. $\dipole^{\mathrm{conv.}}$ SMILES"
hmax_RMSE_dipole = "max. RMSE($\dipole^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_dipole = "max. RMSE($\dipole^{\mathrm{conv.}}$) SMILES"
hmax_RMSE_gap = "max. RMSE($\gap^{\mathrm{conv.}}$), a.u."
hsmax_RMSE_gap = "max. RMSE($\gap^{\mathrm{conv.}}$) SMILES"


def egc_extrema_dictionnary(relevant_egcs):
    lf = latex_table_format(present_negative=True)
    output = {hnum_mol: lf(len(relevant_egcs))}
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
        hmax_dipole,
        hsmax_dipole,
        hmax_RMSE_dipole,
        hsmax_RMSE_dipole,
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
            cur_h = "$\phantom{-}$" + gap_constr
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
