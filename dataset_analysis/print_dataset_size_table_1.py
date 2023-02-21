from bmapqml.utils import loadpkl
from bmapqml.chemxpl.ext_graph_compound import ExtGraphCompound
from bmapqml.chemxpl.valence_treatment import str2ChemGraph
from bmapqml.chemxpl.random_walk import egc_valid_wrt_change_params
import pandas as pd
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit
from sortedcontainers import SortedList
from bmapqml.data import NUCLEAR_CHARGE
from bmapqml.dataset_processing.qm9_format_specs import read_SMILES

"""
The script produces condensed version of the output of minfunc_creation/morfeus_xTB_solvation/*_morfeus_xTB_analyse_wconstr_*.py scripts
in form of a LaTeX table.
Also produces lists of SMILES strings of dataset molecules compiant with QM9* and EGP*, as well as strong or weak gap constraints.
"""

# Datasets considered.

reference_datasets = ["QM9", "EGP"]

# Specifications of QM9* and EGP*

forbidden_bonds = {
    "QM9": [
        (7, 7),
        (8, 8),
        (9, 9),
        (7, 8),
        (7, 9),
        (8, 9),
    ],
    "EGP": [
        (7, 7),
        (8, 8),
        (9, 9),
        (7, 8),
        (7, 9),
        (8, 9),
        (7, 17),
        (8, 17),
        (9, 17),
        (17, 17),
        (7, 35),
        (8, 35),
        (9, 35),
        (17, 35),
        (35, 35),
    ],
}

not_protonated = {"QM9": [8, 9], "EGP": [5, 8, 9, 14, 15, 16, 17, 35]}

nhatoms_range = {"QM9": [1, 9], "EGP": [1, 15]}

possible_elements = {
    "QM9": ["C", "N", "O", "F"],
    "EGP": ["B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br"],
}

# Gap constraints considered.

gap_constraints = [None, "weak", "strong"]

gap_constraint_vals = {"weak": 0.08895587351640835, "strong": 0.17734766509696615}

gap_name = "HOMO_LUMO_gap"


def compliant_egc_list(egcs, gap_constraint=None, chemspace_constraint=None):
    if gap_constraint is None:
        unconv_counter = 0
        egcs_converged = []
    if chemspace_constraint is None:
        chemspace_constraint_dict = {}
    else:
        chemspace_constraint_dict = {
            "forbidden_bonds": forbidden_bonds[chemspace_constraint],
            "not_protonated": not_protonated[chemspace_constraint],
            "nhatoms_range": nhatoms_range[chemspace_constraint],
            "possible_elements": possible_elements[chemspace_constraint],
        }

    output = []
    for egc in egcs:
        if egc is None:
            continue
        if egc_valid_wrt_change_params(egc, **chemspace_constraint_dict):
            gap_val = egc.additional_data["mean"][gap_name]
            if gap_val is None:
                if gap_constraint is None:
                    unconv_counter += 1
                else:
                    continue
            if gap_constraint is None:
                egcs_converged.append(egc)
            else:
                if gap_val < gap_constraint_vals[gap_constraint]:
                    continue
            output.append(egc)
    if gap_constraint is None:
        print("Unconverged calculation counter:", unconv_counter)
        print("Elements for converged:")
        print_present_heavy_elements(egcs_converged)
        print("Number of converged calculations:", len(egcs_converged))
    return output


def print_present_heavy_elements(egc_list):
    possible_ncharges = SortedList()
    for egc in egc_list:
        if egc is None:
            continue
        for ha in egc.chemgraph.hatoms:
            ncharge = ha.ncharge
            if ncharge not in possible_ncharges:
                possible_ncharges.add(ncharge)
    possible_elements = []
    for el, nc in NUCLEAR_CHARGE.items():
        if nc in possible_ncharges:
            possible_elements.append(el)
    print(*possible_elements)


def print_canon_SMILES(egc_list, output_file_name):
    output = open(output_file_name, "w")
    for egc in egc_list:
        if egc is None:
            SMILES = str(egc)
        else:
            SMILES = chemgraph_to_canonical_rdkit(egc.chemgraph, SMILES_only=True)
        print(SMILES, file=output)
    output.close()


def cg_str(entry):
    for i in ["chemgraph", "chemgraph_str"]:
        if i in entry:
            return entry[i]
    return None


def ref_dataset_egcs(reference_dataset):
    # The data pkl files are generated by scripts in minfunc_creation/morfeus_xTB_solvation/*_morfeus_xTB_data_gen.py
    data_pkl_name = reference_dataset + "_morfeus_xTB_data_water.pkl"
    entry_list = loadpkl(data_pkl_name)
    egcs = []
    for entry in entry_list:
        s = cg_str(entry)
        if s is not None:
            cg = str2ChemGraph(s)
            egc = ExtGraphCompound(chemgraph=cg, additional_data=entry)
        else:
            egc = None
        egcs.append(egc)
    return egcs


def dataset_table_row(reference_dataset):
    tot_egc_list = ref_dataset_egcs(reference_dataset)

    table_row = [len(tot_egc_list)]
    SMILES_file_prefix = "SMILES_" + reference_dataset + "_"

    print("Reference dataset:", reference_dataset)
    print_present_heavy_elements(tot_egc_list)

    print_canon_SMILES(tot_egc_list, SMILES_file_prefix + "full.txt")
    for gap_constraint in gap_constraints:
        cur_egc_list = compliant_egc_list(
            tot_egc_list,
            gap_constraint=gap_constraint,
            chemspace_constraint=reference_dataset,
        )
        print_canon_SMILES(
            cur_egc_list,
            SMILES_file_prefix + "gap_constr_" + str(gap_constraint) + ".txt",
        )
        table_row.append(len(cur_egc_list))
        print("Elements for gap constraint", gap_constraint)
        print_present_heavy_elements(cur_egc_list)
    return table_row


def main():
    headers = ["dataset", "tot. size", "comp. size", "weak constr.", "strong constr."]

    table = {}
    for h in headers:
        table[h] = []

    for reference_dataset in reference_datasets:
        table_row = [reference_dataset] + dataset_table_row(reference_dataset)
        for h, el in zip(headers, table_row):
            table[h].append(el)

    df = pd.DataFrame(table)
    df.to_latex("dataset_sizes.tex", index=False)


if __name__ == "__main__":
    main()
