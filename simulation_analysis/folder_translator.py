import sys, os

gap_constraint_values = {"weak": 0.09126989358754387, "strong": 0.17735152497325582}

best_ref_vals = {
    "QM9": {
        "weak": {
            "dipole": 5.2624387362231335,
            "solvation_energy": -0.03604551425739544,
        },
        "strong": {
            "dipole": 3.399629855492714,
            "solvation_energy": -0.021567673361052482,
        },
    },
    "EGP": {
        "weak": {"dipole": 5.216205876716626, "solvation_energy": -0.03604697067138929},
        "strong": {
            "dipole": 3.142098394030278,
            "solvation_energy": -0.029593595306102352,
        },
    },
}


def folder_name_param_dict(folder_name):
    bname = os.path.basename(folder_name)
    bname_split = bname.split("_")
    last_id = -1
    seed = int(bname_split[last_id])
    last_id -= 1
    if bname_split[last_id] == "energy":
        quant = "solvation_energy"
        last_id -= 2
    else:
        quant = "dipole"
        last_id -= 1
    gap_constraint = bname_split[last_id]
    last_id -= 1
    bias_strength = bname_split[last_id]
    last_id -= 1
    if "egp" in bname_split[:last_id]:
        dataset = "EGP"
    else:
        dataset = "QM9"
    return {
        "dataset": dataset,
        "seed": seed,
        "gap_constraint": gap_constraint,
        "quant": quant,
        "bias_strength": bias_strength,
        "best_ref_val": best_ref_vals[dataset][gap_constraint][quant],
        "gap_constr_val": gap_constraint_values[gap_constraint],
    }


if __name__ == "__main__":
    folder_name = sys.argv[1]
    print(folder_name_param_dict(folder_name))
