import glob, os


def val_in_xyz(xyz_file, quant_name):
    i = open(xyz_file, "r")
    lines = i.readlines()
    if len(lines) == 1:
        s = lines[0].split()[2]
        weird = True
    else:
        weird = False
        l = lines[1].split()
        for q in l:
            q_spl = q.split("=")
            if q_spl[0] == quant_name:
                s = q_spl[1]
                break
    return float(s), weird


def extract_hist_size(data_dir):
    out = glob.glob(data_dir + "/*.stdout_*")[0]
    f = open(out, "r")
    l = f.readlines()
    for s in l[::-1]:
        spl = s.split()
        if spl[0] == "HIST":
            return int(spl[4])


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


def data_dir_best_candidate_val(dirname, quant):
    return val_in_xyz(dirname + "/best_candidates/best_candidate_0.xyz", quant)


def find_seed_with_best_candidate(dirname):
    seed_dirs = glob.glob(dirname + "/data_*")
    quant = folder_name_param_dict(seed_dirs[0])["final_minfunc_val"]
    min_val = None
    for seed_dir in seed_dirs:
        cur_val = data_dir_best_candidate_val(seed_dir, quant)
        if (min_val is None) or (min_val < cur_val):
            min_val = cur_val
            min_seed_dir = seed_dir
    return min_seed_dir
