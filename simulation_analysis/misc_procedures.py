import glob, os
import numpy as np
import subprocess


def all_xyz_vals(xyz_file):
    i = open(xyz_file, "r")
    lines = i.readlines()
    if len(lines) == 1:
        return None
    l = lines[1].split()
    output = {}
    for q in l:
        q_spl = q.split("=")
        quant_name = q_spl[0]
        output[quant_name] = q.split(quant_name + "=")[1]
    return output


def val_in_xyz(xyz_file, quant_name, quant_type=float):
    return quant_type(all_xyz_vals(xyz_file)[quant_name])


def extract_hist_size(data_dir):
    stdouts = glob.glob(data_dir + "/*.stdout_*")
    if len(stdouts) == 0:
        out = data_dir + "/joblog.txt"
    else:
        out = stdouts[0]
    f = open(out, "r")
    l = f.readlines()
    f.close()
    for s in l[::-1]:
        spl = s.split()
        if spl[0] == "HIST":
            return int(spl[4])


gap_constraint_values = {"weak": 0.08895587351640835, "strong": 0.17734766509696615}

# Reference data is presented in tuples as value : std (divide by 4 to get RMSE) : ChemGraph string representation.

best_ref_vals = {
    "QM9": {
        "weak": {
            "dipole": (
                5.2624256900482225,
                8.190772053238254e-05,
                "8@7:7#2@8:7@5@6@8:7@6@7:7@7@8:6#3:6#1:6:6",
            ),
            "solvation_energy": (
                -0.036118152765830525,
                9.940944965402545e-05,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "atomization_energy": (
                -5.489771247476437,
                4.704659650874689e-05,
                "6#3@7:6#3@7:6#3@8:6#3@8:6#3@8:6#2@6@7:6#2@8:6#1:6",
            ),
        },
        "strong": {
            "dipole": (
                3.3568648712900764,
                0.06473921973977188,
                "8@7:8@8:7#1@7@8:7#1@5@7:7#1@6@8:6#2@6:6#2:6:6",
            ),
            "solvation_energy": (
                -0.021467490751841414,
                0.0005614388897913243,
                "8@7:8@8:7#2@7:7#2@8:7#1@5@7:6#2@6:6#2@8:6:6",
            ),
            "atomization_energy": (
                -5.489771247476437,
                4.704659650874689e-05,
                "6#3@7:6#3@7:6#3@8:6#3@8:6#3@8:6#2@6@7:6#2@8:6#1:6",
            ),
        },
    },
    "EGP": {
        "weak": {
            "dipole": (
                5.217191799780868,
                0.0022688571698954913,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "solvation_energy": (
                -0.03608624382949954,
                8.482093682841212e-05,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "atomization_energy": (
                -8.130926853176323,
                0.0015182550242717613,
                "6#3@8:6#3@8:6#3@9:6#3@9:6#3@10:6#3@10:6#3@11:6#3@11:6#1@12:6#1@12:6#1@13:6#1@13:6@13:6",
            ),
        },
        "strong": {
            "dipole": (
                3.1670498606166886,
                0.07476372139614328,
                "16@2@4@5@8:16@3@6@7@9:9:9:8:8:8:8:6#2@9:6#2",
            ),
            "solvation_energy": (
                -0.02976335641089712,
                0.0002888209390425601,
                "16@2@3@6@11:16@4@5@7@12:8:8:8:8:7#2:7#2:6#2@10@11:6#2@10@12:6#2:6#2:6#2",
            ),
            "atomization_energy": (
                -8.065834964640352,
                0.0007362348666310831,
                "7#2@12:7#2@13:6#2@4@12:6#2@5@13:6#2@6:6#2@7:6#2@8:6#2@9:6#2@10:6#2@11:6#2@11:6#2:6#2:6#2",
            ),
        },
    },
}

max_ref_std = {
    "QM9": {
        "weak": {
            "dipole": (
                1.6617385375337195,
                1.4871285627241906,
                "8@4@7:8@5:7#1@5@8:7#1@7:6#2@6:6#1:6#1@8:6@8:6",
            ),
            "solvation_energy": (
                -0.02114448445313673,
                0.0027531017138107917,
                "8@6@7:8@5:7#2@8:7#1@6@8:7#1@6:6#1@7:6:6@8:6",
            ),
            "atomization_energy": (
                -4.585098551070718,
                0.08087940612902654,
                "6#3@1:6#2@5:6#2@5@8:6#2@6@7:6#2@6@7:6#1@8:6#1@8:6#1@8:6",
            ),
        },
        "strong": {
            "dipole": (
                1.4202701525407528,
                0.5232115776074905,
                "8@4@5:7@2@6@7:6#2@8:6#2@7@8:6#1@5@6:6#1@7:6#1@8:6#1:6#1",
            ),
            "solvation_energy": (
                -0.00445322014735519,
                0.0010038515425404924,
                "8@4@5:6#2@2@6:6#2@6:6#2@4@5:6#1@6:6#1@6:6",
            ),
            "atomization_energy": (
                -4.585098551070718,
                0.08087940612902654,
                "6#3@1:6#2@5:6#2@5@8:6#2@6@7:6#2@6@7:6#1@8:6#1@8:6#1@8:6",
            ),
        },
    },
    "EGP": {
        "weak": {
            "dipole": (
                2.0119263516172827,
                1.0926926114878324,
                "8@7@8:7#2@7:7#1@6@8:7#1@6:7@5@8:6#1@7:6#1:6:6",
            ),
            "solvation_energy": (
                -0.017652397664566492,
                0.0010603139018516197,
                "8@7:7#2@7:7#1@5@8:7#1@8:7@6@7@8:6#1@6:6#1:6:6",
            ),
            "atomization_energy": (
                -1.9775680617527625,
                0.05185659142279606,
                "6#2@1@3:6#2@3:6#2@3:6",
            ),
        },
        "strong": {
            "dipole": (
                0.44818039816643773,
                0.5003299395549212,
                "7@8:7@9:6#2@3@6:6#2@7:6#2@5@6:6#2@7:6#1@8:6#1@9:6:6",
            ),
            "solvation_energy": (
                -0.0038722831843677006,
                0.0005331684106806952,
                "17@2:8@3@7:6#2@7:6#2@4:6#2@5:6#2@6:6#2@7:6#1",
            ),
            "atomization_energy": (
                -1.9775680617527625,
                0.05185659142279606,
                "6#2@1@3:6#2@3:6#2@3:6",
            ),
        },
    },
}

quant_STD = {
    "QM9": {
        "weak": {
            "dipole": 0.7202490150375122,
            "solvation_energy": 0.003800219983701787,
        },
        "strong": {
            "dipole": 0.5330517652144172,
            "solvation_energy": 0.0027043460958190825,
        },
    },
    "EGP": {
        "weak": {
            "dipole": 0.740092558470457,
            "solvation_energy": 0.0035497259204171077,
        },
        "strong": {
            "dipole": 0.5868571986460336,
            "solvation_energy": 0.003148923583866458,
        },
    },
}

# Both QM9 and EGP have methane as the compound with the largest gap.
best_gap_val = 0.6193166670127856


def folder_name_param_dict(folder_name):
    bname = os.path.basename(folder_name)
    bname_split = bname.split("_")
    last_id = -1
    seed = int(bname_split[last_id])
    last_id -= 1
    if bname_split[last_id] == "energy":
        if bname_split[last_id - 1] == "solvation":
            quant = "solvation_energy"
        else:
            quant = "atomization_energy"
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
        "best_ref_val": best_ref_vals[dataset][gap_constraint][quant][0],
        "gap_constr_val": gap_constraint_values[gap_constraint],
    }


def data_dir_best_candidate_val(dirname, quant):
    return val_in_xyz(dirname + "/best_candidates/best_candidate_0.xyz", quant)


def find_seed_with_best_candidate(dirname):
    seed_dirs = glob.glob(dirname + "/data_*")
    min_val = None
    min_seed_dir = None
    for seed_dir in seed_dirs:
        cur_val = data_dir_best_candidate_val(seed_dir, "final_minfunc_val")
        if (min_val is None) or (min_val > cur_val):
            min_val = cur_val
            min_seed_dir = seed_dir
    return min_seed_dir


# Procedures for analyzing the graph enumeration results.


def MC_enumeration_folder_name_param_dict(folder_name):
    true_folder_name = os.path.basename(folder_name)
    if len(true_folder_name) == 0:
        true_folder_name = os.path.basename(os.path.dirname(folder_name))
    naming_fields = true_folder_name.split("_")
    # The ordering can be looked up in molopt/examples/chemxpl/bias_potential_checks/MC_graph_enumeration_submit.sh
    seed = int(naming_fields[-1])
    implicit_constraint = naming_fields[-2] == "TRUE"
    genetic_used = naming_fields[-3] == "TRUE"
    bias_type = naming_fields[-4]
    nhatoms = int(naming_fields[-5])
    return {
        "seed": seed,
        "implicit_constraint": implicit_constraint,
        "genetic_used": genetic_used,
        "bias_type": bias_type,
        "nhatoms": nhatoms,
    }


import pandas as pd


def table_to_data(table, headers):
    data = {}
    for h_id, h in enumerate(headers):
        if h is None:
            continue
        data[h] = []
        for row in table:
            data[h].append(row[h_id])
    return data


def table_to_data_transpose(table, headers):
    all_rows = [headers] + table
    data = {}
    for row in all_rows:
        for h_id, (h, el) in enumerate(zip(headers, row)):
            if h is None:
                continue
            if h_id == 0:
                used_header = el
                data[used_header] = []
            else:
                data[used_header].append(el)
    return data


def table_to_csv(table, headers, filename):
    """
    Prints tables I normally use to a proper cvs.
    """
    data = table_to_data(table, headers)
    df = pd.DataFrame(data, columns=data.keys(), index=data[headers[0]])
    df.to_csv(filename, index=False)


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


def table_to_latex(
    table, latex_headers, filename, transpose=False, data=None, **to_latex_kwargs
):
    """
    Print tables I normally use into LaTeX.
    """
    if data is None:
        if transpose:
            data = table_to_data_transpose(table, latex_headers)
        else:
            data = table_to_data(table, latex_headers)

    present_negative = False
    for k in data.keys():
        for v in data[k]:
            if isinstance(v, np.float):
                if v < 0:
                    present_negative = True
    lf = latex_table_format(present_negative=present_negative)

    converted_data = {}
    for k, vals in data.items():
        arr = []
        true_k = None
        for v in vals:
            arr.append(lf(v))
            if (true_k is None) and present_negative:
                if isinstance(v, np.float):
                    true_k = "$\phantom{-}$" + k
        if true_k is None:
            true_k = k
        converted_data[true_k] = arr

    df = pd.DataFrame(
        converted_data,
        columns=converted_data.keys(),
        index=converted_data[latex_headers[0]],
    )

    df.to_latex(filename, index=False, escape=False, **to_latex_kwargs)


def brute_table_write(table_in, filename, del_str=None):
    max_colwidth = 0
    for k, l in table_in.items():
        max_colwidth = max(max_colwidth, len(k))
        for s in l:
            max_colwidth = max(max_colwidth, len(s))
    df = pd.DataFrame(table_in, columns=table_in.keys())
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            filename,
            index=False,
            escape=False,
            multicolumn=True,
            multicolumn_format="c",
            multirow=True,
        )
    if del_str is not None:
        subprocess.run(["sed", "-i", "s/&[ ]*" + del_str + "/  /g", filename])
    subprocess.run(["sed", "-i", "s/\$\\\\phantom{(}/\\\\phantom{(}\$/g", filename])
