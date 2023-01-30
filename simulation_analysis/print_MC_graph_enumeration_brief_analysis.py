import os, sys, glob
from bmapqml.utils import mkdir
import numpy as np
from misc_procedures import MC_enumeration_folder_name_param_dict, table_to_csv
from tabulate import tabulate

bias_types = ["none", "weak", "stronger"]


def use_genetic_available_values(nhatoms):
    if nhatoms == 3:
        return [False]
    else:
        return [False, True]


def needed_subfolders(all_folder_dir, implicit_constraint, nhatoms):
    """
    Returns two-level dictionnary: outer level is bias strength (different rows), inner - whether genetic was used (different columns).
    """
    output = {}
    for subfolder in glob.glob(all_folder_dir + "/data_*"):
        param_dict = MC_enumeration_folder_name_param_dict(subfolder)
        if (nhatoms != param_dict["nhatoms"]) or (
            implicit_constraint != param_dict["implicit_constraint"]
        ):
            continue
        bias_type = param_dict["bias_type"]
        if bias_type not in output:
            output[bias_type] = {}
        subdict = output[bias_type]
        genetic_used = param_dict["genetic_used"]
        if genetic_used not in subdict:
            subdict[genetic_used] = []
        subdict[genetic_used].append(subfolder)
    return output


def hist_size(subfolder):
    joblog = subfolder + "/joblog.txt"
    if not os.path.isfile(joblog):
        raise Exception("joblog.txt not in subfolder")
    inp = open(joblog, "r")
    max_step = None
    max_size = None
    for l in inp.readlines():
        lspl = l.split()
        if lspl[:2] == ["HIST", "SIZE:"]:
            cur_size = int(lspl[4])
            if (max_size is None) or (cur_size > max_size):
                max_size = cur_size
                max_step = int(lspl[2])
    inp.close()
    return max_step, max_size


def print_param_influence_summary(all_folder_dir, implicit_constraint, nhatoms):
    headers = {True: "w. genetic", False: "no genetic"}
    funcs = {"mean": np.mean, "std": np.std}

    subfolders = needed_subfolders(all_folder_dir, implicit_constraint, nhatoms)
    use_genetic_values = use_genetic_available_values(nhatoms)
    headers = ["bias type"]
    for quant in ["max step", "max size"]:
        for genetic_used in use_genetic_values:
            for t in ["mean", "std"]:
                headers.append(t + " " + quant + " " + headers[genetic_used] + " ")
    table = []
    for bias_type in bias_types:
        row_data = {}
        for genetic_used in use_genetic_values:
            row_data[genetic_used] = {}
            diff_seed_subfolders = subfolders[bias_type][genetic_used]
            row_data[genetic_used]["max step"] = []
            row_data[genetic_used]["max size"] = []
            for (max_step, max_size) in [
                hist_size(subfolder) for subfolder in diff_seed_subfolders
            ]:
                row_data[genetic_used]["max step"].append(max_step)
                row_data[genetic_used]["max size"].append(max_size)
        row = [bias_type]
        for quant in ["max step", "max size"]:
            for genetic_used in use_genetic_values:
                for t in ["mean", "std"]:
                    row.append(funcs[t](row_data[genetic_used][quant]))
        table.append(row)

    true_all_folder_dir = os.path.basename(all_folder_dir)
    if len(true_all_folder_dir) == 0:
        true_all_folder_dir = os.path.basename(os.path.dirname(all_folder_dir))

    summary_dir = "summary_" + true_all_folder_dir

    mkdir(summary_dir)
    summary_filename = (
        summary_dir
        + "/job_summary_nhatoms_"
        + str(nhatoms)
        + "_implicit_constraint_"
        + str(implicit_constraint)
    )
    txt_output = open(summary_filename + ".txt", "w")
    print(tabulate(table, headers=["# " + s for s in headers]), file=txt_output)
    txt_output.close()
    table_to_csv(table, headers, summary_filename + ".csv")


def main():
    all_folder_dir = sys.argv[1]

    nhatoms_values = [3, 4, 5, 6]

    implicit_constraint_values = [False, True]

    for implicit_constraint in implicit_constraint_values:
        for nhatoms in nhatoms_values:
            print_param_influence_summary(all_folder_dir, implicit_constraint, nhatoms)


if __name__ == "__main__":
    main()
