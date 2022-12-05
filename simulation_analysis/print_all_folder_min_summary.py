import os, sys, glob, subprocess
from bmapqml.utils import mkdir, loadpkl
import numpy as np


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


parent_folder = sys.argv[1]

summary_dir = os.getcwd() + "/summary_" + os.path.basename(parent_folder)

mkdir(summary_dir)

os.chdir(summary_dir)

biases = ["none", "weak", "stronger"]

gap_constraints = ["weak", "strong"]

quantities = ["dipole", "solvation_energy"]

output = open(summary_dir + "/all_output.txt", "w")

for quantity in quantities:
    for gap_constraint in gap_constraints:
        quant_coeff = None
        for bias in biases:
            weird_counter = 0
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            folder = parent_folder + bulk_name
            data_dirs = glob.glob(folder + "/data_*")
            min_data_dir = None
            min_val = None
            vals = []
            eval_nums = []
            for data_dir in data_dirs:
                if quant_coeff is None:
                    quant_pkl = glob.glob(data_dir + "/*_water_*.pkl")
                    if len(quant_pkl) != 1:
                        raise Exception
                    quant_func = loadpkl(quant_pkl[0])
                    quant_coeff = 1.0 / quant_func.coefficients[0]
                    print(
                        "STD for ",
                        quantity,
                        "constraint ",
                        gap_constraint,
                        ":",
                        abs(quant_coeff),
                        file=output,
                    )
                cur_val, weird = val_in_xyz(
                    data_dir + "/best_candidate_0.xyz", quantity
                )
                if weird:
                    weird_counter += 1
                    subprocess.run(
                        [
                            "cp",
                            data_dir + "/best_candidate_0.png",
                            summary_dir
                            + "/weird_"
                            + str(weird_counter)
                            + bulk_name
                            + ".png",
                        ]
                    )
                vals.append(cur_val * quant_coeff)
                eval_nums.append(extract_hist_size(data_dir))
                if (min_val is None) or (min_val > cur_val):
                    min_val = cur_val
                    min_data_dir = data_dir
            for ending in ["xyz", "png"]:
                subprocess.run(
                    [
                        "cp",
                        min_data_dir + "/best_candidate_0." + ending,
                        summary_dir + "/best" + bulk_name + "." + ending,
                    ]
                )
            print(
                "MIN FOUND FOR",
                quantity,
                "constraint ",
                gap_constraint,
                "bias ",
                bias,
                ":",
                min_val * quant_coeff,
                file=output,
            )
            print("MEAN AND STD:", np.mean(vals), np.std(vals), file=output)
            print(
                "Number of evaluations:",
                np.mean(eval_nums),
                np.std(eval_nums),
                file=output,
            )

output.close()
