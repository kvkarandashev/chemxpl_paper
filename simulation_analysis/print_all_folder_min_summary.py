import os, sys, glob, subprocess
from bmapqml.utils import mkdir, loadpkl
import numpy as np


def val_in_xyz(xyz_file, quant_name):
    i = open(xyz_file, "r")
    l = i.readlines()[1].split()
    for q in l:
        q_spl = q.split("=")
        if q_spl[0] == quant_name:
            return float(q_spl[1])
    raise Exception


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
            bulk_name = "_" + bias + "_" + gap_constraint + "_" + quantity
            folder = parent_folder + bulk_name
            data_dirs = glob.glob(folder + "/data_*")
            min_data_dir = None
            min_val = None
            vals = []
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
                cur_val = val_in_xyz(data_dir + "/best_candidate_0.xyz", quantity)
                vals.append(cur_val * quant_coeff)
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

output.close()
