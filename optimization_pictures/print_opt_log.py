from bmapqml.utils import loadpkl
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit
import sys  # , os, copy

pkl_file = sys.argv[1]

output_file = sys.argv[2]

nprocs = 8

print("Considered restart file:", pkl_file)

try:
    cur_data = loadpkl(pkl_file)
except:
    print("INVALID RESTART FILE")

num_global_MC_steps = cur_data["global_MC_step_counter"]

if num_global_MC_steps != 50000:
    print("SIMULATION INCOMPLETE")
    quit()

min_func_name = cur_data["min_function_name"]

minimized_function = cur_data["min_function"]


xTB_res_label = "xTB_res"

histogram = cur_data["histogram"]

# WARNING: I realized too late that first_global_MC_step_encounter and first_MC_step_encounter the values were switched in the code used to generate some restart files.
def first_step_tuple(tp):
    return (tp.first_MC_step_encounter, tp.first_global_MC_step_encounter)


def first_global_MC_step_encounter(tp):
    return min(first_step_tuple(tp))


# Compile list of best candidates found at each step.
best_candidate_val_at_step = [None for _ in range(num_global_MC_steps + 1)]
best_candidate_at_step = [None for _ in range(num_global_MC_steps + 1)]

for entry in histogram:
    cur_val = entry.calculated_data[min_func_name]
    if cur_val is None:
        continue
    discover_step = first_global_MC_step_encounter(entry)
    replace = best_candidate_val_at_step[discover_step] is None
    if not replace:
        replace = best_candidate_val_at_step[discover_step] > cur_val
    if replace:
        best_candidate_val_at_step[discover_step] = cur_val
        best_candidate_at_step[discover_step] = chemgraph_to_canonical_rdkit(
            entry.chemgraph(), SMILES_only=True
        )

# Print the current best candidate at each step.
output = open(output_file, "w")

running_val = None
running_SMILES = None

coeff = minimized_function.coefficients[0] ** (-1)

for global_MC_step, (cur_val, cur_SMILES) in enumerate(
    zip(best_candidate_val_at_step, best_candidate_at_step)
):
    if (cur_val is not None) and ((running_val is None) or (cur_val < running_val)):
        running_val = cur_val
        running_SMILES = cur_SMILES
    print(global_MC_step, running_val * coeff, running_SMILES, file=output)

output.close()
