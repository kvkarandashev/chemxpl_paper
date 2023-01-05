from bmapqml.utils import write_xyz_file
from bmapqml.chemxpl.minimized_functions.morfeus_quantity_estimates import (
    morfeus_coord_info_from_tp,
)
from bmapqml.chemxpl.random_walk import CandidateCompound
from bmapqml.utils import loadpkl, mkdir
from bmapqml.chemxpl.rdkit_draw_utils import draw_chemgraph_to_file
from sortedcontainers import SortedList
from joblib import Parallel, delayed
import sys, os

pkl_file = sys.argv[1]

num_init_best_candidates = 1000

print("Considered restart file:", pkl_file)

try:
    cur_data = loadpkl(pkl_file)
except:
    print("INVALID RESTART FILE")

if cur_data["global_MC_step_counter"] != 50000:
    print("SIMULATION INCOMPLETE")
    quit()

min_func_name = cur_data["min_function_name"]
min_func_modified = cur_data["min_function"]

print("Original minimized function:", min_func_modified)

min_func_modified.xTB_related_kwargs["num_attempts"] = 16

min_func_modified.xTB_related_kwargs["num_conformers"] = 32

min_func_modified_label = "final_minfunc_val"

opt_quant = min_func_modified.quantities[0]

xTB_res_label = "xTB_res"

print("Imporved minimized function:", min_func_modified)

histogram = cur_data["histogram"]
# Get the candidates according to the underconverged minimized function.
best_candidates = SortedList()

for entry in histogram:
    cur_val = entry.calculated_data[min_func_name]
    if cur_val is not None:
        cur_cc = CandidateCompound(entry, entry.calculated_data[min_func_name])
        if len(best_candidates) == num_init_best_candidates:
            if cur_cc > best_candidates[-1]:
                continue
        best_candidates.add(cur_cc)
        if len(best_candidates) == num_init_best_candidates + 1:
            del best_candidates[-1]

# For each of those best candidates we also check how many trajectory points had been considered up to the point that
# a best handidate had been discovered.

for i, bc in enumerate(best_candidates):
    tp = bc.tp
    num_MC_moves = tp.first_global_MC_step_encounter
    # WARNING: I realized too late that first_global_MC_step_encounter and first_MC_step_encounter the values were switched in the code used to generate the restart files.
    tp.calculated_data["req_num_tps"] = sum(
        (tp.first_global_MC_step_encounter <= num_MC_moves) for tp in histogram
    )
    tp.calculated_data["minfunc_sorted_id"] = i
    tp.calculated_data["first_MC_step_encounter"] = num_MC_moves
    tp.calculated_data["first_global_MC_step_encounter"] = tp.first_MC_step_encounter
    del tp.calculated_data[xTB_res_label]

# Calculate the overconverged minimized function values.
def tp_wconv_vals(cand):
    tp = cand.tp
    tp.calculated_data[min_func_modified_label] = min_func_modified(tp)
    for quant_type in ["mean", "std"]:
        tp.calculated_data["overconv_quant_" + quant_type] = tp.calculated_data[
            xTB_res_label
        ][quant_type][opt_quant]
    return cand


resorted_best_candidates = Parallel(n_jobs=40, backend="multiprocessing")(
    delayed(tp_wconv_vals)(bc) for bc in best_candidates
)

i = 0
while i != len(resorted_best_candidates):
    if resorted_best_candidates[i].tp.calculated_data[min_func_modified_label] is None:
        del resorted_best_candidates[i]
    else:
        i += 1

resorted_best_candidates.sort(
    key=lambda x: x.tp.calculated_data[min_func_modified_label]
)

new_dir = os.path.dirname(pkl_file) + "/best_candidates"

mkdir(new_dir)

file_prefix = new_dir + "/best_candidate_"

for cand_id, cand in enumerate(resorted_best_candidates):
    file_basename = file_prefix + str(cand_id)
    xyz_name = file_basename + ".xyz"
    image_name = file_basename + ".png"

    extra_string = ""
    tp = cand.tp
    draw_chemgraph_to_file(
        tp.egc.chemgraph, image_name, size=(600, 400), bw_palette=False
    )
    for val_name, val in tp.calculated_data.items():
        if isinstance(val, float) or isinstance(val, int):
            extra_string += val_name + "=" + str(val) + " "

    coord_info = morfeus_coord_info_from_tp(tp, num_attempts=128)
    coords = coord_info["coordinates"]
    if coords is None:
        xyz_output = open(xyz_name, "w")
        print(cand_id, tp, cand.func_val, file=xyz_output)
        xyz_output.close()
    else:
        write_xyz_file(
            coords,
            xyz_name,
            nuclear_charges=coord_info["nuclear_charges"],
            extra_string=extra_string[:-1],
        )
