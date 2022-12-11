from bmapqml.chemxpl.plotting import Analyze
from folder_translator import folder_name_param_dict
import os

all_simulations = [x[0] for x in os.walk("/data/jan/konstantin", topdown=False)]

#order the paths by folder hierarchy
all_simulations = sorted(all_simulations, key=lambda x: len(x.split("/")), reverse=True)

# pdb.set_trace()
for result_path in all_simulations:
    print(result_path)
    # result_path = "/data/jan/konstantin/xTB_dipsolv_opt_1_weak_strong_dipole/data_xTB_dipsolv_opt_1_weak_strong_dipole_5"
    sim_info = folder_name_param_dict(result_path)
    sim_name = "/data/jan/konstantin_plots/{}_{}_{}_{}_{}".format(
        sim_info["dataset"],
        sim_info["quant"],
        sim_info["gap_constraint"],
        sim_info["bias_strength"],
        sim_info["seed"],
    )
    best_ref_val, gap_constr_val = sim_info["best_ref_val"], sim_info["gap_constr_val"]

    ana = Analyze(
        "{}/restart_file*".format(result_path), quanity=sim_info["quant"], verbose=True
    )
    ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.parse_results()
    PARETO_CORRECTED = ana.pareto_correct(GLOBAL_HISTOGRAM)
    ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val)


# set
