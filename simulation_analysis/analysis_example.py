from bmapqml.chemxpl.plotting import Analyze
from folder_translator import folder_name_param_dict
import os
import pdb

if __name__ == '__main__':

    DATASET = "QM9"
    PLOT_ALL = False

    if PLOT_ALL:
        all_simulations = [x[0] for x in os.walk("/data/jan/konstantin/{}".format(DATASET), topdown=False)]



        #order the paths by folder hierarchy
        all_simulations = sorted(all_simulations, key=lambda x: len(x.split("/")))
        all_simulations = all_simulations[::-1]

    else:
        
        #read from file best_seeds.txt line by line
        all_simulations = [] 
        with open("best_seeds_{}.txt".format(DATASET), "r") as f:
            for line in f:
                line = line.strip()
                if line != "None":
                    all_simulations.append(line)


    for result_path in all_simulations:
        print(result_path)
        if "atomization" not in result_path:
            sim_info = folder_name_param_dict(result_path)
            sim_name = "/data/jan/konstantin_plots/{}_{}_{}_{}_{}".format(
                sim_info["dataset"],
                sim_info["quant"],
                sim_info["gap_constraint"],
                sim_info["bias_strength"],
                sim_info["seed"],
            )
            best_ref_val, gap_constr_val = sim_info["best_ref_val"], sim_info["gap_constr_val"]
            print(best_ref_val, gap_constr_val)
            ana = Analyze(
                "{}/restart_file*".format(result_path), quantity=sim_info["quant"], verbose=True,full_traj=False
            )
            ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.parse_results()
            pdb.set_trace()
            PARETO_CORRECTED = ana.pareto_correct(GLOBAL_HISTOGRAM)
            
            ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val, coloring="encounter")
            ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val, coloring="density")
            
            