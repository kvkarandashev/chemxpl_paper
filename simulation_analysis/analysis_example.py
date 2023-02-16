from bmapqml.chemxpl.plotting import Analyze,Chem_Div
from folder_translator import folder_name_param_dict
import os
import pdb

import pandas as pd

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

    #open file in append mode and create if not exits

    n_steps_log =  open("/data/jan/konstantin_plots/log/steps.txt", "w")



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
                "{}/restart_file*".format(result_path), quantity=sim_info["quant"], verbose=True,full_traj=True
            )
            
            

            ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.parse_results()
            
            global_MC_step_counter = ana.global_MC_step_counter
            print(global_MC_step_counter)
            n_steps_log.write("{}\t{}\n".format(sim_name.split("/")[-1],global_MC_step_counter))
            PARETO_CORRECTED = ana.pareto(mode="trajectory")
            PARETO_CORRECTED.to_csv("/data/jan/konstantin_plots/log/{}.csv".format(sim_name.split("/")[-1]))
            ALL_TRAJECTORIES.to_csv("/data/jan/konstantin_plots/log/traj_{}.csv".format(sim_name.split("/")[-1]))

            #time_ordered_smiles = ana.time_ordered_smiles
            #dump time ordered smiles to np file
            #np.save("/data/jan/konstantin_plots/log/smiles_{}.npy".format(sim_name.split("/")[-1]), time_ordered_smiles)
            #ana_chem = Chem_Div(traj = time_ordered_smiles, subsample = 1000, verbose = True)
            #ana_chem.compute_representations()
            #ana_chem.compute_diversity()
            #plt.plot(ana_chem.N, ana_chem.diversity)
            #plt.show()
            #pdb.set_trace()
            ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val,dataset=DATASET ,coloring="encounter")
            ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val,dataset=DATASET ,coloring="density")
            
    n_steps_log.close()