from bmapqml.chemxpl.plotting import Analyze,Chem_Div
from folder_translator import folder_name_param_dict
import os
import pdb
import pickle

if __name__ == '__main__':
    load_res = False

    ind = 0
    DATASETS = ["EGP","QM9"]
    for DATASET in DATASETS:
        PLOT_ALL = False

        if PLOT_ALL:
            all_simulations = [x[0] for x in os.walk("/store/common/pareto/data/{}".format(DATASET), topdown=False)]



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

        n_steps_log =  open(f"/store/common/pareto/{DATASET}/log/steps.txt", "w")



        for result_path in all_simulations:
            print(result_path)
            if "atomization" not in result_path:
                sim_info = folder_name_param_dict(result_path)
                sim_name = "/store/common/pareto/{}_{}_{}_{}_{}".format(
                    sim_info["dataset"],
                    sim_info["quant"],
                    sim_info["gap_constraint"],
                    sim_info["bias_strength"],
                    sim_info["seed"],
                )
                best_ref_val, gap_constr_val = sim_info["best_ref_val"], sim_info["gap_constr_val"]
                print(best_ref_val, gap_constr_val)
                

                
                if load_res:
                    with open(f"/store/common/pareto/plot_data/RESULTS_{ind}.pkl", "rb") as f:
                        ana = pickle.load(f)
                
                    ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.ALL_HISTOGRAMS, ana.GLOBAL_HISTOGRAM, ana.ALL_TRAJECTORIES
                else:
                    ana = Analyze("{}/restart_file*".format(result_path), quantity=sim_info["quant"], verbose=True,full_traj=True)
                    ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.parse_results()
                    with open(f"/store/common/pareto/plot_data/RESULTS_{ind}.pkl", "wb") as f:
                        pickle.dump(ana, f)


                
                global_MC_step_counter = ana.global_MC_step_counter
                print(global_MC_step_counter)
                n_steps_log.write("{}\t{}\n".format(sim_name.split("/")[-1],global_MC_step_counter))
                PARETO_CORRECTED = ana.pareto(mode="trajectory")
                PARETO_CORRECTED.to_csv("/store/common/pareto/{}/log/{}.csv".format(DATASET,sim_name.split("/")[-1]))
                ALL_TRAJECTORIES.to_csv("/store/common/pareto/{}/log/traj_{}.csv".format(DATASET,sim_name.split("/")[-1]))

                ana.plot_pareto(sim_name, hline=gap_constr_val, vline=best_ref_val,dataset=DATASET ,coloring="density",plot_quantity = sim_info["quant"],labels=True )
                ind += 1
                
        n_steps_log.close()