from bmapqml.chemxpl.plotting import Analyze
from bmapqml.chemxpl.plotting import make_pretty
from folder_translator import folder_name_param_dict
import os
import umap
from bmapqml.chemxpl import rdkit_descriptors
from rdkit import Chem
import pandas as pd
import pdb
import matplotlib.pyplot as plt
import numpy as np

def compute_representations(MOLS, nBits):
    """
    Compute the representations of all unique smiles in the random walk.
    """

    X = rdkit_descriptors.get_all_FP(MOLS, nBits=nBits, fp_type="MorganFingerprint")
    return X
    
def compute_UMAP(MOLS, nBits=1024):
    #convert smiles to rdkit mol objects
    MOLS = [Chem.MolFromSmiles(x) for x in MOLS]
    #add hydrogens
    MOLS = [Chem.AddHs(x) for x in MOLS]

    """
    Compute PCA
    """

    X = compute_representations(MOLS, nBits=nBits)
    reducer = umap.UMAP(random_state=42)
    reducer.fit(X)
    X_2d = reducer.transform(X)
    return X_2d, reducer

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

    qm9_data = pd.read_csv('qm9.csv')


    SMILES_QM9, qm9_ind, GAP, MU = qm9_data['canon_smiles'].values,qm9_data['GDB17_ID'].values, qm9_data['gap'].values, qm9_data['mu'].values
    len_qm9 = len(SMILES_QM9)
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
            

            plt.close("all")
            fs = 24

            plt.rc("font", size=fs)
            plt.rc("axes", titlesize=fs)
            plt.rc("axes", labelsize=fs)
            plt.rc("xtick", labelsize=fs)
            plt.rc("ytick", labelsize=fs)
            plt.rc("legend", fontsize=fs)
            plt.rc("figure", titlesize=fs)

            fig1, ax1 = plt.subplots(figsize=(10,7))

            #ALL_HISTOGRAMS, GLOBAL_HISTOGRAM, ALL_TRAJECTORIES = ana.parse_results()
            #histogram = ALL_HISTOGRAMS[0]
            #histogram.to_csv("histograms.csv")
            histogram = pd.read_csv("histograms.csv")
            COMBINED_SMILES = np.concatenate((  SMILES_QM9,histogram["SMILES"].values))
            #compute UMAP
            #pdb.set_trace()
            X_2d, reducer = compute_UMAP(COMBINED_SMILES)
            ax1.set_xlabel("UMAP 1", fontsize=20)
            ax1.set_ylabel("UMAP 2", fontsize=20)
            #save figure
            
            ax1.scatter(X_2d[:,0][:len_qm9][::10], X_2d[:,1][:len_qm9][::10], c="grey", s=50, marker="x",alpha=0.05,  label="QM9")
            ax1.scatter(X_2d[:,0][len_qm9:][::10], X_2d[:,1][len_qm9:][::10], c="r", s=50, marker="o",alpha=0.05, label="MOSAiCS")
            ax1.legend()
            ax1 = make_pretty(ax1)
            plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

            plt.savefig(f"UMAP.pdf", dpi=300)
            plt.savefig(f"UMAP.png", dpi=300)
            plt.savefig(f"UMAP.svg")

            plt.close("all")
            exit()