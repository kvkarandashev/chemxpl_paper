#!/bin/bash

function run_for () {
    dataset=$1
    if [ "$dataset" == "EGP" ]
    then
        folder=minimization_runs_xTB_dipole_solvation_cheap_leruli/xTB_dipsolv_opt_egp_cheap_3
    else
        folder=minimization_runs_xTB_dipole_solvation_cheap/xTB_dipsolv_opt_cheap_3
    fi
    python print_all_folder_min_summary_8.py /store/konst/chemxpl_related/$folder $dataset
}

for dataset in QM9 EGP
do
    run_for $dataset
done
