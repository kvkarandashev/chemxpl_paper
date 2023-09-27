#!/bin/bash

script_version=1

function run_for () {
    dataset=$1
    python print_opt_log_$script_version.py /store/konst/chemxpl_related/$folder $dataset
}

output_dir_name=running_minima

mkdir -p $output_dir_name

folder=/store/konst/chemxpl_related/minimization_runs_xTB_dipole_solvation_cheap_post_sub/qm9_randomized_init_none_strong_solvation_energy/

for seed in $(seq 1 8)
do
    output_file=$output_dir_name/running_min_qm9_rand_init_$seed.txt
    python ../print_opt_log.py $folder/data_*_$seed/restart_file_$seed.pkl $output_file
done
