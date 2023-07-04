#!/bin/bash

script_version=1

function run_for () {
    dataset=$1
    python print_opt_log_$script_version.py /store/konst/chemxpl_related/$folder $dataset
}

output_dir_name=running_minima

mkdir -p $output_dir_name

for dataset in QM9 EGP
do
    if [ "$dataset" == "EGP" ]
    then
        folder=minimization_runs_xTB_dipole_solvation_cheap_leruli/xTB_dipsolv_opt_egp_cheap_3
    else
        folder=minimization_runs_xTB_dipole_solvation_cheap/xTB_dipsolv_opt_cheap_3
    fi
    for quant in solvation_energy dipole
    do
        for gap_constr in weak strong
        do
            for bias in none weak stronger
            do
                for seed in $(seq 1 8)
                do
                    output_file=$output_dir_name/running_min_${dataset}_${quant}_${gap_constr}_${bias}_$seed.txt
                    end_str=_${bias}_${gap_constr}_${quant}
                    python print_opt_log.py /store/konst/chemxpl_related/$folder*$end_str/data_*${end_str}_$seed/restart_file_$seed.pkl $output_file
                done
            done
        done
    done
done

