#!/bin/bash

script_version="post_sub_13"

function run_for () {
    dataset=$1
    if [ "$dataset" == "EGP" ]
    then
        folder=minimization_runs_xTB_dipole_solvation_cheap_leruli/xTB_dipsolv_opt_egp_cheap_3
    else
        folder=minimization_runs_xTB_dipole_solvation_cheap/xTB_dipsolv_opt_cheap_3
    fi
    python print_all_folder_min_summary_$script_version.py /store/konst/chemxpl_related/$folder $dataset
}

for dataset in QM9 EGP
do
    run_for $dataset
done

for dataset in QM9 EGP
do
    sed "s/\\midrule/\\midrule \\\multicolumn{4}{c}{"$dataset"*}\\\\\\\\ \\\midrule/g" summary_$dataset/cheap_quant_noise_$dataset.tex
done | awk '{
    if ($1=="\\bottomrule" && q==0) {n=-1}
    if (n==-1) {
        if ($1=="\\midrule") {n=0; q=1}
    }
    if (n==0) {print $0};
}' > cheap_quant_noise_table.tex
