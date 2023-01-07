#!/bin/bash

export OMP_NUM_THREADS=1

restart_files=( $(ls */data_*/restart_file*.pkl | tr '\n' ' ') )

for rf in ${restart_files[@]}
do
    if [ ! -f $(dirname $rf)"/best_candidates/best_candidate_0.xyz" ]
    then
        log_name="$(dirname $rf)/analyze_$(basename $rf | cut -d'.' -f1).log"
        python $(dirname $0)/print_best_candidates.py $rf > $log_name
    fi
done
