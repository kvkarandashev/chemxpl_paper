#!/bin/bash

export OMP_NUM_THREADS=1

for rf in */data_*/restart_file_*.pkl
do
    if [ ! -f $(dirname $rf)"/best_candidates/best_candidate_0.xyz" ]
    then
        python $(dirname $0)/print_best_candidates.py $rf
    fi
done
