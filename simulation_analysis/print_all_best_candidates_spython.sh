#!/bin/bash

MYDIR=$(pwd)

restart_files=( $(ls $MYDIR/*/data_*/restart_file*.pkl | tr '\n' ' ') )

for rf in ${restart_files[@]}
do
    rdir=$(dirname $rf)
    joblog=$rdir/joblog.txt
    if [ ! -f $joblog ]
    then
        if [ "$(ls $rdir/*.stdout_* | wc -l)" == "0" ]
        then
            continue
        fi
        joblog=$(ls $rdir/*.stdout_*)
    fi
    nsteps=$(awk '{if ($1 == "HIST") {if ($3 > n) {n = $3}}} END {print n}' $joblog)
    if [ "$nsteps" != 50000 ]
    then
        continue
    fi
    if [ ! -f $rdir"/best_candidates/best_candidate_0.xyz" ]
    then
        cd $rdir
#        log_name="$(dirname $rf)/analyze_$(basename $rf | cut -d'.' -f1).log"
        spython --CPUs=16 --OMP_NUM_THREADS=1 $(dirname $0)/print_best_candidates.py analyze_$(basename $(dirname $rf))  $rf #> $log_name
        cd $MYDIR
    fi
done
