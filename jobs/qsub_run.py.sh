#!/bin/bash

count=0
while [ $count -le 231 ]; do  
    if [ $count -le 79 ]; then
        qsub -cwd -q citron.q@sge234 -S /bin/bash wrapper_run.py.sh $count $1 &
    elif [ $count -le 159 ]; then
        qsub -cwd -q citron.q@sge233 -S /bin/bash wrapper_run.py.sh $count $1 &
    else  #only submitting 72 jobs to sge1045
        qsub -cwd -q citron.q@sge1045 -S /bin/bash wrapper_run.py.sh $count $1 &
    fi
    count=$((count + 1))
done
wait

wait
