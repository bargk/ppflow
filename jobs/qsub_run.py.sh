#!/bin/bash


count=0
while [ $count -le 159 ]; do
    if [ $count -le 79 ]; then
        qsub -cwd -q citron.q@sge234 -S /bin/bash wrapper_run.py.sh $count $1 &
    else
        qsub -cwd -q citron.q@sge233 -S /bin/bash wrapper_run.py.sh $count $1 &
    fi
    count=$((count + 1))
done
wait

