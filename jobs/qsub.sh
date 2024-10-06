#!/bin/bash
path_to_python="/usr/bin/python3"

count=0
while [ $count -le 159 ]; do
    if [ $count -le 79 ]; then
        qsub -cwd -q citron.q@sge234 -S $path_to_python run.py $count $1 &
    else
        qsub -cwd -q citron.q@sge233 -S $path_to_python run.py $count $1 &
    fi
    count=$((count + 1))
done
wait

