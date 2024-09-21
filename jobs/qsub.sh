#!/bin/bash
path_to_python="/usr/bin/python3"


count=0
while [ $count -le 83 ]; do
    qsub -cwd -q citron.q@sge234 -S $path_to_python run.py $count $1
    count=$((count + 1))
done
