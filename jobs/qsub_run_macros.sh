#!/bin/bash
qsub -cwd -q citron.q@sge234 run_macros.sh 1
qsub -cwd -q citron.q@sge234 run_macros.sh 0