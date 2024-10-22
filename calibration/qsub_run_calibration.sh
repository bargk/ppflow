#!/bin/bash
qsub -cwd -q citron.q@sge234 run_calibration.sh 1
qsub -cwd -q citron.q@sge234 run_calibration.sh 1.5
qsub -cwd -q citron.q@sge234 run_calibration.sh 2
