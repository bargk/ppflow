#!/bin/bash
source ~/.bashrc
kinit -R
path_to_python="/gpfs0/citron/users/bargl/MyPythonInstall/install/bin/python"
$path_to_python run.py $1 $2
