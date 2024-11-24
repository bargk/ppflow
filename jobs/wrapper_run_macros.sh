#!/bin/bash

source ~/.bashrc 
macro_dir="/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/"

for i in "$macro_dir"/S01* "$macro_dir"/S02* "$macro_dir"/S03* "$macro_dir"/S04* "$macro_dir"/S05* "$macro_dir"/S06* "$macro_dir"/S07* "$macro_dir"/S08*
do
    root -b "$i""($1)"
done