#!/bin/bash
source ~/.bashrc
number=1.5
# for i in {1..18}; do #in iteration 13-15
#   # Increase scale by 0.1 every 3 iterations
#   if (( (i - 1) % 3 == 0 && i != 1 )); then
#     scale=$(echo "$scale + 0.1" | bc)
#   fi

#   # Execute the ROOT commands
#   root -b -q "data_gen.cxx(${scale},${i})"
#   root -b -q "Lagrange_method.cxx(${i})"
# done
for i in {1..10}; do #in iteration 13-15
  # Execute the ROOT commands
  root -b -q "data_gen.cxx(${number},${i})"
  root -b -q "Lagrange_method.cxx(${i})"
done
