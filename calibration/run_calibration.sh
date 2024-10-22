#!/bin/bash

root -b data_gen.cxx"($1)"
root -b Lagrange_method.cxx"($1, 0)"
root -b LG_graphs.cxx"($1, 0)"
root -b Lagrange_method.cxx"($1, 1)"
root -b LG_graphs.cxx"($1, 1)"