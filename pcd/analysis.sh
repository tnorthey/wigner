#!/bin/bash

ij=$1

paste ../distances/r"$ij".dat inv_chi2_sample.dat > r"$ij"_data.dat
python3 digitize.py r"$ij"_data.dat

mv out.dat r"$ij"_acc.dat
