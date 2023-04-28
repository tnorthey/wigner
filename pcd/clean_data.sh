#!/bin/bash

paste ../distances/r56.dat inv_chi2_sample.dat > tmp
sed '/e-02/d' tmp > tmp2
sed '/e-01/d' tmp2 > cleaned_data.dat
# or just 
#sed '/e-02/d' tmp > cleaned_data.dat
rm tmp tmp2
