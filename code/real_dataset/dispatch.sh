#!/bin/bash

fnLst='slice_localization_data slice_localization_data_f615 slice_localization_data_f615_dep CASP_f991 CASP_f991_dep temp temp_f894 temp_f894_dep'
for fn in $fnLst
do
    echo $fn
    Rscript main.R $fn
done
