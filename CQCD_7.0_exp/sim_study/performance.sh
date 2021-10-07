#!/bin/bash
for nd in 500_100 500_1000 5000_100 5000_1000 25000_100 25000_1000
do
    for mtd in bridge_FB_CQCD Fisher_scoring forward_select Yi2011
    do
        echo "Rscript --vanilla performance.R ${mtd}_${nd}.RData data_35000.RData table.csv"
        Rscript --vanilla performance.R ${mtd}_${nd}.RData data_35000.RData table.csv
        echo "Rscript --vanilla performance.R ${mtd}_dep_${nd}.RData data_dep_35000.RData table.csv"
        Rscript --vanilla performance.R ${mtd}_dep_${nd}.RData data_dep_35000.RData table.csv
    done
done 
