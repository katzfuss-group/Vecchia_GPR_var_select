# Vecchia Gaussian Process Regression and Variable Selection

To reproduce the figures and results in Cao, Guinness, Genton, Katzfuss (2022), please enter the corresponding directory and run the `.R` or `.Rmd` file. 

For `GP_approx_cmp`, first run `Vecchia_FIC_PIC_compare.Rmd' until the 'Plot' Section, then run `gpflow_SVGP.ipynb`, and finally run the 'Plot' Section of `Vecchia_FIC_PIC_compare.Rmd' because the FITC method is implemented in the 'GPytorch' python package.

For `real_dataset`, the raw data should be first downloaded (e.g., from UCI) to the `data` folder. Then, run the corresponding section of `dataset_prep.Rmd` to pre-process the data. Finally, change the `datasetName` variable in the `main.Rmd` to run the real-data application study.


## Reference
Cao, Guinness, Genton, and Katzfuss (2022). Scalable Gaussian-process regression and variable selection using Vecchia approximations. [*arXiv:2108.04211*](https://arxiv.org/abs/2108.04211).
