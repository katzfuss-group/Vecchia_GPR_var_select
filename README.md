# Vecchia Gaussian Process Regression and Variable Selection

## Prerequisite

We have made extensions to the `GpGp` R package, which have not been pulled into the main GpGp repo. To install our version of `GpGp` package with the necessary extensions, please clone and install the `squared_relevance` branch of [this repository](https://github.com/katzfuss-group/GpGp.git).

## Numerical results

To reproduce the figures and results in Cao, Guinness, Genton, Katzfuss (2022), please enter the corresponding directory and run the `.R` or `.Rmd` file. 

For `GP_approx_cmp`, first run `Vecchia_FIC_PIC_compare.Rmd` until the `Plot` Section, then run `gpflow_SVGP.ipynb`, and finally run the `Plot` Section of `Vecchia_FIC_PIC_compare.Rmd` because the FITC method is implemented in the `GPytorch` python package.

For `real_dataset`, the raw data should be first downloaded (e.g., from UCI) to the `data` folder. Then, run the corresponding section of `dataset_prep.Rmd` to pre-process the data. Finally, change the `datasetName` variable in the `main.Rmd` to run the real-data application study.

For the `temperature_data`, leading spaces should first be removed with a text editor.

For the sparse additive model implemented in `SAM` R package, the input is jittered slightly. Otherwise, some datasets, e.g., the slicing dataset, reported error.

## Example

### Preparation

Clean the environment.
```{r clean environment}
rm(list = ls())
```

Install a modified version of the ''GpGp'' R package.
```{r install modified GpGp}
devtools::install_github("https://github.com/katzfuss-group/GpGp.git", 
                         ref = "squared_relevance")
```

Load needed libraries and R files,
```{r libraries and self-defined funcs}
source("func/misc.R", chdir = T)
source("func/FBselect.R", chdir = T)
source("func/stop_cond.R", chdir = T)
source("func/locs_gen.R", chdir = T)
source("func/y_gen.R", chdir = T)
source("func/Vecc_path.R", chdir = T)
source("func/penalty.R", chdir = T)
```

Generate `locs`. `locs` should be a matrix whose rows are covariate (feature) vectors of the responses. `n` is the number of reponses and `d` is the number of covariates.
```{r gen locs}
set.seed(123)
n <- 6400
d <- 1e3
locs <- locs_gen_idp(n = n, d = d)
```

Generate `y`. `y` should be a vector representing the responses.  `sigmaSq0`, `tauSq0`, and `rSq0` are the true parameter values in the Matern kernel with smoothness 2.5.
```{r gen y}
set.seed(123)
sigmaSq0 <- 1
tauSq0 <- 0.05^2
rSq0 <- c(10, 5, 2, 1, 0.5, rep(0, d - 5))^2
covFn <- "matern25_scaledim_sqrelevance"
y <- MVN_gen(locs = locs, parms = c(sigmaSq0, rSq0, tauSq0), covFn = covFn)
y <- y - mean(y)
```

Define the ratio of in-sample data to total data.
```{r in-sample ratio}
pIn <- 0.75
spltObj <- splt_locs_y(locs, y, pIn)
locsTrn <- spltObj$locs1
locsTst <- spltObj$locs2
yTrn <- spltObj$y1
yTst <- spltObj$y2
```

### VREG

Define penalty function, the iterative adapative bridge penalty. `dpen_fun` and `ddpen_fun` are first and second order information.
```{r pen func}
pen_fun <- pen_avg_brdg
dpen_fun <- dpen_avg_brdg
ddpen_fun <- ddpen_avg_brdg
```

Optimization parameters. `minPosi` is a small positive number to prevent `var`, `nugget` and the sum of `sr` reaching zero. `k` is the number of covariates selected each time. `miniQCCD` and `miniGrad` are batch sizes used in QCCD and variable selection. `nAvg` is the kappa in the penalty function. `m` is the conditioning set size in Vecchia approximation. `mini` is whether to use mini-batch subsampling or not. `taper` is whether to use the iterative adaptive bridge penalty or the bridge penalty. `lambVec` stores the lambda values for model estimation on the regularization path. `stop_con_path` is the stopping conditions based on 1% improvement of the OOS score and the inclusion of new covariate.
```{r some opt parms, eval = T}
var0 <- 1
tau0 <- 1e-4
sr0 <- 1e-2
minPosi <- 1e-8
m <- 100
k <- 3
miniQCCD <- 128
miniGrad <- 128
nAvg <- 2
lambVec <- 2^(rev(-6 : floor(log2(n))))
convQCCD <- 1e-4
convCCD <- 1e-4
cAmij <- 1e-4
maxIterQCCD <- 200
maxIterCCD <- 40
covFn <- "matern25_scaledim_sqrelevance"
mini <- T
taper <- T
silent <- F
stop_con_path <- function(scrVec, idxSet, thetaSet){
  min(stop_con_NEW_path(scrVec, idxSet, thetaSet), 
      stop_con_OOS1_path(scrVec, idxSet, thetaSet))
}
```

Call the function `Vecc_path` for the ''VGPR'' algorithm. `pIn` here specifies the ratio of the validation dataset to total training dataset, specifically, `1 - pIn`. `OOS_rmse` is the function for computing the OOS RMSE score. `stop_con_OOS1_fb` is the stopping conditions based on 1% improvement of the OOS score.
```{r Vecc_path parms est and var select}
bgnTime <- Sys.time()
VeccObj <-
  Vecc_path(
    var0 = var0,
    tau0 = tau0,
    sr0 = sr0,
    locs = locsTrn,
    y = yTrn,
    m = m,
    k = k,
    lambVec = lambVec,
    pen_fun = pen_fun,
    dpen_fun = dpen_fun,
    ddpen_fun = ddpen_fun,
    OOS_score = OOS_rmse,
    stop_con_path = stop_con_path,
    stop_con_fb = stop_con_OOS1_fb,
    miniQCCD = miniQCCD,
    miniGrad = miniGrad,
    pIn = pIn,
    spltSeed = 123,
    convQCCD = convQCCD,
    convCCD = convCCD,
    cAmij = cAmij,
    maxIterQCCD = maxIterQCCD,
    maxIterCCD = maxIterCCD,
    covFn = covFn,
    minPosi = minPosi,
    mini = mini,
    taper = taper,
    silent = silent,
    nAvg = nAvg
  )
endTime <- Sys.time()
timeUsed <- as.numeric(endTime - bgnTime, units = "secs")
cat("Vecchia solution path used", timeUsed, "seconds\n")
```

Evaluate posterior prediction performance.
```{r Vecc_path performance, eval = T}
VeccScr <- OOS_rmse(VeccObj$theta, locsTrn, locsTst, yTrn, yTst, m, covFn)
cat("Vecchia solution path selected: ", VeccObj$idx, "covariates\n")
cat("Vecchia solution path estimated theta:", 
    VeccObj$theta[c(1, VeccObj$idx + 1, d + 2)], "\n")
cat("Vecchia solution path RMSE score:", VeccScr, "\n")
cat("Vecchia solution path FPR:", 
    length(setdiff(VeccObj$idx, 1 : 5)) / length(VeccObj$idx), "\n")
```

## Reference
Cao, Guinness, Genton, and Katzfuss (2022). Scalable Gaussian-process regression and variable selection using Vecchia approximations. 
