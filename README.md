
This repository contains all code for running simulations and analyses
for the paper “The out-of-sample $R^2$: estimation and inference” by
Stijn Hawinkel, Willem Waegeman and Steven Maere. A short demo on how to
calculate the $R^2$ for ordinary least squares (OLS) is given below.

``` r
n = 100 #Sample size
x = rnorm(n) #Regressor
y = rnorm(n, x*0.75) #Outcome variable
df = list("x" = cbind(x), "y" = y) #Make dataframe
```

``` r
for(i in list.files("R", full.names = TRUE)){source(i)};rm(i) #Source all R-code
```

``` r
library(parallel) #Load the necessary packages
```

``` r
cvSplits = 200 #Number of cross-validation splits
nFolds = 10 #Number of CV folds
nCores = 4 #Number of cores to use in multithreading
nested_cv_result = nestedCV(df, nFolds, cvSplits = cvSplits, nCores = nCores)
nested_cv_result[[1]]["Bates", "MSEhat"] #The estimated MSE
```

    ## [1] 0.9475463

``` r
nested_cv_result[[1]]["Bates", "SE"] #Corresponding SE
```

    ## [1] 0.1510203

``` r
#Estimate correlation between MSE and MST using bootstrap
bootReps = 100
bootEsts = bootCV(df, nFolds, cvSplits = cvSplits, bootReps = bootReps, nCores = nCores)
matBootEsts = simplify2array(bootEsts)[1,,] #Convert to matrix
rhoEstimate = cor(matBootEsts["MSE",], matBootEsts["margVar",])
rhoEstimate #Considerable correlation
```

    ## [1] 0.7174728

``` r
MST = var(y)*(length(y)+1)/length(y) #The MST
seObj = RsquaredSE(MSE = nested_cv_result[[1]]["Bates", "MSEhat"], 
                   SEMSE = nested_cv_result[[1]]["Bates", "SE"],
                   margVar = MST, 
                    rhoBoot = rhoEstimate, 
                   n = length(y))
round(seObj, 4) #The resulting R2 and standard error
```

    ##     R2   R2SE 
    ## 0.2655 0.0841
