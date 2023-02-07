#'Jackknife the R2 for BCa
jackKnifeCV = function(dat, outerFolds, cvSplits){
    vapply(seq_along(dat$y), FUN.VALUE = double(2), function(i){
        dat$x = dat$x[-i,,drop = FALSE];dat$y = dat$y[-i]
        MSE = mean(unlist(sapply(integer(cvSplits), function(foo){
            simpleCV(dat, nOuterFolds)
        })))
        c("MSE" = MSE, "margVar" = var(dat$y))
    })
}
