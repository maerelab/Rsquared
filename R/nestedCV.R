#' Nested CV per feature by Bates2021
nestedCV = function(dat, nOuterFolds, cvSplits, nInnerFolds = nOuterFolds - 1, nCores = 1,loc = NULL, verbose = FALSE, saveFolder = NULL,
                    toBeSaved = FALSE, mc.preschedule = FALSE,...){
    out = mclapply(mc.cores = nCores, mc.preschedule = mc.preschedule, seq_len(ncol(dat$x)), function(j){
        if(verbose && (j %% 50)==0) cat("Feature", j, "\t")
        if(!toBeSaved ||  (!is.null(saveFolder) && !file.exists(saveFile <- paste0(saveFolder, "/cvSplitReps", j, ".RData")))){
            #Repeat splitting into folds
            cvSplitReps = lapply(seq_len(cvSplits), function(cvs){
                folds = sample(rep(unFolds <- seq_len(nOuterFolds), length.out = nrow(dat$x)))
                lapply(unFolds, function(uf){
                    idTrain = folds!=uf
                    yTrain = dat$y[idTrain]; xTrain <- dat$x[idTrain,j];locTrain = loc[idTrain,]
                    predTest = predLin(xTrain, yTrain, dat$x[!idTrain, j], loc = loc[idTrain,], ...)
                    eOut = (predTest-dat$y[!idTrain])^2
                    eBarOut = mean(eOut)
                    b = var(eOut)/length(eOut)
                    inFolds = sample(rep(unFoldsIn <- seq_len(nInnerFolds), length.out = sum(idTrain)))
                    errHatTilde = mean(vapply(FUN.VALUE = double(1), unFoldsIn, function(inf){
                        idTrainIn = inFolds!=inf
                        predTestIn = predLin(xTrain[idTrainIn], yTrain[idTrainIn], xTrain[!idTrainIn], loc = locTrain[idTrainIn,],...)
                        mean((predTestIn-yTrain[!idTrainIn])^2)
                    }))
                    a = (errHatTilde-eBarOut)^2
                    list("a" = a, "b" = b, "errHatTilde" = errHatTilde, "eOut" = eOut, "margVar" = var(dat$y[!idTrain]))
                })
            })
            seNested = getSEsNested(cvSplitReps, nOuterFolds = nOuterFolds, n = nrow(dat$x))
            if(toBeSaved){save(seNested, file = saveFile)}
        } else if(toBeSaved){
            load(saveFile)
        }
        seNested
    })
    names(out) = colnames(dat$x)
    return(out)
}
#Build confidence intervals
getSEsNested = function(cvSplitReps, nOuterFolds, n){
    ErrNCV = mean(sapply(cvSplitReps, function(y) sapply(y, function(x) x[["errHatTilde"]])))
    MSEhat = mean(sapply(cvSplitReps, function(y) sapply(y, function(x) x[["a"]]))) -
        mean(sapply(cvSplitReps, function(y) sapply(y, function(x) x[["b"]])))
    errOuter0 = lapply(cvSplitReps, function(y) lapply(y, function(x) x[["eOut"]]))
    mseOuter = sapply(errOuter0, function(w) sapply(w, mean))
    errOuter = unlist(errOuter0)
    margVars = sapply(cvSplitReps, function(y) sapply(y, function(x) x[["margVar"]]))
    SEest = sqrt(max(0, nOuterFolds/(nOuterFolds-1)*MSEhat))
    naiveRMSE = sd(errOuter)/sqrt(n)
    maxMSE = naiveRMSE * sqrt(nOuterFolds)
    if(is.na(SEest) || (SEest < naiveRMSE)){ #See below equation (17), prevent implausible values
        SEest = naiveRMSE
    } else if (SEest > maxMSE){
        SEest = maxMSE
    }
    #Correct the bias
    ErrCV = mean(errOuter)
    Bias = (1+(nOuterFolds-2)/nOuterFolds)*(ErrNCV-ErrCV)
    ErrNCVBC = ErrNCV - Bias#Bias correction
    cbind("MSEhat" = c("Naive" = ErrCV, "Bates" = ErrNCVBC),
          "SE" = c("Naive" = naiveRMSE, "Bates" = SEest),
          "CovAndCor" = c("Covariance" = cov(c(mseOuter), c(margVars)), "Correlation" = cor(c(mseOuter), c(margVars))))
}
nestedCVglmnet = function(dat, nOuterFolds, cvSplits, nInnerFolds = nOuterFolds - 1, nCores = 1, alpha, loc = NULL, ...){
    #Repeat splitting into folds
    cvSplitReps = mclapply(mc.cores = nCores, seq_len(cvSplits), function(cvs){
        cat("cvSplit", cvs, "\t")
        folds = sample(rep(unFolds <- seq_len(nOuterFolds), length.out = nrow(dat$x)))
        tmp = lapply(unFolds, function(uf){
            #cat("Fold", uf, "\t")
            idTrain = folds!=uf
            yTrain = dat$y[idTrain]; xTrain <- dat$x[idTrain,]
            predTest = predGlmnet(xTrain, yTrain, dat$x[!idTrain, ], nfolds = max(4,nInnerFolds-1),
                                  alpha = alpha, loc = loc[idTrain, ], ...)
            eOut = (predTest-dat$y[!idTrain])^2
            eBarOut = mean(eOut)
            b = var(eOut)/length(eOut) #Also like this in nestedcv package
            inFolds = sample(rep(unFoldsIn <- seq_len(nInnerFolds), length.out = sum(idTrain)))
            errHatTilde = mean(unlist(sapply(unFoldsIn, function(inf){
                idTrainIn = inFolds!=inf
                predTestIn = predGlmnet(xTrain[idTrainIn,], yTrain[idTrainIn], xTrain[!idTrainIn,],
                                     nfolds = max(4,nInnerFolds-1), alpha = alpha, loc = loc[idTrain,][idTrainIn, ], ...)
                (predTestIn-yTrain[!idTrainIn])^2
            })))
            a = (errHatTilde-eBarOut)^2
            list("a" = a, "b" = b, "errHatTilde" = errHatTilde, "eOut" = eOut, "margVar" = var(dat$y[!idTrain]))
        })
    })
    getSEsNested(cvSplitReps, nOuterFolds = nOuterFolds, n = nrow(dat$x))
}
