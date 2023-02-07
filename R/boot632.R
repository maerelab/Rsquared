#The .632 bootstrap
boot632 = function(dat, id, id0, meanFit){
            predTest = predLin(dat$x[id,], dat$y[id], dat$x[-id,])
            ErrOutOfSample = mean((predTest-dat$y[-id])^2) #Out of sample error
            ErrInSample = mean((predLin(dat$x, dat$y, dat$x)-dat$y)^2) #In sample error
            bootRes = c("ErrOutOfSample" = ErrOutOfSample, "ErrInSample" = ErrInSample)
            list("bootRes" = bootRes, "margVar" = var(dat$y))
        }
process632 = function(x, oobObj){
    procOb = process632In(x$obsBoot, oobObj = oobObj)
    procBoot = sapply(x$bootBoot, process632In)
    procBootParam = sapply(x$bootBootParam, process632In)
    procJack = sapply(x$jackBoot, process632In)
    margVarBoot = sapply(x$bootBoot, function(y) y[[1]]$margVar)
    margVarBootParam = sapply(x$bootBootParam, function(y) y[[1]]$margVar)
    margVarJack = sapply(x$jackBoot, function(y) y[[1]]$margVar)
    list('procOb' = procOb, 'procBoot' = procBoot, "procBootParam" = procBootParam, "procJack" = procJack,
         "margVarBoot" = margVarBoot, "margVarBootParam" = margVarBootParam, "margVarJack" = margVarJack)
}
process632In = function(x, oobObj = NULL){
    bootEsts = matrix(sapply(x, function(y) sum(y$bootRes*expvec)), ncol = length(x))
    MSEhat = rowMeans(bootEsts)
    if(!is.null(oobObj)){
        SEhatNaive = sd(bootEsts)
        SEhat = MSEhat/oobObj[, "MSEhat"]*oobObj[, "SEhat"]
    } else {
        SEhat = SEhatNaive = NULL
        }
    rbind("MSEhat" = MSEhat, "SEhat" = SEhat, "SEhatNaive" = SEhatNaive)
}
#Glmnet
boot632Glmnet = function(dat, id, ...){
    predTest = predGlmnet(dat$x[id,], dat$y[id], dat$x[-id,], ...)
    ErrOutOfSample = mean((predTest-dat$y[-id])^2) #Out of sample error
    ErrInSample = mean((predGlmnet(dat$x, dat$y, dat$x, ...)-dat$y)^2) #In sample error
    sum(c("ErrOutOfSample" = ErrOutOfSample, "ErrInSample" = ErrInSample)*expvec)
}
process632Glmnet = function(x, oobObj){
    procOb = process632InGlmnet(x$obsBoot, oobObj = oobObj)
    procBoot = sapply(x$bootBoot, process632InGlmnet)
    procBootParam = sapply(x$bootBootParam, process632InGlmnet)
    procJack = sapply(x$jackBoot, function(z){
        sum(z$bootRes*expvec)
    })
    margVarBoot = sapply(x$bootBoot, function(y) y$margVar)
    margVarBootParam = sapply(x$bootBootParam, function(y) y$margVar)
    margVarJack = sapply(x$jackBoot, function(y) y$margVar)
    list('procOb' = procOb, 'procBoot' = procBoot, "procBootParam" = procBootParam, "procJack" = procJack,
         "margVarBoot" = margVarBoot, "margVarBootParam" = margVarBootParam, "margVarJack" = margVarJack)
}
process632InGlmnet = function(x, oobObj = NULL){
    MSEhat = mean(tmp <- unlist(x$bootRes))
    if(!is.null(oobObj)){
        SEhatNaive = sd(tmp)
        SEhat = MSEhat/oobObj[, "MSEhat"]*oobObj[, "SEhat"]
    } else {
        SEhat = SEhatNaive = NULL
    }
    rbind("MSEhat" = MSEhat, "SEhat" = SEhat, "SEhatNaive" = SEhatNaive)
}