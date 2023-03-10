#The oob bootstrap (smooths leave-one-out CV)
bootOob = function(dat, id, id0, meanFit){
    id2 = id0[-id]
    Eis = matrix(0, length(id), p <- NCOL(dat$x))
    Nis = vapply(id0, FUN.VALUE = integer(1), function(x) sum(x==id))
    Eis[id2,] = vapply(FUN.VALUE = double(length(id2)), seq_len(p), function(j){
        predTest = predLin(dat$x[id,j], dat$y[id], dat$x[id2,j])
        (predTest-dat$y[id2])^2
    })
    cbind(Eis, "Nis" = Nis)
}
processOob = function(oobObj){
    procOb = processOobIn(oobObj$obsBoot) #Process observed and bootstrap
    procBoot = lapply(oobObj$bootBoot, processOobIn)
    list('procOb' = procOb, 'procBoot' = procBoot)
}
processOobIn = function(x){
    Nmat = sapply(x, function(y) y[, "Nis"])
    n = nrow(Nmat)
    Imat = Nmat==0
    rI = rowSums(Imat)
    IQmat = vapply(FUN.VALUE = Nmat, seq_len(ncol(x[[1]])-1), function(j){
        vapply(FUN.VALUE = double(n), x, function(y){y[,j]})*Imat
    })
    Eis = apply(IQmat, c(1,3), sum)/rI
    errEsts = colMeans(Eis)
    # Following Efron1997, equation (40)
    qMat = colMeans(IQmat)
    Dis = (2+1/(n-1))*(Eis-errEsts)/n + ((Nmat-rowMeans(Nmat)) %*% qMat)/rI
    seEsts = sqrt(colSums(Dis^2))
    cbind("MSEhat" = errEsts, "SEhat" = seEsts, "SEhatNaive" = apply(Eis, 2, sd)/sqrt(n))
}
bootOobGlmnet = function(dat, id, id0, alpha){
    id2 = id0[-id]
    Eis = double(length(id))
    Nis = vapply(id0, FUN.VALUE = integer(1), function(x) sum(x==id))
    Eis[id2] = {
        predTest = predGlmnet(dat$x[id,], dat$y[id], dat$x[id2,], alpha = alpha)
        (predTest-dat$y[id2])^2
    }
    cbind(Eis, "Nis" = Nis)
}
processOobGlmnet = function(oobObj){
    procOb = processOobInGlmnet(oobObj$obsBoot) #Process observed and bootstrap
    procBoot = lapply(oobObj$bootBoot, processOobIn)
    list('procOb' = procOb, 'procBoot' = procBoot)
}
processOobInGlmnet = function(x){
    Nmat = sapply(x$bootRes, function(y) y[, "Nis"])
    n = nrow(Nmat)
    Imat = Nmat==0
    rI = rowSums(Imat)
    IQmat = vapply(FUN.VALUE = Nmat, seq_len(ncol(x$bootRes[[1]])-1), function(j){
        vapply(FUN.VALUE = double(n), x$bootRes, function(y){y[,j]})*Imat
    })
    Eis = apply(IQmat, c(1,3), sum)/rI
    errEsts = colMeans(Eis)
    # Following Efron1997, equation (40)
    qMat = colMeans(IQmat)
    Dis = (2+1/(n-1))*(Eis-errEsts)/n + ((Nmat-rowMeans(Nmat)) %*% qMat)/rI
    seEsts = sqrt(colSums(Dis^2))
    cbind("MSEhat" = errEsts, "SEhat" = seEsts, "SEhatNaive" = apply(Eis, 2, sd)/sqrt(n))
}
