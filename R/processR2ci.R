processR2ci = function(xxx, oracleSd, oracleRho, quants = c(sigLevel/2, 1-sigLevel/2)){
    with(xxx, {
        rhosMat = cbind("bootRhos" = vapply(FUN.VALUE = double(1), seq_along(seEsts), function(x) {
            out = sapply(bootEsts, function(w){w[x, ]})
            cor(out[1,], out[2,])
        }),
        "bootRhosParam" = vapply(FUN.VALUE = double(1), seq_along(seEsts), function(x) {
            out = sapply(bootEstsParam, function(w){w[x, ]})
            cor(out[1,], out[2,])
        }),
        "cvRhos" = vapply(FUN.VALUE = double(1), seEsts, function(x) x[2,"CovAndCor"]),
        "oracleRho" = oracleRho)
        MSEs = vapply(FUN.VALUE = double(1), seEsts, function(x) x["Bates", "MSEhat"])
        bootEstsR2 = sapply(bootEsts, function(w){1-w[, "MSE"]/w[, "margVar"]})
        bootEstsParamR2 = sapply(bootEstsParam, function(w){1-w[, "MSE"]/w[, "margVar"]})
        sesMat = cbind(apply(rhosMat, 2, function(rho){
            vapply(FUN.VALUE = double(1), seq_along(seEsts), function(j){
                RsquaredSE(MSE = seEsts[[j]]["Bates", "MSEhat"], n = n, margVar = margVar,
                           SEMSE = seEsts[[j]]["Bates", "SE"], rhoBoot = rho[j])[2]
            })
        }),
        "bootSEs" = apply(bootEstsR2, 1, sd),
        "bootParamSEs" = apply(bootEstsParamR2, 1, sd),
        "oracle" = oracleSd
        )
        bootQuantiles = apply(bootEstsR2, 1, quantile, quants);bootQuantilesParam = apply(bootEstsR2, 1, quantile, quants)
        cbind("MSEhat" = MSEs, "margVar" = margVar, sesMat,
              "bootQuantiles" = t(bootQuantiles), "bootQuantilesParam" = t(bootQuantilesParam))
    })
}