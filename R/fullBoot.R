#Estimate MSE purely through bootstrap
fullBoot = function(dat, bootReps, bootMethod, paramBoot = FALSE, bootRepsOuter = NULL, jackknife = FALSE){
   id0 = seq_len(length(dat$y))
   obsBoot = makeBootList(dat, id0, bootMethod, Fit = NULL, bootReps = bootReps)
   Fit = getMeanSigma(dat)[[1]]
   bootBoot = if(!is.null(bootRepsOuter)) lapply(integer(bootRepsOuter), function(i){
       datBoot = bootstrap(dat$x, dat$y, id = seq_along(dat$y), paramBoot = FALSE,
                       Mean = Fit$Mean, Sigma = Fit$Sigma)
       makeBootList(datBoot, id0, bootMethod, bootReps = bootReps, Fit)
   }) else NULL
   bootBootParam = if(!is.null(bootRepsOuter)) lapply(integer(bootRepsOuter), function(i){
       datBoot = bootstrap(dat$x, dat$y, id = seq_along(dat$y), paramBoot = TRUE,
                           Mean = Fit$Mean, Sigma = Fit$Sigma)
       makeBootList(datBoot, id0, bootMethod,bootReps = bootReps, Fit)
   }) else NULL
   jackBoot = if(jackknife) lapply(id0, function(i){
       dat$y = dat$y[-i];dat$x = dat$x[-i,,drop = FALSE]
       makeBootList(dat, id0[-length(id0)], bootMethod, bootReps = bootReps)
   }) else NULL
   list("obsBoot" = obsBoot, "bootBoot" = bootBoot,
        "bootBootParam" = bootBootParam, "jackBoot" = jackBoot)
}
simpleBoot = function(dat, id0, bootMethod, meanFit, bootReps){
    lapply(integer(bootReps), function(ii){
        id = sample(id0, replace = TRUE)
        switch(bootMethod,
               "632" = boot632(dat, id, id0, meanFit),
               "oob" = bootOob(dat, id, id0, meanFit))
    })
}
makeBootList = function(datBoot, id0, bootMethod, bootReps, Fit){
        simpleBoot(datBoot, id0, bootMethod, Fit, bootReps = bootReps)
}
fullBootGlmnet = function(dat, bootReps, bootMethod, bootRepsOuter = NULL, jackknife = FALSE, Fit, alpha){
    id0 = seq_len(length(dat$y))
    obsBoot = makeBootListGlmnet(dat, id0, bootMethod, bootReps, alpha = alpha)
    bootBoot = if(!is.null(bootRepsOuter)) lapply(integer(bootRepsOuter), function(i){
        datBoot = bootstrap(dat$x, dat$y, id = seq_along(dat$y), paramBoot = FALSE)
        makeBootListGlmnet(datBoot, id0, bootMethod, bootReps, alpha = alpha)
    }) else NULL
    bootBootParam =  if(!is.null(bootRepsOuter)) lapply(integer(bootRepsOuter), function(i){
        datBoot = bootstrap(dat$x, dat$y, id = seq_along(dat$y), paramBoot = TRUE, Mean = as.vector(Fit$Mean), Sigma = Fit$Sigma)
        makeBootListGlmnet(datBoot, id0, bootMethod, bootReps, alpha = alpha)
    }) else NULL
    jackBoot = if(jackknife) lapply(id0, function(i){
        dat$y = dat$y[-i];dat$x = dat$x[-i, , drop = FALSE]
        list("bootRes" = simpleBootGlmnet(dat, id0[-length(id0)], bootMethod, alpha = alpha),
             "margVar" = var(dat$y))
    }) else NULL
    list("obsBoot" = obsBoot, "bootBoot" = bootBoot, "bootBootParam" = bootBootParam, "jackBoot" = jackBoot)
}
simpleBootGlmnet = function(dat, id0, bootMethod, alpha){
    id = sample(id0, replace = TRUE)
    switch(bootMethod,
           "632" = boot632Glmnet(dat, id, alpha = alpha),
           "oob" = bootOobGlmnet(dat, id, id0, alpha = alpha))
}
makeBootListGlmnet = function(datBoot, id0, bootMethod, bootReps, alpha){
    bootRes = lapply(integer(bootReps), function(i){
        simpleBootGlmnet(datBoot, id0, bootMethod, alpha = alpha)
    })
    margVar = var(datBoot$y)
    list("bootRes" = bootRes, "margVar" = margVar)
}
