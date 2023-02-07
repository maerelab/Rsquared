buildBCa = function(obsPars, bootPars, jackPars, sigLevel){
    z0 = qnorm(mean(obsPars > bootPars, na.rm = TRUE))
    thetaDot = mean(jackPars, na.rm = TRUE)
    aHat = sum(na.rm = TRUE, (thetaDot-jackPars)^3)/(6*sum(na.rm = TRUE, (thetaDot-jackPars)^2)^(3/2))
    if(is.infinite(z0) || is.infinite(aHat))
        return(c(NA, NA))
    alpha1 = pnorm(z0+(z0+qnorm(sigLevel/2))/(1-aHat*(z0+qnorm(sigLevel/2))))
    alpha2 = pnorm(z0+(z0+qnorm(1-sigLevel/2))/(1-aHat*(z0+qnorm(1-sigLevel/2))))
    boot:::norm.inter(bootPars, c(alpha1, alpha2))[, 2]
}
