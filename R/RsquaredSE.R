#' Calculate R² and find a standard error
#'
#' @param MSE The MSE estimated through cross-validation or bootstrap
#' @param margVar The marginal variance, multiplied by (n+1)/n
#' @param SEMSE The standard error on the MSE estimate
#' @param n the sample size
#' @param rhoBoot,covBoot The correlation or covariance between MSE and MST estimates
RsquaredSE = function(MSE, margVar, SEMSE, n, rhoBoot = NULL, covBoot = NULL){
    if(!(is.null(rhoBoot) || is.null(covBoot))){
        stop("Supply either bootstrap correlation or covariance")
    }
    if(is.na(rhoBoot) || is.na(SEMSE)){
        return(c("R2" = 1-MSE/margVar, "R2SE" = NA))
    }
    Grad = c(-1/margVar, MSE/margVar^2)
    SEmargVar = sqrt(2/(n-1))*margVar
    covSSEmarg = if(!is.null(covBoot)) covBoot else rhoBoot*SEMSE*SEmargVar
    covMat = matrix(c(SEMSE^2, covSSEmarg, covSSEmarg, SEmargVar^2), 2, 2)
    if(!isPD(covMat)){
        covMat = nearPD(covMat)$mat
    }
    c("R2" = 1-MSE/margVar, "R2SE" = as.vector(sqrt(Grad %*% covMat %*% Grad)))
}
isPD = function(mat, tol = 1e-6){
    ev = eigen(mat, symmetric = TRUE)$values
    all(ev >= -tol * abs(ev[1L]))
}