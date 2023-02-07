#Fit linear model and make prediction
predLin = function(x, y, xnew, loc = NULL, ...){
    if(!is.null(loc)){
        dat = data.frame("a" = y, "b" = x, loc)
        glsFit = try(silent = TRUE, gls(model = formula("a ~ b"), data = dat, na.action = na.exclude, ...))
    }
    if(is.null(loc) || inherits(glsFit, "try-error")) {
        predOLS(x, y, xnew)
    } else {
        glsFit$coef[1] + xnew * glsFit$coef[2]
    }
}
predGlmnet = function(x,y, xnew, pengls = FALSE,alpha, loc, ...){
    glmNetFit = if(pengls){
        cv.pengls(data = data.frame(x, "y" = y, loc),
                  xNames = colnames(x), outVar = "y", alpha = alpha, ...)
        } else {
            cv.glmnet(x = x, y = y, alpha = alpha, ...)
    }
    predict(glmNetFit, newx =  xnew)
}
predOLS = function(x, y, xnew){
    lmFit = lm.fit(x = cbind(1, x), y = y)
    lmFit$coef[1] + xnew * lmFit$coef[2]
}