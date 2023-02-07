#Extract MSE from simulation object
getMSE = function(x){
    sapply(x, function(y) y$MSE)
}
getR2 = function(x){
    sapply(x, function(y){1-y$MSE/y$margVar})
}