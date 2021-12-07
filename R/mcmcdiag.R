## test convergence ###
#Compute the Gelman-Rubin and Geweke statistics
mcmcdiag <- function(x, nparams){
  #x is an mcmc object
  #nparams is number of params in the object
  Gr <- vector(length = nparams)
  GK <- vector(length = nparams)
  k <- 1
  a <- character(0)
  while (k <= nparams){
    m <- x[1:length(x)][, k]
    gr <- gelman.diag(m)
    print(paste("Feature", k, sep = " "))
    print("Gelman-Rubin")
    print(gr)
    if (gr[[1]][1] <= 2) {
      Gr[k] <- "passed" } else{
        Gr[k] <- "failed" }
    gk <- geweke.diag(m,frac1 = 0.1,frac2 = 0.5)
    suspectGK <- names(which(2 * pnorm(-abs(gk[[1]]$z)) < 0.08))
    if (identical(a, suspectGK)) { GK[k] <- "passed"
    } else if (suspectGK == "var1") {GK[k] <- "failed" }
    k <- k + 1}
  return(list(Gr, GK))
}
