#' Composite kernel association test for SNP-set analysis in pharmacogenetics (PGx) studies.
#' @param G - genotype matrix.
#' @param Tr - treatment vector, 0 indicates placebo, 1 indicates treatment.
#' @param X - non-genetic covariates data matrix.
#' @param y - response vector. Currently continuous and binary responses are supported. Survival response will be added soon.
#' @param trait - response indicator. trait = "continuous" or "binary".
#' @param ker - kernel. ker = "linear", "IBS", "Inter" (interaction kernel) and "RBF" (radial basis function kernel).
#' @param grids - grids of the candidate weights.
#' @param n_a - the number of intervals for manual integration (when integrate function fails). Default n_a = 1000.
#' @param method - method for getting density of A (see details in the reference). Default method is Liu's method.
#' @param subdiv - parameter of Davies' method. Default value is 1E6.
#' @return pvals - p-values of each individual association test.
#' @return finalp - final p-value of the CKAT test.
#' @examples
#' nsamples = 500; nsnps = 10
#' X = rnorm(nsamples,0,1)
#' Tr = sample(0:1,nsamples,replace=TRUE)
#' G = matrix(rbinom(nsamples*nsnps, 1, 0.05), nrow = nsamples, ncol = nsnps)
#' GxT = G*Tr
#' Y0 = 0.5*X + Tr + rnorm(nsamples)
#' CKAT(G, Tr, X, Y0, grids=c(0,0.5,1))
#' @export
#' @importFrom stats dbeta lm glm integrate qchisq dchisq pchisq resid sigma
#' @importFrom CompQuadForm davies liu farebrother

CKAT<- function(G, Tr, X,  y, trait="continuous", ker="linear",
                        grids=c(0,0.5,1),n_a=1000,method="liu", subdiv=10^6){
  #order original data by treatment
  dat = data.frame(y=y,X=X,Tr=Tr,G=G)
  dat = dat[order(-dat$Tr),]
  y = dat[, grepl("^y", names(dat))]
  X = dat[, grepl("^X", names(dat))]
  Tr = dat[, grepl("^Tr", names(dat))]
  G = as.matrix(dat[, grepl("^G", names(dat))])

  n1 = sum(Tr)
  n = length(Tr)
  maf <-  apply(G,2,mean, na.rm=T)/2
  weights = dbeta(maf, 1,25)

  if(ker=="IBS"){
    K_G = IBSKern(G, weight = weights^2);
    K_GxTr = matrix(0,n,n); K_GxTr[1:n1,1:n1] = K_G[1:n1,1:n1]
  }else if(ker=="linear"){
    G.w <- t(t(G) * (weights))
    K_G = G.w %*% t(G.w)
    K_GxTr = matrix(0,n,n); K_GxTr[1:n1,1:n1] = K_G[1:n1,1:n1]
  }else if(ker=="Inter"){
    G.w <- t(t(G) * (weights))
    K_G = InterKern(G.w)
    K_GxTr = matrix(0,n,n); K_GxTr[1:n1,1:n1] = K_G[1:n1,1:n1]
  }
  # else if(ker=="RBF"){
  #   G.w <- t(t(G) * (weights))
  #   K_G = kernelMatrix(rbfdot(sigma=rbfpara), G.w)
  #   K_GxTr = matrix(0,n,n); K_GxTr[1:n1,1:n1] = K_G[1:n1,1:n1]
  # }

  K_G = K_G/n
  K_GxTr = K_GxTr/n
  Ks = list()
  wgrid = 0
  for (w in grids){
    wgrid = wgrid + 1
    K = (1-w)*K_G +w*K_GxTr
    Ks[[wgrid]] = K
  }

  nullCovar = as.matrix(data.frame(Tr,X))
  results = Get_CKAT_P(y=y,Ks=Ks,grids=grids,X=nullCovar,ytype=trait,n_a = n_a, method=method, subdiv=subdiv)
  return(list(pvals=results$pvals, finalp=results$finalp))
}
