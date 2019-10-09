#' @importFrom stats dbeta lm glm integrate qchisq dchisq pchisq
#' @importFrom CompQuadForm davies liu farebrother
#'
Get_CKAT_P <- function(y,  Ks, grids, X = NULL, ytype="continuous", n_a = 1000,method="liu", subdiv=10^4){

  n = length(y)
  if(ytype=="continuous"){
    if (is.null(X)) {
      p = 1
      X1 = matrix(1, nrow=n)
      mod = glm(y~1, family = "gaussian")
    } else {
      p = ncol(X)
      X1 = cbind(1, X)
      mod = glm(y~X, family = "gaussian")
    }

    s2 = sigma(mod)^2
    D0 = diag(n)
    M = resid(mod)/sqrt(s2)
    P0= (D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1))
    S = sapply(Ks, getIndivP, M, P0, s2=1)
    pvals = unlist(S[7,])
    w = min(pvals)
  }else if(ytype=="binary"){
    if (is.null(X)) {
      p = 1
      X1 = matrix(1, nrow=n)
      mod = glm(y~1, family = "binomial")
    } else {
      p = ncol(X)
      X1 = cbind(1, X)
      mod = glm(y~X, family = "binomial")
    }

    mu = mod$fitted.values
    D = diag(mu*(1-mu))
    D.5 = 1/sqrt((mu*(1-mu)))
    res = y - mu
    DX1 = mu * (1-mu) * X1
    gg  = X1 %*% solve(t(X1) %*% (DX1)) %*% t(DX1)
    P0  = D - diag(D)*gg
    S = sapply(Ks, getIndivP_hm, res, mu, D.5, P0)
    pvals = unlist(S[1,])
    w = min(pvals)
  }
  if(w==0){
    finalp = tempint = 0
  }else{
    nlist = length(Ks)
    K_G = Ks[[1]]
    n1 = sum(X[,1])
    K1 = K_G[1:n1,1:n1];
    K2 = K_G[(n1+1):n,(n1+1):n];
    K12 = K_G[1:n1,(n1+1):n];

    eig1val = Get_Lambda(K1)
    eig2val = Get_Lambda(K2)

    eigval_Ks = sapply(Ks, Get_Lambda)
    q_Ks = sapply(eigval_Ks,function(x)Get_Q_quantile(w, x, method="liu"))
    grids[grids==1] = 1-1e-16
    get_c = Vectorize(function(x)min((q_Ks-x)/(1-grids)))

    eig12val = svd(K12)$d
    eig12val = eig12val[eig12val>1e-10]
    mu_B = sum(eig2val)
    v_B = 2*sum(eig2val^2)
    v_C = sum(eig12val^2)


    if(method=="davies"){
      if(w>1e-2){
        F_B = Vectorize(function(x)davies(sqrt(v_B)/(sqrt(v_B+4*v_C))*(get_c(x) - mu_B) + mu_B,eig2val, acc = 1e-4)$Qq)
      }else{
        F_B = Vectorize(function(x)davies(sqrt(v_B)/(sqrt(v_B+4*v_C))*(get_c(x) - mu_B) + mu_B,eig2val, acc = 1e-7)$Qq)
      }
    }else if(method=="farebrother"){
      F_B = Vectorize(function(x)farebrother(sqrt(v_B)/(sqrt(v_B+4*v_C))*(get_c(x) - mu_B) + mu_B,eig2val)$Qq)
    }else if(method=="liu"){
      F_B = Vectorize(function(x)liu(sqrt(v_B)/(sqrt(v_B+4*v_C))*(get_c(x) - mu_B) + mu_B,eig2val))
    }


    upper1 = Get_Q_quantile(1e-10, eig1val)
    lower1 = 0#

    intgr1 = try(integrate(function(x)Get_Q_density(x, eig1val,method="liu")*F_B(x), lower1, upper1,subdivisions = subdiv)$value,silent=TRUE)
    tempint = intgr1
    if(w<1e-2){
      intgr1 = intgr1/(F_B(upper1)-F_B(lower1)) + abs(1-try(integrate(function(x)Get_Q_density(x, eig1val,method="liu"), lower1, upper1,subdivisions = subdiv)$value,silent=TRUE))
    }
    intgr = intgr1
    if(class(intgr)=="try-error"){
      a = seq(lower1,upper1,length.out=n_a)
      delta_a = a[2] - a[1]
      f_A = Get_Q_density(a, eig1val,method="liu")
      dnsty = sum(f_A)*delta_a
      intgr1 = sum(f_A*F_B(a))*delta_a
      intgr = abs(1-dnsty) + intgr1
    }
    finalp = intgr
  }
  return(list(pvals=pvals,pval_original=tempint,finalp=finalp))
}

#Get_Liu_Params_Mod_Lambda function is adapted from R SKAT package https://github.com/cran/SKAT/

Get_Upper <- function(p, lambda){
  param = Get_Liu_Params_Mod_Lambda(lambda)
  muQ <- param$muQ
  muX <- param$muX
  varQ<- param$sigmaQ^2
  varX<- param$sigmaX^2
  delta = 0*param$d
  df<-param$l
  upper = 5*((qchisq(1-1e-10,df=df) - df)/sqrt(2*df) *sqrt(varQ) + muQ)
  lower = 0
  q = mean(c(upper,lower))
  val = dchisq((q-muQ)/sqrt(varQ)*sqrt(varX)+muX, df, delta) * sqrt(varX)/sqrt(varQ)
  flag = 1
  while(abs(val-p)/p>1e-2&flag<100){
    if(val<p){
      upper = q
      q = mean(c(upper,lower))
      val = dchisq((q-muQ)/sqrt(varQ)*sqrt(varX)+muX, df, delta) * sqrt(varX)/sqrt(varQ)
    }else{
      lower = q
      q = mean(c(upper,lower))
      val = dchisq((q-muQ)/sqrt(varQ)*sqrt(varX)+muX, df, delta) * sqrt(varX)/sqrt(varQ)
    }
    flag = flag + 1
  }
  return(q)
}

Get_Q_density <- function(q, lambda, method="liu"){
  if(method=="liu"){
    param = Get_Liu_Params_Mod_Lambda(lambda)
    muQ <- param$muQ
    muX <- param$muX
    varQ<- param$sigmaQ^2
    varX<- param$sigmaX^2
    delta = 0*param$d
    df<-param$l
    dchisq((q-muQ)/sqrt(varQ)*sqrt(varX)+muX, df, delta) * sqrt(varX)/sqrt(varQ)
  }else if(method=="farebrother"){
    sapply(q, function(x) farebrother(x, lambda)$dnsty)
  }else if(method=="mixed"){
    param = Get_Liu_Params_Mod_Lambda(lambda)
    muQ <- param$muQ
    muX <- param$muX
    varQ<- param$sigmaQ^2
    varX<- param$sigmaX^2
    delta = 0*param$d
    df<-param$l
    id1 = which(q>muQ)
    id2 = which(q<=muQ)
    dnsty1 = dchisq((q-muQ)/sqrt(varQ)*sqrt(varX)+muX, df, delta) * sqrt(varX)/sqrt(varQ)
    dnsty2 = sapply(q, function(x) farebrother(x, lambda)$dnsty)
    pmax(dnsty1,dnsty2)
  }
}


Get_Q_quantile<-function(p, lambda,method="davies"){
  # p is the p-value, i.e. right-tail probability
  param = Get_Liu_Params_Mod_Lambda(lambda)
  muQ <-param$muQ
  varQ<-param$sigmaQ^2
  df<-param$l
  q.org<-qchisq(1-p,df=df)
  q.liu = (q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
  if(method=="liu"){
    return(q.liu)
  }else if(method=="davies"){
    lower = 0.2*q.liu
    upper = 5*q.liu
    p.davies = davies(mean(c(lower,upper)), lambda, acc = 1e-8)$Qq
    flag = 1
    while(abs(p-p.davies)/p>1e-8&flag<100){
      if(p.davies>p){
        lower = mean(c(lower,upper))
      }else{
        upper = mean(c(lower,upper))
      }
      p.davies = davies(mean(c(lower,upper)), lambda, acc = 1e-8)$Qq
      flag = flag + 1
    }
    return(mean(c(lower,upper)))
  }
}


#########################
#MiRKAT Helper Functions#
#########################
# Adapt from R MiRKAT package https://github.com/cran/MiRKAT
Get_Var_Elements =function(m4,u1,u2){
  temp1 = u1^2 * u2^2
  a1    = sum(m4 * temp1)
  a2    = sum(u1^2) * sum(u2^2) - sum(temp1)
  a3    = sum(u1*u2)^2 - sum(temp1)
  return(a1+a2+2*a3)
}

getIndivP_hm = function(K, res, mu, D0, P0){
  Q   = t(res)%*% K %*% res
  # K1  = 1/D0 * P01 %*% K %*% t(1/D0 * P01 )
  K1 = P0 %*% (D0*t(D0*K)) %*% P0
  eK  = eigen(K1, symmetric = T)
  # Instead of matching the first two moments, match to the fourth moment
  # Code adapted from SKAT package
  lambda = eK$values[eK$values > 1e-10]
  U   = as.matrix(eK$vectors[,eK$values > 1e-10])
  p.m = length(lambda)
  m4  = (3*mu^2-3*mu +1)/(mu*(1-mu))

  zeta =rep(0,p.m)
  var_i=rep(0,p.m)
  varQ = 0

  for(i in 1:p.m){   # The diagonals
    temp.M1 = sum(U[,i]^2)^2 - sum(U[,i]^4)
    zeta[i] = sum(m4 * U[,i]^4) + 3* temp.M1 # because ( \sum .)^4, not ^2
    var_i[i]= zeta[i] - 1
  }

  if(p.m == 1){
    Cov_Mat = matrix(zeta* lambda^2, ncol=1,nrow=1)
  } else if(p.m > 1){
    Cov_Mat = diag(zeta* lambda^2)
    for(i in 1:(p.m-1)){
      for(j in (i+1):p.m){
        Cov_Mat[i,j] = Get_Var_Elements(m4,U[,i],U[,j])
        Cov_Mat[i,j] = Cov_Mat[i,j]* lambda[i]* lambda[j]
      }
    }
  }
  Cov_Mat       = Cov_Mat + t(Cov_Mat)
  diag(Cov_Mat) = diag(Cov_Mat)/2

  varQ = sum(Cov_Mat) - sum(lambda)^2
  muQ  = sum(lambda)
  lambda.new = lambda * sqrt(var_i)/sqrt(2)
  df         =  sum(lambda.new^2)^2/sum(lambda.new^4)
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p_corrected= 1 - pchisq(Q_corrected ,df = df)

  p_corrected = ifelse(p_corrected <0, 0, p_corrected)
  return(list(p_hm= p_corrected, Q = Q, muQ = muQ, varQ = varQ, df = df))
}

getIndivP <- function(K, res, P0, s2) {
  ## gets all of the parameters necessary for each kernel
  S = t(res)%*%K%*%res/s2
  #S = t(res)%*%K%*%res
  W = P0%*%K%*%P0
  scale = sum(diag(W))
  ee = eigen(W, symmetric = T)
  w = which(ee$values>1e-10)
  pvals = davies(S,ee$values[w], acc = 1e-4)$Qq
  if(pvals<1e-3){
    pvals = davies(S,ee$values[w], acc = 1e-7)$Qq
  }
  return(list(S=S, W = W,scale = scale, m = length(w), evals = ee$values[w], evecs = ee$vectors[,w], pvals = pvals))
}

#######################
#SKAT Helper Functions#
#######################
# Adapt from R SKAT package https://github.com/cran/SKAT/
Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation
  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get_Lambda <- function(K){

  out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)

  lambda1<-out.s$values
  IDX1<-which(lambda1 >= 1e-10)

  # eigenvalue bigger than sum(eigenvalues)/1000
  IDX2<-which(lambda1 > sum(lambda1[IDX1])/1000)

  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda<-lambda1[IDX2]
  return(lambda)

}

##################
#Kernel Functions#
##################

# IBS Kernel
# Adapt from https://msu.edu/~qlu/doc/function.r
IBSKern=function(geno, weight = 1){
  g=geno
  p=dim(g)[2] #number of SNPs
  n=dim(g)[1]   #number of subjects
  wt = weight

  k_g=matrix(NA,ncol=n,nrow=n)

  gtemp1 = gtemp2 = geno;
  gtemp1[geno==2] = 1;gtemp2[geno==1] = 0;gtemp2[geno==2] = 1;
  gtemp = cbind(gtemp1,gtemp2);

  Inner = gtemp%*%diag(wt,nrow=2*p,ncol=2*p)%*%t(gtemp);
  X2 = matrix(diag(Inner),nrow=n,ncol=n);
  Y2 = matrix(diag(Inner),nrow=n,ncol=n,byrow=T);
  Bound = sum(matrix(wt,nrow=1,ncol=p)*2);
  Dis = (X2+Y2-2*Inner);
  k_g = (Bound-Dis);
  return(k_g)
}

# Interactive Kernal
InterKern=function(x){
  x=as.matrix(x)
  x=t(x)
  Kernel=((1+crossprod(x))^2-crossprod(x^2))/2-0.5
  return(Kernel)
}
