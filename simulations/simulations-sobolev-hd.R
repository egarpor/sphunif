#Simulations for "High-dimensional Sobolev tests on hyperspheres" 

#libraries
library(sphunif)
library(goftest)
library(mvtnorm)
library(PearsonDS)
library(MonteCarlo)
library(mnt)
library(MVT)
library(Directional)
library(gofgamma)
library(MASS)
library(xtable)

#required functions

cv.T<-function(alpha){
  require(PearsonDS)
  w=c(1,-3 + pi^2/3,-pi^2 + 10,-35 + 10/3*pi^2 + 1/45*pi^4) #exact values of the sum of powers of the weights.
  u2 <- 2*w[2]
  kappa_j <- sapply(1:4, function(j)
    2^(j - 1) * factorial(j - 1) * (u2^j + w[j]))
  kum=kappa_j
  mom=c(kum[1:2],kum[3]*kum[2]^(-3/2),3+kum[4]*kum[2]^(-2))
  return(qpearson(1-alpha,moments=mom))
}

BHEP_SH<-function (data, a = 1) 
{
  n = dim(data)[1]
  d = dim(data)[2]
  if (is.null(n)) {
    if (is.vector(data)) {
      n = length(data)
      SUMME1 = 0
      SUMME2 = 0
      for (j in 1:n) {
        SUMME2 = SUMME2 + exp(-a^2 * data[j]^2/(2 * (1 + 
                                                       a^2)))
        for (k in 1:n) {
          SUMME1 = SUMME1 + exp(-a^2 * (data[j] - data[k])^2/2)
        }
      }
      ret = 1/n * SUMME1 - (2/sqrt(1 + a^2)) * SUMME2 + 
        n/sqrt(1 + 2 * a^2)
      return(ret)
    }
    else {
      warning("Wrong dimensions of data!")
    }
  }
  else if (!is.numeric(data)) {
    stop("The data contains non-numeric entries! Check the description for more Information.")
  }
  else if (!a > 0) {
    warning("tuning parameter a>0 needed!")
  } else {
    Djk = data %*% t(data)
    Rquad = diag(Djk)
    Dj = matrix(Rquad, n, n)
    Y = exp((-a^2/2) * (Dj - 2 * Djk + t(Dj)))
    Y2 = exp((-a^2/(2 * (1 + a^2))) * Rquad)
    ret = sum(Y)/n - 2 * ((1 + a^2)^(-d/2)) * sum(Y2) + ((1 + 2 * a^2)^(-d/2)) * n
    return(ret)
  }
}



#-------------------------------------------------------------------------------------
#Simulation HD Sobolev - Normal -> Table 1
#-------------------------------------------------------------------------------------



# Hybrid statistic
stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3), nu = 1,
                     type = c("norm", "t")[1]) {
  
  ## Radiii and projections
  
  # Squared radii
  radii_2 <- rowSums(X^2)
  
  # Projections
  projs <- X / sqrt(radii_2)
  
  ## Projections statistic
  
  # Statistic for the projections
  dim(projs) <- c(dim(projs), 1)
  projs_stat <- unif_stat(data = projs, type = "Sobolev",
                          Sobolev_vk2 = Sobolev_vk2)$Sobolev
  
  # Center to have Tn as in the paper
  p <- ncol(X)
  dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
  mean_projs_stat <- sum(Sobolev_vk2 * dpk)
  projs_stat <- projs_stat - mean_projs_stat
  
  # Scale to have Tn / sigman
  sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
  projs_stat <- projs_stat / sd_projs_stat
  
  # Square statistic to have a limiting chi-square
  projs_stat <- projs_stat^2
  
  ## Radii statistic
  
  if (type == "norm") {
    
    # Standardize squared radii
    radii_2 <- (radii_2 - p) / sqrt(2 * p)
    
    # Radii statistic
    radii_stat <- goftest::ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
    
  } else if (type == "t") {
    
    stop("Preliminar version, does not work yet (does not respect significance
         --- issues with the radii_2 distribution).")
    
    # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
    radii_2 <- radii_2 / p
    
    # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
    # radii_2 <- nu * radii_2
    
    # Standardize squared radii
    # radii_2 <- (radii_2 - p) / sqrt(2 * p)
    
    # Radii statistic
    radii_stat <- goftest::ad.test(x = radii_2, null = "pf",
                                   df1 = nu, df2 = p)$p.value
    
  } else {
    
    stop("Invalid type. Choose 'norm' or 't'.")
    
  }
  
  ## Sum of statistics
  
  #return(radii_stat)
  return(u2 * projs_stat + radii_stat)
  
}



#Simulations of Quantiles (parallel computation)

H0.quan<-function(samplesize=100, dimension=100)
{
  require(MASS)
  require(sphunif)
  require(goftest)
  
  BHEP_SH<-function (data, a = 1) 
  {
    n = dim(data)[1]
    d = dim(data)[2]
    if (is.null(n)) {
      if (is.vector(data)) {
        n = length(data)
        
        SUMME1 = 0
        SUMME2 = 0
        for (j in 1:n) {
          SUMME2 = SUMME2 + exp(-a^2 * data[j]^2/(2 * (1 + a^2)))
          for (k in 1:n) {
            SUMME1 = SUMME1 + exp(-a^2 * (data[j] - data[k])^2/2)
          }
        }
        ret = 1/n * SUMME1 - (2/sqrt(1 + a^2)) * SUMME2 + 
          n/sqrt(1 + 2 * a^2)
        return(ret)
      }
      else {
        warning("Wrong dimensions of data!")
      }
    } else if (!is.numeric(data)) {
      stop("The data contains non-numeric entries! Check the description for more Information.")
    }
    else if (!a > 0) {
      warning("tuning parameter a>0 needed!")
    } else {
      #data = standard(data)
      Djk = data %*% t(data)
      Rquad = diag(Djk)
      Dj = matrix(Rquad, n, n)
      Y = exp((-a^2/2) * (Dj - 2 * Djk + t(Dj)))
      Y2 = exp((-a^2/(2 * (1 + a^2))) * Rquad)
      ret = sum(Y)/n - 2 * ((1 + a^2)^(-d/2)) * sum(Y2) + ((1 + 
                                                              2 * a^2)^(-d/2)) * n
      return(ret)
    }
  }
  
  RES=array(0,dim=c(length(dimension),length(samplesize)))
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      X=mvrnorm(n,rep(0,d),diag(1,d))
      j.n=j.n+1
      RES[j.d,j.n]=BHEP_SH(X, a = 1/sqrt(d))
    }
  }
  return(list("RES"=RES))
}

set.seed(0815)
param_list=list("samplesize"=c(100,200), "dimension"=c(100,200,300))
time_START<-Sys.time()
s.hyb.quan<-MonteCarlo(func=H0.quan, nrep=10000, param_list=param_list, ncpus=14)
summary(s.hyb.quan)
time_END<-Sys.time()
difftime(time_END, time_START)

par(mfrow=c(3,3))
q.s.hyb=matrix(0,2,3)
j.n=0
for (n in c("samplesize=100","samplesize=200"))
{ j.n=j.n+1
j.d=0
for (d in c("dimension=100","dimension=200","dimension=300"))
{
  j.d=j.d+1
  q.s.hyb[j.n,j.d]=sort(as.numeric(s.hyb.quan$results$RES[n,d,]))[0.95*10000]
  #plot.ecdf(s.hyb.quan$results$RES2[n,d,])
}
}

# Empirical powers of the test (parallel computation)
H1.s.hyb<-function(Choice=c(1,2),samplesize=c(100,200), dimension=c(100,200))
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(Directional)
  require(rotasym)
  require(mvtnorm)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3)) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    
    # Projections
    projs <- X / sqrt(radii_2)
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(x,df=p))$p.value
    
    ## Sum of statistics
    
    #return(radii_stat)
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
    
  }
  
  BHEP_SH<-function (data, a = 1) 
  {
    n = dim(data)[1]
    d = dim(data)[2]
    if (is.null(n)) {
      if (is.vector(data)) {
        n = length(data)
        
        SUMME1 = 0
        SUMME2 = 0
        for (j in 1:n) {
          SUMME2 = SUMME2 + exp(-a^2 * data[j]^2/(2 * (1 + a^2)))
          for (k in 1:n) {
            SUMME1 = SUMME1 + exp(-a^2 * (data[j] - data[k])^2/2)
          }
        }
        ret = 1/n * SUMME1 - (2/sqrt(1 + a^2)) * SUMME2 + 
          n/sqrt(1 + 2 * a^2)
        return(ret)
      }
      else {
        warning("Wrong dimensions of data!")
      }
    } else if (!is.numeric(data)) {
      stop("The data contains non-numeric entries! Check the description for more Information.")
    }
    else if (!a > 0) {
      warning("tuning parameter a>0 needed!")
    } else {
      #data = standard(data)
      Djk = data %*% t(data)
      Rquad = diag(Djk)
      Dj = matrix(Rquad, n, n)
      Y = exp((-a^2/2) * (Dj - 2 * Djk + t(Dj)))
      Y2 = exp((-a^2/(2 * (1 + a^2))) * Rquad)
      ret = sum(Y)/n - 2 * ((1 + a^2)^(-d/2)) * sum(Y2) + ((1 + 
                                                              2 * a^2)^(-d/2)) * n
      return(ret)
    }
  }
  
  #distributions
  # Simulate non-normal data due to dependency between the radial part and the
  # projections, both marginally correct
  r_non_normal <- function(n, p, mu = c(1, rep(0, p - 1)), rho = 0) {
    
    # Simulate dependency between the radial part and the projections using a
    # Gaussian copula
    stopifnot(p >= 2)
    stopifnot(abs(rho) <= 1)
    rhos <- rho^c(1:p)
    sigma <- rbind(c(1, rhos), cbind(rhos, diag(1, nrow = p, ncol = p)))
    x <- mvrnorm(n = n, mu = rep(0, p + 1), Sigma = sigma)
    
    # Simulate the radial part
    radii2 <- qchisq(pnorm(x[, 1]), df = p)
    
    # Simulate the projections
    projs <- x[, -1, drop = FALSE]
    projs <- projs / sqrt(rowSums(projs^2))
    
    # Return normal-like observations
    return(sqrt(radii2) * projs)
    
  }
  
  # Simulate non-spherically symmetric distributions such that the deviation is
  # on the non-uniform distributions of the projections (type = "projs_vMF") or
  # the deviation is on the conditional radii (type = "radii_cond") given the
  # projections
  r_non_sph_sym <- function(n, p, mu = c(1, rep(0, p - 1)), kappa = 0,
                            type = c("projs_vMF", "radii_cond")) {
    
    stopifnot(p >= 2)
    if (type == "projs_vMF") {
      
      projs <- r_vMF(n = n, mu = mu, kappa = kappa)
      radii <- sqrt(rchisq(n = n, df = p))
      
    } else if (type == "radii_cond") {
      
      projs <- r_unif_sphere(n = n, p = p)
      radii <- rnorm(n, mean = 1 + kappa * (1 + projs %*% mu), sd = 1)
      
    }
    return(radii * projs)
    
  }
  
  #choose sample
  
  RES=array(0,dim=c(length(dimension),length(samplesize)))
  RES2=array(0,dim=c(length(dimension),length(samplesize)))
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      sample0=mvrnorm(n,rep(0,d),diag(1,d))
      
      proba=0.1
      hilf=c(rep(0.5,floor(proba*d)),rep(1.5,d-floor(proba*d)))
      sample1=mvrnorm(n,rep(0,d),(d*hilf/sum(hilf))*diag(1,d))
      proba=0.25
      hilf=c(rep(0.5,floor(proba*d)),rep(1.5,d-floor(proba*d)))
      sample2=mvrnorm(n,rep(0,d),(d*hilf/sum(hilf))*diag(1,d))  
      proba=0.5
      hilf=c(rep(0.5,floor(proba*d)),rep(1.5,d-floor(proba*d)))
      sample3=mvrnorm(n,rep(0,d),(d*hilf/sum(hilf))*diag(1,d))
      
      sample4=r_non_sph_sym(n = n, p = d, kappa = 1, type = "projs_vMF")
      sample5=r_non_sph_sym(n = n, p = d, kappa = 4, type = "projs_vMF")
      
      sample6=rmvt(n,sigma=0.75*diag(1,d),df=10)
      sample7=rmvt(n,sigma=0.8*diag(1,d),df=10)
      sample8=rmvt(n,sigma=0.9*diag(1,d),df=10)
      sample9=rmvt(n,sigma=1.1*diag(1,d),df=10)
      sample10=rmvt(n,sigma=diag(1,d),df=10)
      sample11=rmvt(n,sigma=diag(1,d),df=30)
      sample12=rmvt(n,sigma=diag(1,d),df=100)
      sample13=rmvt(n,sigma=diag(1,d),df=500)
      sample14=rmvt(n,sigma=diag(1,d),df=1000)
      
      
      sample=switch(Choice,sample0,sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10,sample11,sample12,sample13,sample14)
      j.n=j.n+1
      RES[j.d,j.n]=stat_hyb(sample)$p.value<0.05
      RES2[j.d,j.n]=BHEP_SH(sample,a=1/sqrt(d))>q.s.hyb[j.n,j.d]
    }
  }
  return(list("RES"=RES,"RES2"=RES2))
}

set.seed(0815)
param_list=list("Choice"=1:15,"samplesize"=c(100,200), "dimension"=c(100,200,300))
time_START<-Sys.time()
s.hyb.norm<-MonteCarlo(func=H1.s.hyb, nrep=5000, param_list=param_list, ncpus=14)
summary(s.hyb.norm)
time_END<-Sys.time()
difftime(time_END, time_START)
save(s.hyb.norm, file = "norm.RData")

#Arranging the results in a table
Erg=matrix(0,nrow=15,ncol=12)
j.n=0
for (n in c("samplesize=100","samplesize=200"))
{ j.n=j.n+1
RESULT=matrix(0,nrow=15,ncol=6)
for (a in 1:15)
{
  RESULT[a,1]=mean(as.numeric(s.hyb.norm$results$RES[a,j.n,1,]))
  RESULT[a,2]=mean(as.numeric(s.hyb.norm$results$RES[a,j.n,2,]))
  RESULT[a,3]=mean(as.numeric(s.hyb.norm$results$RES[a,j.n,3,]))
  RESULT[a,4]=mean(as.numeric(s.hyb.norm$results$RES2[a,j.n,1,]))
  RESULT[a,5]=mean(as.numeric(s.hyb.norm$results$RES2[a,j.n,2,]))
  RESULT[a,6]=mean(as.numeric(s.hyb.norm$results$RES2[a,j.n,3,]))
}
if(j.n==1) {Erg=RESULT} else {Erg=cbind(Erg,RESULT)}
}

xtable(t(t(round(Erg,2)*100)),digits=0,include.rownames = FALSE)



#-------------------------------------------------------------------------------------
#Simulation HD Sobolev - t simple hypothesis -> Table 2 upper part
#-------------------------------------------------------------------------------------



# Hybrid statistic
stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),
                     type = c("norm", "t")[2],nu) {
  
  ## Radiii and projections
  
  # Squared radii
  radii_2 <- rowSums(X^2)
  
  # Projections
  projs <- X / sqrt(radii_2)
  
  ## Projections statistic
  
  # Statistic for the projections
  dim(projs) <- c(dim(projs), 1)
  projs_stat <- unif_stat(data = projs, type = "Sobolev",
                          Sobolev_vk2 = Sobolev_vk2)$Sobolev
  
  # Center to have Tn as in the paper
  p <- ncol(X)
  dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
  mean_projs_stat <- sum(Sobolev_vk2 * dpk)
  projs_stat <- projs_stat - mean_projs_stat
  
  # Scale to have Tn / sigman
  sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
  projs_stat <- projs_stat / sd_projs_stat
  
  # Square statistic to have a limiting chi-square
  projs_stat <- 1-pchisq(projs_stat^2,df=1)
  
  ## Radii statistic
  
  if (type == "norm") {
    
    # Standardize squared radii
    radii_2 <- (radii_2 - p) / sqrt(2 * p)
    
    # Radii statistic
    radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
    
  } else if (type == "t") {
    
    #stop("Preliminar version, does not work yet (does not respect significance
    #   --- issues with the radii_2 distribution).")
    
    # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
    radii_2 <- radii_2 / p
    
    # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
    # radii_2 <- nu * radii_2
    
    # Standardize squared radii
    # radii_2 <- (radii_2 - p) / sqrt(2 * p)
    
    # Radii statistic
    
    radii_stat <- ad.test(x = radii_2, null = function(x) pf(x,df1=p,df2=nu))$p.value
    
  } else {
    
    stop("Invalid type. Choose 'norm' or 't'.")
    
  }
  
  ## Fisher's method for combination of independent statistics
  
  #return(radii_stat)
  return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
  
}

#check limit distribution chisq with 4 df

H0.quan<-function(samplesize=c(20,50,100), dimension=c(2,3,5))
{
  require(MASS)
  require(sphunif)
  require(goftest)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),
                       type = c("norm", "t")[2],nu) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    
    # Projections
    projs <- X / sqrt(radii_2)
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    
    if (type == "norm") {
      
      # Standardize squared radii
      radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
      
    } else if (type == "t") {
      
      #stop("Preliminar version, does not work yet (does not respect significance
      #   --- issues with the radii_2 distribution).")
      
      # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
      radii_2 <- radii_2 / p
      
      # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
      # radii_2 <- nu * radii_2
      
      # Standardize squared radii
      # radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      
      radii_stat <- ad.test(x = radii_2, null = function(x) pf(x,df1=p,df2=nu))$p.value
      
    } else {
      
      stop("Invalid type. Choose 'norm' or 't'.")
      
    }
    
    ## Fisher's method for combination of independent statistics
    
    #return(radii_stat)
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
    
  }
  
  true.nu=5
  
  RES=array(0,dim=c(length(dimension),length(samplesize)))
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      #X=mvrnorm(n,rep(0,d),diag(1,d))
      X=rmvt(n,sigma=diag(1,d),df=true.nu)
      j.n=j.n+1
      RES[j.d,j.n]=stat_hyb(X, type = "t",nu=true.nu)$statistic
      #RES[j.d,j.n]=BHEP_SH(X, a = 1)
    }
  }
  return(list("RES"=RES))
}

set.seed(0815)
param_list=list("samplesize"=c(100,200,500,1000), "dimension"=c(100,200,300))
time_START<-Sys.time()
s.hyb.quan<-MonteCarlo(func=H0.quan, nrep=5000, param_list=param_list, ncpus=14)
summary(s.hyb.quan)
time_END<-Sys.time()
difftime(time_END, time_START)

#Empirical powers of the test (parallel computation)

H1.s.hyb<-function(Choice=c(1,2),samplesize=c(100,200), dimension=c(100,200))
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(Directional)
  require(rotasym)
  require(mvtnorm)
  require(sn)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),
                       type = c("norm", "t")[2],nu) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    
    # Projections
    projs <- X / sqrt(radii_2)
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    
    if (type == "norm") {
      
      # Standardize squared radii
      radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
      
    } else if (type == "t") {
      
      #stop("Preliminar version, does not work yet (does not respect significance
      #   --- issues with the radii_2 distribution).")
      
      # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
      radii_2 <- radii_2 / p
      
      # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
      # radii_2 <- nu * radii_2
      
      # Standardize squared radii
      # radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      
      radii_stat <- ad.test(x = radii_2, null = function(x) pf(x,df1=p,df2=nu))$p.value
      
    } else {
      
      stop("Invalid type. Choose 'norm' or 't'.")
      
    }
    
    ## Fisher's method for combination of independent statistics
    
    #return(radii_stat)
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
    
  }
  
  #distributions
  # Simulate non-normal data due to dependency between the radial part and the
  # projections, both marginally correct
  r_non_normal <- function(n, p, mu = c(1, rep(0, p - 1)), rho = 0) {
    
    # Simulate dependency between the radial part and the projections using a
    # Gaussian copula
    stopifnot(p >= 2)
    stopifnot(abs(rho) <= 1)
    rhos <- rho^c(1:p)
    sigma <- rbind(c(1, rhos), cbind(rhos, diag(1, nrow = p, ncol = p)))
    x <- mvrnorm(n = n, mu = rep(0, p + 1), Sigma = sigma)
    
    # Simulate the radial part
    radii2 <- qchisq(pnorm(x[, 1]), df = p)
    
    # Simulate the projections
    projs <- x[, -1, drop = FALSE]
    projs <- projs / sqrt(rowSums(projs^2))
    
    # Return normal-like observations
    return(sqrt(radii2) * projs)
    
  }
  
  # Simulate non-spherically symmetric distributions such that the deviation is
  # on the non-uniform distributions of the projections (type = "projs_vMF") or
  # the deviation is on the conditional radii (type = "radii_cond") given the
  # projections
  r_non_sph_sym <- function(n, p, mu = c(1, rep(0, p - 1)), kappa = 0,
                            type = c("projs_vMF", "radii_cond")) {
    
    stopifnot(p >= 2)
    if (type == "projs_vMF") {
      
      projs <- r_vMF(n = n, mu = mu, kappa = kappa)
      radii <- sqrt(rchisq(n = n, df = p))
      
    } else if (type == "radii_cond") {
      
      projs <- r_unif_sphere(n = n, p = p)
      radii <- rnorm(n, mean = 1 + kappa * (1 + projs %*% mu), sd = 1)
      
    }
    return(radii * projs)
    
  }
  
  #choose sample
  
  RES=array(0,dim=c(length(dimension),length(samplesize)))
  
  true.nu=5
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
       sample0=rmvt(n,sigma=diag(1,d),df=true.nu)
       sample1=rmvt(n,sigma=diag(1,d),df=true.nu+1)
       sample2=rmvt(n,sigma=diag(1,d),df=true.nu+2)
       sample3=rmvt(n,sigma=diag(1,d),df=true.nu+3)
       sample4=rmvt(n,sigma=diag(1,d),df=true.nu+4)
       sample5=rmvt(n,sigma=0.9*diag(1,d),df=true.nu)
       sample6=rvmf(n,mu=c(1,rep(0,d-1)),k=5)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       sample7=rvmf(n,mu=c(1,rep(0,d-1)),k=10)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       cov_matrix <- diag(0.75,d)+matrix(0.25,d,d)
       mean_vec <- rep(0,d)      # Mean vector
       alpha0 <- rep(0,d)
       alpha1 <- c(1,rep(0,d-1))
       sample8=rmst(n, xi = mean_vec, Omega = cov_matrix, alpha = alpha0, nu=true.nu)
       sample9=rmst(n, xi = mean_vec, Omega = diag(1,d), alpha = alpha1, nu=true.nu)
       sample10=rmst(n, xi = mean_vec, Omega = cov_matrix, alpha = alpha0, nu=true.nu+3)
       sample11=rmst(n, xi = mean_vec, Omega = diag(1,d), alpha = alpha1, nu=true.nu+3)
 
      sample=switch(Choice,sample0,sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10,sample11)
      j.n=j.n+1
      RES[j.d,j.n]=stat_hyb(sample, type = "t",nu=true.nu)$p.value<0.05
    }
  }
  return(list("RES"=RES))
}

set.seed(0815)
param_list=list("Choice"=1:12,"samplesize"=c(100,200,500,1000), "dimension"=c(100,200,300))
time_START<-Sys.time()
s.hyb.alt.t<-MonteCarlo(func=H1.s.hyb, nrep=5000, param_list=param_list, ncpus=50)
summary(s.hyb.alt.t)
time_END<-Sys.time()
difftime(time_END, time_START)
save(s.hyb.alt.t, file = "t_SH.RData")

#Arranging the results in a table

j.n=0
for (n in c("samplesize=100","samplesize=200","samplesize=500","samplesize=1000"))
{ j.n=j.n+1
RESULT=matrix(0,nrow=12,ncol=3)
for (a in 1:12)
{
  RESULT[a,1]=mean(as.numeric(s.hyb.alt.t$results$RES[a,j.n,1,]))
  RESULT[a,2]=mean(as.numeric(s.hyb.alt.t$results$RES[a,j.n,2,]))
  RESULT[a,3]=mean(as.numeric(s.hyb.alt.t$results$RES[a,j.n,3,]))
}
if(j.n==1) {Erg=RESULT} else {Erg=cbind(Erg,RESULT)}
}
xtable(t(t(round(Erg,2)*100)),digits=0,include.rownames = FALSE)



#-------------------------------------------------------------------------------------
#Simulation HD Sobolev - t composite hypothesis -> Table 2 lower part
#-------------------------------------------------------------------------------------



# Hybrid stat according to Fisher
stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),boot=500) {
  
  ## Radiii and projections
  
  # Squared radii
  radii_2 <- rowSums(X^2)
  radii<-sqrt(radii_2)
  
  # Projections
  projs <- X / radii
  
  ## Projections statistic
  
  # Statistic for the projections
  dim(projs) <- c(dim(projs), 1)
  projs_stat <- unif_stat(data = projs, type = "Sobolev",
                          Sobolev_vk2 = Sobolev_vk2)$Sobolev
  
  # Center to have Tn as in the paper
  p <- ncol(X)
  dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
  mean_projs_stat <- sum(Sobolev_vk2 * dpk)
  projs_stat <- projs_stat - mean_projs_stat
  
  # Scale to have Tn / sigman
  sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
  projs_stat <- projs_stat / sd_projs_stat
  
  # Square statistic to have a limiting chi-square
  projs_stat <- 1-pchisq(projs_stat^2,df=1)
  
  ## Radii statistic
  
  mle.F.sp<-function(data,df1){
    if (df1<=0) {stop("df1 must be positive")} else {
      diff.log.pdf <- function(x,d_1,d_2){((d_1 * x + d_2) * digamma(d_1 / 2 + d_2 / 2) + (-d_1 * x - d_2) * log(d_1 * x + d_2) + (-d_1 * x - d_2) * digamma(d_2 / 2) + (d_1 * x + d_2) * log(d_2) + d_1 * (x - 1)) / (2 * d_1 * x + 2 * d_2)}
      hilf<-function(nu,data,df1){return(sum(diff.log.pdf(data,df1,nu)))}
      return(uniroot(hilf,df1=df1,data=data,lower=0.01,upper=1000)$root)}
  }
  
  para=mle.F.sp(radii_2/p,df1=p)
  
  radii_stat <- ad.test(radii_2/p,null = function(x) pf(x,df1=p,df2=para),estimated=TRUE)$p.value
  
  #bootstrap version
  #y <- rep(0, boot)
  #for (j in 1:boot)
  #{
  #  BS <- rf(n, df1=p,df2=para)
  #  BS_estimator <- mle.F.sp(BS,df1=p)
  #  y[j] = ad.test(BS,null = function(x) pf(x,df1=p,df2=BS_estimator),estimated=T)
  #}
  #radii_stat=mean((y>radii_stat))
  
  return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
}

#Empirical powers of the test (parallel computation)

H1.s.hyb<-function(Choice=1,samplesize=100, dimension=100)
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(Directional)
  require(rotasym)
  require(gofgamma)
  require(sn)
  require(mvtnorm)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),boot=500) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    radii<-sqrt(radii_2)
    
    # Projections
    projs <- X / radii
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    
    mle.F.sp<-function(data,df1){
      if (df1<=0) {stop("df1 must be positive")} else {
        diff.log.pdf <- function(x,d_1,d_2){((d_1 * x + d_2) * digamma(d_1 / 2 + d_2 / 2) + (-d_1 * x - d_2) * log(d_1 * x + d_2) + (-d_1 * x - d_2) * digamma(d_2 / 2) + (d_1 * x + d_2) * log(d_2) + d_1 * (x - 1)) / (2 * d_1 * x + 2 * d_2)}
        hilf<-function(nu,data,df1){return(sum(diff.log.pdf(data,df1,nu)))}
        return(uniroot(hilf,df1=df1,data=data,lower=0.01,upper=1000)$root)}
    }
    
    para=mle.F.sp(radii_2/p,df1=p)
    
    radii_stat <- as.numeric(ad.test(radii_2/p,null = function(x) pf(x,df1=p,df2=para))$statistic)
    
    #bootstrap version
    y <- rep(0, boot)
    for (j in 1:boot)
    {
      BS <- rf(n, df1=p,df2=para)
      BS_estimator <- mle.F.sp(BS,df1=p)
      y[j] = as.numeric(ad.test(BS,null = function(x) pf(x,df1=p,df2=BS_estimator))$statistic)
    }
    radii_stat=mean((y>radii_stat))
    
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
  }
  
  #distributions
  # Simulate non-normal data due to dependency between the radial part and the
  # projections, both marginally correct
  r_non_normal <- function(n, p, mu = c(1, rep(0, p - 1)), rho = 0) {
    
    # Simulate dependency between the radial part and the projections using a
    # Gaussian copula
    stopifnot(p >= 2)
    stopifnot(abs(rho) <= 1)
    rhos <- rho^c(1:p)
    sigma <- rbind(c(1, rhos), cbind(rhos, diag(1, nrow = p, ncol = p)))
    x <- mvrnorm(n = n, mu = rep(0, p + 1), Sigma = sigma)
    
    # Simulate the radial part
    radii2 <- qchisq(pnorm(x[, 1]), df = p)
    
    # Simulate the projections
    projs <- x[, -1, drop = FALSE]
    projs <- projs / sqrt(rowSums(projs^2))
    
    # Return normal-like observations
    return(sqrt(radii2) * projs)
    
  }
  
  # Simulate non-spherically symmetric distributions such that the deviation is
  # on the non-uniform distributions of the projections (type = "projs_vMF") or
  # the deviation is on the conditional radii (type = "radii_cond") given the
  # projections
  r_non_sph_sym <- function(n, p, mu = c(1, rep(0, p - 1)), kappa = 0,
                            type = c("projs_vMF", "radii_cond")) {
    
    stopifnot(p >= 2)
    if (type == "projs_vMF") {
      
      projs <- r_vMF(n = n, mu = mu, kappa = kappa)
      radii <- sqrt(rchisq(n = n, df = p))
      
    } else if (type == "radii_cond") {
      
      projs <- r_unif_sphere(n = n, p = p)
      radii <- rnorm(n, mean = 1 + kappa * (1 + projs %*% mu), sd = 1)
      
    }
    return(radii * projs)
    
  }
  
  #choose sample
  
  RES=array(0,dim=c(length(dimension),length(samplesize)))
  
  true.nu=5
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      cov_matrix <- diag(0.75,d)+matrix(0.25,d,d)
      mean_vec <- rep(0,d)      # Mean vector
      alpha0 <- rep(0,d)
      alpha1 <- c(1,rep(0,d-1))
      
       sample0=rmvt(n,sigma=diag(1,d),df=true.nu)
       sample1=rmvt(n,sigma=0.8*diag(1,d),df=true.nu)
       sample2=rmvt(n,sigma=0.9*diag(1,d),df=true.nu)
       sample3=rmvt(n,sigma=1.1*diag(1,d),df=true.nu)
       sample4=rmvt(n,sigma=1.25*diag(1,d),df=true.nu)
       sample5=rvmf(n,mu=c(1,rep(0,d-1)),k=1)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       sample6=rvmf(n,mu=c(1,rep(0,d-1)),k=3)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       sample7=rvmf(n,mu=c(1,rep(0,d-1)),k=5)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       sample8=rvmf(n,mu=c(1,rep(0,d-1)),k=10)*sqrt(d*matrix(rf(n,df1=d,df2=true.nu),nrow=n,ncol=d))
       sample9=rmst(n, xi = mean_vec, Omega = diag(1,d), alpha = alpha1, nu=true.nu)
       sample10=rmst(n, xi = mean_vec, Omega = cov_matrix, alpha = alpha0, nu=true.nu)
       
      sample=switch(Choice,sample0,sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10)
      j.n=j.n+1
      RES[j.d,j.n]=stat_hyb(sample)$p.value<0.05
    }
  }
  return(list("RES"=RES))
}

set.seed(0815)
param_list=list("Choice"=1:11,"samplesize"=c(100,200,500,1000), "dimension"=c(100,200,300))
time_START<-Sys.time()
s.hyb.alt.t.2<-MonteCarlo(func=H1.s.hyb, nrep=5000, param_list=param_list, ncpus=60)
summary(s.hyb.alt.t.2)
time_END<-Sys.time()
difftime(time_START, time_END)
save(s.hyb.alt.t.2, file = "t_CH.RData")

#Get results: For every sample size we have a new table
j.n=0
for (n in c("samplesize=100","samplesize=200","samplesize=500","samplesize=1000"))
{ j.n=j.n+1
RESULT=matrix(0,nrow=11,ncol=3)
for (a in 1:11)
{
  RESULT[a,1]=mean(as.numeric(s.hyb.alt.t.2$results$RES[a,j.n,1,]))
  RESULT[a,2]=mean(as.numeric(s.hyb.alt.t.2$results$RES[a,j.n,2,]))
  RESULT[a,3]=mean(as.numeric(s.hyb.alt.t.2$results$RES[a,j.n,3,]))
}
if(j.n==1) {Erg=RESULT} else {Erg=cbind(Erg,RESULT)}
}
xtable(t(t(round(Erg,2)*100)),digits=0,include.rownames = FALSE) #Note that the last two rows have to be switched in Table 2


#-------------------------------------------------------------------------------------
#Simulation HD Sobolev - stable hypothesis -> Table 3 
#-------------------------------------------------------------------------------------

#check limit distribution chisq with 4 df

H0.quan<-function(samplesize=c(20,50,100), dimension=c(2,3,5))
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(stabledist)
  
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),
                       type = c("norm", "t","stable")[3],nu=5,alpha0=1) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    
    # Projections
    projs <- X / sqrt(radii_2)
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    
    if (type == "norm") {
      
      # Standardize squared radii
      radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
      
    } else if (type == "t") {
      
      #stop("Preliminar version, does not work yet (does not respect significance
      #   --- issues with the radii_2 distribution).")
      
      # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
      radii_2 <- radii_2 / p
      
      # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
      # radii_2 <- nu * radii_2
      
      # Standardize squared radii
      # radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      
      radii_stat <- ad.test(x = radii_2, null = function(x) pf(x,df1=p,df2=nu))$p.value
      
    } else if(type == "stable") {
      
      r2dist<-function(y,d,alpha0){
        hilf<-function(x,y,d,alpha){return(stabledist::pstable(y/x,alpha=alpha/2,beta=1,gamma=2*(cos(pi*alpha/4))^(2/alpha),delta=0,pm=1)*dchisq(x,df=d))}
        hilf.erg=rep(0,length(y))
        for (i in 1:length(y))
        {
          hilf.erg[i]=integrate(hilf,y=y[i],d=d,alpha=alpha0,lower=0,upper=Inf)$value
        }
        return(hilf.erg)
      }
      
      radii_stat <- ad.test(x = radii_2, null = function(x) r2dist(x,d=p,alpha0=alpha0))$p.value
      
    } else {
      
      stop("Invalid type. Choose 'norm', 't' or 'stable'.")
      
    }
    
    ## Fisher's method for combination of independent statistics
    
    #return(radii_stat)
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4)))
    
  }
  
  true.alpha0=1.2
  
  Erg=array(0,dim=c(length(dimension),length(samplesize)))
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      X=mvrnorm(n,rep(0,d),diag(1,d))
      #X=rmvt(n,sigma=diag(1,d),df=true.nu)
      A=stabledist::rstable(n,true.alpha0/2,1,2*(cos(pi*true.alpha0/4))^(2/true.alpha0),0,pm=1)
      X=matrix(sqrt(A),ncol=d,nrow=n)*X
      j.n=j.n+1
      Erg[j.d,j.n]=stat_hyb(X, type = "stable",alpha0=true.alpha0)$statistic
      #Erg[j.d,j.n]=BHEP_SH(X, a = 1)
    }
  }
  return(list("Erg"=Erg))
}
library(MonteCarlo)
set.seed(0815)
param_list=list("samplesize"=c(50,100), "dimension"=c(50,100))
ZeitAnfang<-Sys.time()
s.hyb.quan<-MonteCarlo(func=H0.quan, nrep=100, param_list=param_list, ncpus=10)
summary(s.hyb.quan)
ZeitEnde<-Sys.time()
difftime(ZeitEnde, ZeitAnfang)

par(mfrow=c(2,2))
j.n=0
for (n in c("samplesize=50","samplesize=100"))
{ j.n=j.n+1
j.d=0
for (d in c("dimension=50","dimension=100"))
{
  plot.ecdf(as.numeric(s.hyb.quan$results$Erg[n,d,]))
  curve(pchisq(x,df=4),from=0,to=25,add=T,col="red")
  #j.d=j.d+2
}
}

#power study

H1.s.hyb<-function(Choice=c(1,2),samplesize=c(100,200), dimension=c(100,200))
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(Directional)
  require(rotasym)
  require(stabledist)
  require(MVT)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),
                       type = c("norm", "t","stable")[3],nu=5,alpha0=1) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    
    # Projections
    projs <- X / sqrt(radii_2)
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    
    if (type == "norm") {
      
      # Standardize squared radii
      radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      radii_stat <- ad.test(x = radii_2, null = function(x) pchisq(sqrt(2*p)*x+p,df=p))$statistic
      
    } else if (type == "t") {
      
      #stop("Preliminar version, does not work yet (does not respect significance
      #   --- issues with the radii_2 distribution).")
      
      # X'X / p ~ F(nu, p) (https://en.wikipedia.org/wiki/Multivariate_t-distribution#Radial_Distribution)
      radii_2 <- radii_2 / p
      
      # nu * F(nu, p) -> chi^2(p) as p -> Inf (https://en.wikipedia.org/wiki/F-distribution#Properties_and_related_distributions)
      # radii_2 <- nu * radii_2
      
      # Standardize squared radii
      # radii_2 <- (radii_2 - p) / sqrt(2 * p)
      
      # Radii statistic
      
      radii_stat <- ad.test(x = radii_2, null = function(x) pf(x,df1=p,df2=nu))$p.value
      
    } else if(type == "stable") {
      
      r2dist<-function(y,d,alpha0){
        hilf<-function(x,y,d,alpha){return(stabledist::pstable(y/x,alpha=alpha/2,beta=1,gamma=2*(cos(pi*alpha/4))^(2/alpha),delta=0,pm=1)*dchisq(x,df=d))}
        hilf.erg=rep(0,length(y))
        for (i in 1:length(y))
        {
          hilf.erg[i]=integrate(hilf,y=y[i],d=d,alpha=alpha0,lower=0,upper=Inf)$value
        }
        return(hilf.erg)
      }
      
      radii_stat <- ad.test(x = radii_2, null = function(x) r2dist(x,d=p,alpha0=alpha0))$p.value
      
    } else {
      
      stop("Invalid type. Choose 'norm', 't' or 'stable'.")
      
    }
    
    ## Fisher's method for combination of independent statistics
    
    #return(radii_stat)
    return(list(statistic=-2 * (log(projs_stat) + log(radii_stat)),p.value=1-pchisq(-2 * (log(projs_stat) + log(radii_stat)),df=4),proj=projs_stat,rad=radii_stat))
    
  }
  
  #distributions
  # Simulate non-normal data due to dependency between the radial part and the
  # projections, both marginally correct
  r_non_normal <- function(n, p, mu = c(1, rep(0, p - 1)), rho = 0) {
    
    # Simulate dependency between the radial part and the projections using a
    # Gaussian copula
    stopifnot(p >= 2)
    stopifnot(abs(rho) <= 1)
    rhos <- rho^c(1:p)
    sigma <- rbind(c(1, rhos), cbind(rhos, diag(1, nrow = p, ncol = p)))
    x <- mvrnorm(n = n, mu = rep(0, p + 1), Sigma = sigma)
    
    # Simulate the radial part
    radii2 <- qchisq(pnorm(x[, 1]), df = p)
    
    # Simulate the projections
    projs <- x[, -1, drop = FALSE]
    projs <- projs / sqrt(rowSums(projs^2))
    
    # Return normal-like observations
    return(sqrt(radii2) * projs)
    
  }
  
  # Simulate non-spherically symmetric distributions such that the deviation is
  # on the non-uniform distributions of the projections (type = "projs_vMF") or
  # the deviation is on the conditional radii (type = "radii_cond") given the
  # projections
  r_non_sph_sym <- function(n, p, mu = c(1, rep(0, p - 1)), kappa = 0,
                            type = c("projs_vMF", "radii_cond")) {
    
    stopifnot(p >= 2)
    if (type == "projs_vMF") {
      
      projs <- r_vMF(n = n, mu = mu, kappa = kappa)
      radii <- sqrt(rchisq(n = n, df = p))
      
    } else if (type == "radii_cond") {
      
      projs <- r_unif_sphere(n = n, p = p)
      radii <- rnorm(n, mean = 1 + kappa * (1 + projs %*% mu), sd = 1)
      
    }
    return(radii * projs)
    
  }
  
  #choose sample
  
  Erg=array(0,dim=c(length(dimension),length(samplesize)))
  
  true.alpha0=1
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      X=mvrnorm(n,rep(0,d),diag(1,d))
      A=rstable(n,true.alpha0/2,1,2*(cos(pi*true.alpha0/4))^(2/true.alpha0),0,pm=1)
      sample0=matrix(sqrt(A),ncol=d,nrow=n)*X
      
      sample1=rmvt(n,sigma=diag(1,d),df=true.alpha0)
      sample2=rmvt(n,sigma=diag(1,d),df=true.alpha0+0.25)
      sample3=rmvt(n,sigma=diag(1,d),df=true.alpha0+0.5)
      alpha1=0.8
      B=rstable(n,alpha1/2,1,2*(cos(pi*alpha1/4))^(2/alpha1),0,pm=1)
      sample8=matrix(sqrt(B),ncol=d,nrow=n)*r_non_normal(n,d,rho=0)
      sample9=matrix(sqrt(B),ncol=d,nrow=n)*r_non_normal(n,d,rho=0.25)
      sample10=matrix(sqrt(B),ncol=d,nrow=n)*r_non_normal(n,d,rho=0.5)
      
      sample=switch(Choice,sample0,sample1,sample2,sample3,sample8,sample9,sample10)
      j.n=j.n+1
      Erg[j.d,j.n]=stat_hyb(sample, type = "stable",alpha0= true.alpha0)$p.value<0.05
    }
  }
  return(list("Erg"=Erg))
}

set.seed(0815)
param_list=list("Choice"=1:7,"samplesize"=c(50,100), "dimension"=c(50,100))
ZeitAnfang<-Sys.time()
s.hyb.alt1<-MonteCarlo(func=H1.s.hyb, nrep=5000, param_list=param_list, ncpus=14)
summary(s.hyb.alt1)
ZeitEnde<-Sys.time()
difftime(ZeitEnde, ZeitAnfang)

# Get results in one table
Erg=matrix(0,nrow=7,ncol=4)
j.n=0
for (n in c("samplesize=50","samplesize=100"))
{ j.n=j.n+1
Ergebnis=matrix(0,nrow=7,ncol=2)
for (a in 1:7)
{
  Ergebnis[a,1]=mean(as.numeric(s.hyb.alt1$results$Erg[a,j.n,1,]))
  Ergebnis[a,2]=mean(as.numeric(s.hyb.alt1$results$Erg[a,j.n,2,]))
}
if(j.n==1) {Erg=Ergebnis} else {Erg=cbind(Erg,Ergebnis)}
}
xtable::xtable(t(t(round(Erg,2)*100)),digits=0,include.rownames = FALSE)



#-------------------------------------------------------------------------------------
#Simulation HD Sobolev - gamma radii composite hypothesis -> Table 4
#-------------------------------------------------------------------------------------



H1.s.hyb<-function(Choice=1,samplesize=100, dimension=100)
{
  require(MASS)
  require(sphunif)
  require(goftest)
  require(Directional)
  require(rotasym)
  require(gofgamma)
  
  stat_hyb <- function(X, Sobolev_vk2 = c(1, 0), u2 = 2 * (pi^2 / 3 - 3),boot=500,alpha=0.05) {
    
    ## Radiii and projections
    
    # Squared radii
    radii_2 <- rowSums(X^2)
    radii<-sqrt(radii_2)
    
    # Projections
    projs <- X / radii
    
    ## Projections statistic
    
    # Statistic for the projections
    dim(projs) <- c(dim(projs), 1)
    projs_stat <- unif_stat(data = projs, type = "Sobolev",
                            Sobolev_vk2 = Sobolev_vk2)$Sobolev
    
    # Center to have Tn as in the paper
    p <- ncol(X)
    dpk <- d_p_k(p = p, k = seq_along(Sobolev_vk2))
    mean_projs_stat <- sum(Sobolev_vk2 * dpk)
    projs_stat <- projs_stat - mean_projs_stat
    
    # Scale to have Tn / sigman
    sd_projs_stat <- sqrt(2 * sum(Sobolev_vk2^2 * dpk))
    projs_stat <- projs_stat / sd_projs_stat
    
    # Square statistic to have a limiting chi-square
    projs_stat <- 1-pchisq(projs_stat^2,df=1)
    
    ## Radii statistic
    para=gamma_est(radii)
    radii_stat <- AD(radii/para[2], para[1])
    Test.value = u2 * projs_stat + radii_stat
    
    y <- rep(0, boot)
    for (j in 1:boot)
    {
      BS <- rgamma(n, para[1], 1)
      BS_k_estimator <- gamma_est(BS)[1]
      BS_lambda_estimator <- mean(BS)/BS_k_estimator
      x_BS = BS/BS_lambda_estimator
      y[j] = AD(x_BS, BS_k_estimator)
    }
    BS.p.value=mean((y>radii_stat))
    return(list(statistic=-2 * (log(projs_stat) + log(BS.p.value)),p.value=1-pchisq(-2 * (log(projs_stat) + log(BS.p.value)),df=4)))
    
  }
  
  
  
  #distributions
  # Simulate non-normal data due to dependency between the radial part and the
  # projections, both marginally correct
  r_non_normal <- function(n, p, mu = c(1, rep(0, p - 1)), rho = 0) {
    
    # Simulate dependency between the radial part and the projections using a
    # Gaussian copula
    stopifnot(p >= 2)
    stopifnot(abs(rho) <= 1)
    rhos <- rho^c(1:p)
    sigma <- rbind(c(1, rhos), cbind(rhos, diag(1, nrow = p, ncol = p)))
    x <- mvrnorm(n = n, mu = rep(0, p + 1), Sigma = sigma)
    
    # Simulate the radial part
    radii2 <- qchisq(pnorm(x[, 1]), df = p)
    
    # Simulate the projections
    projs <- x[, -1, drop = FALSE]
    projs <- projs / sqrt(rowSums(projs^2))
    
    # Return normal-like observations
    return(sqrt(radii2) * projs)
    
  }
  
  # Simulate non-spherically symmetric distributions such that the deviation is
  # on the non-uniform distributions of the projections (type = "projs_vMF") or
  # the deviation is on the conditional radii (type = "radii_cond") given the
  # projections
  r_non_sph_sym <- function(n, p, mu = c(1, rep(0, p - 1)), kappa = 0,
                            type = c("projs_vMF", "radii_cond")) {
    
    stopifnot(p >= 2)
    if (type == "projs_vMF") {
      
      projs <- r_vMF(n = n, mu = mu, kappa = kappa)
      radii <- sqrt(rchisq(n = n, df = p))
      
    } else if (type == "radii_cond") {
      
      projs <- r_unif_sphere(n = n, p = p)
      radii <- rnorm(n, mean = 1 + kappa * (1 + projs %*% mu), sd = 1)
      
    }
    return(radii * projs)
    
  }
  
  #choose sample
  
  Erg=array(0,dim=c(length(dimension),length(samplesize)))
  
  true.nu=5
  
  j.d=0
  for (d in dimension)
  {
    j.d=j.d+1
    j.n=0
    for (n in samplesize)
    { 
      sample0=rvmf(n,mu=c(1,rep(0,d-1)),k=0)*matrix(rgamma(n,2,5),nrow=n,ncol=d)
      sample1=rvmf(n,mu=c(1,rep(0,d-1)),k=0.25)*matrix(rchisq(n,df=2),nrow=n,ncol=d)
      sample2=rvmf(n,mu=c(1,rep(0,d-1)),k=0.25)*matrix(abs(rcauchy(n,2,5)),nrow=n,ncol=d)
      sample3=rvmf(n,mu=c(1,rep(0,d-1)),k=0.5)*matrix(abs(rt(n,2)),nrow=n,ncol=d)
      sample4=rvmf(n,mu=c(1,rep(0,d-1)),k=10)*matrix(rchisq(n,df=20),nrow=n,ncol=d)
      sample5=rvmf(n,mu=c(1,rep(0,d-1)),k=5)*matrix(rchisq(n,df=d),nrow=n,ncol=d)
      sample6=rvmf(n,mu=c(1,rep(0,d-1)),k=2)*matrix(rgamma(n,2,5),nrow=n,ncol=d)
      sample7=rvmf(n,mu=c(1,rep(0,d-1)),k=5)*matrix(rgamma(n,2,5),nrow=n,ncol=d)
      sample8=rvmf(n,mu=c(1,rep(0,d-1)),k=10)*matrix(rgamma(n,2,5),nrow=n,ncol=d)
      sample9=rvmf(n,mu=c(1,rep(0,d-1)),k=20)*matrix(rgamma(n,2,5),nrow=n,ncol=d)
      
      
      sample=switch(Choice,sample0,sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9)
      j.n=j.n+1
      Erg[j.d,j.n]=stat_hyb(sample,boot=500)$p.value<0.05
    }
  }
  return(list("Erg"=Erg))
}

set.seed(0815)
param_list=list("Choice"=1:10,"samplesize"=c(100,200), "dimension"=c(100,200,300,1000))
ZeitAnfang<-Sys.time()
s.hyb.alt3<-MonteCarlo(func=H1.s.hyb, nrep=5000, param_list=param_list, ncpus=14)
summary(s.hyb.alt3)
ZeitEnde<-Sys.time()
difftime(ZeitEnde, ZeitAnfang)
save(s.hyb.alt3, file = "gamma.RData")

#Get results in a single table
j.n=0
Ergebnis=matrix(0,nrow=10,ncol=8)
for (n in c("samplesize=100","samplesize=200"))
{ j.n=j.n+1
for (a in 1:10)
{
  if (j.n==1) {
    Ergebnis[a,1]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,1,]))
    Ergebnis[a,2]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,2,]))
    Ergebnis[a,3]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,3,]))
    Ergebnis[a,4]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,4,]))
  } else { 
    Ergebnis[a,5]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,1,]))
    Ergebnis[a,6]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,2,]))
    Ergebnis[a,7]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,3,]))
    Ergebnis[a,8]=mean(as.numeric(s.hyb.alt3$results$Erg[a,j.n,4,]))
  }
}
}
xtable::xtable(t(t(round(Ergebnis,2)*100)),digits=0,include.rownames = FALSE)

