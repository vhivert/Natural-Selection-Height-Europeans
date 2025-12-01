#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# num.pops=13
# nSNP=1000

# Parse arguments (ensure correct order)
path      <- args[1]
num.pops      <- as.integer(args[2])
nSNP      <- as.integer(args[3])
h2_eq     <- as.numeric(args[4])
spousal.corr <- as.numeric(args[5])
n.reps <- as.numeric(args[6])
num.cores     <- as.integer(args[7])

library(parallel)

cat("Arguments received:\n")
cat(paste("path:", path, "\nh2_eq:", h2_eq, "\nnpop:", num.pops, "\nnSNP:", nSNP, "\nr:", spousal.corr, "\nnrep:", n.reps, "\ncores:", num.cores, "\n"))


setwd(path)
simAM <- function(N=10000,M=1000,h2_eq=0.5,r=0.2,ng=10,npop=13){
  ## Theory
  varY_eq <- 1
  rho     <- r * h2_eq
  varE    <- (1 - h2_eq) * varY_eq
  sdE     <- sqrt(varE)
  varA_eq <- h2_eq * varY_eq
  varA_0  <- varA_eq * (1 - (1-0.5/M) * rho) # Eq. 7.19a Lynch & Walsh
  varY_0  <- varA_0 + varE
  h2_0    <- varA_0 / (varA_0 + varE)
  Vb      <- varA_0 / M
  sdVb    <- sqrt(Vb)
  alpha   <- (varA_eq/varA_0-1)/(2*M-1) # rho/(2*M*(1-rho)+rho)
  theta   <- varA_eq/varA_0
  
  ## Tools
  nmates<- 100
  rsq   <- r * r
  s1_rsq<- sqrt(1-rsq) 
  p     <- runif(n=M,min=0.01,max=0.99)
  p0    <- runif(n=M,min=0.01,max=0.99)
  dp    <- 2 * p
  sh    <- sqrt(dp*(1-p))
  BETA  <- abs( rnorm(n=M,mean=0,sd=sdVb/sh) )
  m_g   <- sum(dp*BETA)
  FREQ1 <- matrix(NA,nrow=M,ncol=npop)
  FREQ0 <- matrix(NA,nrow=M,ncol=npop)
  
  colnames(FREQ1) <- colnames(FREQ0) <- paste0("FREQ_pop",1:npop)

  for(pop in 1:npop){
    # cat(paste0("Simulating pop#",pop,".\n"))
    ## Simulate base population
    X  <- do.call("cbind",lapply(1:M,function(j) rbinom(N,size=2,prob=p[j])))
    X0 <- do.call("cbind",lapply(1:M,function(j) rbinom(N,size=2,prob=p0[j])))
    g  <- c(X %*% BETA)
    e  <- rnorm(N,0,sd=sdE)
    Y  <- (g-m_g) + e
    ## Start AM generations
    for(igen in 1:ng){
      sex   <- rbinom(n=N,size=1,prob=0.5)
      indm  <- which(sex==1)
      indf  <- which(sex==0)
      Ys    <- (Y-mean(Y))/sd(Y)
      idm1  <- sample(indm,N,replace = TRUE)
      mates <- do.call("rbind",lapply(1:N,function(i) sample(indf,nmates,replace = TRUE)))
      idm2  <- do.call("rbind",lapply(1:N,function(i){
        x  <- Ys[idm1[i]]; y <- Ys[mates[i,]]
        pc <- exp( 0.5*(2.*r*x*y-rsq*(x*x + y*y))/(1-rsq) )/s1_rsq
        pc <- pc / sum(pc)
        return(sample(mates[i,],1,prob = pc))
      }))
      rp <- cor(Ys[idm1],Ys[idm2])
      # cat(paste0("\tr_emp = ",rp,".\n"))

      ## Offspring
      X1 <- do.call("cbind",lapply(1:M,function(j) rbinom(N,1,prob=0.5*X[idm1,j])))
      X2 <- do.call("cbind",lapply(1:M,function(j) rbinom(N,1,prob=0.5*X[idm2,j])))
      X  <- X1+X2
      g  <- c(X %*% BETA)
      e  <- rnorm(N,0,sd=sdE)
      Y  <- (g-m_g) + e
      ## Neutral variants
      X10 <- do.call("cbind",lapply(1:M,function(j) rbinom(N,1,prob=0.5*X0[idm1,j])))
      X20 <- do.call("cbind",lapply(1:M,function(j) rbinom(N,1,prob=0.5*X0[idm2,j])))
      X0  <- X10 + X20
    }
    
    ## Calculate allele frequencies
    FREQ1[,pop] <- colMeans(X)/2
    FREQ0[,pop] <- colMeans(X0)/2
  }
  
  RESULTS <- cbind.data.frame(ANC_FREQ=c(p,p0),BETA=c(BETA,rep(0,M)),rbind(FREQ1,FREQ0))
  
  return(RESULTS)
}

Qx_func <- function(simu, num.pops=10, nSNP=1000){
  neutral=simu[c(nSNP+1):nrow(simu),]
  neutral.pop=(neutral[,-c(1,2)]) #An L x num.pop matrix
  epsilon=rowMeans(neutral.pop)
  var.cov.snps <- epsilon * ( 1 - epsilon)
  scaled.snp.by.pop  <- neutral.pop / c ( sqrt ( var.cov.snps ))
 
  uncentered.cov.mat <- cov ( scaled.snp.by.pop)
  
  S=diag(1/var.cov.snps)
  T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops - 1 ) ) , times = num.pops - 1 ) , ncol = num.pops , nrow = num.pops - 1 )
  #uncentered.cov.mat.2=(1/999)*t(neutral.pop)%*%S%*%as.matrix(neutral.pop)
  #F=T.mat%*%t(neutral.pop)%*%S%*%as.matrix(neutral.pop)%*%t(T.mat)
  
  freqs = t(simu[1:nSNP,-c(1,2)])
  effects = simu$BETA[1:nSNP]
  epsilon.freqs=colMeans(freqs)
  var <- sum ( effects^2  * epsilon.freqs * ( 1 - epsilon.freqs ) )
  
  # mean center and drop a pop
  F.mat <- T.mat %*% ( uncentered.cov.mat ) %*% t ( T.mat )
  F.inv <- solve ( F.mat )
  C.inv <- solve ( t ( chol ( F.mat ) ) )
  
  # Qx and components
  contribs <- t ( t ( freqs ) * effects )
  std.contribs <- C.inv %*% T.mat %*% contribs
  adj.cov.mat <- t ( std.contribs ) %*% std.contribs
  Fst.comp <- sum ( diag ( adj.cov.mat ) ) / var
  LD.comp <- sum ( adj.cov.mat [ row ( adj.cov.mat ) != col ( adj.cov.mat ) ] ) / var
  Qx <- Fst.comp + LD.comp 
  return(c(Fst.comp,LD.comp,Qx))
}

out=paste0("Simu_AM_",num.pops,"Pops_10Gens_",nSNP,"SNPs_",h2_eq,"h2eq_r",spousal.corr,".txt")

results <- do.call(rbind, mclapply(1:n.reps, function(i) {
  simu <- simAM(N = 1000, M = nSNP, h2_eq = h2_eq, r = spousal.corr, ng = 10, npop = num.pops)
  Qx_func(simu, num.pops = num.pops, nSNP = nSNP)
}, mc.cores = num.cores))

# Label columns
colnames(results) <- c("Fst_component", "LD_component", "Qx")

write.table(results,out,quote = F,sep = "\t",row.names = F,col.names = T)
##################################################
