## This is code for the sub-functions contained within the PCSALM AECM wrapper function, including:
## Initialization methods, Convergence criteria, Expected value calculations
## Parameter estimation equations, Factor Analyzer Decomposition functions,
## Model selection functions, other fixes required for the operation of
## the AECM algorithm, and functions used to easily sort the model results

# Required libraries ####
library(Bessel)

dsal <- function(p,x,mu,alpha,inv.sig,det.sig){
  d1 <- apply(x,1,function(v){v-mu})
  d2 <- apply(x,1,function(v){log((v-mu)%*%inv.sig%*%(v-mu))})
  nu <- (2-p)/2
  ## Fraction 1
  e1 <- log(2)-(p/2)*log((2*pi))-0.5*log(det.sig)
  e2 <- t(d1)%*%inv.sig%*%alpha
  ## Fraction 2
  e3 <- log(2+c(alpha%*%inv.sig%*%alpha))
  e4 <- (d2-e3)*(nu/2)
  ## Bessel Function
  e5 <- exp(0.5*(e3+d2))
  e6 <- log(besselK(x=e5,nu=nu,expon.scaled=TRUE))-e5
  l.dist <- e1+e2+e4+e6
  dist <- exp(l.dist)
  return(dist)
}

###################################
# SAL k-means Initialization

init.SAL <- function(x,n,p,G){
  index <- kmeans(x=x,centers=G,iter.max=200,nstart=50)$cluster
  wt <- matrix(0,G,n)
  i.parm <- list(mu=matrix(NA,p,G),sigma=array(NA,dim=c(p,p,G)),inv.sig=array(NA,dim=c(p,p,G)),det.sig=numeric(G));
  for(g in 1:G){
    wt[g,which(index == g)] <- 1;
    temp <- cov.wt(x=x,wt=wt[g,],center=TRUE,method="ML")
    i.parm$mu[,g] <- temp$center;
    i.parm$sigma[,,g] <- temp$cov;
    #print(i.parm$sigma[,,g])
    i.parm$inv.sig[,,g] <- solve(i.parm$sigma[,,g]);
    i.parm$det.sig[g] <- det(i.parm$sigma[,,g])
  }
  i.parm$pi.g <- rowMeans(wt)
  return(i.parm)
}

##########################################
# Aitkens Convergence

aitkens <- function(lval,i){
  if(i > 2){
    a <- (lval[i] - lval[i-1])/(lval[i-1]-lval[i-2]);
    inf.l <- lval[i-1] + (1/(1-a))*(lval[i]-lval[i-1]);
    return(inf.l)
  } else {
    inf.l <- lval[i] + 0.02;
    return(inf.l)
  }
}

conv.val <- function(a1,loglik,i,type){
  if(type==1){v <- a1[i]-loglik[i-1]}
  if(type==2){v <- a1[i]-loglik[i]}
  if(type==3){v <- a1[i]-a1[i-1]}
  return(v)
}

####################################
# Expectations: E1[W], E2[1/W]

e1e2 <- function(p,x,mu,alpha,inv.sig){
  a <- 2 + c(alpha%*%inv.sig%*%alpha)
  b <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%(v-mu)})
  t1 <- exp((log(a)+log(b))/2)
  t2 <- log(besselK(x=t1,nu=(2-p)/2 + 1,expon.scaled=TRUE)) - t1
  t3 <- log(besselK(x=t1,nu=(2-p)/2, expon.scaled=TRUE)) - t1
  t4 <- t2 - t3
  t5 <- c(log(b)-log(a))/2
  t6 <- sign((2-p)/2)*exp(log(2)+log(abs((2-p)/2)) - log(b))
  E1 <- exp(t5+t4)
  E2  <- exp(t4-t5) - t6
  out <- list();
  out$E1 <- E1
  out$E2 <- E2
  return(out)
}

##############################################
# INF Likelihood Check

check.val <- function(x,mu){
  v1 <- apply(x,1,function(v){ifelse(sqrt(sum((v-mu)^2))>1e-10 , 0, 1)})
  return(sum(v1))
}

###############################################
# M-Step - SAL EM

mstep <- function(x,p,n.g,zig,E1,E2,mu.star){
  v1 <- zig*E1
  v2 <- zig*E2
  v3 <- zig%*%x
  v4 <- v2%*%x
  d <- sum(v1)*sum(v2)- n.g^2
  mu <- c((sum(v1)*v4 - n.g*v3)/d)
  cv <- check.val(x,mu)
  if(cv == 0){
    alpha <- c((sum(v2)*v3 - n.g*v4)/d)
  }
  else{
    mu <- mu.star
    alpha <- as.vector(zig%*%t(t(x)-mu)/sum(v1))
  }
  #Scale Matrix
  S <- (1/n.g)*matrix(c(colSums(v2*t(apply(x,1,function(v){(v-mu)%o%(v-mu)})))),p,p)
  R <- (-1/n.g)*colSums(zig*t(apply(x,1,function(v){(v-mu)})))%o%alpha
  A <- (1/n.g)*sum(v1)*alpha%o%alpha
  sigma <- S + R + t(R) + A
  inv.sig <- solve(sigma)
  det.sig <- det(sigma)
  out <- list();
  out$mu <- mu;
  out$alpha <- alpha;
  out$sigma <- sigma;
  out$inv.sig <- inv.sig;
  out$det.sig <- det.sig;
  return(out)
}

#########################################
###### Section 2: SAL EM  ##############
#######################################


SAL.EM <- function(x,G,max.it=1000,print.res=FALSE){
  n <- nrow(x)
  p <- ncol(x)
  zig <- E1 <- E2 <- clust <- matrix(0,G,n)
  n.g <- pi.g <- det.sig <- numeric(G)
  sigma <- inv.sig <- array(NA,dim=c(p,p,G))
  mu <- matrix(NA,p,G)
  Mix.dist <- matrix(NA,n,G)
  loglik <- NULL
  a1 <- 0
  a2 <- 1
  i <- 1
  while(abs(a1- a2) > 0.01){
    #cat("SAL Iteration",i,"\n")
    ## Initilization ##
    if(i == 1){
      init.parm <- init.SAL(x=x,n=n,p=p,G=G)
      mu <- init.parm$mu
      sigma <- init.parm$sigma
      inv.sig <- init.parm$inv.sig
      det.sig <- init.parm$det.sig
      pi.g <- init.parm$pi.g
      alpha <- matrix(0,p,G)
    }
    ###### M-step #######
    else{
      n.g <- rowSums(zig)
      pi.g <- rowMeans(zig)
      for (g in 1:G){
        mstep1 <- mstep(x=x,p=p,n.g=n.g[g],zig=zig[g,],E1=E1[g,],E2=E2[g,],mu.star=mu[,g])
        mu[,g] <- mstep1$mu
        alpha[,g] <- mstep1$alpha
        sigma[,,g] <- mstep1$sigma
        inv.sig[,,g] <- mstep1$inv.sig
        det.sig[g] <- mstep1$det.sig
      }                                                                                              
    }
    #### E - Step #####
    for(g in 1:G){
      Mix.dist[,g] <- pi.g[g]*dsal(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g],det.sig=det.sig[g])
    }
    ## Log-Likelihood ##
    loglik[i] <- sum(log(rowSums(Mix.dist)))
    if(i > 2){ if(loglik[i] < loglik[i-1]) stop(paste("Log-Likelihood Decreased on Iteration ",i,sep="")) }
    a1 <- aitkens(lval=loglik,i=i)
    if(i > 1) {a2 <- loglik[i-1]}
    if(i==max.it)break
    ## Expectations ##
    for (g in 1:G){
      zig[g,] <- Mix.dist[,g]/rowSums(Mix.dist)
      e1e2 <- e1e2(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g])
      E1[g,] <- e1e2$E1
      E2[g,] <- e1e2$E2
    }
    i <- i+1
  }
  if(print.res==TRUE) cat("Number of SAL Iterations:",i,"\n")
  for(g in 1:G){
    clust[g,which(round(zig[g,])==1)] <- g
  }
  v <- (G-1) + 2*G*p + G*(p*(p-1))/2
  bic <- 2*loglik[length(loglik)] - v*log(n)
  constraint <- sum(log(apply(zig,2,max)))
  icl <- bic + constraint
  em.out <- list();
  em.out$mu <- mu;
  em.out$alpha <- alpha;
  em.out$sigma <- sigma;
  em.out$n.g <- n.g;
  em.out$zig <- zig;
  em.out$inv.sig <- inv.sig;
  em.out$det.sig <- det.sig;
  em.out$pi.g <- pi.g;
  em.out$icl <- icl
  em.out$Soln <- colSums(clust); 
  return(em.out)
}


################################################
####  Section 3: PCSALM Subfunctions ##########
##############################################


#####################################
# A1: CM Step 1

cm1.step <- function(x,p,n.g,zig,vig,E1,E2,E1.tilde,E2.tilde,eta,mu.star){
  v1 <- zig*(vig*E2+((1-vig)/eta)*E2.tilde)
  v2 <- zig*(vig*E1+(1-vig)*E1.tilde)
  v4 <- zig*(vig+((1-vig)/sqrt(eta)))
  A <- sum(v1)
  B <- sum(v2)
  D <- sum(v4)
  mu <- c((B*(v1%*%x)-D*(v4%*%x))/(B*A-D^2))
  cv <- check.val(x,mu)
  if(cv == 0){
    alpha <- c((A*(v4%*%x)-D*(v1%*%x))/(B*A-D^2))
  }
  else{
    mu <- mu.star
    alpha <- c(v4%*%t(t(x)-mu)/B)
  }
  out <- list();
  out$mu <- mu;
  out$alpha <- alpha;
  return(out)
}

#####################################
# A1: CM Step 2

cm2.step <- function(x,p,zig,vig,E2.tilde,mu,alpha,inv.sig,eta.cap){
  v1 <- zig*(1-vig)
  v2 <- zig*(1-vig)*E2.tilde
  v3 <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%(v-mu)})
  v4 <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%alpha})
  a <- p*sum(v1)
  b <- c(v1%*%v4)
  c <- c(v2%*%v3)
  sq.eta <- (-b+sqrt(b^2+4*a*c))/(2*a)
  ###############################
  eta <- max(1,sq.eta^2)
  eta <- min(eta,eta.cap)
  return(eta) 
}

cm2.const1 <- function(x,p,zig,vig,E2.tilde,mu,alpha,inv.sig){
  v1 <- zig*(1-vig)
  v2 <- zig*(1-vig)*E2.tilde
  v3 <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%(v-mu)})
  v4 <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%alpha})
  a <- p*sum(v1)
  b <- c(v1%*%v4)
  c <- c(v2%*%v3)
  return(c(a,b,c))
}

cm2.const2 <- function(coef,eta.cap){
  a <- coef[1]
  b <- coef[2]
  c <- coef[3]
  sq.eta <- (-b+sqrt(b^2+4*a*c))/(2*a)
  eta <- max(1,sq.eta^2)
  eta <- min(eta,eta.cap)
  return(eta) 
}

##########################################
# FA Decomp Functions

woodbury <- function(q,Lambda,tLambda,Psi){
  i.Psi <- diag(diag(1/Psi))
  id <- solve(diag(q)+tLambda%*%i.Psi%*%Lambda)
  inv.sig <- i.Psi - i.Psi%*%Lambda%*%id%*%tLambda%*%i.Psi
  d1 <- diag(q)-tLambda%*%inv.sig%*%Lambda
  det.sig <- det(Psi)/det(d1)
  out <- list(inv.sig=inv.sig,det.sig=det.sig)
  return(out)
}

Sg <- function(x,p,n.g,zig,vig,E1,E2,E1.tilde,E2.tilde,eta,mu,alpha){
  v1 <- zig*(vig*E2+((1-vig)/eta)*E2.tilde)
  v2 <- zig*(vig*E1+(1-vig)*E1.tilde)
  v4 <- zig*(vig+((1-vig)/sqrt(eta)))
  B <- sum(v2)
  #Covariance
  m1 <- (1/n.g)*matrix(c(colSums(v1*t(apply(x,1,function(v){(v-mu)%o%(v-mu)})))),p,p)
  m2 <- (-1/n.g)*colSums(v4*t(apply(x,1,function(v){(v-mu)})))%o%alpha
  m3 <- (1/n.g)*B*alpha%o%alpha
  S <- m1+m2+t(m2)+m3
  return(S)
}

s.tilde <- function(p,G,pi.g,S){
  S.tilde <- matrix(0,p,p)
  for(g in 1:G){
    S.tilde <- S.tilde + pi.g[g]*S[,,g]
  }
  return(S.tilde)
}

##############################################
# Model Selection

free.covpar <- function(G,p,q,model){
  if(model == "CCCC") out.p <- (p*q-q*(q-1))/2 + 1
  if(model == "CCCU") out.p <- (p*q-q*(q-1))/2 + p	
  if(model == "CCUC") out.p <- (p*q-q*(q-1))/2 + G
  if(model == "CCUU") out.p <- (p*q-q*(q-1))/2 + (G+(p-1))
  if(model == "CUCU") out.p <- (p*q-q*(q-1))/2 + (1+G*(p-1))
  if(model == "CUUU") out.p <- (p*q-q*(q-1))/2 + G*p
  if(model == "UCCC") out.p <- G*(p*q-q*(q-1))/2 + 1
  if(model == "UCCU") out.p <- G*(p*q-q*(q-1))/2 + p
  if(model == "UCUC") out.p <- G*(p*q-q*(q-1))/2 + G
  if(model == "UCUU") out.p <- G*(p*q-q*(q-1))/2 + (G+(p-1))
  if(model == "UUCU") out.p <- G*(p*q-q*(q-1))/2 + (1+G*(p-1))
  if(model == "UUUU") out.p <- G*(p*q-q*(q-1))/2 + G*p
  return(out.p)
}

model.select <- function(like,G,p,n,q,model,zig,vig,soln){
  k <- free.covpar(G=G,p=p,q=q,model=model)
  v <- (G-1) + 2*G*(p+1) + k
  c <- length(like)
  bic <- 2*like[c] - v*log(n)
  aic <- 2*like[c] - v*2
  aic3 <- 2*like[c] - v*3
  constraint <- sum(log(apply(zig,2,max)));
  icl <- bic + constraint
  constraint2 <- sum(log(apply(vig,2,max)));
  n.good <- length(which(soln!=0))
  icl.mod <- 2*like[c] - v*log(n) + constraint + constraint2
  out <- list();
  out$bic <- bic
  out$aic <- aic
  out$aic3 <- aic3
  out$icl <- icl
  out$icl.mod <- icl.mod
  return(out)
}

#####################################################
# Labeling for Solution

sol.labels <- function(zig,zig.star,vig,vig.star,n,G){
  n.check <- length(which(rowSums(zig) < 5))
  if(n.check!=0){ zig <- zig.star; 
  #print("Locked n.g")
  }
  clust <- good <- matrix(0,G,n)
  for(g in 1:G){
    clust[g,which(round(zig[g,])==1)] <- g
    for(N in 1:n){
      if(clust[g,N] == g){
        good[g,N] <- round(vig[g,N])
      }
      else{good[g,N] <- 0}
    }
  }
  out <- list();
  out$clust <- clust
  out$good <- good
  return(out)
}

sol.labels2 <- function(zig,zig.star,vig,vig.star,n,G){
  n.check <- length(which(rowSums(zig) < 5))
  if(n.check!=0){ zig <- zig.star; 
  #print("Locked n.g")
  }
  clust <- apply(zig,2,function(x){which(x==max(x))})
  good <- rep(0,n)
  for(g in 1:G){
    good[which(clust == g)] <- round(vig[g,which(clust == g)])
  }
  out <- list();
  out$clust <- clust
  out$good <- good
  return(out)
}




################################################
###  Section 4: Factor Analyzer Decomp ########
##############################################


###############################################
# FA Decomp initilizations

init.Parm <- function(p,q,G,S){
  SE <- eigen(S)
  #print(eigen(S))
  if(length(which(SE$values < 0)) > 0)stop(paste("S matrix is not positive definite"))
  tLambda <- sqrt(SE$values)[1:q]*t(SE$vectors)[1:q,]
  Lambda <- t(tLambda)
  if(q == 1){ pd <- tLambda%*%Lambda }
  else{ pd <- Lambda%*%tLambda }
  Psi <- diag(diag(S-pd))
  w <- det(Psi)^(1/p)
  Delta <- Psi/w
  out <- list()
  if(q == 1){
    out$Lambda <- c(Lambda)
    out$tLambda <- c(tLambda)
    out$pd <- pd
  }
  else{
    out$Lambda <- Lambda
    out$tLambda <- tLambda
    out$pd <- pd
  }
  out$w <- w
  out$Delta <- Delta
  out$Psi <- Psi
  return(out)
}

##############################################
## CCCC / CCCU / UCCC / UCCU / UCUC / UUUU ##

decomp.Parm <- function(p,q,Lambda,tLambda,S,inv.sig){
  Beta <- tLambda%*%inv.sig
  tBeta <- t(Beta)
  Theta <- diag(q)-(Beta%*%Lambda)+(Beta%*%S%*%tBeta)
  inv.Theta <- solve(Theta) 
  Lambda <- S%*%tBeta%*%inv.Theta
  tLambda <- t(Lambda)
  v1 <- diag(S-Lambda%*%Beta%*%S)
  tr <- (1/p)*sum(v1)
  dg <- diag(v1)
  out <- list()
  out$Lambda <- Lambda
  out$tLambda <- tLambda
  out$tr <- tr*diag(p)
  out$dg <- dg
  return(out)
}

##################
## UCUU / UUCU ##

decomp2.Parm <- function(p,q,Lambda,tLambda,inv.Delta,S,inv.sig){
  Beta <- tLambda%*%inv.sig
  tBeta <- t(Beta)
  Theta <- diag(q)-(Beta%*%Lambda)+(Beta%*%S%*%tBeta)
  inv.Theta <- solve(Theta)
  Lambda <- S%*%tBeta%*%inv.Theta
  tLambda <- t(Lambda)
  A <- S-Lambda%*%Beta%*%S
  v2 <- diag(inv.Delta%*%A)
  w <- (1/p)*sum(v2)
  out <- list()
  out$Lambda <- Lambda
  out$tLambda <- tLambda
  out$dg <- diag(diag(A))
  out$w <- w
  return(out)
}

################################
## CCUU / CCUC / CUCU / CUUU ##

# Lambda - CCUU / CCUC #

cc.lambda <- function(p,q,n.g,w,Lambda,tLambda,S,inv.sig){
  Beta <- tLambda%*%inv.sig
  tBeta <- t(Beta)
  Theta <- diag(q)-(Beta%*%Lambda)+(Beta%*%S%*%tBeta)
  sum1 <- (n.g/w)*S%*%tBeta
  sum2 <- (n.g/w)*Theta
  out <- list()
  out$Beta <- Beta
  out$Theta <- Theta
  out$sum1 <- sum1
  out$sum2 <- sum2
  return(out)
}

# Lambda CUCU / CUUU #

cu.lambda <- function(p,q,n.g,Delta,Lambda,tLambda,S,inv.sig){
  Beta <- tLambda%*%inv.sig
  tBeta <- t(Beta)
  Theta <- diag(q)-(Beta%*%Lambda)+(Beta%*%S%*%tBeta)
  inv.delta <- diag(1/diag(Delta))
  R <- n.g*inv.delta%*%S%*%tBeta
  sum1 <- array(NA,dim=c(q,q,p))
  for(j in 1:p){
    sum1[,,j] <- n.g*inv.delta[j,j]*Theta
  }
  out <- list()
  out$Beta <- Beta
  out$Theta <- Theta
  out$R <- R
  out$sum1 <- sum1
  return(out)
}

# Psi - CCUC #

ccuc.psi <- function(p,Lambda,tLambda,S,Beta,Theta){
  A <- S-2*Lambda%*%(Beta%*%S)
  A1 <- A + Lambda%*%(Theta%*%tLambda)
  tr <- (1/p)*sum(diag(A1))
  return(tr)
}

# Psi - CCUU / CUCU / CUUU #

not.ccuc.psi <- function(p,Delta,Lambda,tLambda,S,Beta,Theta){
  A <- S-2*Lambda%*%(Beta%*%S)+Lambda%*%(Theta%*%tLambda)
  v1 <- diag(solve(Delta)%*%A)
  w <- (1/p)*sum(v1)
  dg <- diag(diag(A))
  out <- list()
  out$w <- w
  out$dg <- dg
  return(out)
}

#########################################
###  Section 5: Sorting Results ########
#######################################

aecm.output <- function(aecm,id){
  output <- list();
  if(is.list(aecm)==F){
    output$converge <- F
    output$Soln <- NA
    output$loglik <- NA
    output$aic <- NA
    output$aic3 <- NA
    output$bic <- NA
    output$icl <- NA
    output$icl.mod <- NA
    output$penalty <- NA
    output$rho <- NA
    output$tab <- NA
    output$ari <- NA
    output$time <- NA
    output$id <- id
  }
  if(is.list(aecm)==T){
    output$converge <- T
    output$Soln <- aecm$Soln
    output$loglik <- aecm$loglik
    output$aic <- aecm$aic
    output$aic3 <- aecm$aic3
    output$bic <- aecm$bic
    output$icl <- aecm$icl
    output$icl.mod <- aecm$icl.mod
    output$penalty <- aecm$penalty
    output$rho <- aecm$rho
    output$tab <- table(aecm$Soln,label)
    output$ari <- ARI(aecm$Soln,label)
    output$time <- aecm$time
    output$id <- id
  }
  return(output)
}


sim.select <- function(list,n=192){
  bic <- icl <- icl.mod <- ari <- time <- conv <- NULL
  for(i in 1:n){
    bic[i] <- c(list[[i]]$bic)
    icl[i] <- c(list[[i]]$icl)
    icl.mod[i] <- c(list[[i]]$icl.mod)
    ari[i] <- c(list[[i]]$ari)
    time[i] <- c(list[[i]]$time)
    conv[i] <- c(list[[i]]$converge)
  }
  out <- list();
  out$converge <- which(conv==T)
  out$max.bic <- which(bic==max(bic,na.rm=T))
  out$max.icl <- which(icl==max(icl,na.rm=T))
  out$max.icl.mod <- which(icl.mod==max(icl.mod,na.rm=T))
  out$max.ari <- which(ari==max(ari,na.rm=T))
  out$max.time <- which(time==min(time,na.rm=T))
  return(out)
}
