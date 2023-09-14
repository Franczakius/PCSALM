## This code is the PCSALM Wrapper function
## It contains sub-functions found in: PCSALM | Required Functions.R

# Required Functions ####
library(timeR)
library(MixGHD)

PCSALM <- function(x,q,G,model,init,max.it=2000,c.crit=1e-2,c.type=1,eta.cap=1000,eta.C=F,print.res=FALSE){
  timer1 <- createTimer(F,precision = "ms")
  timer1$start("event1")
  n <- nrow(x)
  p <- ncol(x)
  if( p<=q ){stop(print("Number of factors must be less than the number of variables: q < p"))}
  if( missing(model)){model="UUUU"}
  if(model==1){model <- "CCCC"}
  if(model==2){model <- "CCCU"}
  if(model==3){model <- "CCUC"}
  if(model==4){model <- "CCUU"}
  if(model==5){model <- "CUCU"}
  if(model==6){model <- "CUUU"}
  if(model==7){model <- "UCCC"}
  if(model==8){model <- "UCUC"}
  if(model==9){model <- "UCCU"}
  if(model==10){model <- "UCUU"}
  if(model==11){model <- "UUCU"}
  if(model==12){model <- "UUUU"}
  if( model!="UCCC" & model!="UCCU" & model!="UCUC" & model!="UCUU" & model!="UUCU" & model!="UUUU" & model!="CCCC" & model!="CCCU" & model!="CCUC" & model!="CCUU" & model !="CUCU" & model!="CUUU"){
    stop(print("Ineligble model: Factor Analyzer Decomp. DNE"))
  }
  zig <- vig <- E1 <- E2 <- E1.tilde <- E2.tilde <- matrix(0,G,n)
  S <- array(NA,dim=c(p,p,G))
  mu <- matrix(NA,p,G)
  Mix.dist <- Group.dist <- Good.dist <- Bad.dist <- matrix(NA,n,G)
  loglik <- a1 <- n.g <- pi.g <- det.sig <- w <- NULL
  conv.value <- 1
  i <- 1
  if(model=="CCCC"||model=="CCCU"){
    inv.sig <- matrix(0,p,p)
    det.sig <- 0
    while(conv.value > c.crit){
      #cat("PCSALM Iteration",i,"\n")
      ## Initialization ##
      if(i == 1){
        if(missing(init)){
          init.parm <- SAL.EM(x=x,G=G,max.it=200,print.res=print.res)
        }
        else{init.parm <- init}
        pi.g <- init.parm$pi.g
        mu <- init.parm$mu
        n.g <- init.parm$n.g
        zig.star <- init.parm$zig
        for(g in 1:G){
          inv.sig <- inv.sig + pi.g[g]*init.parm$inv.sig[,,g]
          det.sig <- det.sig + pi.g[g]*init.parm$det.sig[g]
        }
        alpha <- matrix(0,p,G)
        rho <- rep(0.999,G)
        eta <- rep(1.001,G)
      }
      ## M - Steps ##
      else {
        ## Alternation 1: CM1 and CM2 steps ##
        n.g <- rowSums(zig)
        pi.g <- rowMeans(zig)
        rho <- rowSums(zig*vig)/n.g
        eta.coef <- numeric(3)
        ## Check Group Size ##
        n.check <- length(which(n.g < 5))
        #print(n.check)
        if(n.check==0){
          zig.star <- zig
        }
        if(n.check!=0){
          zig <- zig.star
          n.g <- rowSums(zig)
        }
        for (g in 1:G){
          #rho[g] <- max(c(0.5,rho[g]))
          cmstep1 <- cm1.step(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu.star=mu[,g])
          mu[,g] <- cmstep1$mu
          alpha[,g] <- cmstep1$alpha
          if(eta.C == F){ 
            eta[g] <- cm2.step(x=x,p=p,zig=zig[g,],vig=vig[g,],E2.tilde=E2.tilde[g,],mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig,eta.cap=eta.cap)
            S[,,g] <- Sg(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu=mu[,g],alpha=alpha[,g])
          }
          if(eta.C == T){ eta.coef <- eta.coef + pi.g[g]*cm2.const1(x=x,p=p,zig=zig[g,],vig=vig[g,],E2.tilde=E2.tilde[g,],mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig) }
        }
        if(eta.C == T){ 
          eta1 <- cm2.const2(coef=eta.coef,eta.cap=eta.cap)
          for(g in 1:G){
            eta[g] <- eta1
            S[,,g] <- Sg(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu=mu[,g],alpha=alpha[,g])
          }
        }
        ## Alternation 2: M step ##
        S.tilde <- s.tilde(p=p,G=G,pi.g=pi.g,S=S)
        if(i == 2){
          init.parm <- init.Parm(p=p,q=q,G=G,S=S.tilde)
          Lambda <- init.parm$Lambda
          tLambda <- init.parm$tLambda
          inv.sig <- solve(init.parm$pd+init.parm$Psi)
        }
        dc.parm <- decomp.Parm(p=p,q=q,Lambda=Lambda,tLambda=tLambda,S=S.tilde,inv.sig=inv.sig)
        Lambda <- dc.parm$Lambda
        tLambda <- dc.parm$tLambda
        if(model=="CCCU"){ Psi <- dc.parm$dg }
        if(model=="CCCC"){ Psi <- dc.parm$tr }
        wood <- woodbury(q=q,Lambda=Lambda,tLambda=tLambda,Psi=Psi)
        inv.sig <- wood$inv.sig
        det.sig <- wood$det.sig
      }
      ## E - Step ##
      for (g in 1:G){
        Good.dist[,g] <- rho[g]*dsal(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig,det.sig=det.sig)
        Good.dist[,g] <- sapply(Good.dist[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
        Bad.dist[,g] <- dsal(p=p,x=x,mu=mu[,g],alpha=(sqrt(eta[g])*alpha[,g]),inv.sig=((1/eta[g])*inv.sig),det.sig=((eta[g]^p)*det.sig))
        Bad.dist[,g] <- sapply(Bad.dist[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
        Group.dist[,g] <- Good.dist[,g]+(1-rho[g])*Bad.dist[,g]
        Mix.dist[,g] <- pi.g[g]*Group.dist[,g]
      }
      ## Log-likelihood ##
      loglik[i] <- sum(log(rowSums(Mix.dist)))
      if(i > 2){ if(loglik[i] < loglik[i-1]) stop(paste("Log-Likelihood Decreased on Iteration ",i,sep=""))}
      a1[i] <- aitkens(lval=loglik,i=i)
      if(i > 2){conv.value <- conv.val(a1=a1,loglik=loglik,i=i,type=c.type)}
      ## Expectations ##
      for (g in 1:G){
        zig[g,] <- Mix.dist[,g]/rowSums(Mix.dist)
        vig[g,] <- Good.dist[,g]/Group.dist[,g]
        e1e2 <- e1e2(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig)
        E1[g,] <- e1e2$E1
        E2[g,] <- e1e2$E2
        e1e2.tilde <- e1e2(p=p,x=x,mu=mu[,g],alpha=sqrt(eta[g])*alpha[,g],inv.sig=(1/eta[g])*inv.sig)
        E1.tilde[g,] <- e1e2.tilde$E1
        E2.tilde[g,] <- e1e2.tilde$E2
      }
      if(i==max.it)break
      i <- i+1
    }
    sigma <- solve(inv.sig)
  }
  else{
    if(model=="UCCC"||model=="UCCU"||model=="UCUC"||model=="UCUU"||model=="UUCU"||model=="UUUU"){
      init.Lambda <- array(NA,dim=c(p,q,G))
      init.tLambda <- array(NA,dim=c(q,p,G))
      init.Psi <- init.Delta <- init.pd <- array(NA,dim=c(p,p,G))
      init.w <- numeric(G)
      tLambda <- array(NA,dim=c(q,p,G))
      Lambda <- array(NA,dim=c(p,q,G))
      inv.sig <- sigma <- Psi <- array(NA,dim=c(p,p,G))
      if(model=="UCUU"){
        w <- numeric(G)
        dg <- array(NA,dim=c(p,p,G))
      }
      if(model=="UUCU"){
        w <- numeric(G)
        k <- numeric(G)
        Delta <- inv.Delta <- array(NA,dim=c(p,p,G))
      }
    }
    if(model=="CCUC"||model=="CCUU"||model=="CUCU"||model=="CUUU"){
      init.Lambda <- array(NA,dim=c(p,q,G))
      init.tLambda <- array(NA,dim=c(q,p,G))
      init.Psi <- init.Delta <- init.pd <- array(NA,dim=c(p,p,G))
      init.w <- numeric(G)
      tBeta <- array(NA,dim=c(p,q,G))
      Beta <- array(NA,dim=c(q,p,G))
      Theta <- tTheta <- array(NA,dim=c(q,q,G))
      inv.sig <- sigma <- Psi <- array(NA,dim=c(p,p,G))
      if(model=="CUCU"||model=="CUUU"){
        k <- numeric(G)
        Lambda <- matrix(NA,p,q)
        Delta <- dg <- array(NA,dim=c(p,p,G))
      }
    }
    while(conv.value > c.crit){
      #cat("PCSALM Iteration",i,"\n")
      ## Initilization ##
      if(i == 1){
        if(missing(init)){
          init.parm <- SAL.EM(x=x,G=G,max.it=200)
        }
        else{init.parm <- init}
        pi.g <- init.parm$pi.g
        mu <- init.parm$mu
        ng.star <- init.parm$n.g
        zig.star <- init.parm$zig
        inv.sig <- init.parm$inv.sig
        det.sig <- init.parm$det.sig
        alpha <- matrix(0,p,G)
        rho <- rep(0.999,G)
        eta <- rep(1.001,G)
      }
      ## M - Steps ##
      else {
        ## Alternation 1: CM1 and CM2 Steps ##
        n.g <- rowSums(zig)
        pi.g <- rowMeans(zig)
        rho <- rowSums(zig*vig)/n.g
        eta.coef <- numeric(3)
        ## Check Group Size ##
        n.check <- length(which(n.g < 5))
        if(n.check==0){
          zig.star <- zig
          ng.star <- n.g
        }
        if(n.check!=0){
          zig <- zig.star
          n.g <- ng.star
        }
        for (g in 1:G){
          #rho[g] <- max(c(0.5,rho[g]))
          cmstep1 <- cm1.step(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu.star=mu[,g])
          mu[,g] <- cmstep1$mu
          alpha[,g] <- cmstep1$alpha
          if(eta.C == F){ 
            eta[g] <- cm2.step(x=x,p=p,zig=zig[g,],vig=vig[g,],E2.tilde=E2.tilde[g,],mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g],eta.cap=eta.cap)
            S[,,g] <- Sg(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu=mu[,g],alpha=alpha[,g])
          }
          if(eta.C == T){ eta.coef <- eta.coef + pi.g[g]*cm2.const1(x=x,p=p,zig=zig[g,],vig=vig[g,],E2.tilde=E2.tilde[g,],mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g]) }
        }
        if(eta.C == T){ 
          eta1 <- cm2.const2(coef=eta.coef,eta.cap=eta.cap)
          for(g in 1:G){
            eta[g] <- eta1
            S[,,g] <- Sg(x=x,p=p,n.g=n.g[g],zig=zig[g,],vig=vig[g,],E1=E1[g,],E2=E2[g,],E1.tilde=E1.tilde[g,],E2.tilde=E2.tilde[g,],eta=eta[g],mu=mu[,g],alpha=alpha[,g])
          }
        }
        ## Alternation 2: M-Step ##
        if(i==2){
          for(g in 1:G){
            init.parm <- init.Parm(p=p,q=q,G=G,S=S[,,g])
            init.Lambda[,,g] <- init.parm$Lambda
            init.tLambda[,,g] <- init.parm$tLambda
            init.pd[,,g] <- init.parm$pd
            init.Psi[,,g] <- init.parm$Psi
            init.Delta[,,g] <- init.parm$Delta
            init.w[g] <- init.parm$w
          }
        }
        ## FA Decomposition ##
        if(model=="UCCC"||model=="UCCU"||model=="UCUC"||model=="UUUU"){
          if(model=="UCCU"||model=="UCCC"){ Psi <- matrix(0,p,p)}
          for(g in 1:G){
            if(i == 2){
              Lambda[,,g] <- init.Lambda[,,g]
              tLambda[,,g] <- init.tLambda[,,g]
              inv.sig[,,g] <- solve(init.pd[,,g]+init.Psi[,,g])
            }
            dc.parm <- decomp.Parm(p=p,q=q,Lambda=Lambda[,,g],tLambda=tLambda[,,g],S=S[,,g],inv.sig=inv.sig[,,g])            
            Lambda[,,g] <- dc.parm$Lambda
            tLambda[,,g] <- dc.parm$tLambda
            if(model=="UUUU"){ Psi[,,g] <- dc.parm$dg }
            if(model=="UCUC"){ Psi[,,g] <- dc.parm$tr }
            if(model=="UCCU"){ Psi <- Psi + pi.g[g]*dc.parm$dg }
            if(model=="UCCC"){ Psi <- Psi + pi.g[g]*dc.parm$tr }
          }
        }
        if(model=="UCUU"){
          if(i==2){
            Delta <- matrix(0,p,p)
            for(g in 1:G){
              Delta <- Delta + pi.g[g]*init.Delta[,,g]
            }
            for(g in 1:G){
              Lambda[,,g] <- init.Lambda[,,g]
              tLambda[,,g] <- init.tLambda[,,g]
              Psi[,,g] <- init.w[g]*Delta
              inv.sig[,,g] <- solve(init.pd[,,g]+Psi[,,g])
            }
          }
          inv.Delta <- solve(Delta)
          dg.sum <- matrix(0,p,p)
          for(g in 1:G){
            dc.parm <- decomp2.Parm(p=p,q=q,Lambda=Lambda[,,g],tLambda=tLambda[,,g],inv.Delta=inv.Delta,S=S[,,g],inv.sig=inv.sig[,,g])            
            Lambda[,,g] <- dc.parm$Lambda
            tLambda[,,g] <- dc.parm$tLambda
            w[g] <- dc.parm$w
            dg.sum <- dg.sum + (n.g[g]/w[g])*dc.parm$dg
          }
          k <- prod(diag(dg.sum)^(1/p))
          Delta <- (1/k)*dg.sum
          for(g in 1:G){ Psi[,,g] <- w[g]*Delta }  
        }
        if(model=="UUCU"){
          if(i == 2){
            w.sum <- 0 
            for(g in 1:G){
              w.sum <- w.sum + pi.g[g]*init.w[g]
            }
            for(g in 1:G){
              Lambda[,,g] <- init.Lambda[,,g]
              tLambda[,,g] <- init.tLambda[,,g]
              Delta[,,g] <- init.Delta[,,g]
              Psi[,,g] <- w.sum*Delta[,,g]
              inv.sig[,,g] <- solve(init.pd[,,g]+Psi[,,g])
            }
          }
          w.sum <- 0 
          for(g in 1:G){
            inv.Delta[,,g] <- solve(Delta[,,g])
            dc.parm <- decomp2.Parm(p=p,q=q,Lambda=Lambda[,,g],tLambda=tLambda[,,g],inv.Delta=inv.Delta[,,g],S=S[,,g],inv.sig=inv.sig[,,g])            
            Lambda[,,g] <- dc.parm$Lambda
            tLambda[,,g] <- dc.parm$tLambda
            w.sum <- w.sum + pi.g[g]*dc.parm$w
            k[g] <- prod(diag(dc.parm$dg)^(1/p))
            Delta[,,g] <- (1/k[g])*dc.parm$dg
          }
          for(g in 1:G){ Psi[,,g] <- w.sum*Delta[,,g] }
        }
        if(model=="CCUC"||model=="CCUU"){
          if(i == 2){
            Delta <- matrix(0,p,p)
            Lambda <- matrix(0,p,q)
            tLambda <- matrix(0,q,p)
            for(g in 1:G){
              Lambda <- Lambda + pi.g[g]*init.Lambda[,,g]
              tLambda <- tLambda + pi.g[g]*init.tLambda[,,g]
              Delta <- Delta + pi.g[g]*init.Delta[,,g]
            }
            for(g in 1:G){
              Psi[,,g] <- Delta*init.w[g]
              w[g] <- init.w[g]
            }
            if(q==1){
              Lambda <- c(Lambda)
              tLambda <- c(tLambda)
              pd <- Lambda%o%tLambda
            }
            else{ pd <- Lambda%*%tLambda}
            for(g in 1:G){
              inv.sig[,,g] <- solve(pd+Psi[,,g])
            }
          }
          sum1 <- matrix(0,p,q)
          sum2 <- matrix(0,q,q)
          for(g in 1:G){
            dc.parm1 <- cc.lambda(p=p,q=q,n.g=n.g[g],w=w[g],Lambda=Lambda,tLambda=tLambda,S=S[,,g],inv.sig=inv.sig[,,g])
            sum1 <- sum1 + dc.parm1$sum1
            sum2 <- sum2 + dc.parm1$sum2
            Beta[,,g] <- dc.parm1$Beta
            Theta[,,g] <- dc.parm1$Theta
          }
          Lambda <- sum1%*%solve(sum2)
          tLambda <- t(Lambda)
          if(model=="CCUC"){
            for(g in 1:G){
              w[g] <- ccuc.psi(p=p,Lambda=Lambda,tLambda=tLambda,S=S[,,g],Beta=Beta[,,g],Theta=Theta[,,g])
              Psi[,,g] <- w[g]*diag(p)
            }
          }
          if(model=="CCUU"){
            dg.sum <- matrix(0,p,p)
            for(g in 1:G){
              dc.parm2 <- not.ccuc.psi(p=p,Delta=Delta,Lambda=Lambda,tLambda=tLambda,S=S[,,g],Beta=Beta[,,g],Theta=Theta[,,g])
              w[g] <- dc.parm2$w
              dg.sum <- dg.sum + (n.g[g]/w[g])*dc.parm2$dg 
            }
            k <- prod(diag(dg.sum)^(1/p))
            Delta <- (1/k)*(dg.sum)
            for(g in 1:G){
              Psi[,,g] <- w[g]*Delta
            }
          }
        }
        if(model=="CUCU"||model=="CUUU"){
          if(i==2){
            Lambda <- matrix(0,p,q)
            tLambda <- matrix(0,q,p)
            for(g in 1:G){
              Lambda <- Lambda + pi.g[g]*init.Lambda[,,g]
              tLambda <- tLambda + pi.g[g]*init.tLambda[,,g]
              Psi[,,g] <- init.Psi[,,g]
              w[g] <- init.w[g]
              Delta[,,g] <- init.Delta[,,g]
            }
            if(q==1){
              Lambda <- c(Lambda)
              tLambda <- c(tLambda)
              pd <- Lambda%o%tLambda
            }
            else{ pd <- Lambda%*%tLambda}
            for(g in 1:G){
              inv.sig[,,g] <- solve(pd+Psi[,,g])
            }
          }
          if(model=="CUCU"){
            R <- matrix(0,p,q)
            sum1 <- array(0,dim=c(q,q,p))
            for(g in 1:G){
              dc.parm <- cu.lambda(p=p,q=q,n.g=n.g[g],Delta=Delta[,,g],Lambda=Lambda,tLambda=tLambda,S=S[,,g],inv.sig=inv.sig[,,g])
              R <- R + dc.parm$R
              for(j in 1:p){
                sum1[,,j] <- sum1[,,j] + dc.parm$sum1[,,j]
              }
              Beta[,,g] <-dc.parm$Beta
              Theta[,,g] <- dc.parm$Theta
            }
            for(j in 1:p){
              if(q > 1){Lambda[j,] <- R[j,]%*%solve(sum1[,,j])}
              if(q==1){Lambda[j] <- R[j,]%*%solve(sum1[,,j])}
            }
            tLambda <- t(Lambda)
            w.sum <- 0
            for(g in 1:G){
              dc.parm2 <- not.ccuc.psi(p=p,Delta=Delta[,,g],Lambda=Lambda,tLambda=tLambda,S=S[,,g],Beta=Beta[,,g],Theta=Theta[,,g])
              w[g] <- dc.parm2$w
              w.sum <- w.sum + pi.g[g]*w[g]
              k[g] <- prod(diag(dc.parm2$dg)^(1/p))
              Delta[,,g] <- (1/k[g])*dc.parm2$dg
            }
            for(g in 1:G){
              #k[g] <- (n.g[g]/w.sum)*k[g]
              #Delta[,,g] <- (n.g[g]/(w.sum*k[g]))*dg[,,g]
              Psi[,,g] <- w.sum*Delta[,,g]
            }
          }
          if(model=="CUUU"){
            R <- matrix(0,p,q)
            sum1 <- array(0,dim=c(q,q,p))
            for(g in 1:G){
              dc.parm <- cu.lambda(p=p,q=q,n.g=n.g[g],Delta=Delta[,,g],Lambda=Lambda,tLambda=tLambda,S=S[,,g],inv.sig=inv.sig[,,g])
              R <- R + dc.parm$R
              for(j in 1:p){
                sum1[,,j] <- sum1[,,j] + dc.parm$sum1[,,j]
              }
              Beta[,,g] <-dc.parm$Beta
              Theta[,,g] <- dc.parm$Theta
            }
            for(j in 1:p){
              if(q > 1){Lambda[j,] <- R[j,]%*%solve(sum1[,,j])}
              if(q==1){Lambda[j] <- R[j,]%*%solve(sum1[,,j])}
            }
            tLambda <- t(Lambda)
            for(g in 1:G){
              dc.parm2 <- not.ccuc.psi(p=p,Delta=Delta[,,g],Lambda=Lambda,tLambda=tLambda,S=S[,,g],Beta=Beta[,,g],Theta=Theta[,,g])
              w[g] <- dc.parm2$w
              k[g] <- prod(diag(dc.parm$dg))^(1/p)
              Delta[,,g] <- (1/k[g])*dc.parm2$dg
              Psi[,,g] <- dc.parm2$dg
            }
          }
        }
        ## Inverse Sigma via Woodbury ##
        for(g in 1:G){
          if(model=="UCUC"||model=="UCUU"||model=="UUCU"||model=="UUUU"){
            wood <- woodbury(q=q,Lambda=Lambda[,,g],tLambda=tLambda[,,g],Psi=Psi[,,g])
          }
          if(model=="UCCC"||model=="UCCU"){
            wood <- woodbury(q=q,Lambda=Lambda[,,g],tLambda=tLambda[,,g],Psi=Psi)
          }
          if(model=="CCUC"||model=="CCUU"||model=="CUCU"||model=="CUUU"){
            wood <- woodbury(q=q,Lambda=Lambda,tLambda=tLambda,Psi=Psi[,,g])
          }
          inv.sig[,,g] <- wood$inv.sig
          det.sig[g] <- wood$det.sig
        }
      }
      ## E - Step ##
      for (g in 1:G){
        Good.dist[,g] <- rho[g]*dsal(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g],det.sig=det.sig[g])
        Good.dist[,g] <- sapply(Good.dist[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
        Bad.dist[,g] <- (1-rho[g])*dsal(p=p,x=x,mu=mu[,g],alpha=sqrt(eta[g])*alpha[,g],inv.sig=(1/eta[g])*inv.sig[,,g],det.sig=(eta[g]^p)*det.sig[g])
        Bad.dist[,g] <- sapply(Bad.dist[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
        Group.dist[,g] <- Good.dist[,g]+Bad.dist[,g]
        Mix.dist[,g] <- pi.g[g]*Group.dist[,g]
      }
      ## Log-Likelihood ##
      loglik[i] <- sum(log(rowSums(Mix.dist)))
      if(i > 2){ if( loglik[i] < loglik[i-1] ) stop(paste("Log-Likelihood Decreased on Iteration ",i,sep="")) }
      a1[i] <- aitkens(lval=loglik,i=i)
      if(i > 2){conv.value <- conv.val(a1=a1,loglik=loglik,i=i,type=c.type)}
      ## Expectations ##
      for (g in 1:G){
        zig[g,] <- Mix.dist[,g]/rowSums(Mix.dist)
        vig[g,] <- Good.dist[,g]/Group.dist[,g]
        e1e2 <- e1e2(p=p,x=x,mu=mu[,g],alpha=alpha[,g],inv.sig=inv.sig[,,g])
        E1[g,] <- e1e2$E1
        E2[g,] <- e1e2$E2
        e1e2.tilde <- e1e2(p=p,x=x,mu=mu[,g],alpha=sqrt(eta[g])*alpha[,g],inv.sig=(1/eta[g])*inv.sig[,,g])
        E1.tilde[g,] <- e1e2.tilde$E1
        E2.tilde[g,] <- e1e2.tilde$E2
      }
      # Iteration Count
      if(i==max.it)break
      i <- i+1
    }
    for(g in 1:G){
      sigma[,,g] <- solve(inv.sig[,,g])
    }
  }
  # Labeling
  label <- sol.labels2(zig=zig,zig.star=zig.star,vig=vig,vig.star=vig.star,n=n,G=G)
  clust <- label$clust
  good <- label$good
  if(i==max.it) i <- i+1
  if(print.res==TRUE) cat("PCSALM - ",model,"Iterations:",i-1,"\n")
  em.out <- list();
  em.out$pi <- pi.g; 
  em.out$rho <- rho;
  em.out$mu <- mu; 
  em.out$alpha <- alpha;
  em.out$sigma <- sigma;
  em.out$inv.sig <- inv.sig;
  em.out$det.sig <- det.sig;
  em.out$Lambda <- Lambda;
  em.out$Psi <- Psi;
  em.out$eta <- eta;
  em.out$clust <- clust;
  em.out$Good <- good;
  em.out$Soln <- clust*good;
  select <- model.select(like=loglik,G=G,n=n,p=p,q=q,model=model,zig=zig,vig=vig,soln=good)
  em.out$loglik <- loglik[i-1]
  em.out$bic <- select$bic;
  em.out$aic <- select$aic;
  em.out$aic3 <- select$aic3;
  em.out$icl <- select$icl;
  em.out$icl.mod <- select$icl.mod;
  em.out$penalty <- abs(select$icl-select$icl.mod);
  em.out$iter <- i-1
  timer1$stop("event1")
  em.out$time <- getTimer(timer1)$timeElapsed
  return(em.out)
}