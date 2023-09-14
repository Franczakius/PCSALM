rm(list=ls())

# Required functions ####
source('PCSALM | Required Functions.R')
source('PCSALM | AECM.R')

# AIS Data ####
## Data cleaning ####

data(ais,package="DAAG")
X <- as.matrix(ais[,1:11])
head(X)
label <- ais[,12]

## Model building ####

ais.output1 <- ais.init <- list();
for(G in 1:4){
  try(init <- ais.init[[G]] <- SAL.EM(X,G,max.it=250,print.res=F),silent=T)
  for(q in 1:5){
    print("Factors, Group")
    print(c(q,G))
    for(m in 1:12){
      aecm <- NA
      id <- c(q,G,m,1)
      h <- m + (q-1)*12 + (G-1)*60
      try(aecm <- PCSALM(X,q,G,m,init=init,print.res=F),silent=T)
      ais.output1[[h]] <- aecm.output(aecm,id)
    }
  }
}

## Model Selection / Assessment ####

ari <- bic <- icl <- micl <- conv <- NULL
for(h in 1:240){
  ari[h] <- ais.output1[[h]]$ari
  icl[h] <- ais.output1[[h]]$icl
  micl[h] <- ais.output1[[h]]$icl.mod
  bic[h] <- ais.output1[[h]]$bic
  conv[h] <- ais.output1[[h]]$converge  
}
which(icl==max(icl,na.rm=T))
which(bic==max(bic,na.rm=T))
which(micl==max(micl,na.rm=T))
which(ari==max(ari,na.rm=T))
AIS_Result <- ais.output1[[166]]
AIS_Result

# Italian Olive Oil Data ####
## Data ####

data(olive,package="pgmm")
X <- as.matrix(olive[,-c(1,2)])
label <- olive[,1]

## Model Building ####

olive.output1 <- olive.init <- list();
for(G in 1:4){
  try(init <- olive.init[[G]] <- SAL.EM(X,G,max.it=250,print.res=F),silent=T)
  for(q in 1:5){
    print("Factors, Group")
    print(c(q,G))
    for(m in 1:12){
      aecm <- NA
      id <- c(q,G,m,1)
      h <- m + (q-1)*12 + (G-1)*60
      try(aecm <- PCSALM(X,q,G,m,init=init,print.res=F),silent=T)
      olive.output1[[h]] <- aecm.output(aecm,id)
    }
  }
}

## Model Selection / Assessment ####

ari <- bic <- icl <- micl <- conv <- NULL
for(h in 1:240){
  ari[h] <- olive.output1[[h]]$ari
  icl[h] <- olive.output1[[h]]$icl
  micl[h] <- olive.output1[[h]]$icl.mod
  bic[h] <- olive.output1[[h]]$bic
  conv[h] <- olive.output1[[h]]$converge  
}
which(icl==max(icl,na.rm=T))
which(bic==max(bic,na.rm=T))
which(micl==max(micl,na.rm=T))
which(ari==max(ari,na.rm=T))
Olive_Result <- olive.output1[[179]]
Olive_Result
