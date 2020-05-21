#Data path and markers are specified
path <- "/media/dokada/HD-LBVU3/MDS/test2/"
dir.create(path)
library(flowCore)
data_path <- "/media/dokada/HD-LBVU3/MDS/nakamurasan/Tcell_gated"
markers <- c("AmCyan-A","PerCP-Cy5-5-A","FITC-A","APC-A","PE-A","Pacific Blue-A")
genes <- c("CD4","CD8","CD45RA","CD45RO","CD25","CCR7")
files <- sort(list.files(data_path))[1:20]



########################
#STEP1:Decide the initial values of EM with one person
library(mixtools)
set.seed(1000)
file1 <- files[1]  
fcs <- read.FCS(paste0(data_path,"/",file1),transformation=FALSE)
expr <- fcs@exprs[,markers]
mu.ini.mat <- NULL         
lambda.ini.mat <- NULL　　　　　　
sigma.ini.mat <- NULL　　
for(i in 1:ncol(expr)){
  model <- normalmixEM(x=expr[,i], k=2,maxit=10000)
  mu.ini.mat <- cbind(mu.ini.mat,model$mu)
  lambda.ini.mat <- cbind(lambda.ini.mat,model$lambda)
  sigma.ini.mat <- cbind(sigma.ini.mat,model$sigma)
}
write.csv(mu.ini.mat,file=paste0(path,"mu_ini_mat.csv"),quote=F)
write.csv(lambda.ini.mat,file=paste0(path,"lambda_ini_mat.csv"),quote=F)
write.csv(sigma.ini.mat,file=paste0(path,"sigma_ini_mat.csv"),quote=F)


#############################################
#STEP2:Calculate cutoff
#the fnction calculating cutoff for one marker with EM
mycutoff.n2 <- function(x,mu.ini,sigma.ini,lambda.ini,fixed_range=NULL){
  model <- normalmixEM(x=x, k=2,maxit=100000,mu=mu.ini,sigma=sigma.ini,lambda=lambda.ini) 
  f2 <- function(x){　                       #lamda=The final mixing proportions
      abs(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) - model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2]))
  }
  bunpu1_idx <- which.min(model$mu)
  bunpu2_idx <- setdiff(c(1,2),bunpu1_idx)
  if(is.null(fixed_range)){
    grid <- seq(model$mu[bunpu1_idx] - 1*model$sigma[bunpu1_idx],model$mu[bunpu2_idx] + 1*model$sigma[bunpu2_idx],length=10000)
  }else{
    grid <- seq(fixed_range,model$mu[bunpu2_idx] + 1*model$sigma[bunpu2_idx],length=10000)
  }
  cutoff <- grid[which.min(f2(grid))]
  return(cutoff)
}


cell.percent1031 <- function(expr,mu.ini.mat,sigma.ini.mat,lambda.ini.mat,file_name,resampling=F){
  genes <- colnames(expr)
  if(resampling==F){
    cutoff.vec  <- rep(NA,ncol(expr))
    for(i in 1:ncol(expr)){
       fixed_range <- NULL
       cutoff.vec[i] <- mycutoff.n2(expr[,i],mu.ini=mu.ini.mat[,i],sigma.ini=sigma.ini.mat[,i],lambda.ini=lambda.ini.mat[,i],fixed_range=NULL)
    }
  }else{
    cutoffs_tmp <- matrix(NA,iter,6)
    for(j in 1:iter){
      cells <- sample(nrow(expr),each_sample)
　　　 for(i in 1:ncol(expr)){
         cutoff.vec[i] <- mycutoff.n2(expr[,i],mu.ini=mu.ini.mat[,i],sigma.ini=sigma.ini.mat[,i],lambda.ini=lambda.ini.mat[,i],image=F)
      }
    　cutoffs_tmp[j,] <-cutoff.vec
    }
    cutoff.vec <- apply(cutoffs_tmp,2,median)
  }
  return(cutoff.vec)
}

#計算
set.seed(1000)
mu.ini.mat.t <- read.csv(paste0(path,"mu_ini_mat.csv"),header=T,row.names=1)
lambda.ini.mat.t <- read.csv(paste0(path,"lambda_ini_mat.csv"),header=T,row.names=1)
sigma.ini.mat.t <- read.csv(paste0(path,"sigma_ini_mat.csv"),header=T,row.names=1)
markers <- c("AmCyan-A","PerCP-Cy5-5-A","FITC-A","APC-A","PE-A","Pacific Blue-A")
genes <- c("CD4","CD8","CD45RA","CD45RO","CD25","CCR7")
kekka.cutoff <- NULL
for(i in 1:length(files)){
  file1 <- files[i]
  fcs <- read.FCS(paste0(data_path,"/",file1),transformation=FALSE)
  expr <- fcs@exprs[,markers]
  colnames(expr) <- genes
  result <- cell.percent1031(expr=expr,mu.ini.mat=mu.ini.mat.t,sigma.ini.mat=sigma.ini.mat.t,lambda.ini.mat=lambda.ini.mat.t)
  kekka.cutoff <- rbind(kekka.cutoff,result)
  cat(i,"\n")
}
colnames(kekka.cutoff) <- genes
write.csv(kekka.cutoff,file=paste0(path,"kct_tmp.csv"),quote=F) 

###STEP3:Correction of the cutoff value
library(MASS)
cutoff.norm.des <- function(y){
  cr <- c(NA,NA)  
　　fit <- fitdistr(y, densfun = "normal")
　　param <- fit$estimate
  cr[1] <- qnorm(0.01, mean=param["mean"], sd=param["sd"],lower.tail=T)
  cr[2] <- qnorm(0.01, mean=param["mean"], sd=param["sd"],lower.tail=F)
  result.list <- list(param=param, cr=cr)
  return(result.list)
}
kekka.cutoff <- read.csv(paste0(path,"kct_tmp.csv"),header=T)
kc　<- kekka.cutoff[,c("CD4","CD8","CD45RA","CD45RO","CD25","CCR7")]
param.mat <- t(apply(kc,2,function(x){cutoff.norm.des(x)$param})) 

#for  CD25, fit mixture normal distribution to cutoff distribution
library(mixtools)
model <- normalmixEM(x=kc[,"CD25"], k=2,maxit=100000) 
f2 <- function(x){　                       #lamda=The final mixing proportions
    abs(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) - model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2]))
}
bunpu1_idx <- which.min(model$mu)
bunpu2_idx <- setdiff(c(1,2),bunpu1_idx)
grid <- seq(model$mu[bunpu1_idx] - 1*model$sigma[bunpu1_idx],model$mu[bunpu2_idx] + 1*model$sigma[bunpu2_idx],length=10000)
param.mat["CD25",] <- c(model$mu[bunpu2_idx],model$sigma[bunpu2_idx])
cutoff.range <- t(apply(param.mat,1,function(x){c(qnorm(0.01, mean=x[1], sd=x[2],lower.tail=T),qnorm(0.01, mean=x[1], sd=x[2],lower.tail=F))}))
colnames(cutoff.range) <- c("lower","upper")
rownames(cutoff.range) <-genes


kct_comp <- kc   
outlier.number <- rep(NA,length(markers))
base.cutoff.vector <- rep(NA,length(markers))　
for(i in 1:length(markers)){
  tmp <- union(which(kc[,i] < cutoff.range[i,"lower"]),which(kc[,i] > cutoff.range[i,"upper"]))  
　　base.cutoff <- param.mat[i,"mean"]
  kct_comp[tmp,i] <- base.cutoff
  outlier.number[i] <- length(tmp)
  base.cutoff.vector[i] <- base.cutoff
}
write.csv(outlier.number,file=paste0(path,"outlier.number.csv"),quote=F)
write.csv(kct_comp,file=paste0(path,"kct_comp.csv"),quote=F)


###STEP4 : Gating and cell population matrix
cell.percent3b <- function(expr,cutoff.vec){

  cutoff.vec <- as.numeric(cutoff.vec)

  pn.mat <- sapply(1:6,function(i){ifelse(expr[,i] < cutoff.vec[i],0,1)})

  clcnt.name <- NULL
  clcnt <- NULL
  for(c1 in c(0,1)){
    for(c2 in c(0,1)){
      for(c3 in c(0,1)){
        for(c4 in c(0,1)){
          for(c5 in c(0,1)){
            for(c6 in c(0,1)){
              tmp <- intersect(which(pn.mat[,1]==c1),which(pn.mat[,2]==c2))
              tmp <- intersect(tmp,which(pn.mat[,3]==c3))
              tmp <- intersect(tmp,which(pn.mat[,4]==c4))
              tmp <- intersect(tmp,which(pn.mat[,5]==c5))
              tmp <- intersect(tmp,which(pn.mat[,6]==c6))
              clcnt <- c(clcnt,length(tmp))
              clcnt.name <- c(clcnt.name,paste0(c1,c2,c3,c4,c5,c6))
            }
          }
        }
      }
    }
  }
  clpro <- clcnt/nrow(pn.mat)
  names(clpro) <- clcnt.name
  return(clpro)
}


library(flowCore)
kc.comp <- read.csv(paste0(path,"kct_comp.csv"),header=1,row.names=1)
kekka.proportion <- NULL
for(i in 1:length(files)){
  file1 <- files[i]
  fcs <- read.FCS(paste0(data_path,"/",file1),transformation=FALSE)
  expr <- fcs@exprs[,markers]
  
  result <- cell.percent3b(expr,kc.comp[i,])
  kekka.proportion <- rbind(kekka.proportion,result)
  cat(i,"\n")
}
write.csv(kekka.proportion,file=paste0(path,"kpt.csv"),quote=F,row.names=F)
