# Développement du packkage de l'approche STATIS qui est une méthode de classification multibloc

#------------- Les librarys necessaires

library(dplyr)
library(diceR)
library(MatrixCorrelation)
library(FactoMineR)
library(MASS)
library(cluster) 
library(NMF)
library(visdat)

##Définiton de la fonction trace et de la fonction de RV

Trace <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}

# Calcul de RV pour deux tableaux donnés
RV <- function(X1, X2, center = TRUE){
  
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(center){
    X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
    X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
  }
  #if(cr){
  #X_1<- scale(X_1,center=T,scale=T)#centrage et réduction des données d'entrée
  #X_2<- scale(X_2,center=T,scale=T)   
  #}
  AA <- tcrossprod(X1)
  BB <- tcrossprod(X2)
  
  RV <- Trace(AA%*%BB) / (Trace(AA%*%AA)*Trace(BB%*%BB))^0.5
  RV
}


#lon_desire <- length(res.imp_dep$res.imp)

#H_comp<- vector(mode = "list", length = lon_desire)

cluster_statis <- function(data_){
  modele <- kmeans(na.omit(data_),2)
  cluster15=modele$cluster
  cluster15=multi.split(cluster15, split.char = "/", mnames = NULL)
  return(cluster15)
}

Calcul_H <- function(X){
  for (i in seq_along(X$res.imp)) {
  
    H_comp[[i]]=cluster_statis(as.data.frame(X$res.imp[i]))
  
  }
  return(H_comp)
}


## Calcul de W obtenue par H*t(H)

calcul_W <- function(X){
  W1=data.matrix(X, 
                 rownames.force = NA)%*%t(data.matrix(X, rownames.force = NA))
  return(W1)
  
}

## Calcul des poids W_i

Calcul_Wi <- function(X){
  
  lon_desire <- length(X$res.imp)
  
  W_comp<- vector(mode = "list", length = lon_desire)
  for (i in seq_along(X$res.imp)) {
    W_comp[[i]] <- calcul_W(cluster_statis(as.data.frame(X$res.imp[i]))) 
  }
  return(W_comp)

}

## Calcul de RV entre les tableaux pour une liste de tableau donnée

calcul_RV <- function(X){
  
  lon_desire <- length(X$res.imp)#Nombre de tableaux imputés
  
  R_comp<- vector(mode = "list", length = lon_desire)
  
  R_ <- rep(x=NA, times=lon_desire)#initialisation de la liste
  
  for(i in seq_along(X$res.imp)){
    
    for(j in seq_along(X$res.imp)){#Calcul de RV(i,j)
      
      R_[j] <- RV(as.data.frame(X$res.imp[i]),as.data.frame(X$res.imp[j]))
      
    }
    R_comp[[i]]=R_#Stockage de R_ij
  }
  return(R_comp)
}

# Calcul des poids à partir de la matrice S

Calcul_poids <- function(X){
  S <- as.matrix(
    array(as.numeric(unlist(calcul_RV(X))),
          dim = c(length(calcul_RV(X)), length(calcul_RV(X))))
    )
  #Normalisation des poids
  SVD_ <- svd(S); poids <- SVD_$u[,1]/sum(SVD_$u[,1])
  return(poids)
}


calcul_Wtilde <- function(X){
  
  lon_desire <- length(X$res.imp)
  
  W_comp_pon<- vector(mode = "list", length = lon_desire)
  poids <- Calcul_poids(X)
  for (i in seq_along(X$res.imp)) {
    W_comp_pon[[i]]=poids[i]*as.matrix(Calcul_Wi(X))[[i]]
    }
  #Calcul de W_tild=som_prod des H et poids
  #W_comp_pon <- as.matrix(W_comp_pon)
W_tild <- Reduce('+', as.list(as.matrix(W_comp_pon)))
return(W_tild)
}

# Calcul de W non pondéré

W_ <- function(X){
  
  lon_desire <- length(X$res.imp)
  
  W_comp_pon<- vector(mode = "list", length = lon_desire)
  for (i in seq_along(X$res.imp)) {
    W_comp_pon[[i]]=1/length(lon_desire)*as.matrix(Calcul_Wi(X))[[i]]
  }
  #Calcul de W_tild=som_prod des H et poids
  #W_comp_pon <- as.matrix(W_comp_pon)
  W_tild <- Reduce('+', as.list(as.matrix(W_comp_pon)))
  return(W_tild)
}

#W_stat <- calcul_Wtilde(res.imp_dep)
# 
# #Résultats du consensus
# 
# km <- kmeans(W_tild, centers = 2, nstart=25)
# ss <- silhouette(km$cluster, dist(W_tild))
# c_sh_statis=mean(ss[, 3])
# c_s_statis=c_sh_statis

