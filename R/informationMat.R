#' Calculate Information Matrix
#'
#' Calculates the information matrix for a given adjacency matrix and design
#' @param A An adjacency Matrix
#' @param X design matrix
#' @param setToZero we may insist that we delete the row corresponding to the last treatment effect
#' for identifiability, equivalent to setting that to be zero. setToZero=1 means we do this.
#' @keywords Calculate an information matrix
#' @export


informationMat<-function(A,X,setToZero=1){
#A<-as.matrix(A)
#X<-as.matrix(X)
#s<-dim(X);
#p<-s[2];
F<-cbind((rep(1,dim(X)[1])),X[,1:(dim(X)[2]-setToZero)],A%*%X)
infMatrix<-t(F)%*%F
}

informationMatFull<-function(A,X){ # For compatibility
infMatrix<-informationMat(A,X,setToZero=0)
}

informationMatViral<-function(A,X,thetaPrior,K){ # InfMat for viral parameter
  infMatrix<-informationMatFull(A,X)
  sigmaSq<-1
  n<-dim(X)[1]
  p<-length(thetaPrior)-2
  XX<-cbind(rep(1,n),X)
 infMatrix<-infMatrix[1:(p+1),1:(p+1)]
  infMatrix<-cbind(infMatrix,t(XX)%*%A%*%K%*%XX%*%thetaPrior[1:(p+1)])
 #cat(dim(t(X)%*%A%*%K%*%X%*%thetaPrior[2:(p+1)]))
  #infMatrix<-cbind(infMatrix,c(t(XX)%*%A%*%K%*%XX%*%thetaPrior[1:(p+1)],t(X)%*%A%*%A%*%K%*%X%*%thetaPrior[2:(p+1)])) #Almost works
  viralInf<-t(A%*%K%*%XX%*%thetaPrior[1:(p+1)])%*%A%*%K%*%XX%*%thetaPrior[1:(p+1)]+(sigmaSq/2)*sum(diag((t(K)%*%t(A)+A%*%K)*(t(K)%*%t(A)+A%*%K)))
  infMatrix<-rbind(infMatrix,c(t(infMatrix[1:(p+1),p+2]),viralInf))
 }
