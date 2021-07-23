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

informationMatViral<-function(A,X,thetaPrior,Kh){ # InfMat for viral parameter
  infMatrix<-informationMatFull(A,X)
    sigmaSq<-1
  n<-dim(X)[1]
  ph<-length(thetaPrior)-2
  XX<-cbind(rep(1,n),X, A%*%X)
  AKh<-A%*%Kh # To avoid calculating so many times
 infMatrix<-cbind(infMatrix,t(XX)%*%AKh%*%XX%*%thetaPrior[1:(ph+1)])
  viralInf<-t(AKh%*%XX%*%thetaPrior[1:(ph+1)])%*%AKh%*%XX%*%thetaPrior[1:(ph+1)]+(sigmaSq/2)*sum(diag((t(AKh)+AKh)*(t(AKh)+AKh))) # The AR parameter
#  infMatrix<-cbind(infMatrix,t(XX)%*%A%*%Kh%*%XX%*%thetaPrior[1:(ph+1)])
 # viralInf<-t(A%*%Kh%*%XX%*%thetaPrior[1:(ph+1)])%*%A%*%Kh%*%XX%*%thetaPrior[1:(ph+1)]+(sigmaSq/2)*sum(diag((t(Kh)%*%t(A)+A%*%Kh)*(t(Kh)%*%t(A)+A%*%Kh))) # The AR parameter
  infMatrix<-rbind(infMatrix,c(t(infMatrix[1:(ph+1),ph+2]),viralInf))
 }
