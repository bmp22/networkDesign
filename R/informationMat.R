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
A<-as.matrix(A)
X<-as.matrix(X)
#s<-dim(X);
#p<-s[2];
F<-cbind((rep(1,dim(X)[1])),X[,1:(dim(X)[2]-setToZero)],A%*%X)
infMatrix<-t(F)%*%F
}

informationMatFull<-function(A,X){ # For compatibility
infMatrix<-informationMat(A,X,setToZero=0)
}

