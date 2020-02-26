#' Calculate the optimality criteria function for a given network design
#'
#' Calculate the optimality criteria (A_s optimality) for a given network design
#' @param infMatrix An information Matrix
#' @param p number of treatments
#' @export
#'

evaluateNetworkDesign<-function(infMatrix,p){

    ATrial<-0
    ATrialNet<-0
    invInfMatrix<-solve(infMatrix)#Needs to calculate twice so might as well store
    #TODO: Review this matrix algebra in R
    for (count1 in (2:p)){
      for (count2 in (count1+1):(p+1)){
        cVec<-rep(0,2*p)
        cVec[count1]<-1
        cVec[count2]<- -1
        cVec[p+1]<-0
       ATrial<-ATrial+(cVec%*%invInfMatrix%*%(cVec))[1,1]

        cVec2<-rep(0,2*p)
        cVec2[count1+p-1]<- 1
        cVec2[count2+p-1]<- -1
       ATrialNet<-ATrialNet+(cVec2%*%invInfMatrix%*%(cVec2))[1,1]
        }
    }

    ATrial<-ATrial/(p*(p-1)/2)

    ATrialNet<-ATrialNet/(p*(p-1)/2)
    return(list("ATrial"=ATrial,"ATrialNet"=ATrialNet))
}
