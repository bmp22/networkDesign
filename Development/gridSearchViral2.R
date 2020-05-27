#' Search over all possible designs to find the global optimal
#'
#' Evaluates all possible designs by searching over a grid to find a global optimal. Can be slow for high A>20 or high p>2
#' @param A An adjacency matrix
#' @param p The number of treatments in the experiment,
#' @param isomorphisms Whether to ignore designs that are isomorphic to designs evaluated earlier
#' @param ignoreLastNode Whether we set the last node as zero. Can be useful to satisfy constraints in double blocking structures
#' @param algorithm Which algorithm to use for finding which designs to evaluate; currently 'sequential', 'random' and 'ce' are implemented
#' @param specialNodes Are there any nodes which are fixed that we don't want to apply treatments to.
#' @keywords For a given adjacency matrix, find the optimal design with p treatemtnes
#' @export
#'
gridSearchViral<-function(A,p,thetaPrior,isomorphisms=FALSE,specialNodes=NULL, ignoreLastNode=FALSE, algorithm="sequential"){

  ### Initiate parameters
  n<-dim(A)[1]   # How many experimental units.
  b<-length(specialNodes)  # How many blocks
  AOptVal<-Inf
  AOptDesign<-NULL
  numEval<-0

  ### Parameters for CE algorithms
  trialCE<-Inf
  numTrialsCE<-0
  successfulCE<-1
  numCERandStarts<-100 # How many random starts we do: TODO- put in parameters

### Work out isomorphims if we are using them
  z<- NULL
  if (isomorphisms==TRUE){
    ex1Igraph<-igraph::graph_from_adjacency_matrix(A)
    z<-isomorphisms(ex1Igraph,ex1Igraph)
    numIso<-length(z)
  }
  else{
    z[[1]]<-rep(1:n)
    numIso<-1
  }
  zReduced<-lapply(z, function(z) utils::head(z, n-b))

  mult<-NULL
  for (i in (1:(n-b))){
    mult[i]<-p^((n-b)-i) #Todo: change for really big p
  }

  lex.cmp <- function(vec1, vec2) {
    index <- which.min(vec1 == vec2)  # find the first diff
    sign(vec1[index] - vec2[index])   # assumes numeric
  }


  testValidDesign<-function(td){

    x<-(match(c(1:p),td))
    if (anyNA(x)){return(FALSE)}
    if ((all(sort(x)==x))==FALSE){return(FALSE)} # Check whether design has elements 1:p with first
    # occurrence of k occurring before first occurrence of k+1
    tdm<-sum(td*mult)

    for (i in (1:numIso)){
      if(tdm>sum(td[zReduced[[i]]]*mult)){return(FALSE)}
    }
    return(TRUE)
  }


  testDesign<-c(rep(1,n-b),specialNodes)

  ### Which rows are we intersted in for the block design- treatment rows, and block rows, assuming last one to be zero for each
  #neworder<-c(1:(p),(p+b+p+1):(p+b+p+b))
  #if (algorithm=="CE"){testDesign<-c(1,nextDesign(testDesign[1:(n-b-1)],p,algorithm="random"))}
  if (algorithm=="CE"){testDesign<-c(1,sample(((1:(n-1))%%p+1),size=(n-1),replace=FALSE))} #Balances start

  if (ignoreLastNode==FALSE){
    neworder<-c(1:(p),(p+b+p+2):(p+b+p+b+1)) # Works for block designs
  }
  else
  {
    neworder<-c(1:(p),(p+b+p+2):(p+b+p+b)) # Works for LS designs
  }
  if (is.null(specialNodes)){
    neworder<-c(1:p,(p+2):(2*p+1))
  }
  #neworder<-c(1,2,3,13,14,16,17)# For 3x3 latin square
  #neworder<-c(1,2,3,4,17,18,19,21,22,23)# For 3x3 latin square
  #cat(neworder)
  stopCriterion<-FALSE
  while (stopCriterion==FALSE){

    # if  (testLexSmallestFunction(testDesign[1:(n-b)])){
    if  (testValidDesign(testDesign[1:(n-b)])){
      #  print(testDesign)
      # update GUI console
      # flush.console()
#      cat("\r",testDesign)
      testWeight<-matrix(rep(0,(n*(p+b))),nrow=n)
      # testWeight
      for (j in 1:(n)){
        testWeight[j,testDesign[j]]<-1
      }
      # for (j in (n+1):(n+b)){
      #   testWeight[j,specialNodes[j-n+b]]<-1
      # }
      # print(testWeight)
      #Test if we can find a new optimum

      ### New bit for simple viral model
     # thetaPrior<-c(rep(1,p+1),0.5)
      gammaPar<-thetaPrior[p+2]
      X<-testWeight
      n<-dim(X)[1]
      KInv<-solve((diag(n)-gammaPar*A))
      F<-cbind((KInv%*%cbind(rep(1,n),X)),KInv%*%A%*%KInv%*%cbind(rep(1,n),X)%*%thetaPrior[1:(p+1)])
      infMatrix<-(t(F)%*%F)
     # infMatrix<-informationMatFull(A,testWeight)
      # Gives a standard information matrix with rows corresponding to information about mean and then 2m parameters
      neworder=c(1:(p),(p+2))
     infMatrix<-infMatrix[neworder,neworder] # Inf matrix for main and block effects
#cat(infMatrix)
      if ((determinant(infMatrix)$modulus)>-10){
        # designEval$ATrial<- evaluateNetworkDesign(infMatrix,p)
        numEval<-numEval+1;
        invInfMatrix<-solve(infMatrix,tol=1e-21)
        #invInfMatrix<-invInfMatrix[neworder,neworder] # Inf matrix for main and block effects
        ATrial<-invInfMatrix[p+1,p+1]


        #        ATrialNet<-0

    #     for (count1 in (2:p)){
    #       for (count2 in (count1+1):(p+1)){
    #       cVec<-rep(0,p)
    #          cVec[count1]<-1
    #         if (count2<(p+1)){cVec[count2]<- -1}
    #         ATrial<-ATrial+(cVec%*%invInfMatrix%*%(cVec))[1,1]
    #
    #         #  ATrial<-sum(diag(invInfMatrix)[2:p])
    #  #        cVec2<-rep(0,2*p)
    #  #        cVec2[count1+p-1]<- 1
    #  #        cVec2[count2+p-1]<- -1
    #  #       ATrialNet<-ATrialNet+(cVec2%*%invInfMatrix%*%(cVec2))[1,1]
    #         #
    #       }
    #     }
    #
    #     ATrial<-ATrial/(p*(p-1)/2)
    # #    ATrialNet<-ATrialNet/(p*(p-1)/2)

        if (ATrial<AOptVal){
          AOptVal<-ATrial
          AOptDesign<-testWeight
          successfulCE<-trialCE
          numTrialsCE<-0}

     #    if (ATrialNet<AOptNetVal){
     #    AOptNetVal<-ATrialNet
      #   AOptNetDesign<-testWeight
      #   }

      }
    }# End evaluate if canonical ordering

    #cat(testDesign,"\n")
    numTrialsCE<-numTrialsCE+1
    if ((algorithm=="CE") && (numTrialsCE%%(p+1)==0)){
      testDesign<-c(nextDesign(testDesign[1:(n-b)],p,algorithm=algorithm,par=trialCE),specialNodes) # Put back to start of cycle
      trialCE<-(trialCE+1)
      if (trialCE>(n-b)){trialCE<-1}
    }# If we've cycled through all designs for ExpUnit, advance to next.

    testDesign<-c(nextDesign(testDesign[1:(n-b)],p,algorithm=algorithm,par=trialCE),specialNodes)

    if (length(testDesign)!=n){stopCriterion<-TRUE}
    if ((algorithm=="random") && (numEval == 100 )){stopCriterion<-TRUE}
    if ((algorithm=="CE") && (numTrialsCE>(n-b)*p)){
      if (numCERandStarts>0){
        testDesign<-c(1,nextDesign(testDesign[1:(n-b-1)],p,algorithm="random"),specialNodes)
        numCERandStarts<-numCERandStarts-1
      }
      else{
        stopCriterion<-TRUE
      }
    }
  }
#  cat("\n")
  return(list("AOptVal"=AOptVal,"AOptDesign"=AOptDesign))
         #"AOptNetVal"=AOptNetVal,"AOptNetDesign"=AOptDesign,"numEval"=numEval))

}
