# TODO- rewrite CE algorithm
#' Search over all possible designs to find the global optimal
#'
#' Evaluates all possible designs by searching over a grid to find a global optimal. Can be slow for high A>20 or high p>2
#' @param A An adjacency matrix
#' @param p The number of treatments in the experiment
#' @param isomorphisms Whether to ignore designs that are isomorphic to designs evaluated earlier
#' @param ignoreLastNode Whether we set the last node as zero. Can be useful to satisfy constraints in double blocking structures
#' @param algorithm Which algorithm to use for finding which designs to evaluate; currently 'sequential', 'random' and 'ce' are implemented
#' @param specialNodes Are there any nodes which are fixed that we don't want to apply treatments to.
#' @param networkEffects Should the optimal design for the network effects be found (default is to find direct effects)
#' @keywords For a given adjacency matrix, find the optimal design with p treatments.
#' @export
#'
gridSearch<-function(A,p,isomorphisms=FALSE, blockList=NULL, ignoreLastNode=FALSE, algorithm="sequential",networkEffects=FALSE){
  ### Initiate parameters
  n<-dim(A)[1] # How many experimental units
  b<-length(blockList) # How many blocks
  AOptVal<-Inf #
  AOptDesign<-NULL
  numEval<-0 # How many designs we have evaluated

  if (b>0) {specialNodes<-c((p+1):(p+b)) # For compatibility- delete later
  }
  else
  {specialNodes<-NULL}

  ## Make a copy of A, and expand, then add block structure to network
  B<-matrix(0,ncol=(n+b),nrow=(n+b))
  B[1:n,1:n]<-A
  if (b>0){ #Assign nodes linked to blocks
  for (i in 1:b){
    B[blockList[[i]],n+i]<-1
    B[n+i,blockList[[i]]]<-1
  }
  }
  A<-B

  #### Set parameters for coordinate exchange algorithm
 # numTrialsCE<-0
  lastCEWin<-0
  numCERandStarts<-100 # How many random starts we do: TODO- put in parameters

  ### Work out isomorphisms of the network if we are using them
  z<- NULL
  if (isomorphisms==TRUE){
    ex1Igraph<-igraph::graph_from_adjacency_matrix(A)
    z<-isomorphisms(ex1Igraph,ex1Igraph) # Can this be faster in some other routine
    numIso<-length(z)
  }
  else{
    z[[1]]<-rep(1:n)
    numIso<-1
  }
  zReduced<-lapply(z, function(z) utils::head(z, n))

  ### For use in testing if a design is valid
  mult<-NULL
  for (i in (1:(n))){
    mult[i]<-p^((n)-i) #Todo: change for really big p
  }

# Thest whether a design is valid-
  testValidDesign<-function(td){
    x<-(match(c(1:p),td))
    if (anyNA(x)){return(FALSE)}
    if ((all(sort(x)==x))==FALSE){return(FALSE)} # Check whether design has elements 1:p with first
    # occurrence of k occurring before first occurrence of k+1
    tdm<-sum(td*mult)

    for (i in (1:numIso)){
      if(tdm>sum(td[zReduced[[i]]]*mult)){return(FALSE)}  # Check whether design is  first in lexicographic order
    }
    return(TRUE)
  }


  # Pick a starting design. Initiate the design so that 1s on all real nodes, and each block
  # gets a special treatment
  if (b>0){testDesign<-c(rep(1,n),c((p+1):(p+b)))} else {testDesign<-c(rep(1,n))}

  # For the CE algorithm, pck a random starting design instead
  if (algorithm=="CE"){testDesign<-c(1,sample(((1:(n-1))%%p+1),size=(n-1),replace=FALSE))
  AOptTest<-testDesign
  if (b>0){testDesign<-c(testDesign,c((p+1):(p+b)))}
  }

  # For block designs, we generally set the treatment effect of the last treatment to zero, and
  # set all the direct effects of the blocks to zero
  # For doubly blocked designs, such as Latin Squares, we might want to also set some block effect
  # such as the last to zero
  if (ignoreLastNode==FALSE){
    setToZero<-c((p+1):(2*p+b+1))
  }
  else
  {
    setToZero<-c((p+1):(2*p+b+1),2*p+2*b+1)
  }

  # For non-blocked designs, set the last treatment to be zero for identifiability.
  if (is.null(blockList)){ # For non-blocked designs
    setToZero<-p+1
    optVars<-c(2:(p+1))
    if (networkEffects==TRUE){optVars=c((p+2):(2*p+1))} # Do we want treatment effects, or network effects.
  }

### We now need to rewrite optVars to take account of those we set to zero,
### i.e. if we want the second and 4th variable, and 3rd is set to zero,
### We want the second and 3rd of the non-zero variables
  binVars<-rep(0,2*p+2*b+1)
  binVars[optVars]<-1
  binVars<-binVars[-setToZero]
  newOptVars<-(1:length(binVars))[as.logical(binVars)]

  ### Now the main loop- keep going until we stop at some point
  stopCriterion<-FALSE
  while (stopCriterion==FALSE){
    if  (testValidDesign(testDesign[1:(n)])){
      testWeight<-matrix(rep(0,((n+b)*(p+b))),nrow=(n+b))
      for (j in 1:(n+b)){
        testWeight[j,testDesign[j]]<-1
      }
      infMatrix<-informationMatFull(A,testWeight)
      # Gives a standard information matrix with rows corresponding to information about mean
      # and then 2m parameters
      infMatrix<-infMatrix[-setToZero,-setToZero] # Inf matrix reduced just for the non-zero effects
      if ((determinant(infMatrix)$modulus)>-10){ # If the matrix is invertible

        invInfMatrix<-solve(infMatrix,tol=1e-21)

        ##ATrial<-mean(diag(invInfMatrix)[newOptVars])  # This is A-optimality criteria, not pairwise.

        ### A rather convoluted way to get the utility, but I can't think of anything faster...
        ATrial<-0
        for (count1 in (1:(length(optVars)-1))){
          for (count2 in (count1+1):(length(optVars))){
            cVec<-rep(0,2*p+2*b+1)
            cVec[optVars[count1]]<-1
            cVec[optVars[count2]]<- -1
            cVec<-cVec[-setToZero]
            ATrial<-ATrial+(cVec%*%invInfMatrix%*%(cVec))[1,1]
        }
        }
        ATrial<-ATrial/(p*(p-1)/2)


        if (ATrial<AOptVal){ # If we've found an optimal, store it and update the CE algorithm parameters
          AOptVal<-ATrial
          AOptDesign<-testWeight
          AOptTest<-testDesign
          lastCEWin<-numEval
          }
      }
    }# End evaluate if a valid design

  #  numTrialsCE<-numTrialsCE+1 # How many trials have we had - counter for CE algorithm.
numEval<-numEval+1

    ### Now pick the next design
    if ((algorithm=="sequential") | (algorithm=="random")){
     testDesign<-c(nextDesign(testDesign[1:n],p,algorithm=algorithm,par=numEval))
     }

     if (algorithm=="CE") {
        testDesign<-nextDesign(AOptTest[1:n],p,par=numEval)
    }

    # Add on again the special block nodes which are fixed
    if (b>0) {testDesign<-c(testDesign,c((p+1):(p+b)))}

    # Stop if we've reached a stop criterion
    if (length(testDesign)!=(n+b)){stopCriterion<-TRUE}
    if ((algorithm=="random") && (numEval == 100 )){stopCriterion<-TRUE}
    if ((algorithm=="CE") && ((numEval-lastCEWin)>(n)*(p-1))){
      if (numCERandStarts>0){
        testDesign<-c(1,nextDesign(testDesign[1:(n-1)],p,algorithm="random"),specialNodes)
        numCERandStarts<-numCERandStarts-1
      }
      else{
        stopCriterion<-TRUE
      }
    }
  }

  return(list("AOptVal"=AOptVal,"AOptDesign"=AOptDesign,"numEval"=numEval))
  ##TODO: Consider replacing with (as.vector(AOptDesign%*%c(1:n)[1:n]))
}

### A utility function for CE. TODO: Incorporate in nextDesign
nextCEDesign<-function(startDesign,p,CECount){
  if (CECount<=0){return(startDesign)}
  n<-length(startDesign)
  startDesign<-startDesign-rep(1,n) # Take from 1->n notation to 0->n-1 for nicer modular arithmetic
  varToChange<-(((CECount-1)%/%(p-1))%%n)+1  #
  amountToChange<-(CECount-1)%%(p-1)+1
  startDesign[varToChange]<-(startDesign[varToChange]+amountToChange)%%p

  return(startDesign+rep(1,n))
}
