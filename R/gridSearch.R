#' Find the optimal design for a given network.
#'
#' Finds the optimal design for a given network. Various algorithms are implemented to search a network to find the optimal design
#' It can also be used to find optimal designs for experiments that contain blocking.
#' Can be slow for high A>20 or high p>2 .
#'
#' @param A An adjacency matrix
#' @param p The number of treatments in the experiment
#' @param isoSearch Whether to ignore designs that are isomorphic to designs evaluated earlier
#' @param ignoreLastNode Whether we set the last node as zero. Set to true to satisfy constraints in double blocking structures
#' @param algorithm Which algorithm to use for finding which designs to evaluate; currently 'sequential', 'random' , 'ce' and 'exchange' are implemented
#' @param blockList A list, each element of which is a collection of nodes which together form a block.
#' @param indirect  Set to FALSE is we do not model indirect effects. Default is TRUE, unless a weightPrior is set in which case is set to false
#' @param networkEffects Set to true if the optimal design for the network effects be found (default is to find direct effects)
#' @param viralOpt If we have a autoregressive model, do we want to estimate the autoregressive parameter(TRUE) or just estimate the difference in treatment effects
#' @param weightPrior A prior distribution on the parameters. Each row corresponds to a vector of unknown parameters.
#' @param weights Weights given to each element of weightPrior
#' @keywords For a given network with adjacency matrix A, find the optimal design with p treatments.
#' @export
#'
gridSearch<-function(A,p,isoSearch=FALSE, blockList=NULL, ignoreLastNode=FALSE, algorithm="sequential",indirect=NULL, networkEffects=FALSE, weightPrior=NULL, viralOpt=FALSE,weights=1){

  ### Initiate parameters
  if (!is.matrix(A)){stop("A must be a matrix")}
  n<-dim(A)[1] # How many experimental units
  b<-length(blockList) # How many blocks
  AOptVal<-Inf #
  AOptDesign<-NULL
  numEval<-0 # How many designs we have evaluated
  numFunEval <- 0 # How many times we evaluated the information matrix

  ### If we have block design, add nodes to the network to represent effect of blocks
  if (b>0) {
    specialNodes<-c((p+1):(p+b))
    } else   {specialNodes<-NULL}
  # Make a copy of A, and expand, then add block structure to network
  B<-matrix(0,ncol=(n+b),nrow=(n+b))
  B[1:n,1:n]<-A
  if (b>0){ #Assign nodes linked to blocks
  for (i in 1:b){
    B[blockList[[i]],n+i]<-1
    B[n+i,blockList[[i]]]<-1
  }
  }
  A<-B


  ### Do we want network effects in the model
  if (is.null(indirect)){
    if (is.null(weightPrior)){indirect<-TRUE} else {indirect<-FALSE}
  }else{
    if ((indirect==FALSE)&(networkEffects==TRUE)){stop("To estimate network effects must include them in the model")}
    }

  #### Set parameters for coordinate exchange algorithm
  lastCEWin<-0
  numCERandStarts<-10000 # How many random starts we do: TODO- put in parameters


  ### Work out isomorphisms of the network if we are using them
  z<- NULL
  if (isoSearch==TRUE){
    ex1Igraph<-igraph::graph_from_adjacency_matrix(A)
    z<-igraph::isomorphisms(ex1Igraph,ex1Igraph) # Can this be faster in some other routine
    numIso<-length(z)
  }
  else{
    z[[1]]<-rep(1:(n+b))
    numIso<-1
  }
  zUnlist<-matrix(unlist(z),ncol=(n+b),byrow = TRUE)



 ### Test whether a design is valid- whether there is symmetry in labels or symmetry in networks
  mult2<-(2^(c((n+b):1)))
 testValidDesign<-function(td){
  if (!is.null(weightPrior)){return(TRUE)} # if we have a prior on treatments, any design can be optimal. TODO- could be improved
  x<-(match(c(1:p),td))
  if (anyNA(x)){return(FALSE)}
  if (p==2){if (td[1]==2){return(FALSE)}}
    else{ if ((all(sort(x)==x))==FALSE){return(FALSE)} # Check whether design has elements 1:p with first
        }
    for (i in (1:numIso)){
    if (sum(sign(td[zUnlist[i,]]-td)*mult2)<0){return(FALSE)}
    }
  return(TRUE)
}

  ### Pick a starting design. Initiate the design so that 1s on all real nodes, and each block
  ### gets a special treatment
  if (b>0){testDesign<-c(rep(1,n),c((p+1):(p+b)))} else {testDesign<-c(rep(1,n))}

  # For the CE or exchange algorithm, pick a random starting design instead
  if ((algorithm=="CE")|(algorithm=="exchange")){testDesign<-c(1,sample(((1:(n-1))%%p+1),size=(n-1),replace=FALSE))
  AOptTest<-testDesign
  if (b>0){testDesign<-c(testDesign,c((p+1):(p+b)))}
  }
  oldDesign<-testDesign # For exchange algorithm

  ##### This section now defines what we want to estimate, and what we assume to be zero in the model.
  ##### TODO: Could probably be simplified
  optVars<-c(2:(p+1)) # By default estimate the treatment effects

  # For block designs, we generally set the treatment effect of the last treatment to zero, and
  # set all the direct effects of the blocks to zero
  # For doubly blocked designs, such as Latin Squares, we might want to also set some block effect
  # such as the last to zero
  if (ignoreLastNode==FALSE){ # For Blocked designs
    setToZero<-c((p+1):(2*p+b+1))
  }else{
    setToZero<-c((p+1):(2*p+b+1),2*p+2*b+1) # For doubly blocked designs, last block is zero
  }
  # For non-blocked designs, set the last treatment to be zero for identifiability.
  if (is.null(blockList)){ # For non-blocked designs
    setToZero<-p+1
  }

  if (indirect==FALSE){setToZero<-c(setToZero,(p+b+2):(p+b+p+1))} # If we don't want indirect effects in our model, set them to zero as well

  if (networkEffects==TRUE){optVars=c((p+b+2):(2*p+b+1))} # If we want to estimate treatment effects, or network effects.

  if (!is.null(weightPrior)){ # If we have an autoregressive component
    sigmaSq<-1 # Assume variance is 1. TODO: COuld be put in parameters
   if (indirect==FALSE){
     setToZero<-c((p+1),(p+b+2):(2*(p+b)+1))#ignore all network nodes plus last treatment
   }else{
     setToZero<-c((p+1))
     if (!is.null(blockList)){
     if (ignoreLastNode==FALSE){ # For Blocked designs
       setToZero<-c((p+1):(2*p+b+1))
     }else{
       setToZero<-c((p+1):(2*p+b+1),2*p+2*b+1) # For doubly blocked designs, last block is zero
     }
}
   }
    numParam<-2*p+2*b+2
    if (viralOpt==TRUE){ # Do we want to estimate the viral parameter itself...
      optVars<-c(2*p+2*b+2) # Override what we wanted to estimate
    }#else{
      #optVars<-c(2:(p+1))} #... or just the treatment effects
  } else{
    numParam<-2*p+2*b+1  }


### We now need to rewrite optVars to take account of those we set to zero,
### i.e. if we want the second and 4th variable, and 3rd is set to zero,
### We want the second and 3rd of the non-zero variables
  binVars<-rep(0,numParam)
  binVars[optVars]<-1
  binVars<-binVars[-setToZero]
  newOptVars<-(1:length(binVars))[as.logical(binVars)]


### Now the main loop- keep going until we stop at some stopCriterion
  stopCriterion<-FALSE
  while (stopCriterion==FALSE){
    if  (testValidDesign(testDesign[1:(n+b)])){
      ATrial<-0
      for (i in 1: length(weights)){
        APartTrial<-0
      testWeight<-matrix(rep(0,((n+b)*(p+b))),nrow=(n+b))
      for (j in 1:(n+b)){
        testWeight[j,testDesign[j]]<-1
      }
      if (is.null(weightPrior)){
        # if we do not have a prior, i.e. no autoregressive design,
        # calculate an information matrix without ar parameter
      infMatrix<-informationMatFull(A,testWeight)
      # Gives a standard information matrix with rows corresponding to
      #information about mean and then 2m parameters
      }else{ # calculate information matrix with a viral parameter
      if (length(weights)>1){thetaPrior<-weightPrior[i,]}else{thetaPrior<-weightPrior} # Must be a clever way to do this in R
        K<-solve((diag(n+b)-thetaPrior[2*p+2*b+2]*t(A)))
        infMatrix<-informationMatViral(A=A,X=testWeight,thetaPrior=thetaPrior,Kh=K) #
      }


      infMatrix<-infMatrix[-setToZero,-setToZero] # Inf matrix reduced just for the non-zero effects

      if ((determinant(infMatrix)$modulus)>-10){ # If the matrix is invertible
        ## Various different ways to get inverse of information matrix, but as it is positive definite the Cholesky Decomposition
        ## Cholesky decomposition seems about 10% faster.
        #invInfMatrix<-solve(infMatrix,tol=1e-21)
        #invInfMatrix<-pd.solve(infMatrix)
        invInfMatrix<-chol2inv(chol(infMatrix))

         ### A rather convoluted way to get the utility, but I can't think of anything faster...
        if (length(optVars)>1){
          for (count1 in (1:(length(optVars)-1))){
            for (count2 in (count1+1):(length(optVars))){
              #cVec<-rep(0,2*p+2*b+1)
              cVec<-rep(0,numParam)
              cVec[optVars[count1]]<-1
              cVec[optVars[count2]]<- -1
              cVec<-cVec[-setToZero]
              APartTrial<-APartTrial+(cVec%*%invInfMatrix%*%(cVec))[1,1]
            }
          }

          APartTrial<-APartTrial/(p*(p-1)/2)
        }else{
          APartTrial<-invInfMatrix[newOptVars,newOptVars]

        }
      }
      ATrial<-ATrial+weights[i]*APartTrial
    }


    if ((ATrial>0)&(ATrial<AOptVal)){ # If we've found an optimal, store it and update the CE algorithm parameters
      AOptVal<-ATrial
      AOptDesign<-testWeight
      AOptTest<-testDesign
      lastCEWin<-numEval
    }
    numFunEval<-numFunEval+1
  }# End evaluate if a valid design
  numEval<-numEval+1 # Counter for CE Algorithm

  ### Now pick the next design
   if ((algorithm=="sequential") | (algorithm=="random")){
       testDesign<-c(nextDesign(testDesign[1:n],p,algorithm=algorithm,par=numEval))
     }

     if (algorithm=="CE") {
        testDesign<-nextDesign(AOptTest[1:n],p,par=numEval)
    }

    if (algorithm=="exchange"){
      if ((numEval==(lastCEWin+1))|runif(1)<0.01){ #If we have found an improvement
      oldDesign<-testDesign # Update the design
      if (runif(1)<0.3){
        testDesign[sample(2:n,1)]<-sample(1:p,1)
      }else{
        a<-sample(2:n,2)
        testDesign[a]<-testDesign[rev(a)]
      }
      }else{
        testDesign<-oldDesign # revert to the previous design
        if (runif(1)<0.001){ # Start again with small probability to avoid getting stuck. TODO: Allow this parameter to change
        testDesign<-sort(c(1:p,nextDesign(testDesign[1:(n-p)], p ,algorithm="random"),specialNodes))
        }
        }
    #  cat(numEval," ",testDesign,"\n")
    }
    # Add on again the special block nodes which are fixed
    if (b>0) {testDesign<-c(testDesign,c((p+1):(p+b)))}

    ### Stop if we've reached a stop criterion
    if (length(testDesign)!=(n+b)){stopCriterion<-TRUE}
    if (((algorithm=="random")|(algorithm=="exchange")) && (numEval == 10000 )){stopCriterion<-TRUE}
    if ((algorithm=="CE") && ((numEval-lastCEWin)>(n)*(p-1))){
      if (numCERandStarts>0){
        testDesign<-sort(c(1:p,nextDesign(testDesign[1:(n-p)], p ,algorithm="random"),specialNodes))
        numCERandStarts<-numCERandStarts-1
        lastCEWin<-numEval
      }
      else{ stopCriterion<-TRUE }
    }

  }

  return(list("AOptVal"=AOptVal,"AOptDesign"=(as.vector(AOptDesign%*%c(1:(p+b)))[1:n]),"numFunEval"=numFunEval))
}

