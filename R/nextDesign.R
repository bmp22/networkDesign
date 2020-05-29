#' Get the next design sequentially
#'
#' Calculates the next valid design in a set of length n with k possible labels
#' @param design The current design
#' @param k the number of treatments
#' @param algorithm Can be "sequential" in which case all designs evaluated, or "random", in which case random designs evaluated
#' @param par A counting variable for the CE algorithm
#' @export

nextDesign<-function(design,k,algorithm="sequential",par=NULL){
  n<-length(design)
  if (algorithm=="sequential"){
    if(n==1){
      design<-NULL
    }
    else
    {
      if (design[n]<min(max(design[1:n-1])+1,k))
      {
        design[n]<-design[n]+1;
      }
      else
      {
        design<-c(nextDesign(design[1:(n-1)],k),1) # Set the last digit to 1 and recurse

      }
    }
    return(design)
  }
  if (algorithm=="random"){
    nextDesign<-sample(1:k, size=n, replace=TRUE)
  }

  if (algorithm=="CE"){
    CECount<-par
    startDesign<-design
    if (CECount<=0){return(startDesign)}
    n<-length(startDesign)
    startDesign<-startDesign-rep(1,n) # Take from 1->n notation to 0->n-1 for nicer modular arithmetic
    varToChange<-(((CECount-1)%/%(p-1))%%n)+1  #
    amountToChange<-(CECount-1)%%(p-1)+1
    startDesign[varToChange]<-(startDesign[varToChange]+amountToChange)%%p

    return(startDesign+rep(1,n))
  }
  return(nextDesign)
}

