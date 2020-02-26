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
    nextDesign<-design
    nextDesign[par]<-(design[par]+1)
    if(nextDesign[par]>k){nextDesign[par]<-1}

  }
  return(nextDesign)
}

