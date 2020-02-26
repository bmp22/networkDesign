#' For a big network, apply treatments to a subset and search
#'
#' For a big network where we can apply treatments to a subset of them, find the optimal
#' design using a Fedorov exchange type algorithm.
#'
#' @param A adjacency Matrix of the network
#' @param treatments number of treatments to search over
#' @param sampleSize Number of nodes given a treatment, as opposed to null treatment
#' @param numIterations How many iterations to use in the exchange algorithm
#' @param restartProb The probability of starting from a random restart. Set to 1 to evaluate random
#' balanced designs.
#' @export
#'

bigNetworkSearch<-function(A,treatments,sampleSize,restartProb=0,numIterations){

if (dim(A)[1]!=dim(A)[2]){stop("Adjacency Matrix is not square")}
AOptVal<-Inf   #
AOptDesign<-NA
AOptNetVal<-Inf
AOptNetDesign<-NA
runningAOptVal<-NA
runningAOptNetVal<-NA

nextRandomDesign<-function(inSample,sampleDesign,labels,popSize){
  # Calculates the next valid design in a set of length n with
  # k possible labels
  n<-length(inSample);

  firstPoss<-sample(1:n,1); # Pick a random something in the design
  newLabel<-sample(2:labels,1); # pick a random label that is not 1;


  if (stats::runif(1)<0.1){ sampleDesign[firstPoss]<-newLabel}
  else {secondPoss<-inSample[1]
  while(secondPoss%in%inSample){
    secondPoss<-sample(1:popSize,1) # and a second random number in design
  }
  inSample[firstPoss]<-secondPoss
  }
  return(list("inSample"=inSample,"sampleDesign"=sampleDesign))
}



exchangeIterations<-numIterations; #How many iterations of algorithms to use.
p<-treatments+1; # first treatment is special treatment meaning not selected
n<-dim(A)[1]; # maximumInterval examined


a<-sample(1:n);
inSample<-a[1:sampleSize];	# Choose random elements for start
testDesign<-rep((2:p),sampleSize/treatments)

#while (length(testDesign)~=n-1)
for (myCount in 1:exchangeIterations){
  testWeight<-cbind(1,matrix(rep(0,n*(p-1)),nrow=n))
 for (j in 1:sampleSize){
  testWeight[inSample[j],testDesign[j]]<-1;
  testWeight[inSample[j],1]<-0;
}




  #Test if we can find a new optimum
  infMatrix<-informationMat(A,testWeight)
  if (det(infMatrix)>0.01){
    designEval<-evaluateNetworkDesign(infMatrix,p)
    ATrial<-designEval$ATrial
    ATrialNet<-designEval$ATrialNet

    if (ATrial<AOptVal){
      AOptVal<-ATrial
      AOptDesign<-testWeight
    }

    if (ATrialNet<AOptNetVal){
      AOptNetVal<-ATrialNet
      AOptNetDesign<-testWeight
    }

} # End evaluate if determinant non zero
oldDesign<-testDesign
oldSample<-inSample

nextDesign<-nextRandomDesign(inSample,testDesign,p,n)
testDesign<-nextDesign$sampleDesign
inSample<-nextDesign$inSample
if (stats::runif(1)<restartProb) {#Random restarts
  a<-sample(1:n)
  inSample<-a[1:sampleSize]	# Choose random elements for start
  #testDesign=2*ones(1,sampleSize); %Start with everything given first non special treatment
  testDesign<-rep((2:p),sampleSize/treatments)
}



runningAOptVal[myCount]<-AOptVal;
runningAOptNetVal[myCount]<-AOptNetVal;
}
return(list("AOptVal"=AOptVal,"AOptDesign"=AOptDesign,"AOptNetVal"=AOptNetVal,"AOptNetDesign"=AOptNetDesign,"runningAOptVal"=runningAOptVal,"runningAOptNetVal"=runningAOptNetVal))
}
