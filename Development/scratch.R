myAdj<-JRSSExamples$ex1
myAdj<-JRSSExamples$ex2
myAdj<-as_adjacency_matrix(make_lattice(c(4,4)))





##### Fisher information matrix when gamma variable
A<-ex3
n<-dim(A)[1]
p<-2
thetaPrior<-c(1,1,1,0.9)
gammaPar<-thetaPrior[p+2]
X<-gridSearchViral(A,m,gammaPar=gammaPar)$AOptDesign
KInv<-solve((diag(n)-gammaPar*A))
F<-cbind((KInv%*%cbind(rep(1,n),X)),KInv%*%A%*%KInv%*%cbind(rep(1,n),X)%*%thetaPrior[1:(m+1)])
F<-as.matrix(F)
(infMatrix<-(t(F)%*%F))
selRows<-c(1:(m),(m+2))
det(infMatrix[selRows,selRows])^(1/4)

#### If we just want to estimate gamma

for (i in 1:99){out[i]<-gridSearchViral(ex2,p=2,thetaPrior=c(0,1,1,i*0.01))$AOptVal}
plot(out)

### Generate examples for social network models for JRSS A papers

setwd("~/Dropbox/ViralNetworks/SocialNetworks")
networkList<-as.list(JRSSExamples)
 shortname<-"ex"
#or
#networkList<-list(as_adjacency_matrix(make_lattice(c(3,3))),as_adjacency_matrix(make_lattice(c(4,4))))
#shortname<-"lattice"
p<-2

for (j in (1:length(networkList))){
pdf(file=paste(shortname,j,".pdf", sep=""))
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (i in 1:9){
  myAdj<-networkList[[j]]
  diag(myAdj)<-0 # Some of the network diagonals were mistakenly non-zero
  myDes<-gridSearchViral(as.matrix(myAdj),p,thetaPrior = c(0,10*c(1:p),0.1*i),isomorphisms=TRUE)
    g2<-graph_from_adjacency_matrix(myAdj)
  #l <- layout_in_circle(g2)
  ##or
  l <- layout_on_grid(g2)
  E(g2)$arrow.mode<-0
  colrs<-c("red","blue","purple","orange")
  V(g2)$color<-colrs[myDes$AOptDesign%*%c(1:p)]
  plot.igraph(g2,edge.arrow.size=0.4, vertex.size=30,vertex.label.color="white",margin=0,layout=l,main=paste("p=",0.1*i))
}
dev.off()
}

### Generate random spatial fill
#library(ggplot2)
par(mfrow=c(1,1))
set.seed(2020)
n<-15
x<-runif(n)
y<-runif(n)
dat<-data.frame(x,y)
p<-2

d<-matrix(rep(0,n*n),nrow = n)
for (i in 1:(n-1)){
  for (j in ((i+1):n)){
    d[i,j]<-sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
    d[j,i]<-d[i,j]
  }
}

#plot(x,y)
myAdj<-d

desOut<-list()
for (i in 1:9){
myDes<-gridSearchViral(as.matrix(myAdj),p,thetaPrior = c(0,rep(1,p),0.1*i))
desOut[[i]]<-cbind(rho=i*0.1,dat,des=as.factor(myDes$AOptDesign%*%c(1:p)))
}

alldata <-do.call(rbind,desOut)
gr<-ggplot(alldata, aes(x=x, y=y, fill=des)) +
            geom_point(shape=21,size=5) +
            facet_wrap(~rho, ncol=3) +
  theme(panel.background = element_rect(fill = NA, color = "black", linetype="solid", size=2))
gr
## Saved as spatial.pdf

#### Redesign block functionality in gridSearch

gridSearchv2(matrix(rep(0,9*9),nrow = 9),p=3,blockList=list(c(1,2,3),c(4,5,6),c(7,8,9)))

##### Fisher information matrix when gamma variable

## USe example 3 from JRSS A
A<-JRSSExamples$ex3
diag(A)<-0 # No self loops!
n<-dim(A)[1]
p<-2
b<-0


thetaPrior<-c(1,2,3,0.4)
gammaPar<-thetaPrior[p+2]
testWeight<-matrix(rep(0,(n*(p+b))),nrow=n)

testDesign<-c(rep(1,5),rep(2,15))
for (j in 1:(n)){
  testWeight[j,testDesign[j]]<-1
}
X<-cbind(rep(1,n),testWeight)
K<-solve((diag(n)-gammaPar*t(A)))
infMatrix<-informationMatFull(A,testWeight)[1:(p+1),1:(p+1)]
infMatrix<-cbind(infMatrix,t(X)%*%A%*%K%*%X%*%thetaPrior[1:(p+1)])
sigmaSq<-1
viralInf<-t(A%*%K%*%X%*%thetaPrior[1:(p+1)])%*%A%*%K%*%X%*%thetaPrior[1:(p+1)]+(1/sigmaSq)*sum(diag((t(K)%*%t(A)+A%*%K)))
infMatrix<-rbind(infMatrix,c(t(infMatrix[1:(p+1),p+2]),viralInf))
(infMatrix)
det(infMatrix)
selRows<-c(1:(p),(p+2))
reducedInfMatrix<-infMatrix[selRows,selRows]
solve(reducedInfMatrix)

#### Now let's just right a function to calculate the information matrix for a viral prior
viralInf<-function(design,A,beta,sigmaSq,rho){
  n<-dim(A)[1]
  p<-length(beta)-1 #beta has beta_0,beta_p
  b<-0 # no blocks
  testWeight<-matrix(rep(0,(n*(p+b))),nrow=n)
  for (j in 1:(n)){
    testWeight[j,design[j]]<-1
  }
  X<-cbind(rep(1,n),testWeight)
  K<-solve((diag(n)-rho*t(A)))
  infMatrix<-informationMatFull(A,testWeight)[1:(p+1),1:(p+1)]
  infMatrix<-cbind(infMatrix,t(X)%*%A%*%K%*%X%*%beta)
  viralInf<-t(A%*%K%*%X%*%beta[1:(p+1)])%*%A%*%K%*%X%*%beta[1:(p+1)]+(1/sigmaSq)*sum(diag((t(K)%*%t(A)+A%*%K)))
  infMatrix<-rbind(infMatrix,c(t(infMatrix[1:(p+1),p+2]),viralInf))
  return(infMatrix)
}

A<-JRSSExamples$ex2
diag(A)<-0
# optimal designs for ex2 for rho =0.1,0.2,0.3,0.4
des1<-c(1,rep(2,9))
des2<-c(1,rep(2,8),1)
des3<-c(1,1,rep(2,8))
des4<-c(1,2,1,2,2,1,1,2,2,2)

num1<-0
outInf<-matrix(nrow=512,ncol=9)
#outData<-data.frame()
for (i in 1:9){
des<-c(rep(1,9),2)
for (j in 1:512){
foo<-viralInf(design=des,A=A,beta=c(1,2,3),sigmaSq=1,rho=i/10)
p<-2
selRows<-c(1:(p),(p+2))
outInf[j,i]<-solve(foo[selRows,selRows])[2,2]
num1[j]<-sum(des==1)
des<-nextDesign(des,k = 2,algorithm = "sequential")

}

#outData<-rbind(outData,c(des,outInf))
}

mr<-0
outInf<-cbind(outInf,rowSums(outInf),rowSums(outInf[,1:3])) # Last column is weighted average
# corresponding to vague prior
for (j in 1:11) {mr[j]<-which(outInf[,j]==min(outInf[,j]))}

ex3Vis<-data.frame(outInf)
colnames(ex3Vis)<-c(paste("p",c(1:9),sep=""),"bayes1","bayes2")
ex3Vis$col<-"1Non"
ex3Vis$col[mr]<-"2Opt"
ex3Vis$col[mr[c(10,11)]]<-"3Bayes"
ex3Vis$col<-as.factor(ex3Vis$col)
ex3Vis$num1<-as.factor(num1)
#install.packages("GGally")
#library(GGally)

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Dark2",direction=1) + scale_fill_brewer(palette="Dark2",direction=1)
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))

gp<-ggscatmat(
  ex3Vis, columns=c("p2","p4","p6","p8","bayes1","bayes2"), color="col"
  #lower = list(
   # continuous = "smooth",
   # combo = "facetdensity",
  # mapping = aes(colour=col)
  #  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  #)
)
gp

gh<-ggplot( ex3Vis, aes(x=num1,y=bayes1)) +
  geom_boxplot()+
  xlab("Number of Treatments Set to 1")+
  ylab("Optimal function value")

gh

###############
x<-NULL
for (i in 1:999){
x[i]<-solve(viralInf(design=c(1,1,1,1,1,1,1,1,1,2),A=A,beta=c(1,2,3),sigmaSq=1,rho=i/1000)[selRows,selRows])[3,3]
}
plot(c(1:999)/1000,x,xlab="rho",ylab="InfCrit")
