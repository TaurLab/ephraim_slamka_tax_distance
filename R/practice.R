


# how to use get.pairwise, get.dist, view.dist (written in functions.R) ----------------------------

library(yingtools2)
library(tidyverse)
library(grid)
library(gridExtra) 
#load functions.R
source("R/functions.R")

dist1 <- dist(mtcars)

# function to view the distance matrix
g.dist1 <- view.dist(dist1)
g.dist1

# convert the distance matrix to a table of pairwise distances
pairwise <- get.pairwise(dist1)

# convert the pairwise table back to distance matrix
dist2 <- get.dist(pairwise)
g.dist2 <- view.dist(dist2)

# visually compare the distances side by side
grid.arrange(g.dist1,g.dist2)




# applying a function over a list -----------------------------------------
library(foreach)
library(purrr)

data <- 
  list(
    c(1,2,3,4),
    c(2,3,4,5,6,7,8,9),
    c(10,20,30,40)
  )

# apply a function to each vector
lapply(data,mean)
lapply(data,length)

# same as above but simplies the result to a vector
sapply(data,mean)
sapply(data,length)

# custom-written functions can also be used.
custom_fun1 <- function(x) {
  answer <- mean(x) + 11 - median(x)
  return(answer)
}

custom_fun2 <- function(x) {
  answer <- (sum(x) + 1) / length(x)
  return(answer)
}

lapply(data,custom_fun1)

lapply(data,custom_fun2)

lapply(mock.phy, custom_fun1)


<<<<<<< HEAD

#ignore what follows, just an attempt at learning, any mention of these functions was just an attempt at learning-----------
#this isn't finished, I'm just trying to figure out how all of this works.
athing<- function(list){sqrt(sum(list)) }
aathing<- function(list) {n<-length(list)
  pwr<- 1/n
  (sum(list))^pwr}

aabthing<-function(x) {n<-length(x)
  pwr<- 1/n
  this<-((sum(x))^pwr)
  this/n}

aacthing<-function(x) {n<-length(x)
pwr<- 1/n
(sum(x))^pwr}
geomean<- function(set){n<-length(set)
  pwr<-1/n
  (prod(set))^pwr}
dist.taxhorn.aabthing<- get.taxdist(phy, fn=aabthing, method="horn")
g.hc.aabthing<-view.hclust(dist.taxhorn.aabthing, title= "aabthing")
grid.draw(g.hc.aabthing)

#working here 
abthing<-function(x){
 n<-(sum(x))*(prod(x))
 m<-(n^(1/length(x)))/x
 mean(m)/length(m)
 #pwr<- 1/(length(n))
# (prod(n))^pwr  
}
acthing<-function(x){
  n<-(sum(x))*(prod(x))
  m<-(n^(1/length(x)))/x
  mean(m)/mean(x)
  #pwr<- 1/(length(n))
  # (prod(n))^pwr
}
#test these two ^^vv
geomean2<- function(x){
 
  n<-length(x)
pwr<- (1/n)
(prod(x))^pwr}

bthing<- function(x){
n<-x^x
m<-sum(n)
m/length(n)
}
bathing<-function(x){
  n<-x^(1/(x^2))
  m<-sum(n)
  m/length(n)
}
bbthing<-function(x){
  n<-x^(sin(x))
  m<-sum(n)
  m/length(n)
}
bcthing<-function(x){
  n<-x^tan(x)
  m<-sum(n)
  m/length(n)
}
#this is the area I'm working on
#note for future ephraim, if elses are doomed to fail in 
#nonspecific cases. try using a median or mean or something, idk.
#be clever future E. past E believes in you.
bdthing<-function(x){
  
  d<-mean(x)
  n<-(x)^x
  o<-(x)^(1/x)
  w<-mean(n)
  s<-mean(o)
  ifelse(x>d, return(w), return(s))
}
frog<-c(0,2,44,57,0.47,9)
mean(frog[-1])

nsk<-function(x){
  mean(x[-1])
}
wnsk1<-function(x){
  set<-x[-1]
  weights<-length(set):1
  sum(set*weights)/sum(weights)
}
wnsk2<-function(x){
  set<-x[-1]
  #weights<-length(set):1
  weights<-c(2,3,4,5,4,2,1)
  sum(set*weights)/sum(weights)
}
weighted.mean <- function(x) {
  weights <- length(x):1
  sum( x*weights ) / sum(weights)
}

data2<-
  list(
    c(0.001,0.001,0.001),
    c(0.0001,0.0002,0.0003),
    c(0.0001,0.001,0.01),
    c(0.1,0.5,1),
    c(0.00001,1,1),
    c(20,11,10)
    
  )

=======
>>>>>>> c1449da5594d94c609fb8bcd29ffff94a8006736


g.v.manhattan <- view_violin(dist.manhattan,title="manhattan")
g.v.bray <- view_violin(dist.bray,title="bray")
g.v.euclidean <- view_violin(dist.euclidean,title="euclidean")
g.v.horn <- view_violin(dist.horn,title="horn")
g.v.unifrac <- view_violin(dist.unifrac,title="unifrac")
g.v.wunifrac <- view_violin(dist.wunifrac,title="wunifrac")
g.v.taxhorn.mean <- view_violin(dist.taxhorn.mean,title="taxhorn.mean")
g.v.taxhorn.weightedmean <- view_violin(dist.taxhorn.weightedmean,title="taxhorn.weightedmean")
#g.v.taxhorn.nsk<-view_violin(dist.taxhorn.nsk,title="taxhorn.nsk")

#the one below is the final taxHorn
g.v.taxhorn.wnsk1<-view_violin(dist.taxhorn.wnsk1,title="taxhorn")

vlist2 <- list( 
  g.v.bray,
  #g.v.manhattan,
  g.v.euclidean,
  g.v.horn,
  g.v.unifrac,
  #g.v.wunifrac,
  #g.v.taxhorn.mean,
  #g.v.taxhorn.weightedmean,
  #g.v.taxhorn.nsk,
  g.v.taxhorn.wnsk1)


#do.call(grid.arrange,vlist2)

mg2 <- marrangeGrob(vlist2,ncol=3,nrow=2)
ggsave("plots/csefviolin_compare_groups.pdf",mg2,width=20,height=12)
shell.exec(normalizePath("plots/csefviolin_compare_groups.pdf"))


