


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







