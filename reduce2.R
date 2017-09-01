#rm(list=ls())
#setwd("~/Documents/Research/jul14")

# helper functions
toBits <- function (x, nBits = nodes) {tail(rev(as.numeric(intToBits(x))),nBits)}
bitsToInt<-function(x) {packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")}

# ints takes a binary vector and returns the base 10 integers corresponding to the non-zero elements (helper function for pset)
ints = function(a) {
  num_vec = c()
  g = length(a)
  for (i in 1:g) if (a[i] > 0) {num_vec = c(num_vec, 2^(g-i))}
  return(num_vec)
}

# pset takes a vector of size n (where n is the sum of vector components) and a response vector 
# (corresponding to the number of observations of each hyper edge) and increments the response vector Y 
# for each of the vector's subsets of size {1,2,...,(n-1)}. 
# For example, if we have 3 nodes and the input vector is <1,1,1>, then the response vector elements corresponding
# to <0,0,1>, <0,1,0>, <1,0,0>, <0,1,1>, <1,0,1> are incremented by 1. 

pset = function(vec, Y, X) {
  int_vec = ints(vec)
  g = sum(vec)
  a = c()
  for (i in 1:(2^g-2)) { # -2 excludes the case of all zeros and all ones
    a = c(a,int_vec %*% toBits(i, nBits=g))
  }
  Y[a] = Y[a] + Y[bitsToInt(vec)]
  X[a,] = sweep(X[a,], 2, -X[bitsToInt(vec),])
  return(list(Y=Y, X=X))
}

# reduce takes a response vector 'Y', a design matrix 'X', and a vector 'size' that should be in decreasing order 
# (i.e., size = c(5,4)) - this vector indicates the sets of hyperedges that you want to remove from the design matrix/ 
# response vector. For example, the vector <5,4> will remove edges of size 5 and 4, while incrementing the response 
# vector for all subsets of observed edges of size 5, 4. Nodes is the number of nodes in the design matrix. 
# The function returns a list containing a 'reduced' response vector and a 'reduced' design matrix

reduce2 = function(Y, X, size,nodes) { 
  Xo = X # original binary design matrix
  gamma = apply(X, 1, sum)
  for (i in 1:length(size)) {
    Xr = Xo[which(gamma==size[i]),]
    Xr = array(Xr, dim = c(length(Xr)/nodes, nodes))
    for(j in 1:(nrow(Xr))) {
      Y = pset(Xr[j,],Y, X)$Y
      X = pset(Xr[j,],Y, X)$X }
  }
  return(list(Y=Y[-which(gamma %in% size)], X=X[-which(gamma %in% size),]))
}