setwd("")

library(gRim)
library(numDeriv)
library(RBGL)

Rcpp::sourceCpp('code/IPS.cpp')
source('code/helper_functions.R')


source("http://www.biostat.umn.edu/~hodges/RPLMBook/Datasets/10_Slovenian_stomach_cancer/Slovenia_stomach_cancer_data.txt")


### Format the example Slovenia data
N           <- length(O)
A           <- -Q
diag(A)     <- 0
X           <- cbind(rep(1,N), SEc)
colnames(X) <- c("Int", "SEc")


### Fit the copCS model
res <- copCS(Y          = O, 
             X          = X, 
             offset     = log(E), 
             A          = A, 
             conf.level = 0.95, 
             conf.int   = 1,     #Produce a CI? 0==NO 
             boot.iter  = 500)


### Generate an inverse correlation from an adjacency matrix A and a 
### nearest-neighbor correlation rho
Q <- IPS(rho=0.513, A=A)
R <- solve(Q)

R[1:5,1:5] * A[1:5,1:5]
Q[1:5,1:5] * A[1:5,1:5]




