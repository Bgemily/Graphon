# Generate Data -----------------------------------------------------------

# L: n1*r, R: n2*r, M: n1*n2
n1 = 10; r = 5; n2 = 10
missfrac = 0.1

set.seed(1983)
L_true = matrix(rnorm(n1*r), n1, r)
set.seed(831)
R_true = matrix(rnorm(n2*r), n2, r)
M_true = L_true %*% t(R_true)


Projector = function(M_complete, miss_id = imiss)
{
  M_miss = M_complete
  M_miss[miss_id] = 0
  return(M_miss)
}


set.seed(19)
imiss = sample(seq(n1*n2), n1*n2*missfrac, replace = F)
M_obs = Projector(M_true, imiss)


# Debug OptS() --------------------------------------------------------------

OptS = function(M, X, Y, S0, lr, subMaxIter){
  St = S0
  for (t in 1:subMaxIter) {
    St = St - lr * t(X) %*% Projector(X%*%St%*%t(Y)-M) %*% Y
  }
  return(St)
}

M = M_true
library(rsvd)
rsvd = rsvd(M, r) ## This step can make X and Y match (rule out the error caused by random sign)
X = rsvd$u; Y = rsvd$v; S_true = rsvd$d
S0 = matrix(rnorm(r*r), r, r)
OptS(M, X, Y, S0, 1, 100)
