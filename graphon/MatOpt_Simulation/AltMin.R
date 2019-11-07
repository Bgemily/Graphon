# Generate Data -----------------------------------------------------------

# L: n1*r, R: n2*r, M: n1*n2
n1 = 10; r = 5; n2 = 10
missfrac = 0.3

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



# Initialization ----------------------------------------------------------

library(rsvd)
LR_decomp = function(M)
{
  rsvd = rsvd(M, r)
  D = diag(rsvd$d); U = rsvd$u; V = rsvd$v #### Note: assume r is given.
  L0 = U %*% sqrt(D)
  R0 = V %*% sqrt(D)
  return(list(L0=L0, R0=R0))
}



# AltMin ------------------------------------------------------------------

OptR = function(M, L, R0, lr, MaxIter){
  Rt = R0
  for (i in 1:MaxIter) {
    Rt = Rt - lr * t(Projector(L%*%t(Rt)-M)) %*% L
  }
  return(Rt)
}

OptL = function(M, L0, R, lr, MaxIter){
  Lt = L0
  for (i in 1:MaxIter) {
    Lt = Lt - lr * Projector(Lt%*%t(R)-M) %*% R
  }
  return(Lt)
}


AltMin = function(M, L0, R0, lr, MaxIter, subMaxIter){
  Lt = L0; Rt = R0
  for (t in 1:MaxIter) {
    Lt = OptL(M, Lt, Rt, lr, subMaxIter)
    Rt = OptR(M, Lt, Rt, lr, subMaxIter)
  }
  return(list(L_est = Lt, R_est = Rt))
}

stT = Sys.time()

LR_0 = LR_decomp(M_obs/(1-missfrac))
L0 = LR_0$L0; R0 = LR_0$R0 
AltMin_est = AltMin(M_obs, L0, R0, 0.06, 100, 100)
L_est = AltMin_est$L_est; R_est = AltMin_est$R_est
M_est = L_est %*% t(R_est)

endT = Sys.time()


norm(Projector(M_est-M_obs), 'f')/norm(M_obs, 'f')
norm(M_est-M_true, 'f')/norm(M_true, 'f')
endT - stT
