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


# Initialization ----------------------------------------------------------

library(r)
LR_decomp = function(M)
{
  r = r(M, r)
  D = diag(r$d); U = r$u; V = r$v #### Note: assume r is given.
  L0 = U %*% sqrt(D)
  R0 = V %*% sqrt(D)
  return(list(L0=L0, R0=R0))
}


# Projected GD ------------------------------------------------------------

IncohProj = function(X, X0, c, mu){
  n = dim(X)[1]; r  = dim(X)[2]
  for (i in 1:n) {
    if (norm(X[i,], '2') > sqrt(c*mu*r/n)*norm(X0, '2')) {
      X[i,] = X[i,] / norm(X[i,], '2') * sqrt(c*mu*r/n) * norm(X0, '2')
    }
  }
  return(X)
}


ProjectedGD = function(M, L0, R0, c, mu, lr, MaxIter)
{
  Lt = L0; Rt = R0;
  for (t in 1:MaxIter) {
    L = Lt - lr * Projector(Lt%*%t(Rt)-M) %*% Rt
    L = IncohProj(L, L0, c, mu)
    R = Rt - lr * t(Projector(Lt%*%t(Rt)-M)) %*% Lt
    R = IncohProj(R, R0, c, mu)
    Lt = L
    Rt = R
  }
  return(list(Lt = Lt, Rt = Rt))
}


c = 1
mu = 1
lr = 0.01
MaxIter = 10000


stT = Sys.time()

LR_0 = LR_decomp(M_obs/(1-missfrac))
L0 = LR_0$L0; R0 = LR_0$R0 

LR_est_proj = ProjectedGD(M_obs, L0, R0, c, mu, lr, MaxIter)
L_est_proj = LR_est_proj$Lt; R_est_proj = LR_est_proj$Rt
M_est = L_est_proj%*%t(R_est_proj)

endT = Sys.time()

norm(Projector(M_est-M_obs), 'f')/norm(M_obs, 'f')
norm((M_est-M_true), 'f')/norm(M_true, 'f')
endT - stT





