# Simulation of incomplete matrix low-rank decomposition.


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

library(rsvd)
LR_decomp = function(M)
{
  D = diag(rsvd(M, r)$d); U = rsvd(M, r)$u; V = rsvd(M, r)$v #### Note: assume r is given.
  L0 = U %*% sqrt(D)
  R0 = V %*% sqrt(D)
  return(list(L0=L0, R0=R0))
}


# Vanilla GD --------------------------------------------------------------

VanillaGD = function(M, L0, R0, lr, MaxIter)
{
  Lt = L0; Rt = R0;
  for (t in 1:MaxIter) {
    L = Lt - lr * Projector(Lt%*%t(Rt)-M) %*% Rt
    R = Rt - lr * t(Projector(Lt%*%t(Rt)-M)) %*% Lt
    Lt = L
    Rt = R
  }
  return(list(Lt = Lt, Rt = Rt))
}

L0 = LR_decomp(M_obs/(1-missfrac))$L0; R0 = LR_decomp(M_obs/(1-missfrac))$R0 
LR_est = VanillaGD(M_obs, L0, R0, 0.01, 10000)
L_est = LR_est$Lt; R_est = LR_est$Rt

norm(Projector(L_est%*%t(R_est)-M_obs), 'f')/norm(M_obs, 'f')
norm(L_est%*%t(R_est)-M_true, 'f')/norm(M_true, 'f')

# abs((L_est%*%t(R_est)-M_true)/M_true)

