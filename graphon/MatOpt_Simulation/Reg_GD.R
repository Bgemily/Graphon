# Generate Data -----------------------------------------------------------

# L: n1*r, R: n2*r, M: n1*n2
n1 = 10; r = 5; n2 = 10
missfrac = 0.2

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


# Regularized GD ----------------------------------------------------------

RowScalSpecInit = function(M, r, beta1, beta2, p) ## beta1,2: target row norm bounds. p: observe probability.
{
  rsvd = rsvd(M/p, r)
  D = diag(rsvd$d); U = rsvd$u; V = rsvd$v #### Note: assume r is given.
  L0 = U %*% sqrt(D)
  R0 = V %*% sqrt(D)
  for (i in 1:dim(M)[1]) {
    if (norm(L0[i,], '2') > sqrt(2/3)*beta1) {
      L0[i,] = L0[i,]/norm(L0[i,], '2')*sqrt(2/3)*beta1
    }
  }
  for (i in 1:dim(M)[2]) {
    if (norm(R0[i,], '2') > sqrt(2/3)*beta2) {
      R0[i,] = R0[i,]/norm(R0[i,], '2')*sqrt(2/3)*beta2
    }
  }
  return(list(L0=L0, R0=R0))
}


G0Derivative = function(z) I(z>1)*2*(z-1)


ReguGD = function(M, L0, R0, rho, beta1, beta2, betaT, lr, MaxIter){
  Lt = L0; Rt = R0;
  for (i in 1:MaxIter) {
    L = Lt - lr * (Projector(Lt%*%t(Rt)-M)%*%Rt + 
                     rho*(3/beta1^2*Lt)*G0Derivative(3/(2*beta1^2)*apply(Lt, 1, function(x) norm(x, '2')^2)) + 
                     rho*(3/betaT^2*Lt)*G0Derivative(3/(2*betaT^2)*norm(Lt, 'f')^2))
    R = Rt - lr * (t(Projector(Lt%*%t(Rt)-M))%*%Lt + 
                     rho*(3/beta2^2*Rt)*G0Derivative(3/(2*beta2^2)*apply(Rt, 1, function(x) norm(x, '2')^2)) + 
                     rho*(3/betaT^2*Rt)*G0Derivative(3/(2*betaT^2)*norm(Rt, 'f')^2))
    Lt = L
    Rt = R
  }
  return(list(Lt=Lt, Rt=Rt))
}

Ct = 5
Cd = 6500
C1 = 0.4
C2 = 1
MaxSingVal = max(rsvd(M_obs/(1-missfrac), r)$d) ### = C1 * sqrt(norm(M_obs, 'f')^2/(1-missfrac)^r)
mu = C2 * sqrt(n1*n2) / (r*MaxSingVal) * max(abs(M_obs))
MinSingVal = min(rsvd(M_obs/(1-missfrac), r)$d)
kappa = MaxSingVal / MinSingVal
delta = MinSingVal/(Cd*r^1.5*kappa)
delta0 = delta / 6
betaT = sqrt(Ct*r*MaxSingVal)
beta1 = betaT * sqrt(3*mu*r/n1)
beta2 = betaT * sqrt(3*mu*r/n2)
# rho = 8 * (1-missfrac) * delta0^2
rho = 1


stT = Sys.time()

LR_init =  RowScalSpecInit(M_obs, r, beta1, beta2, 1-missfrac)
L0 = LR_init$L0; R0 = LR_init$R0

# L0 = LR_decomp(M_obs/(1-missfrac))$L0; R0 = LR_decomp(M_obs/(1-missfrac))$R0

LR_est_regu = ReguGD(M_obs, L0, R0, rho, beta1, beta2, betaT, 0.01, 10000)
L_est_reg = LR_est_regu$Lt; R_est_reg = LR_est_regu$Rt
M_est = L_est_reg%*%t(R_est_reg)

endT = Sys.time()

norm(Projector(M_est-M_obs), 'f')/norm(M_obs, 'f')
norm(M_est-M_true, 'f')/norm(M_true, 'f')
endT - stT




