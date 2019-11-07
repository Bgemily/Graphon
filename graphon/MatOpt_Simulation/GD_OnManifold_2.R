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



# GD on Manifold ----------------------------------------------------------

FirstTime = 0
OptS = function(Z, X, Y){
  S = t(X) %*% Z %*% Y
  return(S)
}


F_manifold = function(Z, X, Y){
  # S0 = 1/(1-missfrac) * t(X) %*% M %*% Y ######## Improve?
  # S = t(X) %*% Z %*% Y
  return(1/2 * norm(Projector(Z-X%*%t(X) %*% Z %*% Y%*%t(Y)), 'f')^2)
}


OptSpace = function(M, X0, Y0, MaxIter){
  Xt = X0; Yt = Y0; Zt = M
  n1 = dim(M)[1]; n2 = dim(M)[2];
  for (i in 1:MaxIter) {
    St = OptS(Zt, Xt, Yt)
    gradX = (diag(1, n1)-Xt%*%t(Xt)) %*% (Xt%*%St%*%t(Yt)-Zt) %*% Yt %*% t(St)
    gradY = (diag(1, n2)-Yt%*%t(Yt)) %*% t((Xt%*%St%*%t(Yt))-Zt) %*% Xt %*% St
    SVD_gradX = rsvd(-gradX, r); U_gradX = SVD_gradX$u; V_gradX = SVD_gradX$v; d_gradX = SVD_gradX$d
    SVD_gradY = rsvd(-gradY, r); U_gradY = SVD_gradY$u; V_gradY = SVD_gradY$v; d_gradY = SVD_gradY$d
    F_XtYt = F_manifold(Zt, Xt, Yt) 
    grad_norm_sq = norm(gradX, 'f')^2 + norm(gradY, 'f')^2
    
    # lr = 0.001
    # X = Xt %*% V_gradX %*% diag(cos((d_gradX)*lr)) %*% t(V_gradX) + U_gradX %*% diag(sin((d_gradX)*lr)) %*% t(V_gradX)
    # Y = Yt %*% V_gradY %*% diag(cos((d_gradY)*lr)) %*% t(V_gradY) + U_gradY %*% diag(sin((d_gradY)*lr)) %*% t(V_gradY)
    
    m = 1
    alpha = 1
    while (m<=12) {
      lr = alpha / 2^(m-1)
      X = Xt %*% V_gradX %*% diag(cos((d_gradX)*lr)) %*% t(V_gradX) + U_gradX %*% diag(sin((d_gradX)*lr)) %*% t(V_gradX)
      Y = Yt %*% V_gradY %*% diag(cos((d_gradY)*lr)) %*% t(V_gradY) + U_gradY %*% diag(sin((d_gradY)*lr)) %*% t(V_gradY)
      if ((F_manifold(Zt, X, Y) - F_XtYt) < (-0.5*lr*(grad_norm_sq))) {
        break
      }
      m = m+1
    }
    
    
    Xt = X; Yt = Y
    Zt = Xt%*%St%*%t(Yt) + Projector(M-Xt%*%St%*%t(Yt))
  }
  return(list(Xt = Xt, Yt = Yt, Zt = Zt))
}

library(rsvd)
Init_OptSpace = function(M){
  n1 = dim(M)[1]; n2 = dim(M)[2]
  miss_position = matrix(0, n1, n2); miss_position[imiss] = 1
  row_degree = rowSums(miss_position)
  col_degree = colSums(miss_position)
  for (j in 1:n2) {
    if (col_degree[j]>2*length(imiss)/n2) M[,j] = 0
  }
  for (i in 1:n1) {
    if (row_degree[i]>2*length(imiss)/n1) {
      M[i,] = 0
    }
  }
  rsvd_TrM = rsvd(1/(1-missfrac)*M, r)
  return(list(X0 = rsvd_TrM$u, Y0 = rsvd_TrM$v, S0 = diag(rsvd_TrM$d)))
}

init_optspace = Init_OptSpace(M_obs)
X0 = init_optspace$X0; Y0 = init_optspace$Y0; S0 = init_optspace$S0
GD_OnManifold = OptSpace(M_obs, X0, Y0, 3000)
X_est = GD_OnManifold$Xt; Y_est = GD_OnManifold$Yt
# S_est = 1/(1-missfrac) * t(X_est) %*% M_obs %*% (Y_est)
# S_est = OptS(M_obs, X_est, Y_est)
M_est = GD_OnManifold$Zt


norm(Projector(M_est-M_obs), 'f')/norm(M_obs, 'f')
norm(M_est-M_true, 'f')/norm(M_true, 'f')

