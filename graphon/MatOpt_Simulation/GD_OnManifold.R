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



# GD on Manifold ----------------------------------------------------------

FirstTime = 0
OptS = function(M, X, Y, S0, lr, subMaxIter){
  St = S0
  for (t in 1:subMaxIter) {
    St = St - lr * t(X) %*% Projector(X%*%St%*%t(Y)-M) %*% Y
  }
  return(St)
}


F_manifold = function(M, X, Y, S0, lr, subMaxIter){
  # S0 = 1/(1-missfrac) * t(X) %*% M %*% Y ######## Improve?
  S = OptS(M, X, Y, S0, lr, subMaxIter)
  return(1/2 * norm(Projector(M-X%*%S%*%t(Y)), 'f')^2)
}


OptSpace = function(M, X0, Y0, S0, MaxIter, sublr, subMaxIter){
  Xt = X0; Yt = Y0; St = S0
  n1 = dim(M)[1]; n2 = dim(M)[2];
  for (i in 1:MaxIter) {
    # St = 1/(1-missfrac) * t(Xt) %*% M %*% (Yt) ###### Improve?
    St = OptS(M, Xt, Yt, St, sublr, subMaxIter)
    gradX = (diag(1, n1)-Xt%*%t(Xt)) %*% Projector(Xt%*%St%*%t(Yt)-M) %*% Yt %*% t(St)
    gradY = (diag(1, n2)-Yt%*%t(Yt)) %*% t(Projector(Xt%*%St%*%t(Yt))-M) %*% Xt %*% St
    SVD_gradX = rsvd(-gradX, r); U_gradX = SVD_gradX$u; V_gradX = SVD_gradX$v; d_gradX = SVD_gradX$d
    SVD_gradY = rsvd(-gradY, r); U_gradY = SVD_gradY$u; V_gradY = SVD_gradY$v; d_gradY = SVD_gradY$d
    F_XtYt = F_manifold(M, Xt, Yt, St, sublr, subMaxIter)
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
      if ((F_manifold(M, X, Y, St, lr, subMaxIter) - F_XtYt) < (-0.5*lr*(grad_norm_sq))) {
        break
      }
      m = m+1
    }
    Xt = X; Yt = Y
  }
  return(list(Xt = Xt, Yt = Yt))
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


stT = Sys.time()

init_optspace = Init_OptSpace(M_obs)
X0 = init_optspace$X0; Y0 = init_optspace$Y0; S0 = init_optspace$S0
GD_OnManifold = OptSpace(M_obs, X0, Y0, S0, 10000, 1, 50)
X_est = GD_OnManifold$Xt; Y_est = GD_OnManifold$Yt
# S_est = 1/(1-missfrac) * t(X_est) %*% M_obs %*% (Y_est)
S_est = OptS(M_obs, X_est, Y_est, S0, 1, 50)
M_est = X_est %*% S_est %*% t(Y_est)

endT = Sys.time()

norm(Projector(M_est-M_obs), 'f')/norm(M_obs, 'f')
norm(M_est-M_true, 'f')/norm(M_true, 'f')
endT - stT
