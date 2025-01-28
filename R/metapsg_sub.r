
##### Architecture
# inputs : list per grps / then per variable[dim1] x studies [dim2]
# -> beta_t [[1:G]] [1:ng,1,1:s]
# -> var_t [[1:G]] [1:ng,1:ng,1:s]
############################################

sourceCpp("arma_inv.cpp")

library(abind)

#############################################
# Sub Functions to dev
#############################################

# "Proximal" quadratic form
#############################################
meta_regul_solver_g <- function(b, m, x, rho=1, mrho_inv) { # applied on a group g
  # b = beta_t
  # m = V^-1
  # x = z-u
  # rho = augmented lagrangian
  # mrho_inv = (2 V^-1 + rho I)^-1
  K=dim(b)[2]
  result <- matrix(nrow=dim(b)[1],ncol=K)
  for (k in 1:K){
    #if (is.na(mrho_inv[,,k])){
    #  result[,k] <- solve(2*m[,,k] + diag(rho, dim(m[,,k])[1])) %*% (2*m[,,k] %*% b[,k] + rho * x[,k])
    #} else {
    result[,k] <- mrho_inv[,,k] %*% (2*m[,,k] %*% b[,k] + rho * x[,k])
    #}
  }
  return(result)
}

# Calculation of Inverse Variance Matrix
#############################################
f_inv_var <- function(covar, same_var=FALSE){
  K <- dim(covar)[3]
  n <- dim(covar)[1]
  result <- array(NA, dim=c(n,n,K))
  if (same_var==TRUE){
    tempo_val <- solve(covar[,,1])
    for (k in 1:K){result[,,k] <- tempo_val}
  } else {
    for (k in 1:K){result[,,k] <- solve(covar[,,k])}
  }
  return(result)
}

# C++ version
f_cpp_inv_var <- function(covar, same_var=FALSE){
  K <- dim(covar)[3]
  n <- dim(covar)[1]
  result <- array(NA, dim=c(n,n,K))
  if (same_var==TRUE){
    tempo_val <- solve(covar[,,1])
    for (k in 1:K){result[,,k] <- tempo_val}
  } else {
    for (k in 1:K){result[,,k] <- solve(covar[,,k])}
  }
  return(result)
}


# Calculation of the solver Matrix : S = (2 V^-1 + rho I)^-1
#############################################
f_mrho_inv <- function(covar_inv, rho, same_var=FALSE){
  K <- dim(covar_inv)[3]
  n <- dim(covar_inv)[1]
  result <- array(NA,dim=c(n,n,K))
  if (same_var==TRUE){
    tempo_val <- solve(2*covar_inv[,,1] + diag(rho,dim(covar_inv)[1]))
    for (k in 1:K){result[,,k] <- tempo_val}
    #result <- apply(result, MARGIN=3, function(x) x<-tempo_val, simplify=FALSE)
  } else {
    for (k in 1:K){result[,,k] <- solve(2*covar_inv[,,k] + diag(rho,dim(covar_inv)[1]))}
    #result <- apply(covar_inv, MARGIN=3, function(x) x <- solve(2*x + diag(rho,dim(x)[1])), simplify=FALSE)
  }
  #result <- do.call(abind::abind, c(result, along = 3))
  return(result)
}

# C++ version
f_cpp_mrho_inv <- function(covar_inv, rho, same_var=FALSE){
  K <- dim(covar_inv)[3]
  n <- dim(covar_inv)[1]
  result <- array(NA,dim=c(n,n,K))
  if (same_var==TRUE){
    tempo_val <- arma_inv(2*covar_inv[,,1] + diag(rho,dim(covar_inv)[1]))
    for (k in 1:K){result[,,k] <- tempo_val}
    #result <- apply(result, MARGIN=3, function(x) x<-tempo_val, simplify=FALSE)
  } else {
    for (k in 1:K){result[,,k] <- arma_inv(2*covar_inv[,,k] + diag(rho,dim(covar_inv)[1]))}
    #result <- apply(covar_inv, MARGIN=3, function(x) x <- arma_inv(2*x + diag(rho,dim(x)[1])), simplify=FALSE)
  }
  #result <- do.call(abind::abind, c(result, along = 3))
  return(result)
}

# Proximal operator function for L2,1 norm
#############################################
prox_L21 <- function(x, lambda, w_p){
  l21_norm_w <- x^2 |> apply(FUN=sum, MARGIN=1) |> sqrt()
  result <- pmax(1 - lambda * w_p/l21_norm_w, 0) * x
  return(result)
}

# Proximal operator function for G2,1 norm
#############################################
prox_G21 <- function(x, lambda, w_g) {
  result <- max(1 - lambda * w_g / sqrt(sum(x^2)), 0) * x
  return(result)
}

# Proximal operator function for L1 norm
#############################################
prox_L1 <- function(x, lambda, w_v){
  result <- pmax(1 - lambda * w_v/abs(x), 0) * x
  return(result)
}

# Proximal operator function for L1,2 norm
#############################################
prox_L12 <- function(x, lambda, w_s){
  l12_norm_w <- x^2 |> apply(FUN=sum, MARGIN=2) |> sqrt()
  result <- pmax(1 - lambda * w_s/l12_norm_w, 0) |>  {\(z) x * z}()
  return(result)
}

# To define weights: w_g & w_p
###################################

f_weight_p_per_group <- function(beta_t, adaptive=FALSE){ #1D groups for each variable = 1 variable x studies
  if (adaptive){ 
    w_p <- beta_t^2 |> apply(FUN=sum, MARGIN=1) |> sqrt() |> {\(x) 1/x}()
  } else{ # To put variables weights to 1
    w_p <- rep(sqrt(dim(beta_t)[2]),dim(beta_t)[1])
  }
  return(w_p)
}


f_weight_g_per_group <- function(beta_t, adaptive=FALSE){ #2D groups of all variables x studies
  if (adaptive){ 
    w_g <- beta_t |> norm(type="F") |> {\(x) 1/x}()
  } else{ # To put group weights to 1
    w_g <- sqrt(dim(beta_t)[1]*dim(beta_t)[2])
  }
  return(w_g)
}


f_weight_p_per_group_no_scale <- function(beta_t, adaptive=FALSE){ #1D groups for each variable = 1 variable x studies
  if (adaptive){ 
    w_p <- beta_t^2 |> apply(FUN=sum, MARGIN=1) |> sqrt() |> {\(x) 1/x}()
  } else{ # To put variables weights to 1
    w_p <- rep(1,dim(beta_t)[1])
  }
  return(w_p)
}


f_weight_g_per_group_no_scale <- function(beta_t, adaptive=FALSE){ #2D groups of all variables x studies
  if (adaptive){ 
    w_g <- beta_t |> norm(type="F") |> {\(x) 1/x}()
  } else{ # To put group weights to 1
    w_g <- 1
  }
  return(w_g)
}


# SSE
###################################
f_sse_s <- function(beta_h, beta_t){
  result <- (beta_h - beta_t)^2 |> sum()
  return(result)
}

f_sse_v_inv <- function(beta_h, beta_t, inv_var, n_s){
  result <- 0
  for (s in 1:length(n_s)){
    result <- result + t(beta_h[,s] - beta_t[,s]) %*% (inv_var[,,s]/n_s[s]) %*% (beta_h[,s] - beta_t[,s])
  }
  return(result)
}

f_sse_v_inv_us <- function(beta_h, beta_t, inv_var, n_s){
  result <- 0
  for (s in 1:length(n_s)){
    result <- result + t(beta_h[,s] - beta_t[,s]) %*% inv_var[,,s] %*% (beta_h[,s] - beta_t[,s])
  }
  return(result)
}


# MIC
###################################

#f_mic_s <- function(beta_h, beta_t, n_s){
#  result <- 0
#  p <- 0
#  G <- length(beta_h)
#  for (g in 1:G){
#    result <- result + f_sse_s(beta_h[[g]], beta_t[[g]])
#    p <- p + dim(beta_h[[g]])[1]
#  }
#  result <- result + p * sum(log(n_s)/n_s)
#  return(result)
#}

f_mic_s <- function(beta_h, beta_t, n_s){
  result <- 0
  p <- 0
  G <- length(beta_h)
  for (g in 1:G){
    result <- result + f_sse_s(beta_h[[g]], beta_t[[g]])
    for (s in 1:length(n_s)){
      p <- p + length(which(beta_h[[g]][,s]!=0))*log(n_s[s])/n_s[s]
    }
  }
  result <- result + p
  return(result)
}

f_mic_v_inv <- function(beta_h, beta_t, n_s, inv_var){
  result <- 0
  p <- 0
  G <- length(beta_h)
  for (g in 1:G){
    result <- result + f_sse_v_inv(beta_h[[g]], beta_t[[g]], inv_var[[g]], n_s)
    for (s in 1:length(n_s)){
      p <- p + length(which(beta_h[[g]][,s]!=0))*log(n_s[s])/n_s[s]
    }
  }
  result <- result + p
  return(result)
}

f_mic_v_inv_per_var <- function(beta_h, beta_t, n_s, inv_var){
  result <- 0
  p <- 0
  G <- length(beta_h)
  nstud <- length(n_s)
  for (g in 1:G){
    result <- result + f_sse_v_inv(beta_h[[g]], beta_t[[g]], inv_var[[g]], n_s)
    for (s in 1:length(n_s)){
      p <- p + length(which(beta_h[[g]][,s]!=0))*log(n_s[s])/n_s[s]
    }
  }
  result <- result + p/nstud
  return(result)
}

f_mic_v_inv_us <- function(beta_h, beta_t, n_s, inv_var){
  result <- 0
  p <- 0
  G <- length(beta_h)
  for (g in 1:G){
    result <- result + f_sse_v_inv_us(beta_h[[g]], beta_t[[g]], inv_var[[g]], n_s)
    for (s in 1:length(n_s)){
      p <- p + length(which(beta_h[[g]][,s]!=0))*log(n_s[s])/n_s[s]
    }
  }
  result <- result + p
  return(result)
}


# Getting the starting lambda...
###################################

# To find lambda max if only variable x studies structure is considered for regularization
find_lambda_max_in_g_for_vsnorm <- function (sumstat, w_p, coeff=1){
  G <- length(sumstat)
  K <- length(sumstat[[1]]$beta[1,])
  v_lambda_max <- rep(NA, G)
  
  for (g in 1:G){
    v_this_g <- array(NA, dim=dim(sumstat[[g]]$beta))
    for (k in 1:K){
      v_this_g[,k] <- (2/ w_p[[g]]) * sqrt(rowSums(sumstat[[g]]$beta^2)) * abs(sumstat[[g]]$cov_inv[,,k] %*% sumstat[[g]]$beta[,k]) / abs(sumstat[[g]]$beta[,k])
    }
    v_lambda_max[g] <- max(v_this_g, na.rm=TRUE)
  }
  result <- max(v_lambda_max)/coeff
  return(result)
}

find_lambda_max_in_g_for_gnorm <- function (sumstat, w_g, coeff=1){
  G <- length(sumstat)
  K <- length(sumstat[[1]]$beta[1,])
  v_lambda_max <- rep(NA, G)
  
  for (g in 1:G){
    v_this_g <- array(NA, dim=dim(sumstat[[g]]$beta))
    for (k in 1:K){
      v_this_g[,k] <- (2/ w_g[[g]]) * norm(sumstat[[g]]$beta, type="F") * abs(sumstat[[g]]$cov_inv[,,k] %*% sumstat[[g]]$beta[,k]) / abs(sumstat[[g]]$beta[,k])
    }
    v_lambda_max[g] <- max(v_this_g, na.rm=TRUE)
  }
  result <- max(v_lambda_max)/coeff
  return(result)
}

### OLDER FUNCTIONS

# # getting the norm of the convex optimisation problem (step z) using beta_t the same shape as inputs
# f_lambda_max_inputs <- function (sumstat, w_p, alpha){
  # max_v <- 0
  # for (g in 1:length(sumstat)){
    # max_v <- max(max_v, sumstat[[g]]$beta^2 |> apply(FUN=sum, MARGIN=1) |> sqrt |> {.*w_p[[g]]} |> max)
  # }
  # result <- max_v / alpha
  # return(result)
# }

# # Same function but with beta as direct list of elements
# f_lambda_max_list <- function (beta, w_p, alpha){
  # max_v <- 0
  # for (g in 1:length(beta)){
    # max_v <- max(max_v, beta[[g]]^2 |> apply(FUN=sum, MARGIN=1) |> sqrt |> {.*w_p[[g]]} |> max)
  # }
  # result <- max_v / alpha
  # return(result)
# }

# General frobenius norm from data matrix in list (of groups)
##################################################

f_f_norm_mat_in_list <- function(beta_h){
  result <- sapply(FUN=function(x) x<- (x^2 |> sum()), beta_h, simplify=TRUE) |> sum() |> sqrt()
}

# Creation of a grid from lambda_min = eps*lambda_max to lambda_max, in log scale
##################################################
f_grid_loge <- function(lambda_max, eps=1e-03, npts=20){
  a <- log(lambda_max*eps, base = exp(1))
  b <- log(lambda_max, base = exp(1))
  result <- seq(from=a, to=b, length.out = npts) |> exp()
  return(result)
}

f_grid_log10 <- function(lambda_max, eps=1e-03, npts=20){
  a <- log(lambda_max*eps, base = 10)
  b <- log(lambda_max, base = 10)
  result <- seq(from=a, to=b, length.out = npts) |> exp()
  return(result)
}

#  To run the meta psg solver on a grid
##################################################

find_lambda_max <- function(sumstat, n_s, alpha=0.1, w_g, w_p, f_test=1.5, f_selec=1.2){
  result <- 0
  l_max <- min(find_lambda_max_in_g_for_gnorm(sumstat, w_g), find_lambda_max_in_g_for_vsnorm(sumstat, w_p))
  test <- metapsg_gv(sumstat, n_s, alpha=alpha, lambda=l_max, w_p=w_p, w_g=w_g)
  while (test$info[dim(test$info)[1],]$z_norm==0) {
    l_max <- l_max / 2
    test <- metapsg_gv(sumstat, n_s, alpha=alpha, lambda=l_max, w_p=w_p, w_g=w_g)
  }
  test <- metapsg_gv(sumstat, n_s, alpha=alpha, lambda=f_test*l_max, w_p=w_p, w_g=w_g)
  if(test$info[dim(test$info)[1],]$z_norm!=0){result <- f_test * l_max}else{result <- f_selec * l_max}
  return(result)
}
