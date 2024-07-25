#library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("/home1/6_AMLAP/PE/METAPSG/scripts/src/arma_inv.cpp")
source("/home1/6_AMLAP/PE/METAPSG/scripts/metapsg_sub.r")


##### Architecture
# inputs : list per grps / then per variable[dim1] x studies [dim2]
# -> beta_t [[1:G]] [1:ng,1,1:s]
# -> var_t [[1:G]] [1:ng,1:ng,1:s]
############################################


# META_PSG to solve Group (g) + Variable (v) regularization.
#
#  Most suited to problems with k=2
#############################################

metapsg_gv <- function(sumstat, n_s, alpha=0.1, lambda, same_var = FALSE, gwas = FALSE, 
                               rho=1, t_inc=2, t_dec=2, mu=10,
                               w_p, w_g, 
                               max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02) {
  #
  # print('start running gv solver')
  #
  # Shape of inputs:
  # sumstat[[g]]$beta[1:n_g,1:k]
  # sumstat[[g]]$cov_inv[1:n_g,1:n_g,1:k]
  # n_s : vector with number of individuals for each study (same order than beta & cov_inv) 
  G <- length(sumstat) # number of groups
  S <- dim(sumstat[[1]]$cov_inv)[3] # number of studies
  
  lambda1 <- (1-alpha)*lambda
  lambda2 <- alpha*lambda
    
  # Extraction of beta estimates (beta inputs)
  beta_t <- list()
  beta_t <- mapply(function(x) x$beta, sumstat, SIMPLIFY=FALSE)
  
  P <- length(unlist(beta_t))/S
    
  # Calculation of Inverse Variance Matrix
  inv_var_t <- list()
  for (g in 1:G){inv_var_t[[g]] <- sumstat[[g]]$cov_inv}
  
  # Calculation of the solver Matrix : S = (2 V^-1 + rho I)^-1
  mrho_inv <- list()
  for (g in 1:G){
    # print(g)
    mrho_inv[[g]] <- f_cpp_mrho_inv(covar_inv=inv_var_t[[g]], rho=rho, same_var=same_var)
  }
  
  # Starting values for b, z & u
  b <- beta_t
  z <- mapply(function(x) array(0, dim=dim(x)), beta_t, SIMPLIFY=FALSE)
  u <- mapply(function(x) array(0, dim=dim(x)), beta_t, SIMPLIFY=FALSE)
  
  # ADMM Solver
  ########################################
  cc_info <- data.frame(iter=0, r_norm=NA, s_norm=NA, rho=NA, b_norm=NA, z_norm=NA, u_norm=NA)
  kt <- 0
  
  # To force ADMM stop before max_iter
  for (iter in 1:max_iter) {
    
    b_old <- b
    z_old <- z
    u_old <- u
    x <- mapply(`-`, z_old, u_old, SIMPLIFY=FALSE)

    # Step 1 - l2-regularized meta-analysis problem solving (update of B according to r = z_old - u_old)
    #b <- mapply(FUN=meta_regul_solver_g, beta_t, inv_var_t, r, rho, mrho_inv, SIMPLIFY=FALSE)
    for (g in 1:G){
      b[[g]] <- meta_regul_solver_g(b=beta_t[[g]], m=inv_var_t[[g]], x=x[[g]], rho=rho, mrho_inv=mrho_inv[[g]])
    }

    # Step 2 - Convex optimisation problem (update of Z according to t = b + u_old)
    tt <- mapply(`+`, b, u_old, SIMPLIFY=FALSE)
    z <- mapply(prox_L21, tt, lambda2/rho, w_p, SIMPLIFY=FALSE)
    z <- mapply(prox_G21, z, lambda1/rho, w_g, SIMPLIFY=FALSE)
    
    # Step 3 - Dual variable update (update of U)
    for (g in 1:G){
      u[[g]] <- u_old[[g]] + b[[g]] - z[[g]]
    }
    
    # update count iter
    kt <- kt + 1
    
    # Defining convergency criteria
    ##########################################
    cc_r_norm <- (unlist(b) - unlist(z)) |> as.matrix() |> norm(type="F")
    cc_s_norm <- (unlist(z) - unlist(z_old)) |> as.matrix() |> {\(x) rho * x}() |> norm(type="F")
    cc_rho <- rho
    cc_b_norm <- unlist(b) |> as.matrix() |> norm(type="F")
    cc_z_norm <- unlist(z) |> as.matrix() |> norm(type="F")
    cc_u_norm <- unlist(u) |> as.matrix() |> norm(type="F")
    
    cc_eps_pri <- sqrt(P*S) * tol_abs + tol_rel * max(norm(as.matrix(unlist(b)),type="F"), norm(as.matrix(-unlist(z)),type="F"))
    cc_eps_dual <- sqrt(P*S) * tol_abs + tol_rel * norm(as.matrix(unlist(u)),type="F")
    
    cc_info <- rbind(cc_info,c(kt, cc_r_norm, cc_s_norm, cc_rho, cc_b_norm, cc_z_norm, cc_u_norm))
    
    
    if (cc_r_norm < cc_eps_pri && cc_s_norm < cc_eps_dual) {
        break
    }
    
    # Varying Penalty Parameter 
    ##########################################
    if (cc_r_norm > mu * cc_s_norm) {
        rho <- rho * t_inc
        u <- mapply(FUN=function(x) x<-x/t_inc, u, SIMPLIFY=FALSE)
        mrho_inv <- list()
        for (g in 1:G){mrho_inv[[g]] <- f_cpp_mrho_inv(covar_inv=inv_var_t[[g]], rho=rho, same_var=same_var)}
    } else if (cc_s_norm > mu * cc_r_norm) {
        rho <- rho/t_dec
        u <- mapply(FUN=function(x) x<-x*t_dec, u, SIMPLIFY=FALSE)
        mrho_inv <- list()
        for (g in 1:G){mrho_inv[[g]] <- f_cpp_mrho_inv(covar_inv=inv_var_t[[g]], rho=rho, same_var=same_var)}
    }
    
  }

  ### Calculation of the MIC 
    if (gwas){val_mic <- f_mic_s(z, beta_t, n_s)}else{val_mic <- f_mic_v_inv(z, beta_t, n_s, inv_var_t)}
    
  #### Return outputs
    return(list(beta_select=z, info=cc_info[-1,], MIC=val_mic))
}


#  To run metapsg solvers on a 1D grid based on a alpha value
###############################################################

f_selec_mod_metapsg_gv <- function(sumstat, n_s, alpha, adaptive=TRUE, same_var = FALSE, gwas = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-02, npts=20){
  selmod <- data.frame(alpha=NA, lambda=NA, MIC=NA, iter=NA, norm_betah=NA)
  results <- list()
  G <- length(sumstat)
  if(adaptive){
    w_p <- list()
    for (g in 1:G){w_p[[g]] <- f_weight_p_per_group(sumstat[[g]]$beta, adaptive=TRUE)}
    w_g <- list()
    for (g in 1:G){w_g[[g]] <- f_weight_g_per_group(sumstat[[g]]$beta, adaptive=TRUE)}
  }else{
    w_p <- list()
    for (g in 1:G){w_p[[g]] <- f_weight_p_per_group(sumstat[[g]]$beta, adaptive=FALSE)}
    w_g <- list()
    for (g in 1:G){w_g[[g]] <- f_weight_g_per_group(sumstat[[g]]$beta, adaptive=FALSE)} 
  }
  # To build the grid of 
  lambda_max <- find_lambda_max(sumstat=sumstat, n_s=n_s, alpha=alpha, w_g=w_g, w_p=w_p)
  lambda_grid <- f_grid_loge(lambda_max, eps=eps, npts=npts)
  # To run  meta-psg solver on the grid
  for (i in 1:length(lambda_grid)){
    results[[i]] <- list(info_grid=NA, beta_select=NA, info=NA)
    # temporary stock
    tempo <- metapsg_gv(sumstat, n_s, alpha, lambda_grid[i], same_var = FALSE, gwas = gwas, rho, t_inc, t_dec, mu, w_p, w_g, max_iter, tol_abs, tol_rel)
    # information for model selection
    temp_selmod <- data.frame(alpha=NA, lambda=NA, MIC=NA, iter=NA, norm_betah=NA, ng_hit=NA, nv_hit=NA)
    temp_selmod$alpha <- alpha
    temp_selmod$lambda <- lambda_grid[i]
    temp_selmod$MIC <- tempo$MIC
    temp_selmod$iter <- dim(tempo$info)[1]
    temp_selmod$norm_betah <- tempo$info[dim(tempo$info)[1], "z_norm"]
    # To get number of non zero groups
    G <- length(sumstat)
    ng_hit <- 0
    nv_hit <- 0
    for (g in 1:G){
        g_norm <- sqrt(sum(tempo$beta_select[[g]]^2))
        if(g_norm!=0){ng_hit <- ng_hit +1}
        v_norm <- tempo$beta_select[[g]]^2
        v_norm <- apply(v_norm, MARGIN = 1, FUN = sum) |> sqrt()
        nv_hit <- nv_hit + length(which(v_norm!=0))
    }      
    temp_selmod$ng_hit <- ng_hit
    temp_selmod$nv_hit <- nv_hit
    if(i==1){selmod <- temp_selmod}else{selmod <- rbind(selmod, temp_selmod)}
    # information for model selection
    results[[i]]$info_grid <- data.frame(alpha, lambda_grid[i])
    results[[i]]$beta_select <- tempo$beta_select
    results[[i]]$info <- tempo$info 
  }
  res <- list(selmod=selmod, results=results)
  return(res)
}
