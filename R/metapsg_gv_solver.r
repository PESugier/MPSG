#library(tidyverse)



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

#' Fit the MPSG model along a grid of penalty values
#'
#' Fits the penalised multivariate meta-analysis model with group structure
#' (variables nested in groups, effects estimated in several studies) with an
#' ADMM algorithm, for one value of `alpha` and a grid of values of the
#' penalty `lambda`. This is the interface used for the analyses of Sugier
#' et al. (2025).
#'
#' The penalty `lambda` is split into `(1 - alpha) * lambda` for the
#' group-level penalty and `alpha * lambda` for the variable-level penalty.
#' The grid contains `npts` values of `lambda`, equally spaced on the log
#' scale between `eps * lambda_max` and `lambda_max`, where `lambda_max` is
#' computed from the data by an internal heuristic. The function returns the
#' fits for all values of the grid and does not select a model.
#'
#' @param sumstat List with one element per group (for instance a gene). Each
#'   element is a list with two components: `beta`, a numeric matrix of
#'   dimension `n_g x K` of estimated effects (one row per variable of the
#'   group, one column per study), and `cov_inv`, a numeric array of dimension
#'   `n_g x n_g x K` of inverse covariance matrices of these estimates (see
#'   [f_cpp_inv_var()]).
#' @param n_s Numeric vector of length `K`: number of individuals in each
#'   study, in the same order as the columns of `beta`.
#' @param alpha Single numeric value: share of the penalty given to the
#'   variable level (`alpha * lambda`); the group level receives
#'   `(1 - alpha) * lambda`.
#' @param adaptive Logical. If `TRUE` (the default), adaptive weights are
#'   used: the weight of a variable is the inverse of the Euclidean norm of
#'   its estimated effects across studies, and the weight of a group is the
#'   inverse of the Frobenius norm of its matrix `beta`. If `FALSE`, the
#'   weight of a variable is `sqrt(K)` and the weight of a group is
#'   `sqrt(n_g * K)`.
#' @param same_var Logical. Currently ignored: the solver is always called
#'   with `same_var = FALSE`.
#' @param gwas Logical. If `FALSE` (the default), the goodness-of-fit term of
#'   the modified information criterion (MIC) is the sum, over studies, of the
#'   quadratic form of the differences between penalised and input estimates
#'   with the matrix `cov_inv / n_s`. If `TRUE`, it is the plain sum of
#'   squared differences.
#' @param rho Positive numeric: initial value of the ADMM penalty parameter.
#' @param t_inc,t_dec Numeric factors by which `rho` is multiplied (`t_inc`)
#'   or divided (`t_dec`) when the norm of the primal residual is larger than
#'   `mu` times the norm of the dual residual, or the reverse.
#' @param mu Numeric: threshold on the ratio of the primal and dual residual
#'   norms that triggers an update of `rho`.
#' @param max_iter Maximum number of ADMM iterations for each value of
#'   `lambda`.
#' @param tol_abs,tol_rel Absolute and relative tolerances of the ADMM
#'   stopping criterion on the primal and dual residuals.
#' @param eps Ratio between the smallest and the largest value of the grid of
#'   `lambda`.
#' @param npts Number of values in the grid of `lambda`.
#'
#' @return A list with two components:
#'   * `selmod`: a data frame with one row per value of `lambda` (in
#'     increasing order) and columns `alpha`, `lambda`, `MIC` (modified
#'     information criterion), `iter` (number of ADMM iterations),
#'     `norm_betah` (Frobenius norm of the penalised estimates), `ng_hit`
#'     (number of groups with at least one non-zero estimate) and `nv_hit`
#'     (number of variables with at least one non-zero estimate).
#'   * `results`: a list with one element per value of `lambda`, in the same
#'     order, each a list with `info_grid` (data frame with the values of
#'     `alpha` and `lambda`), `beta_select` (list with one `n_g x K` matrix of
#'     penalised estimates per group) and `info` (data frame of the ADMM
#'     iterations: `iter`, `r_norm`, `s_norm`, `rho`, `b_norm`, `z_norm`,
#'     `u_norm`).
#'
#' @references Sugier, P.-E., Asgari, Y., Sedki, M., Truong, T. and Liquet, B.
#'   (2025). Biostatistics, 26(1), kxaf037. \doi{10.1093/biostatistics/kxaf037}
#'
#' @seealso [f_selec_mod_metapsg_gv_parallel()] for several values of `alpha`
#'   and model selection, [f_cpp_inv_var()].
#'
#' @export
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
