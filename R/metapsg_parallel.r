
# Function to create alpha/lambda grid
f_metapsg_gv_parallel_grid <- function(sumstat, n_s, alpha, w_g, w_p, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=20){
    #setup parallel backend to use many processors
    #
	n_alpha <- length(alpha)
    cl_grid <- makeCluster(n_alpha, outfile="") #not to overload your computer
    print(cl_grid)
    registerDoParallel(cl_grid)
	#
    # To build the grid of lambda for each value of alpha
    found_grid <- foreach(i=1:n_alpha, .packages = "MPSG") %dopar% {
        #redefine functions
        lambda_max <- find_lambda_max(sumstat=sumstat, n_s=n_s, alpha=alpha[i], w_g=w_g, w_p=w_p, rho=rho, t_inc=t_inc, t_dec=t_dec, mu=mu, max_iter=max_iter, tol_abs=tol_abs, tol_rel=tol_rel)
        lambda_grid <- f_grid_loge(lambda_max, eps=eps, npts=npts)
		list(alpha=alpha[i], lambdas=lambda_grid)
    }
	stopCluster(cl_grid)
	#
	list_grid<-cbind(alpha=rep(found_grid[[1]]$alpha, npts), lambda=found_grid[[1]]$lambdas)
	for (i in seq_len(n_alpha)[-1]){list_grid <- rbind(list_grid,cbind(alpha=rep(found_grid[[i]]$alpha, npts), lambda=found_grid[[i]]$lambdas))}
	#
	return(list_grid)
}   

#' Fit and select the MPSG model on a grid of alpha and lambda, in parallel
#'
#' Fits the penalised multivariate meta-analysis model with group structure
#' with an ADMM algorithm on a two-dimensional grid of values of `alpha` and
#' `lambda`, using parallel worker processes, and returns the model with the
#' smallest modified information criterion (MIC). This is the interface used
#' for the analyses of Sugier et al. (2025).
#'
#' For each value of `alpha`, a grid of `npts` values of `lambda` is built as
#' in [f_selec_mod_metapsg_gv()] (one worker process per value of `alpha`).
#' The models of the whole grid are then fitted on a cluster of `ncpu` worker
#' processes. The MIC is always computed with the inverse covariance matrices
#' (as with `gwas = FALSE` in [f_selec_mod_metapsg_gv()]). The function prints
#' information on the clusters and the MIC values along the grid.
#'
#' @inheritParams f_selec_mod_metapsg_gv
#' @param alpha Numeric vector of values of `alpha` (share of the penalty
#'   given to the variable level, see [f_selec_mod_metapsg_gv()]).
#' @param eps Ratio between the smallest and the largest value of the grid of
#'   `lambda`, for each value of `alpha`.
#' @param npts Number of values in the grid of `lambda`, for each value of
#'   `alpha`.
#' @param ncpu Number of worker processes used to fit the models of the grid.
#'
#' @return A list with four components describing the model with the smallest
#'   MIC (the first one in case of ties):
#'   * `model`: list with one `n_g x K` matrix of penalised estimates per
#'     group.
#'   * `info_model`: one-row data frame with columns `alpha`, `lambda`, `MIC`,
#'     `iter` (number of ADMM iterations) and `norm_betah` (Frobenius norm of
#'     the penalised estimates).
#'   * `results`: data frame of the ADMM iterations of this model (`iter`,
#'     `r_norm`, `s_norm`, `rho`, `b_norm`, `z_norm`, `u_norm`).
#'   * `MIC`: numeric vector of the MIC of all the models of the grid, ordered
#'     by value of `alpha` (as given) and, within each, by increasing `lambda`.
#'
#' @references Sugier, P.-E., Asgari, Y., Sedki, M., Truong, T. and Liquet, B.
#'   (2025). Biostatistics, 26(1), kxaf037. \doi{10.1093/biostatistics/kxaf037}
#'
#' @seealso [f_selec_mod_metapsg_gv()], [f_cpp_inv_var()].
#'
#' @export
f_selec_mod_metapsg_gv_parallel <- function(sumstat, n_s, alpha, adaptive=TRUE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=20, ncpu){
  if (isTRUE(same_var)) warning("'same_var' is not supported and is ignored.")
  #
  #Definitions step
  #
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
  grid_to_run <- f_metapsg_gv_parallel_grid(sumstat=sumstat,n_s=n_s,alpha=alpha,w_g,w_p,rho=rho,t_inc=t_inc,t_dec=t_dec,mu=mu,max_iter=max_iter,tol_abs=tol_abs,tol_rel=tol_rel,eps=eps,npts=npts)
  grid_to_run <- grid_to_run
  #
  #setup parallel backend to use many processors
  #
  cl <- makeCluster(ncpu, outfile="") #not to overload your computer
  print(cl)
  registerDoParallel(cl)
  #
  # To run  meta-psg solver on the grid
  #
  # by using parallel mode
  #
  results <- foreach(i=1:dim(grid_to_run)[1], .packages = "MPSG") %dopar% {
    #redefine functions
    # temporary stock
    tempo <- metapsg_gv(sumstat, n_s, alpha=grid_to_run[i,1], lambda=grid_to_run[i,2], gwas=FALSE, rho=rho, t_inc=t_inc, t_dec=t_dec, mu=mu, w_p=w_p, w_g=w_g, max_iter=max_iter, tol_abs=tol_abs, tol_rel=tol_rel)
    # information for model selection
    list(info_grid=data.frame(alpha=grid_to_run[i,1], lambda=grid_to_run[i,2], MIC=tempo$MIC, iter=dim(tempo$info)[1], norm_betah=tempo$info[dim(tempo$info)[1], "z_norm"]), beta_select=tempo$beta_select, info=tempo$info)
  }
  #
  stopCluster(cl)
  # To build dataframe about resume info on the models along the grid
  #
  print('para ok')
  print(length(results))
  #
  MIC_values <- rep(NA,dim(grid_to_run)[1])
  for (i in 1:dim(grid_to_run)[1]){MIC_values[i]<-results[[i]]$info_grid$MIC}
  min_index <- which(MIC_values == min(MIC_values), arr.ind = TRUE)[1]
  #
  print("MIC_values")
  print(MIC_values)
  #
  print("min_index")
  print(min_index)
  #
  #
  res <- list(model=results[[min_index]]$beta_select, info_model=results[[min_index]]$info_grid, results=results[[min_index]]$info, MIC=MIC_values)
  return(res)
}

# TEST
# alpha<-c(0.5,0.9)
# grid_to_run <- f_metapsg_gv_parallel_grid(sumstat, n_s, alpha, w_g, w_p, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=5)
# grid_to_run <- as.matrix(grid_to_run)
# out <- f_selec_mod_metapsg_gv_parallel(sumstat, n_s, alpha, adaptive=TRUE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=5, ncpu=10)


