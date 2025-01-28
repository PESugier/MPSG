library(tictoc)
library(doParallel)
library(doSNOW)

source("metapsg_gv_solver.r")

# Function to create alpha/lambda grid
f_metapsg_gv_parallel_grid <- function(sumstat, n_s, alpha, w_g, w_p, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=20){
    #setup parallel backend to use many processors
    #
	n_alpha <- length(alpha)
    cl_grid <- makeCluster(n_alpha, outfile="") #not to overload your computer
    print(cl_grid)
    registerDoParallel(cl_grid)
    registerDoSNOW(cl_grid)
	#
    # To build the grid of lambda for each value of alpha
    found_grid <- foreach(i=1:n_alpha) %dopar% {
        #redefine functions
        source("metapsg_gv_solver.r")
        lambda_max <- find_lambda_max(sumstat=sumstat, n_s=n_s, alpha=alpha[i], w_g=w_g, w_p=w_p)
        lambda_grid <- f_grid_loge(lambda_max, eps=eps, npts=npts)
		list(alpha=alpha[i], lambdas=lambda_grid)
    }
	stopCluster(cl_grid)
	#
	list_grid<-cbind(alpha=rep(found_grid[[1]]$alpha, npts), lambda=found_grid[[1]]$lambdas)
	for (i in 2:n_alpha){list_grid <- rbind(list_grid,cbind(alpha=rep(found_grid[[i]]$alpha, npts), lambda=found_grid[[i]]$lambdas))}
	#
	return(list_grid)
}   

f_selec_mod_metapsg_gv_parallel <- function(sumstat, n_s, alpha, adaptive=TRUE, same_var = FALSE, rho=1, t_inc=2, t_dec=2, mu=10, max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02, eps=1e-03, npts=20, ncpu){
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
  grid_to_run <- f_metapsg_gv_parallel_grid(sumstat=sumstat,n_s=n_s,alpha=alpha,w_g,w_p,same_var=same_var,rho=rho,t_inc=t_inc,t_dec=t_dec,mu=mu,max_iter=max_iter,tol_abs=tol_abs,tol_rel=tol_rel,eps=eps,npts=npts)
  grid_to_run <- grid_to_run
  #
  #setup parallel backend to use many processors
  #
  cl <- makeCluster(ncpu, outfile="") #not to overload your computer
  print(cl)
  registerDoParallel(cl)
  registerDoSNOW(cl)
  #
  # To run  meta-psg solver on the grid
  #
  # by using parallel mode
  #
  results <- foreach(i=1:dim(grid_to_run)[1]) %dopar% {
    #redefine functions
    source("/home1/6_AMLAP/PE/METAPSG/scripts/metapsg_gv_solver.r")
    # temporary stock
    tempo <- metapsg_gv(sumstat, n_s, alpha=grid_to_run[i,1], lambda=grid_to_run[i,2], same_var = FALSE, gwas=FALSE, rho, t_inc, t_dec, mu, w_p, w_g, max_iter, tol_abs, tol_rel)
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


