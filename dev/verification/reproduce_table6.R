# Full reproduction of the application of Sugier et al. (2025), Table 6:
# 331 genes, alpha = 0.1..0.9, 30 lambda values, with and without adaptive
# weights, compared with the results saved in paper/application/results_mpsg/.
#
# Usage (from the repository root, after installing the package):
#   Rscript dev/verification/reproduce_table6.R
# Number of cores: environment variable MPSG_NCPU (default: all cores but one).

library(MPSG)
t_start <- Sys.time()
dir <- "paper/application"
load(file.path(dir, "summary_stat_application_gaussian_prior.RData"))  # SUMSTAT_bglm_gaussian
load(file.path(dir, "info_annot_data_BC_TC.RData"))                    # genes, annoted_gene, annoted_snp, n_s
G <- length(genes)
ncpu <- as.integer(Sys.getenv("MPSG_NCPU", max(1, parallel::detectCores() - 1)))
cat(sprintf("%d genes, %d SNPs, %d cores\n", G, length(annoted_gene), ncpu))

# Inputs built as in paper/application/script_metapsg_BC_TC.r
sumstat <- vector("list", G)
ridged <- character(0)
for (g in seq_len(G)) {
  ids <- which(annoted_gene == genes[g]); mg <- length(ids)
  s <- SUMSTAT_bglm_gaussian[[g]]
  covar <- array(0, c(mg, mg, 2)); covar[, , 1] <- s$sigmah1; covar[, , 2] <- s$sigmah2
  ci <- tryCatch(f_cpp_inv_var(covar = covar, same_var = FALSE), error = function(e) NULL)
  if (is.null(ci)) {
    ridged <- c(ridged, genes[g])
    for (k in 1:2) covar[, , k] <- covar[, , k] + 0.01 * diag(diag(covar[, , k]))
    ci <- f_cpp_inv_var(covar = covar, same_var = FALSE)
  }
  sumstat[[g]] <- list(beta = cbind(s$betah1, s$betah2), cov_inv = ci)
}
cat("Ridge added to:", paste(ridged, collapse = ", "), "\n")

configs <- list(
  list(label = "adaptive weights", adaptive = TRUE,
       ref = file.path(dir, "results_mpsg/RESULTS_MetaPSG_Best_Model_adaptive_weights.RData")),
  list(label = "standard weights", adaptive = FALSE,
       ref = file.path(dir, "results_mpsg/RESULTS_MetaPSG_Best_Model.RData")))

summary_lines <- character(0)
for (cfg in configs) {
  t0 <- Sys.time()
  res <- f_selec_mod_metapsg_gv_parallel(sumstat, n_s, alpha = seq(0.1, 0.9, 0.1), adaptive = cfg$adaptive,
                                         same_var = FALSE, rho = 1, t_inc = 2, t_dec = 2, mu = 10,
                                         max_iter = 1000, tol_abs = 1e-03, tol_rel = 1e-02,
                                         eps = 1e-03, npts = 30, ncpu = ncpu)
  elapsed <- difftime(Sys.time(), t0, units = "mins")

  model <- res$model
  for (g in seq_len(G)) rownames(model[[g]]) <- annoted_snp[annoted_gene == genes[g], ]
  names(model) <- genes
  selected <- model[vapply(model, function(b) norm(b) != 0, logical(1))]

  ref_env <- new.env(); load(cfg$ref, envir = ref_env)
  ref <- ref_env$metapsg_pleiotropic_genes

  same_genes <- identical(names(selected), names(ref))
  common <- intersect(names(selected), names(ref))
  n_sel <- function(b) sum(rowSums(b^2) > 0)
  same_counts <- all(vapply(common, function(g) n_sel(selected[[g]]) == n_sel(ref[[g]]), logical(1)))
  max_diff <- max(vapply(common, function(g) max(abs(unname(selected[[g]]) - unname(ref[[g]]))), numeric(1)))

  out <- c(sprintf("== %s ==", cfg$label),
           sprintf("time: %.1f min", as.numeric(elapsed)),
           sprintf("best model: alpha = %s, lambda = %s, MIC = %s",
                   res$info_model$alpha, signif(res$info_model$lambda, 8), signif(res$info_model$MIC, 8)),
           sprintf("genes selected: %d (reference: %d)", length(selected), length(ref)),
           sprintf("same genes, same order: %s", same_genes),
           sprintf("same number of selected SNPs per gene: %s", same_counts),
           sprintf("max |beta - beta_ref| over common genes: %.3g", max_diff),
           "genes (selected SNPs / reference):",
           paste0("  ", vapply(union(names(selected), names(ref)), function(g)
             sprintf("%-8s %4s / %4s", g,
                     if (g %in% names(selected)) n_sel(selected[[g]]) else "-",
                     if (g %in% names(ref)) n_sel(ref[[g]]) else "-"), character(1))))
  cat(out, sep = "\n")
  summary_lines <- c(summary_lines, out, "")
}

summary_lines <- c(summary_lines, sprintf("total time: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
writeLines(summary_lines, "dev/verification/reproduce_table6_results.txt")
cat("Summary written to dev/verification/reproduce_table6_results.txt\n")
