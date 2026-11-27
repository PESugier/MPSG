# Generates reference results with the ORIGINAL code of the paper (git tag v1.1),
# on the first 15 genes of the application data. Used as a regression test fixture.
#
# Usage (from the repository root):  Rscript dev/verification/make_reference_15genes.R <output.rds>
# Requires: git, Rtools (sourceCpp), packages Rcpp, RcppArmadillo, Matrix, abind.

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else "dev/verification/ref_15genes.rds"
root <- normalizePath(".")
out_file <- file.path(root, out_file)

# 1. Data (loaded before changing directory)
load(file.path(root, "paper/application/summary_stat_application_gaussian_prior.RData"))
load(file.path(root, "paper/application/info_annot_data_BC_TC.RData"))

# 2. Original code extracted from tag v1.1 into a temporary directory
tmp <- file.path(tempdir(), "mpsg_v1.1")
dir.create(tmp, showWarnings = FALSE)
for (f in c("metapsg_gv_solver.r", "metapsg_sub.r", "arma_inv.cpp")) {
  txt <- system2("git", c("-C", shQuote(root), "show", paste0("v1.1:R/", f)), stdout = TRUE)
  if (!is.null(attr(txt, "status"))) stop("git show failed for ", f)
  writeLines(sub("\r$", "", txt), file.path(tmp, f))
}
old_wd <- setwd(tmp)
suppressMessages({ library(Matrix); library(Rcpp); library(RcppArmadillo) })
source("metapsg_gv_solver.r")   # also sources metapsg_sub.r and compiles arma_inv.cpp
setwd(old_wd)

# 3. Inputs in the original format (as in paper/application/script_metapsg_BC_TC.r)
n_genes <- 15
sumstat <- list()
for (g in seq_len(n_genes)) {
  ids <- which(annoted_gene == genes[g]); mg <- length(ids)
  s <- SUMSTAT_bglm_gaussian[[g]]
  covar <- array(0, c(mg, mg, 2)); covar[, , 1] <- s$sigmah1; covar[, , 2] <- s$sigmah2
  ci <- tryCatch(f_cpp_inv_var(covar = covar, same_var = FALSE), error = function(e) NULL)
  if (is.null(ci)) {
    for (k in 1:2) covar[, , k] <- covar[, , k] + 0.01 * diag(diag(covar[, , k]))
    ci <- f_cpp_inv_var(covar = covar, same_var = FALSE)
  }
  beta <- cbind(s$betah1, s$betah2); rownames(beta) <- annoted_snp[ids, ]
  sumstat[[g]] <- list(beta = beta, cov_inv = ci)
}
names(sumstat) <- genes[seq_len(n_genes)]

# 4. Original grid search (sequential function) for 4 configurations
configs <- list()
for (adaptive in c(TRUE, FALSE)) for (alpha in c(0.3, 0.7)) {
  r <- f_selec_mod_metapsg_gv(unname(sumstat), n_s, alpha = alpha, adaptive = adaptive, eps = 1e-3, npts = 8)
  sm <- r$selmod
  configs[[length(configs) + 1]] <- list(
    adaptive = adaptive, alpha = alpha, eps = 1e-3, npts = 8,
    selmod = data.frame(lambda = as.numeric(sm$lambda), MIC = as.numeric(sm$MIC), iter = as.numeric(sm$iter),
                        ng_hit = as.numeric(sm$ng_hit), nv_hit = as.numeric(sm$nv_hit)),
    beta_select = lapply(r$results, function(x) x$beta_select))
  cat(sprintf("adaptive = %s, alpha = %.1f : done\n", adaptive, alpha))
}

saveRDS(list(sumstat = unname(sumstat), n_s = n_s, configs = configs,
             provenance = list(code = "git tag v1.1", R = R.version.string,
                               platform = R.version$platform, date = as.character(Sys.Date()))),
        out_file)
cat("Saved:", out_file, "\n")
