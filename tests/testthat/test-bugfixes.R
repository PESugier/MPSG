# Tests of the bug fixes of step 4 (see dev/DECISIONS.md, "Bugs connus").

ref <- readRDS(test_path("fixtures", "ref_15genes.rds"))
small <- ref$sumstat[1:5]

test_that("an observed effect exactly equal to 0 no longer causes an infinite loop", {
  s <- small
  s[[2]]$beta[1, 1] <- 0
  setTimeLimit(elapsed = 120)
  on.exit(setTimeLimit(), add = TRUE)
  res <- f_selec_mod_metapsg_gv(s, ref$n_s, alpha = 0.5, adaptive = FALSE, eps = 1e-3, npts = 3)
  expect_true(all(is.finite(as.numeric(res$selmod$lambda))))
})

test_that("effects that are all zero give an explicit error", {
  s <- lapply(small, function(g) { g$beta[] <- 0; g })
  expect_error(f_selec_mod_metapsg_gv(s, ref$n_s, alpha = 0.5, adaptive = FALSE, npts = 3),
               "starting value for lambda")
})

test_that("same_var = TRUE is ignored with a warning", {
  expect_warning(f_selec_mod_metapsg_gv(small, ref$n_s, alpha = 0.5, same_var = TRUE, npts = 2),
                 "same_var")
})

test_that("the parallel grid search works with a single value of alpha", {
  skip_on_cran()
  cfg <- ref$configs[[1]]
  res <- f_selec_mod_metapsg_gv_parallel(ref$sumstat, ref$n_s, alpha = cfg$alpha, adaptive = cfg$adaptive,
                                         eps = cfg$eps, npts = cfg$npts, ncpu = 1)
  best <- which.min(cfg$selmod$MIC)
  expect_equal(as.numeric(res$info_model$MIC), cfg$selmod$MIC[best], tolerance = 1e-6)
  expect_equal(unname(unlist(res$model)), unname(unlist(cfg$beta_select[[best]])), tolerance = 1e-6)
})
