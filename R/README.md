# MPSG R codes
In this section, you can find the R codes to run MPSG approach (Meta-analysis model adapted for Pleiotropy Selection with Group structure)

Created by: Pierre-Emmanuel Sugier<br>
Creation date: 25 July 2024<br>
Update: July 2024<br>
<br>

## Table of Contents
- [metapsg_gv_solver.r](#main-function)
- [metapsg_parallel.r](#parallelization)
- [metapsg_sub.r](#sub)
- [RcppExports.r](#Rcpp)

## main function

The R code for the main function metapsg_gv() is 
metapsg_gv_solver.r contains the R code for the main function of the mpsg package: metapsg_gv()

## parallelization
To run mpsg approach by using parallelization on a grid of lambda values generated automatically considering a specified alpha value, you can run the f_selec_mod_metapsg_gv_parallel() function of metapsg_parallel.r

## sub
metapsg_sub.r contains some subfunctions called by metapsg_gv_solver.r code to be able to run the metapsg_gv() function

## Rcpp
RcppExports.r to run external functions for computational optimisation
