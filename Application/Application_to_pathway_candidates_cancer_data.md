``` r
library(mvtnorm)
source("Function_PLS1.R")
```

## Generate some data for PLS1

``` r
set.seed(1)
beta <- c(5,0)
n <- 100
X <- mvtnorm::rmvnorm(n,mean=c(1,1),sigma=matrix(c(3,1,1,2),ncol=2,nrow=2,byrow=T))
Y <- X%*%matrix(beta,ncol=1)+rnorm(n,sd=1)
```

## Loss and gradient path for BSS-PLS1

-   Lambda is 0

``` r
t0 <- rep(0.5,2)
plot.landscape.path.BSS(lambda=0,t0=t0)
```
![](landscape_lam0.png)

-   Lambda is 150

``` r
t0 <- rep(0.5,2)
plot.landscape.path.BSS(lambda=150,t0=t0)
#rgl.snapshot('landscape_lam150.png', fmt = 'png')
```
![](landscape_lam150.png)

-   Lambda is 450

``` r
t0 <- rep(0.5,2)
plot.landscape.path.BSS(lambda=450,t0=t0)
#rgl.snapshot('landscape_lam450.png', fmt = 'png')
```
![](landscape_lam450.png)
