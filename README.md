
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dmlpanel

<!-- badges: start -->

<!-- badges: end -->

DML for High-Dimensional Linear Panel Data Models in R.

## Introduction

This is an implementation of Debiased Machine Learning (DML) for
high-dimensional linear panel data models with common parameters and
random coefficients, as introduced by *“Debiased Machine Learning
Inference for Functionals of Unoberseved Heterogeneity”*, by Argañaraz,
F., and Escanciano, J.C., Argañaraz and Escanciano (2025). Specifically,
this function works with the following model:

``` math
Y_i = W_i\beta_0 + V_i\alpha_i + \varepsilon_i,
```

where $`Y_{i}=(Y_{i1},...,Y_{iT})`$ is a $`T\times1`$ vector of
dependent variables, for a finite and fixed $`T,`$
$`X_{i}=(W_{i},V_{i})`$ are strictly exogenous covariates, for
$`T\times p`$ and $`T\times q`$ matrices $`W_{i}`$ and $`V_{i},`$
respectively, $`\alpha_{i}`$ is a $`q\times1`$ vector of individual
heterogeneity (UH) or random coefficients, and $`\varepsilon_{i}`$ are
idiosyncratic (unobserved) errors. In this setting, $`p`$ is allowed to
be large, in particular, $`p>>nT`$. Thus, $`W_i`$ might contain a large
number of regressors. No assumptions are imposed on the unknown
distribution of $`\alpha_i`$. Note $`\alpha_i`$ might be a vector, but
$`rank(V) < T`$.

`dmlpanel` conduct robust inference for the following functionals:

``` math
\psi_0 = \mathbb{E}\left[C^{\prime}_2\alpha\right],
```

``` math
\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right],
```

``` math
\psi_0 = C^{\prime}_1\beta_0,
```

where $`C_1`$ and $`C_2`$ are known $`p \times p_1`$ and
$`q \times p_2`$ matrices, and $`\Omega`$ is a known $`q \times q`$
matrix. Both $`C_1`$ and $`C_2`$ must be of full-column rank.

Note that `dmlpanel` provides point estimates and standard errors for
$`\psi_0`$, when $`C_1`$ and $`C_2`$ are provided as inputs. Moreover,
the function runs inferences for first moments of $`\alpha_i`$ and
pre-specified coordinates of $`\beta_0`$ (default). If $`C_1`$ and
$`C_2`$ are not provided, `dmlpanel` will do the latter only. The same
applies for
$`\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right]`$, when
$`\Omega`$ is provided. In addition, the function will conduct inference
for variances and covariances among the coordinates in $`\alpha_i`$.
When $`\Omega`$ is not given, `dmlpanel` will do the latter only
(default).

## Some important observations

- This function is in a (VERY!) preliminary stage of development. So, be
  aware of this, if you want to use this function in applications.
- The procedure provides locally robust inference, robust to $`\beta_0`$
  and $`\eta_0`$, where $`\eta_0`$ is the conditional density of
  $`\alpha_i`$, which could depend on $`X`$. Hence, a “fixed effect”
  approach is undertaken.
- Estimations are also robust to other nuisance parameters, that appear
  in the construction of the debiased moments, used in inference. We
  emphatically recommend reading our paper to see how these moments are
  constructed and be aware of other important notions.
- Estimation and inference is always conducted by cross-fitting, using
  $`L`$ folds.
- The high-dimensional parameter $`\beta_0`$ is estimated by an
  extension to multidimensional UH of the Lasso procedure developed by
  *“Inference in High-Dimensional Panel Models With an Application to
  Gun Control”*, by Belloni, A., Chernozhukov, V., Hansen, C., and
  Kozbur, D., Belloni et al. (2016). The regularization parameter (the
  “lambda”) is as recommended by these authors. An important constant
  appears that scales this regularization parameter. This is `re` in the
  function. While this should be greater than one, `dmlpanel` allows for
  an user’s choice.

## Additional Considerations

The user needs to provide `indreg`, a natural number. `dmlpanel` will
report point estimates and standard errors for the common parameters
corresponding to `W_i[,1:indreg]`. Hence, for now, you will need to
arrange your regressors accordingly. Clearly, `indreg` should be
“small”, relative to $`p`$.

Following *“Identifying Distributional Characteristics in Random
Coefficients Panel Data Models”*, by Arellano, M, and Bonhomme, S.
(Arellano and Bonhomme (2012)), to identify
$`\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right]`$, the
conditional variance-covariance matrix of $`\varepsilon_i`$ is
restricted. Let this matrix be $`\Sigma`$, that is
$`\Sigma_i:=\mathbb{V}\left[\varepsilon_i| X_i\right]`$. it is assumed
that

``` math
vec(\Sigma_i)=S_2\omega_{0i},
```

where $`vec`$ is the vectorization operation, $`S_2`$ is a known
selection matrix and $`\omega_{0i}`$ is an $`m`$ dimensional vector of
parameters, possibly depending on $`X_i`$. `dmlpanel` allows the user to
specify $`S_2`$. If this is not the case, iid homoscedastic errors are
assumed (default).

## Installation

You can install the development version of `dmlpanel` from
[GitHub](https://github.com/) with:

``` r
# install pak if not installed
install.packages("pak")
# install dmlpanel from github
pak::pak("argafacu/dmlpanel")
```

## Example

This is a basic example of how to use this function. Suppose that $`V`$
has two columns, that is, we will allow two regressors to have random
coefficients, say $`\alpha_{1i}`$ and $`\alpha_{2i}`$. Let one of these
be a random intercept. We will assume that $`\varepsilon_i`$ are idd and
homoscedastic. The following example will lead to robust inference for
first moments of the entries of $`\alpha_i`$, its variances, and
covariance. In addition, it will report point estimates and standard
errors for each component in $`\beta_0[1:8]`$.

``` r
library(dmlpanel); library(fda)
#> Loading required package: splines
#> Loading required package: fds
#> Loading required package: rainbow
#> Loading required package: MASS
#> Loading required package: pcaPP
#> Loading required package: RCurl
#> Loading required package: deSolve
#> 
#> Attaching package: 'fda'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
# Upload the data ---------------------------------------------------------

data("birthpanel", package = "dmlpanel")
df <- birthpanel
TotT <- 3 #There are three periods in this data


ID <- df$id 
unique_IDs <- unique(ID)
ID2 <- match(ID, unique_IDs)
ID <- ID2 #Repeated id's

y <- df$birthwght #outcome variable 

z <- cbind(df$male, df$age, df$agesq, df$kessner2, df$kessner3, 
           df$novisit, df$visit2, df$visit3)

v <- cbind(df$smoke, matrix(1, nrow = nrow(z), ncol = 1)) #regressors with random coefficents. Note that we allow for a random intercept. 


basis <- create.bspline.basis(rangeval=c(min(df$age),max(df$age)), nbasis=5, norder=4,
                              breaks=NULL, dropind=NULL, quadvals=NULL, values=NULL,
                              basisvalues=NULL)


agesplines <- eval.basis(as.vector(df$age), basis, Lfdobj=0, returnMatrix=FALSE) #Create a flexible basis for age

zexpand <- cbind(z, agesplines)
w <- zexpand
w <- round(w,digits =2) #Regressors with common coefficients

Results <- dmlpanel(y=y,v=v,w=w,id=ID,C1=NULL,C2=NULL, Omega=NULL, 
                    S2=NULL, TotT=TotT,  L = 5, re=0.1, indreg=8)

Results$FM #View inferences for first moments of alpha
#>        Est.   S.e.
#> V.1 -160.82  16.97
#> V.2 3501.61 422.94
Results$SM #View inferences for variances of alpha
#>         Est.      S.e.
#> V.1 108604.8  21146.55
#> V.2 122731.4 133842.33
Results$Cov #View inference for covariance among the components of alpha
#> $v1v2
#> $v1v2$Est.
#> [1] -56924.83
#> 
#> $v1v2$S.e.
#> [1] 14571.58
Results$CP #View inferences for beta_0[1:8]
#>       Est.  S.e.
#> W.1 132.40 23.10
#> W.2  11.78 32.67
#> W.3  -0.22  0.57
#> W.4 -11.11 17.21
#> W.5 -45.30 39.55
#> W.6  -1.07  6.92
#> W.7   3.48 19.74
#> W.8  -0.75 23.65
```

Suppose instead that we want to model the conditional var-cov of errors
differently. We could

``` r
S2 <- matrix(c(
  1, 0, 0, 0, 1, 0, 0, 0, 1,
  0, 0, 0, 0, 1, 0, 0, 0, 2
), nrow = 2, byrow = TRUE)

S2 <- t(S2)

Results <- dmlpanel(y=y,v=v,w=w,id=ID,C1=NULL,C2=NULL, Omega=NULL, 
                    S2=S2, TotT=TotT,  L = 5, re=0.1, indreg=8)

Results$FM #View inferences for first moments of alpha
#>        Est.   S.e.
#> V.1 -161.29  16.96
#> V.2 3409.05 418.97
Results$SM #View inferences for variances of alpha
#>         Est.     S.e.
#> V.1 108907.0 27564.13
#> V.2 161808.4 67554.14
Results$Cov #View inference for covariance among the components of alpha
#> $v1v2
#> $v1v2$Est.
#> [1] -63362.88
#> 
#> $v1v2$S.e.
#> [1] 19921.58
Results$CP #View inferences for beta_0[1:8]
#>       Est.  S.e.
#> W.1 130.04 23.00
#> W.2   6.16 31.90
#> W.3  -0.08  0.56
#> W.4  -6.86 17.36
#> W.5 -45.43 40.42
#> W.6  -1.71  7.20
#> W.7   7.85 19.52
#> W.8   4.79 23.60
```

The idea is the same for user specified matrices $`C_1`$ and $`C_2`$.
For example,

``` r
C2 <- matrix(0,ncol(v),1)
C2[1,1] <- 1
C2[2,1] <- 2

Results <- dmlpanel(y=y,v=v,w=w,id=ID,C1=NULL,C2=C2, Omega=NULL, 
                    S2=NULL, TotT=TotT,  L = 5, re=0.1, indreg=8)

Results$FM #View inferences for E[C_2'alpha] first moments of alpha
#>         Est.   S.e.
#> C2.1 6722.26 843.80
#> V.1  -161.31  16.90
#> V.2  3441.79 420.44
Results$SM #View inferences for variances of alpha
#>          Est.     S.e.
#> V.1 107168.02 21060.29
#> V.2  66442.55 92892.95
Results$Cov #View inference for covariance among the components of alpha
#> $v1v2
#> $v1v2$Est.
#> [1] -54206
#> 
#> $v1v2$S.e.
#> [1] 14721.84
Results$CP #View inferences for beta_0[1:8]
#>       Est.  S.e.
#> W.1 127.78 23.26
#> W.2   5.95 32.39
#> W.3  -0.05  0.58
#> W.4 -20.49 17.67
#> W.5 -53.02 41.24
#> W.6  -2.76  7.10
#> W.7   8.98 19.68
#> W.8   4.43 23.92
```
