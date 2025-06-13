
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
# install devtools if not installed
install.packages("devtools")
# install dmlpanel from github
devtools::install_github("argafacu/dmlpanel")
```

## Example

This is a basic example which shows you how to solve a common problem
such as

``` r
library(dmlpanel)
## basic example code
```
