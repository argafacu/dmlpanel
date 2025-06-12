# dmlpanel
DML for High-Dimensional Linear Panel Data Models in R. 

## Introduction
This is an implementation of Debiased Machine Learning (DML) for high-dimensional linear panel data models with common parameters and random coefficients, as introduced by _"Debiased Machine Learning Inference for Functionals of Unoberseved Heterogeneity"_, by Argañaraz, F., and Escanciano, J.C., Argañaraz and Escanciano (2025). Specifically, this function works with the following model:
```math
Y_i = W_i\beta_0 + V_i\alpha_i + \varepsilon_i,
```
where $Y_{i}=(Y_{i1},...,Y_{iT})$ is a $T\times1$ vector of dependent variables, for a finite and fixed $T,$ $X_{i}=(W_{i},V_{i})$ are strictly exogenous covariates, for $T\times p$ and $T\times q$ matrices $W_{i}$ and $V_{i},$ respectively, $\alpha_{i}$ is a $q\times1$ vector of individual heterogeneity (UH) or random coefficients, and $\varepsilon_{i}$ are idiosyncratic (unobserved) errors. In this setting, $p$ is allowed to be large, in particular, $p>>nT$. Thus, $W_i$ might contain a large number of regressors. No assumptions are imposed on the unknown distribution of $\alpha_i$. Note $\alpha_i$ might be a vector, but $rank(V) < T$. 

`dmlpanel` conduct robust inference for the following functionals: 
```math
\psi_0 = \mathbb{E}\left[C^{\prime}_2\alpha\right],
```

```math
\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right],
```
```math
\psi_0 = C^{\prime}_1\beta_0,
```

where $C_1$ and $C_2$ are known $p \times p_1$ and $q \times p_2$ matrices, and $\Omega$ is a known $q \times q$ matrix. Both $C_1$ and $C_2$ must be of full-column rank. 

Note that `dmlpanel` provides point estimates and standard errors for $\psi_0$, when $C_1$ and $C_2$ are provided as inputs. Moreover, the function runs inferences for first moments of $\alpha_i$ and pre-specified coordinates of $\beta_0$ (default). If  $C_1$ and $C_2$ are not provided, `dmlpanel` will do the latter only.  The same applies for $\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right]$, when $\Omega$ is provided. In addition, the function will conduct inference for variances and covariances among the coordinates in $\alpha_i$. When $\Omega$ is not given, `dmlpanel` will do the latter only (default). 

## Some important observations 

+ This function is in a (VERY!) preliminary stage of development. So, be aware of this, if you want to use this function in applications.
+  The procedure provides locally robust inference, robust to $\beta_0$ and $\eta_0$, where $\eta_0$ is the conditional density of $\alpha_i$, which could depend on $X$. Hence, a "fixed effect" approach is undertaken.
+ Estimations are also robust to other nuisance parameters, that appear in the construction of the debiased moments, used in inference. We emphatically recommend reading our paper to see how these moments are constructed and be aware of other important notions.
+ Estimation and inference is always conducted by cross-fitting, using $L$ folds.
+ The high-dimensional parameter $\beta_0$ is estimated by an extension to multidimensional UH of the Lasso procedure developed by  _"Inference in High-Dimensional Panel Models With an Application to Gun Control"_, by Belloni, A., Chernozhukov, V., Hansen, C., and Kozbur, D., Belloni et al. (2016). The regularization parameter (the "lambda") is as recommended by these authors. An important constant appears that scales this regularization parameter. This is `re` in the function. While this should be greater than one, `dmlpanel` allows for an user's choice.

## Additional Considerations 

The user needs to provide `indreg`, a natural number. `dmlpanel` will report point estimates and standard errors for the common parameters corresponding to `W_i[,1:indreg]`. Hence, for now, you will need to arrange your regressors accordingly. Clearly, `indreg` should be "small", relative to $p$. 

Following _"Identifying Distributional Characteristics in Random Coefficients Panel Data Models"_, by Arellano, M, and Bonhomme, S. (Arellano and Bonhomme (2012)), to identify $\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right]$, the conditional variance-covariance matrix of $\varepsilon_i$ is restricted. Let this matrix be $\Sigma$, that is  $\Sigma_i:=\mathbb{V}\left[\varepsilon_i| X_i\right]$. it is assumed that 
```math
vec(\Sigma_i)=S_2\omega_{0i},
```
where $vec$ is the vectorization operation, $S_2$ is a known selection matrix and $\omega_{0i}$ is an $m$ dimensional vector of parameters, possibly depending on $X_i$. `dmlpanel` allows the user to specify $S_2$. If this is not the case, iid homoscedastic errors are assumed (default). 



Inputs are:
```
  #y: An nT x 1 vector of outcome variables. It must be provided.
  #v: An nT x q matrix of regressors with "random" coefficients. It must be provided.
  #w: An nT x p matrix of regressors with common coefficients. It must be provided.
  #id: An nT x 1 vector of repeated id's. It must be provided.
  #C1: A p x p1 matrix. Defult is NULL.
  #C2: A q x p2 matrix.D efult is NULL.
  #Omega: A q x q matrix.
  #S2: A selection matrix for modeling vec(Sigma) (that is, the var-cov matrix of idiosyncratic errors). The number of columns of this matrix determines the dimension of \omega_0 and thus the degree of parametrization of the variance of errors. If this matrix is not provided, iid errors are assumed (default)
  #ToT: Number of "periods" in the panel (i.e., number of times a unit in the cross-section dimension is observed). It must be provided.
  #L: number of folds (scalar). Default is 5.
  #re: regularization constant (scalar). Default is 1.1.
  #indreg: A natural number indicating inference will be conducted for the first indreg regressors in W.It must be provided.
```

Output is: 

A list with four elements: 
1. FM: A matrix of dimension $(p_1 + q) \times 2$, containing point estimates and standard errors of $\psi_0 = \mathbb{E}\left[C^{\prime}_1\alpha\right]$ and first moments of UH.
2. SM: A matrix of dimension $(1 + q) \times 2$, containing point estimates and standard errors of $\psi_0 = \mathbb{E}\left[\alpha^{\prime}\Omega \alpha\right]$ and variances of UH.
3. Cov: A list containing all possible covairances among the entries in $\alpha_i$ along their standard errors.
4. CP: A matrix of dimension $indreg \times 2$, containg point estimates and standard errors in  `\beta_0[1:indreg]`.
