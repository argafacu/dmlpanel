#' DML for High-Dimensional Linear Panel Data Models in R
#'
#' `dmlpanel()` implements Debiased Machine Learning (DML) for high-dimensional linear panel data models with common parameters and random coefficients, as introduced by "Debiased Machine Learning Inference for Functionals of Unoberseved Heterogeneity", by Argañaraz, F., and Escanciano, J.C., Argañaraz and Escanciano (2025).

#' @param y    An nT x 1 vector of outcome variables. It must be provided.
#' @param v    An nT x q matrix of regressors with "random" coefficients. It must be provided.
#' @param w    An nT x p matrix of regressors with common coefficients. It must be provided.
#' @param id   An nT x 1 vector of repeated id's. It must be provided.
#' @param C1   A p x p1 matrix. Defult is NULL.
#' @param C2   A q x p2 matrix.D efult is NULL.
#' @param Omega A q x q matrix.
#' @param S2   A selection matrix for modeling vec(Sigma) (that is, the var-cov matrix of idiosyncratic errors). The number of columns of this matrix determines the dimension of omega_0 and thus the degree of parametrization of the variance of errors. If this matrix is not provided, iid errors are assumed (default).
#' @param TotT Number of "periods" in the panel (i.e., number of times a unit in the cross-section dimension is observed). It must be provided.
#' @param L    number of folds (scalar). Default is 5.
#' @param re   regularization constant (scalar). Default is 1.1.
#' @param indreg A natural number indicating inference will be conducted for the first indreg regressors in W.It must be provided.
#' @return     A list with parameter estimates and standard errors.
#' @importFrom MASS ginv
#' @export

dmlpanel <- function(y,v,w,id,C1=NULL,C2=NULL, Omega=NULL, S2=NULL, TotT,  L = 5, re=1.1, indreg){


n <- length(unique(id)) # Identify the number of cross-sectional units
foldid <- rep.int(1:L, times = ceiling(n/L))[sample.int(n)] #Randomly split the data
I <- split(1:n, foldid) #Collect folds
seqn <- seq(1:n)

if(!is.null(C2)){

if(qr(C2)$rank != ncol(C2)) {
    stop("Error: The rank of C2 is not equal to its number of columns.")
}

Resultsfm1 <- matrix(NA,ncol(C2),2)


ymat <- vector("list", length = n)
Wmat <- vector("list", length = n)
Qmat <- vector("list", length = n)
Hmat <- vector("list", length = n)
Mmat <- list()
Smat <- list()
Gammamat <- list()


  for(l in 1:length(I)){
    WQWmat <- matrix(0,ncol(w), ncol(w))
    HWmat <- matrix(0,ncol(v), ncol(w))

    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    for(j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
      HWmat <- HWmat + h2%*%w2

      ymat[[j1]] <- y2
      Wmat[[j1]] <- w2
      Qmat[[j1]] <- q2
      Hmat[[j1]] <- h2
    }

    Mmat[[l]] <- WQWmat/length(IDs_not_fold)
    Smat[[l]] <- HWmat/length(IDs_not_fold)


    decompM <- eigen(Mmat[[l]])
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M

    Gammamat[[l]] <- -t(C2)%*%Smat[[l]]%*%matMinv

  }


betalasso <- vector("list", L)

    for(l in 1:length(I)){

      IDs_not_fold <- seqn[!seqn%in%I[[l]]]
      n_cl <- length(IDs_not_fold)

      gamma <- 0.1/(log(max(ncol(w),n_cl*TotT)))

      cn <- (re/(sqrt(n_cl*TotT)))*stats::qnorm(1-gamma/(2*ncol(w))) #Regularization parameter from Belloni et al. (2016)


      covW <- Reduce("+", lapply(IDs_not_fold, function(i) {
        t(Wmat[[i]])%*%Qmat[[i]]%*%Wmat[[i]]
      }))

      WY <- Reduce("+", lapply(IDs_not_fold, function(i) {
        t(Wmat[[i]])%*%Qmat[[i]]%*%ymat[[i]]
      }))

      lowcov <- Reduce("+", lapply(IDs_not_fold, function(i) {
        t(Wmat[[i]][,1:(floor(ncol(w)/40))])%*%Qmat[[i]]%*%Wmat[[i]][,1:(floor(ncol(w)/40))]
      }))

      lowWY <- Reduce("+", lapply(IDs_not_fold, function(i) {
        t(Wmat[[i]][,1:(floor(ncol(w)/40))])%*%Qmat[[i]]%*%ymat[[i]]
      }))



      intbeta <- numeric(ncol(w))

      intbeta[1:ncol(lowcov)] <- solve(lowcov)%*%lowWY


      for(k in 1:ncol(w)){

        sumphiint <- 0

        for (i in IDs_not_fold){

          W_tilde <- Qmat[[i]]%*%Wmat[[i]]
          Y_tilde <- Qmat[[i]]%*%ymat[[i]]

          W_k <- W_tilde[,k]


          sumphiint <-   sumphiint +  sum(outer(W_k, W_k)*outer(Y_tilde[,1], Y_tilde[,1])) #Compute the sum of all possible inner products
        }



        phiintk <- sqrt(sumphiint/(n_cl*TotT))

        Ak <- WY[k]/(n_cl*TotT) - (t(intbeta[-k])%*%covW[-k,k])/(n_cl*TotT)

        Bk <- covW[k,k]/(n_cl*TotT)

        if(Ak < -phiintk*cn){
          intbeta[k] <- (Ak+phiintk*cn)/Bk
        }else{
          if(Ak > phiintk*cn){
            intbeta[k] <- (Ak-phiintk*cn)/Bk
          }else{
            intbeta[k] <- 0
          }
        }

      }



      beta <- intbeta   #Compute an initial beta

      rp <- 1

      while(rp <=10){

        for(k in 1:ncol(w)){
          sumphi <- 0

          for (i in IDs_not_fold){

            W_tilde <- Qmat[[i]]%*%Wmat[[i]]
            Y_tilde <- Qmat[[i]]%*%ymat[[i]]

            eps_res <-   Y_tilde -  W_tilde%*%beta

            W_k <- W_tilde[,k]


            sumphi <-    sumphi  +  sum(outer(W_k, W_k)*outer(eps_res[,1], eps_res[,1]))    #Compute the sum of all possible inner products
          }

          phik <- sqrt(sumphi/(n_cl*TotT))  #Compute refined loadings

          Ak <- WY[k]/(n_cl*TotT) - (t(beta[-k])%*%covW[-k,k])/(n_cl*TotT)

          Bk <- covW[k,k]/(n_cl*TotT)

          if(Ak < -phik*cn){
            beta[k] <- (Ak+phik*cn)/Bk
          }else{
            if(Ak > phik*cn){
              beta[k] <- (Ak-phik*cn)/Bk
            }else{
              beta[k] <- 0
            }
          }

        }
        rp <- rp + 1
      }

      betalasso[[l]] <- beta #Estimated beta by Lasso and cross-fitting
    }


psifm0 <- 0

gfmlr <- matrix(NA,ncol(C2),n)


for(l in 1:length(I)){
  for(j1 in I[[l]]){

    gfmlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifm0

  }

}

meangfmlr <- rowMeans(gfmlr)


psifmest <- vector("list", length(I))

for(l in 1:length(I)){
  IDs_not_fold <- seqn[!seqn%in%I[[l]]]
  inneR <- 0

  for(j1 in IDs_not_fold){
    inneR <- inneR + (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]])
  }

  psifmest[[l]] <- inneR/length(IDs_not_fold) #Estimated psi by cross-fitting

}



gfmestlr <- matrix(NA,ncol(C2),n)

for(l in 1:length(I)){

  for(j1 in I[[l]]){
    gfmestlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifmest[[l]]
  }
}

meangfmestlr <- rowMeans(gfmestlr)

secondmeangfmestlr <- gfmestlr%*%t(gfmestlr)/n

varest <- secondmeangfmestlr -  meangfmestlr%*%t(meangfmestlr) #Variance estimator of LR moment for first moment of UH

seest <- sqrt(diag(varest))/sqrt(n)


Resultsfm1[, 1] <- round(meangfmlr,2) #Matrix of results
Resultsfm1[, 2] <- round(seest,2) #Matrix of results

colnames(Resultsfm1) <- c("Est.", "S.e.")
rownames(Resultsfm1) <- paste0("C2.", 1:ncol(C2))




  psifmestmatcf <- matrix(NA,ncol(v),L)
  C2 <- diag(ncol(v))

  ymat <- vector("list", length = n)
  Wmat <- vector("list", length = n)
  Qmat <- vector("list", length = n)
  Hmat <- vector("list", length = n)
  Mmat <- list()
  Smat <- list()
  Gammamat <- list()


  for(l in 1:length(I)){
    WQWmat <- matrix(0,ncol(w), ncol(w))
    HWmat <- matrix(0,ncol(v), ncol(w))

    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    for(j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
      HWmat <- HWmat + h2%*%w2

      ymat[[j1]] <- y2
      Wmat[[j1]] <- w2
      Qmat[[j1]] <- q2
      Hmat[[j1]] <- h2
    }

    Mmat[[l]] <- WQWmat/length(IDs_not_fold)
    Smat[[l]] <- HWmat/length(IDs_not_fold)


    decompM <- eigen(Mmat[[l]])
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M

    Gammamat[[l]] <- -t(C2)%*%Smat[[l]]%*%matMinv

  }


  psifm0 <- 0
  gfmlr <- matrix(NA,ncol(C2),n)


  for(l in 1:length(I)){
    for(j1 in I[[l]]){

      gfmlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifm0

    }

  }

  meangfmlr <- rowMeans(gfmlr)


  psifmest <- vector("list", length(I))

  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]
    inneR <- 0

    for(j1 in IDs_not_fold){
      inneR <- inneR + (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]])
    }

    psifmest[[l]] <- inneR/length(IDs_not_fold) #Estimated psi by cross-fitting

    psifmestmatcf[,l] <- psifmest[[l]]
  }




  gfmestlr <- matrix(NA,ncol(C2),n)

  for(l in 1:length(I)){

    for(j1 in I[[l]]){
      gfmestlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifmest[[l]]
    }
  }

  meangfmestlr <- rowMeans(gfmestlr)

  secondmeangfmestlr <- gfmestlr%*%t(gfmestlr)/n

  varest <- secondmeangfmestlr -  meangfmestlr%*%t(meangfmestlr) #Variance estimator of LR moment for first moment of UH

  seest <- sqrt(diag(varest))/sqrt(n)

  Resultsfm2 <- matrix(NA,ncol(v),2)
  Resultsfm2[,1] <- round(meangfmlr,2)
  Resultsfm2[,2] <- round(seest,2)

  colnames(Resultsfm2) <- c("Est.", "S.e.")
  rownames(Resultsfm2) <- paste0("V.", 1:ncol(v))

  Resultsfm <- rbind(Resultsfm1,Resultsfm2) #Matrix of results



}else{

  psifmestmatcf <- matrix(NA,ncol(v),L)
  C2 <- diag(ncol(v))
  Resultsfm2 <- matrix(NA,ncol(C2),2)
  ymat <- vector("list", length = n)
  Wmat <- vector("list", length = n)
  Qmat <- vector("list", length = n)
  Hmat <- vector("list", length = n)
  Mmat <- list()
  Smat <- list()
  Gammamat <- list()


  for(l in 1:length(I)){
    WQWmat <- matrix(0,ncol(w), ncol(w))
    HWmat <- matrix(0,ncol(v), ncol(w))

    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    for(j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
      HWmat <- HWmat + h2%*%w2

      ymat[[j1]] <- y2
      Wmat[[j1]] <- w2
      Qmat[[j1]] <- q2
      Hmat[[j1]] <- h2
    }

    Mmat[[l]] <- WQWmat/length(IDs_not_fold)
    Smat[[l]] <- HWmat/length(IDs_not_fold)


    decompM <- eigen(Mmat[[l]])
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M

    Gammamat[[l]] <- -t(C2)%*%Smat[[l]]%*%matMinv

  }


  betalasso <- vector("list", L)

  for(l in 1:length(I)){

    IDs_not_fold <- seqn[!seqn%in%I[[l]]]
    n_cl <- length(IDs_not_fold)

    gamma <- 0.1/(log(max(ncol(w),n_cl*TotT)))

    cn <- (re/(sqrt(n_cl*TotT)))*stats::qnorm(1-gamma/(2*ncol(w))) #Regularization parameter from Belloni et al. (2016)


    covW <- Reduce("+", lapply(IDs_not_fold, function(i) {
      t(Wmat[[i]])%*%Qmat[[i]]%*%Wmat[[i]]
    }))

    WY <- Reduce("+", lapply(IDs_not_fold, function(i) {
      t(Wmat[[i]])%*%Qmat[[i]]%*%ymat[[i]]
    }))

    lowcov <- Reduce("+", lapply(IDs_not_fold, function(i) {
      t(Wmat[[i]][,1:(floor(ncol(w)/40))])%*%Qmat[[i]]%*%Wmat[[i]][,1:(floor(ncol(w)/40))]
    }))

    lowWY <- Reduce("+", lapply(IDs_not_fold, function(i) {
      t(Wmat[[i]][,1:(floor(ncol(w)/40))])%*%Qmat[[i]]%*%ymat[[i]]
    }))



    intbeta <- numeric(ncol(w))

    intbeta[1:ncol(lowcov)] <- solve(lowcov)%*%lowWY


    for(k in 1:ncol(w)){

      sumphiint <- 0

      for (i in IDs_not_fold){

        W_tilde <- Qmat[[i]]%*%Wmat[[i]]
        Y_tilde <- Qmat[[i]]%*%ymat[[i]]

        W_k <- W_tilde[,k]


        sumphiint <-   sumphiint +  sum(outer(W_k, W_k)*outer(Y_tilde[,1], Y_tilde[,1])) #Compute the sum of all possible inner products
      }



      phiintk <- sqrt(sumphiint/(n_cl*TotT))

      Ak <- WY[k]/(n_cl*TotT) - (t(intbeta[-k])%*%covW[-k,k])/(n_cl*TotT)

      Bk <- covW[k,k]/(n_cl*TotT)

      if(Ak < -phiintk*cn){
        intbeta[k] <- (Ak+phiintk*cn)/Bk
      }else{
        if(Ak > phiintk*cn){
          intbeta[k] <- (Ak-phiintk*cn)/Bk
        }else{
          intbeta[k] <- 0
        }
      }

    }



    beta <- intbeta   #Compute an initial beta

    rp <- 1

    while(rp <=10){

      for(k in 1:ncol(w)){
        sumphi <- 0

        for (i in IDs_not_fold){

          W_tilde <- Qmat[[i]]%*%Wmat[[i]]
          Y_tilde <- Qmat[[i]]%*%ymat[[i]]

          eps_res <-   Y_tilde -  W_tilde%*%beta

          W_k <- W_tilde[,k]


          sumphi <-    sumphi  +  sum(outer(W_k, W_k)*outer(eps_res[,1], eps_res[,1]))    #Compute the sum of all possible inner products
        }

        phik <- sqrt(sumphi/(n_cl*TotT))  #Compute refined loadings

        Ak <- WY[k]/(n_cl*TotT) - (t(beta[-k])%*%covW[-k,k])/(n_cl*TotT)

        Bk <- covW[k,k]/(n_cl*TotT)

        if(Ak < -phik*cn){
          beta[k] <- (Ak+phik*cn)/Bk
        }else{
          if(Ak > phik*cn){
            beta[k] <- (Ak-phik*cn)/Bk
          }else{
            beta[k] <- 0
          }
        }

      }
      rp <- rp + 1
    }

    betalasso[[l]] <- beta #Estimated beta by Lasso and cross-fitting
  }


  psifm0 <- 0
  gfmlr <- matrix(NA,ncol(C2),n)


  for(l in 1:length(I)){
    for(j1 in I[[l]]){

      gfmlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifm0

    }

  }

  meangfmlr <- rowMeans(gfmlr)


  psifmest <- vector("list", length(I))

  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]
    inneR <- 0

    for(j1 in IDs_not_fold){
      inneR <- inneR + (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]])
    }

    psifmest[[l]] <- inneR/length(IDs_not_fold) #Estimated psi by cross-fitting
    psifmestmatcf[,l] <- psifmest[[l]]
  }



  gfmestlr <- matrix(NA,ncol(C2),n)

  for(l in 1:length(I)){

    for(j1 in I[[l]]){
      gfmestlr[,j1] <- (t(C2)%*%Hmat[[j1]] - Gammamat[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]])%*%(ymat[[j1]] -Wmat[[j1]]%*%betalasso[[l]]) - psifmest[[l]]
    }
  }

  meangfmestlr <- rowMeans(gfmestlr)

  secondmeangfmestlr <- gfmestlr%*%t(gfmestlr)/n

  varest <- secondmeangfmestlr -  meangfmestlr%*%t(meangfmestlr) #Variance estimator of LR moment for first moment of UH

  seest <- sqrt(diag(varest))/sqrt(n)


    Resultsfm2[,1] <- round(meangfmlr,2) #Matrix of results
    Resultsfm2[,2] <- round(seest,2)     #Matrix of results

    colnames(Resultsfm2) <- c("Est.", "S.e.")
    rownames(Resultsfm2) <- paste0("V.", 1:ncol(v))

    Resultsfm <- Resultsfm2
    }


if(is.null(S2)){

if(!is.null(Omega)){

Resultssm1 <- matrix(NA,1,2)
Resultssm2 <- matrix(NA,ncol(v),2)

Resultscov <- list()

varerrors <- vector("list", L)

  for(l in 1:length(I)){
    sigma2u2 <- 0
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    for (j1 in IDs_not_fold) {
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)


      sigma2u2 <- sigma2u2 + t(y2 - w2%*%betalasso[[l]])%*%q2%*%(y2 - w2%*%betalasso[[l]])/(TotT - ncol(v)) #Variance of errors, iid

    }

    varerrors[[l]] <-  sigma2u2/length(IDs_not_fold) #Estimator of the variance of errors using cross-fitting from Arellano and Bonhomme (2012)

  }


Gammaomegamat <- vector("list",L)
Gammabetamat <- vector("list",L)


for(l in 1:length(I)){
  IDs_not_fold <- seqn[!seqn%in%I[[l]]]


  A2 <-  0
  B2 <- matrix(0,TotT^2,1)


  for (j1 in IDs_not_fold){
    y2 <- y[id == j1]
    w2 <- w[id == j1, , drop = FALSE]
    v2 <- v[id == j1, , drop = FALSE]
    q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
    h2 <- MASS::ginv(v2)

    Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
    H2 <- h2%x%h2

    A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%cbind(c(diag(TotT)))*1
    B2 <- B2 + Q2%*%cbind(c(diag(TotT)))*1

  }

  A2 <- A2/length(IDs_not_fold)
  B2 <- B2/length(IDs_not_fold)


  #Compute inverses
  decompB <- svd(B2)
  U <-  decompB$u
  D <-  decompB$d
  V <-  decompB$v


  diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

  matBinv <-  V%*%diagmatinvB%*%t(U)

  Gammaomegamat[[l]] <- A2%*%matBinv

}


for(l in 1:length(I)){
  IDs_not_fold <- seqn[!seqn%in%I[[l]]]

  M2 <- matrix(0,ncol(w), ncol(w))
  L2 <- matrix(0,1,ncol(w))

  for (j1 in IDs_not_fold){
    y2 <- y[id == j1]
    w2 <- w[id == j1, , drop = FALSE]
    v2 <- v[id == j1, , drop = FALSE]
    q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
    h2 <- MASS::ginv(v2)

    Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
    H2 <- h2%x%h2

    M2 <- M2 + t(w2)%*%q2%*%w2


    L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

  }

  M2 <- M2/length(IDs_not_fold)


  decompM <- eigen(M2)
  diagmat <- Re(decompM$values)
  U <- decompM$vectors

  diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

  matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


  L2 <- -L2/length(IDs_not_fold)

  Gammabetamat[[l]] <- L2%*%matMinv
}




gsmlr <- vector("list",n)


psism0 <- 0


for(l in 1:length(I)){

  for (j1 in I[[l]]){
    y2 <- y[id == j1]
    w2 <- w[id == j1, , drop = FALSE]
    v2 <- v[id == j1, , drop = FALSE]
    q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
    h2 <- MASS::ginv(v2)

    Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
    H2 <- h2%x%h2

    U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])

    gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

  }
}


psiseest  <- vector("list", L)

for(l in 1:length(I)){
  funcg0  <- 0

  IDs_not_fold <- seqn[!seqn%in%I[[l]]]

  for (j1 in  IDs_not_fold){
    y2 <- y[id == j1]
    w2 <- w[id == j1, , drop = FALSE]
    v2 <- v[id == j1, , drop = FALSE]
    q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
    h2 <- MASS::ginv(v2)

    Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
    H2 <- h2%x%h2

    U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


    funcg0 <- funcg0 + (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])

    }

  psiseest[[l]] <- funcg0/length(IDs_not_fold)
}


gsmestlr  <- vector("list",n)

for(l in 1:length(I)){

  for (j1 in I[[l]]){
    y2 <- y[id == j1]
    w2 <- w[id == j1, , drop = FALSE]
    v2 <- v[id == j1, , drop = FALSE]
    q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
    h2 <- MASS::ginv(v2)

    Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
    H2 <- h2%x%h2

    U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


    gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) -   psiseest[[l]]
  }
}


gsmestlr <- unlist(gsmestlr)

varest <- mean((gsmestlr  - mean(gsmestlr))^2)

seest <- sqrt(varest)/sqrt(n)

Resultssm1[1,1] <- round(mean(unlist(gsmlr)),2) #Matrix of results
Resultssm1[1,2] <- round(seest,2)     #Matrix of results

colnames(Resultssm1) <- c("Est.", "S.e.")
rownames(Resultssm1) <- "Omega.1"


for(uh in 1:ncol(v)){
  Omega <- matrix(0, ncol(v), ncol(v))
  Omega[uh,uh] <- 1


  Gammaomegamat <- vector("list",L)
  Gammabetamat <- vector("list",L)


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]


    A2 <-  0
    B2 <- matrix(0,TotT^2,1)


    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%cbind(c(diag(TotT)))*1
      B2 <- B2 + Q2%*%cbind(c(diag(TotT)))*1

    }

    A2 <- A2/length(IDs_not_fold)
    B2 <- B2/length(IDs_not_fold)


    #Compute inverses
    decompB <- svd(B2)
    U <-  decompB$u
    D <-  decompB$d
    V <-  decompB$v


    diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

    matBinv <-  V%*%diagmatinvB%*%t(U)

    Gammaomegamat[[l]] <- A2%*%matBinv
  }


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    M2 <- matrix(0,ncol(w), ncol(w))
    L2 <- matrix(0,1,ncol(w))

    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      M2 <- M2 + t(w2)%*%q2%*%w2


      L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

    }

    M2 <- M2/length(IDs_not_fold)


    decompM <- eigen(M2)
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U) #Compute suitable estimator of the Moore-Penrose inverse of M


    L2 <- -L2/length(IDs_not_fold)

    Gammabetamat[[l]] <- L2%*%matMinv
  }


  gsmlr <- vector("list",n)
  psism0 <- 0


  for(l in 1:length(I)){

    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

    }
  }


  gsmestlr <- vector("list",n)

  for(l in 1:length(I)){

    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
    }
  }

  gfmestlrv <- vector("list",n)
  for(l in 1:length(I)){
    psifmestvect <- psifmest[[l]]
    for(j1 in I[[l]]){
      gfmestlrv[[j1]] <- gfmestlr[uh,j1] + psifmestvect[uh]
    }
  }


    gsmestlr <- unlist(gsmestlr)
    gfmestlrv <- unlist(gfmestlrv)

    psivar0 <- 0

    est <- mean(unlist(gsmlr)) - Resultsfm2[uh,1]^2


  var <- var(gsmestlr) - 2*mean(gfmestlrv)*stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*(stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*var(gfmestlrv))

  seest <- sqrt(var)/sqrt(n)

  Resultssm2[uh,1] <- round(est,2) #Matrix of results
  Resultssm2[uh,2] <- round(seest,2)     #Matrix of results

  }

colnames(Resultssm2) <- c("Est.", "se")
rownames(Resultssm2) <- paste0("V.", 1:ncol(v))


Resultssm <- rbind(Resultssm1,Resultssm2)


#Covariances
pairs <- utils::combn(ncol(v), 2)

for (uh in 1:ncol(pairs)){
  i <- pairs[1, uh]
  j <- pairs[2, uh]
  Omega <- matrix(0, ncol(v), ncol(v))
  Omega[i, j] <- 1/2
  Omega[j, i] <- 1/2


  Gammaomegamat <- vector("list",L)
  Gammabetamat <- vector("list",L)


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]


    A2 <-  0
    B2 <- matrix(0,TotT^2,1)


    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%cbind(c(diag(TotT)))*1
      B2 <- B2 + Q2%*%cbind(c(diag(TotT)))*1

    }

    A2 <- A2/length(IDs_not_fold)
    B2 <- B2/length(IDs_not_fold)


    #Compute inverses
    decompB <- svd(B2)
    U <-  decompB$u
    D <-  decompB$d
    V <-  decompB$v


    diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

    matBinv <-  V%*%diagmatinvB%*%t(U)

    Gammaomegamat[[l]] <- A2%*%matBinv
  }


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    M2 <- matrix(0,ncol(w), ncol(w))
    L2 <- matrix(0,1,ncol(w))

    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      M2 <- M2 + t(w2)%*%q2%*%w2


      L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

    }

    M2 <- M2/length(IDs_not_fold)


    decompM <- eigen(M2)
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U) #Compute suitable estimator of the Moore-Penrose inverse of M


    L2 <- -L2/length(IDs_not_fold)

    Gammabetamat[[l]] <- L2%*%matMinv
  }


  gsmlr <- vector("list",n)
  psism0 <- 0


  for(l in 1:length(I)){

    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

    }
  }


  gsmestlr <- vector("list",n)

  for(l in 1:length(I)){
    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
    }
  }

  gfmestlrv1 <-  gfmestlrv2 <- vector("list",n)
  for(l in 1:length(I)){
    psifmestvect1 <-  psifmestvect2 <- psifmest[[l]]
    for(j1 in I[[l]]){
      gfmestlrv1[[j1]] <- gfmestlr[i,j1] + psifmestvect[i]
      gfmestlrv2[[j1]] <- gfmestlr[j,j1] + psifmestvect[j]
    }
  }



    gsmestlr <- unlist(gsmestlr)
    gfmestlrv1 <- unlist(gfmestlrv1)
    gfmestlrv2 <- unlist(gfmestlrv2)


    psivar0 <- 0

    est <- mean(unlist(gsmlr)) - Resultsfm2[i,1]*Resultsfm2[j,1]


  var <- var(gsmestlr) + (mean(gfmestlrv2)^2)*var(gfmestlrv1) + (mean(gfmestlrv1)^2)*var(gfmestlrv2) - 2*mean(gfmestlrv2)*stats::cov(gsmestlr,gfmestlrv1) -2*mean(gfmestlrv1)*stats::cov(gsmestlr,gfmestlrv2) + 2*mean(gfmestlrv1)*mean(gfmestlrv2)*stats::cov(gfmestlrv1,gfmestlrv2)

  seest <- sqrt(var)/sqrt(n)


  Resultscov[[paste0("v", i, "v", j)]] <- list("Est." = est, "S.e." = seest)

}
}else{

Resultssm <- matrix(NA,ncol(v),2)
Resultscov <- list()

for(uh in 1:ncol(v)){
  Omega <- matrix(0, ncol(v), ncol(v))
  Omega[uh,uh] <- 1

  if(uh==1){

    varerrors <- vector("list", L)

    for(l in 1:length(I)){
      sigma2u2 <- 0
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      for (j1 in IDs_not_fold) {
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        #Variance of errors, iid
        sigma2u2 <- sigma2u2 + t(y2 - w2%*%betalasso[[l]]) %*% q2 %*% (y2 - w2%*%betalasso[[l]]) / (TotT - ncol(v))   #Get an estimator of the variance of errors using cross-fitting from Arellano and Bonhomme (2012)

      }

      varerrors[[l]] <-  sigma2u2/length(IDs_not_fold)
    }
  }


  Gammaomegamat <- vector("list",L)
  Gammabetamat <- vector("list",L)


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]


    A2 <-  0
    B2 <- matrix(0,TotT^2,1)


    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%cbind(c(diag(TotT)))*1
      B2 <- B2 + Q2%*%cbind(c(diag(TotT)))*1

    }

    A2 <- A2/length(IDs_not_fold)
    B2 <- B2/length(IDs_not_fold)


    #Compute inverses
    decompB <- svd(B2)
    U <-  decompB$u
    D <-  decompB$d
    V <-  decompB$v


    diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

    matBinv <-  V%*%diagmatinvB%*%t(U)

    Gammaomegamat[[l]] <- A2%*%matBinv
  }


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    M2 <- matrix(0,ncol(w), ncol(w))
    L2 <- matrix(0,1,ncol(w))

    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      M2 <- M2 + t(w2)%*%q2%*%w2


      L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

    }

    M2 <- M2/length(IDs_not_fold)


    decompM <- eigen(M2)
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U) #Compute suitable estimator of the Moore-Penrose inverse of M


    L2 <- -L2/length(IDs_not_fold)

    Gammabetamat[[l]] <- L2%*%matMinv
  }


gsmlr <- vector("list",n)
psism0 <- 0


  for(l in 1:length(I)){

    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

    }
  }


  gsmestlr <- vector("list",n)

  for(l in 1:length(I)){
    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmestlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
    }
  }

  gfmestlrv <- vector("list",n)
  for(l in 1:length(I)){
    psifmestvect <- psifmest[[l]]
    for(j1 in I[[l]]){
      gfmestlrv[[j1]] <- gfmestlr[uh,j1] + psifmestvect[uh]
  }
  }


  gsmestlr <- unlist(gsmestlr)
  gfmestlrv <- unlist(gfmestlrv)

  psivar0 <- 0

  est <- mean(unlist(gsmlr)) - Resultsfm2[uh,1]^2


  var <- var(gsmestlr) - 2*mean(gfmestlrv)*stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*(stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*var(gfmestlrv))

  seest <- sqrt(var)/sqrt(n)

  Resultssm[uh,1] <- round(est,2)       #Matrix of results
  Resultssm[uh,2] <- round(seest,2)     #Matrix of results

}

colnames(Resultssm) <- c("Est.", "S.e.")
rownames(Resultssm) <- paste0("V.", 1:ncol(v))


#Covariances
pairs <- utils::combn(ncol(v), 2)

for (uh in 1:ncol(pairs)){
  i <- pairs[1, uh]
  j <- pairs[2, uh]
  Omega <- matrix(0, ncol(v), ncol(v))
  Omega[i, j] <- 1/2
  Omega[j, i] <- 1/2


  Gammaomegamat <- vector("list",L)
  Gammabetamat <- vector("list",L)


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]


    A2 <-  0
    B2 <- matrix(0,TotT^2,1)


    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%cbind(c(diag(TotT)))*1
      B2 <- B2 + Q2%*%cbind(c(diag(TotT)))*1

    }

    A2 <- A2/length(IDs_not_fold)
    B2 <- B2/length(IDs_not_fold)


    #Compute inverses
    decompB <- svd(B2)
    U <-  decompB$u
    D <-  decompB$d
    V <-  decompB$v


    diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

    matBinv <-  V%*%diagmatinvB%*%t(U)

    Gammaomegamat[[l]] <- A2%*%matBinv
  }


  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    M2 <- matrix(0,ncol(w), ncol(w))
    L2 <- matrix(0,1,ncol(w))

    for (j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      M2 <- M2 + t(w2)%*%q2%*%w2


      L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

    }

    M2 <- M2/length(IDs_not_fold)


    decompM <- eigen(M2)
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U) #Compute suitable estimator of the Moore-Penrose inverse of M


    L2 <- -L2/length(IDs_not_fold)

    Gammabetamat[[l]] <- L2%*%matMinv
  }


  gsmlr <- vector("list",n)
  psism0 <- 0


  for(l in 1:length(I)){

    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

    }
  }


  gsmestlr <- vector("list",n)

  for(l in 1:length(I)){
    for (j1 in I[[l]]){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
      h2 <- MASS::ginv(v2)

      Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
      H2 <- h2%x%h2

      U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - cbind(c(diag(TotT)))*as.numeric(varerrors[[l]])


      gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
    }
  }

  gfmestlrv1 <-  gfmestlrv2 <- vector("list",n)
  for(l in 1:length(I)){
    psifmestvect1 <-  psifmestvect2 <- psifmest[[l]]
    for(j1 in I[[l]]){
      gfmestlrv1[[j1]] <- gfmestlr[i,j1] + psifmestvect1[i]
      gfmestlrv2[[j1]] <- gfmestlr[j,j1] + psifmestvect2[j]
    }
  }

    gsmestlr <- unlist(gsmestlr)
    gfmestlrv1 <- unlist(gfmestlrv1)
    gfmestlrv2 <- unlist(gfmestlrv2)


    psivar0 <- 0

    est <- mean(unlist(gsmlr)) - Resultsfm2[i,1]*Resultsfm2[j,1]



  var <- var(gsmestlr) + (mean(gfmestlrv2)^2)*var(gfmestlrv1) + (mean(gfmestlrv1)^2)*var(gfmestlrv2) - 2*mean(gfmestlrv2)*stats::cov(gsmestlr,gfmestlrv1) -2*mean(gfmestlrv1)*stats::cov(gsmestlr,gfmestlrv2) + 2*mean(gfmestlrv1)*mean(gfmestlrv2)*stats::cov(gfmestlrv1,gfmestlrv2)

  seest <- sqrt(var)/sqrt(n)

  Resultscov[[paste0("v", i, "v", j)]] <- list("Est." = est, "S.e." = seest)

}

}

}else{


  if(!is.null(Omega)){

    Resultssm1 <- matrix(NA,1,2)
    Resultssm2 <- matrix(NA,ncol(v),2)

    Resultscov <- list()

    varerrors <- vector("list", L)

    for(l in 1:length(I)){
      varutvL <- matrix(0,TotT^2,1)
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      for (j1 in IDs_not_fold) {
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)


        residual <- y2 - w2 %*%betalasso[[l]]- v2 %*%psifmestmatcf[,l]
        part1 <- S2%*%MASS::ginv((diag(TotT)%x%q2)%*%S2)
        part2 <- residual%x%(q2%*%(y2 - w2 %*% betalasso[[l]]))
        varutvL <- varutvL + part1 %*% part2

      }

      varerrors[[l]] <-  varutvL/length(IDs_not_fold) #Variance of errors, time-varying variances in levels (robust) from Arellano and Bonhomme (2012)

    }


    Gammaomegamat <- vector("list",L)
    Gammabetamat <- vector("list",L)

    omegadot <- diag(ncol(S2))

    for(l in 1:length(I)){
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]


      A2 <-  matrix(0,1,ncol(omegadot))
      B2 <- matrix(0,TotT^2,ncol(S2))


      for (j1 in IDs_not_fold){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
        H2 <- h2%x%h2

        A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%S2%*%omegadot
        B2 <- B2 + Q2%*%S2%*%omegadot

      }

      A2 <- A2/length(IDs_not_fold)
      B2 <- B2/length(IDs_not_fold)


      #Compute inverses
      decompB <- svd(B2)
      U <-  decompB$u
      D <-  decompB$d
      V <-  decompB$v


      diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

      matBinv <-  V%*%diagmatinvB%*%t(U)

      Gammaomegamat[[l]] <- A2%*%matBinv

    }


    for(l in 1:length(I)){
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      M2 <- matrix(0,ncol(w), ncol(w))
      L2 <- matrix(0,1,ncol(w))

      for (j1 in IDs_not_fold){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
        H2 <- h2%x%h2

        M2 <- M2 + t(w2)%*%q2%*%w2


        L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

      }

      M2 <- M2/length(IDs_not_fold)


      decompM <- eigen(M2)
      diagmat <- Re(decompM$values)
      U <- decompM$vectors

      diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

      matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


      L2 <- -L2/length(IDs_not_fold)

      Gammabetamat[[l]] <- L2%*%matMinv
    }




    gsmlr <- vector("list",n)


    psism0 <- 0


    for(l in 1:length(I)){

      for (j1 in I[[l]]){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
        H2 <- h2%x%h2

        U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]

        gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

      }
    }


    psiseest  <- vector("list", L)

    for(l in 1:length(I)){
      funcg0  <- 0

      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      for (j1 in  IDs_not_fold){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
        H2 <- h2%x%h2

        U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


        funcg0 <- funcg0 + (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])

      }

      psiseest[[l]] <- funcg0/length(IDs_not_fold)
    }


    gsmestlr  <- vector("list",n)

    for(l in 1:length(I)){

      for (j1 in I[[l]]){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
        H2 <- h2%x%h2

        U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


        gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) -   psiseest[[l]]
      }
    }


    gsmestlr <- unlist(gsmestlr)

    varest <- mean((gsmestlr  - mean(gsmestlr))^2)

    seest <- sqrt(varest)/sqrt(n)

    Resultssm1[1,1] <- round(mean(unlist(gsmlr)),2) #Matrix of results
    Resultssm1[1,2] <- round(seest,2)     #Matrix of results

    colnames(Resultssm1) <- c("Est.", "S.e.")
    rownames(Resultssm1) <- "Omega.1"


    for(uh in 1:ncol(v)){
      Omega <- matrix(0, ncol(v), ncol(v))
      Omega[uh,uh] <- 1


      Gammaomegamat <- vector("list",L)
      Gammabetamat <- vector("list",L)


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]


        A2 <-  matrix(0,1,ncol(omegadot))
        B2 <- matrix(0,TotT^2,ncol(S2))


        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%S2%*%omegadot
          B2 <- B2 + Q2%*%S2%*%omegadot

        }

        A2 <- A2/length(IDs_not_fold)
        B2 <- B2/length(IDs_not_fold)


        #Compute inverses
        decompB <- svd(B2)
        U <-  decompB$u
        D <-  decompB$d
        V <-  decompB$v


        diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

        matBinv <-  V%*%diagmatinvB%*%t(U)

        Gammaomegamat[[l]] <- A2%*%matBinv

      }


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]

        M2 <- matrix(0,ncol(w), ncol(w))
        L2 <- matrix(0,1,ncol(w))

        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          M2 <- M2 + t(w2)%*%q2%*%w2


          L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

        }

        M2 <- M2/length(IDs_not_fold)


        decompM <- eigen(M2)
        diagmat <- Re(decompM$values)
        U <- decompM$vectors

        diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

        matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


        L2 <- -L2/length(IDs_not_fold)

        Gammabetamat[[l]] <- L2%*%matMinv
      }


      gsmlr <- vector("list",n)
      psism0 <- 0


      for(l in 1:length(I)){

        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

        }
      }


      gsmestlr <- vector("list",n)

      for(l in 1:length(I)){

        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
        }
      }

      gfmestlrv <- vector("list",n)
      for(l in 1:length(I)){
        psifmestvect <- psifmest[[l]]
        for(j1 in I[[l]]){
          gfmestlrv[[j1]] <- gfmestlr[uh,j1] + psifmestvect[uh]
        }
      }


      gsmestlr <- unlist(gsmestlr)
      gfmestlrv <- unlist(gfmestlrv)

      psivar0 <- 0

      est <- mean(unlist(gsmlr)) - Resultsfm2[uh,1]^2


      var <- var(gsmestlr) - 2*mean(gfmestlrv)*stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*(stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*var(gfmestlrv))

      seest <- sqrt(var)/sqrt(n)

      Resultssm2[uh,1] <- round(est,2) #Matrix of results
      Resultssm2[uh,2] <- round(seest,2)     #Matrix of results

    }

    colnames(Resultssm2) <- c("Est.", "se")
    rownames(Resultssm2) <- paste0("V.", 1:ncol(v))


    Resultssm <- rbind(Resultssm1,Resultssm2)


    #Covariances
    pairs <- utils::combn(ncol(v), 2)

    for (uh in 1:ncol(pairs)){
      i <- pairs[1, uh]
      j <- pairs[2, uh]
      Omega <- matrix(0, ncol(v), ncol(v))
      Omega[i, j] <- 1/2
      Omega[j, i] <- 1/2


      Gammaomegamat <- vector("list",L)
      Gammabetamat <- vector("list",L)


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]


        A2 <-  matrix(0,1,ncol(omegadot))
        B2 <- matrix(0,TotT^2,ncol(S2))


        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%S2%*%omegadot
          B2 <- B2 + Q2%*%S2%*%omegadot

        }

        A2 <- A2/length(IDs_not_fold)
        B2 <- B2/length(IDs_not_fold)


        #Compute inverses
        decompB <- svd(B2)
        U <-  decompB$u
        D <-  decompB$d
        V <-  decompB$v


        diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

        matBinv <-  V%*%diagmatinvB%*%t(U)

        Gammaomegamat[[l]] <- A2%*%matBinv

      }


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]

        M2 <- matrix(0,ncol(w), ncol(w))
        L2 <- matrix(0,1,ncol(w))

        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          M2 <- M2 + t(w2)%*%q2%*%w2


          L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

        }

        M2 <- M2/length(IDs_not_fold)


        decompM <- eigen(M2)
        diagmat <- Re(decompM$values)
        U <- decompM$vectors

        diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

        matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


        L2 <- -L2/length(IDs_not_fold)

        Gammabetamat[[l]] <- L2%*%matMinv
      }


      gsmlr <- vector("list",n)
      psism0 <- 0


      for(l in 1:length(I)){

        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

        }
      }


      gsmestlr <- vector("list",n)

      for(l in 1:length(I)){
        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
        }
      }

      gfmestlrv1 <-  gfmestlrv2 <- vector("list",n)
      for(l in 1:length(I)){
        psifmestvect1 <-  psifmestvect2 <- psifmest[[l]]
        for(j1 in I[[l]]){
          gfmestlrv1[[j1]] <- gfmestlr[i,j1] + psifmestvect[i]
          gfmestlrv2[[j1]] <- gfmestlr[j,j1] + psifmestvect[j]
        }
      }



      gsmestlr <- unlist(gsmestlr)
      gfmestlrv1 <- unlist(gfmestlrv1)
      gfmestlrv2 <- unlist(gfmestlrv2)


      psivar0 <- 0

      est <- mean(unlist(gsmlr)) - Resultsfm2[i,1]*Resultsfm2[j,1]


      var <- var(gsmestlr) + (mean(gfmestlrv2)^2)*var(gfmestlrv1) + (mean(gfmestlrv1)^2)*var(gfmestlrv2) - 2*mean(gfmestlrv2)*stats::cov(gsmestlr,gfmestlrv1) -2*mean(gfmestlrv1)*stats::cov(gsmestlr,gfmestlrv2) + 2*mean(gfmestlrv1)*mean(gfmestlrv2)*stats::cov(gfmestlrv1,gfmestlrv2)

      seest <- sqrt(var)/sqrt(n)


      Resultscov[[paste0("v", i, "v", j)]] <- list("Est." = est, "S.e." = seest)

    }
  }else{

    Resultssm <- matrix(NA,ncol(v),2)
    Resultscov <- list()

    for(uh in 1:ncol(v)){
      Omega <- matrix(0, ncol(v), ncol(v))
      Omega[uh,uh] <- 1

      if(uh==1){

        varerrors <- vector("list", L)

        for(l in 1:length(I)){
          varutvL <- matrix(0,TotT^2,1)
          IDs_not_fold <- seqn[!seqn%in%I[[l]]]

          for (j1 in IDs_not_fold) {
            y2 <- y[id == j1]
            w2 <- w[id == j1, , drop = FALSE]
            v2 <- v[id == j1, , drop = FALSE]
            q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
            h2 <- MASS::ginv(v2)

            residual <- y2 - w2 %*%betalasso[[l]]- v2 %*%psifmestmatcf[,l]
            part1 <- S2%*%MASS::ginv((diag(TotT)%x%q2)%*%S2)
            part2 <- residual%x%(q2%*%(y2 - w2 %*% betalasso[[l]]))
            varutvL <- varutvL + part1 %*% part2

          }

          varerrors[[l]] <- varutvL /length(IDs_not_fold) #Variance of errors, time-varying variances in levels (robust) from Arellano and Bonhomme (2012)
        }
      }


      Gammaomegamat <- vector("list",L)
      Gammabetamat <- vector("list",L)


      omegadot <- diag(ncol(S2))

      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]


        A2 <-  matrix(0,1,ncol(omegadot))
        B2 <- matrix(0,TotT^2,ncol(S2))


        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%S2%*%omegadot
          B2 <- B2 + Q2%*%S2%*%omegadot

        }

        A2 <- A2/length(IDs_not_fold)
        B2 <- B2/length(IDs_not_fold)


        #Compute inverses
        decompB <- svd(B2)
        U <-  decompB$u
        D <-  decompB$d
        V <-  decompB$v


        diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

        matBinv <-  V%*%diagmatinvB%*%t(U)

        Gammaomegamat[[l]] <- A2%*%matBinv

      }


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]

        M2 <- matrix(0,ncol(w), ncol(w))
        L2 <- matrix(0,1,ncol(w))

        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          M2 <- M2 + t(w2)%*%q2%*%w2


          L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

        }

        M2 <- M2/length(IDs_not_fold)


        decompM <- eigen(M2)
        diagmat <- Re(decompM$values)
        U <- decompM$vectors

        diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

        matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


        L2 <- -L2/length(IDs_not_fold)

        Gammabetamat[[l]] <- L2%*%matMinv
      }


      gsmlr <- vector("list",n)
      psism0 <- 0


      for(l in 1:length(I)){

        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

        }
      }


      gsmestlr <- vector("list",n)

      for(l in 1:length(I)){
        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmestlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
        }
      }

      gfmestlrv <- vector("list",n)
      for(l in 1:length(I)){
        psifmestvect <- psifmest[[l]]
        for(j1 in I[[l]]){
          gfmestlrv[[j1]] <- gfmestlr[uh,j1] + psifmestvect[uh]
        }
      }


      gsmestlr <- unlist(gsmestlr)
      gfmestlrv <- unlist(gfmestlrv)

      psivar0 <- 0

      est <- mean(unlist(gsmlr)) - Resultsfm2[uh,1]^2


      var <- var(gsmestlr) - 2*mean(gfmestlrv)*stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*(stats::cov(gsmestlr, gfmestlrv) - 2*mean(gfmestlrv)*var(gfmestlrv))

      seest <- sqrt(var)/sqrt(n)

      Resultssm[uh,1] <- round(est,2)       #Matrix of results
      Resultssm[uh,2] <- round(seest,2)     #Matrix of results

    }

    colnames(Resultssm) <- c("Est.", "S.e.")
    rownames(Resultssm) <- paste0("V.", 1:ncol(v))


    #Covariances
    pairs <- utils::combn(ncol(v), 2)

    for (uh in 1:ncol(pairs)){
      i <- pairs[1, uh]
      j <- pairs[2, uh]
      Omega <- matrix(0, ncol(v), ncol(v))
      Omega[i, j] <- 1/2
      Omega[j, i] <- 1/2


      Gammaomegamat <- vector("list",L)
      Gammabetamat <- vector("list",L)


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]


        A2 <-  matrix(0,1,ncol(omegadot))
        B2 <- matrix(0,TotT^2,ncol(S2))


        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          A2 <- A2 + t(cbind(c(Omega)))%*%H2%*%S2%*%omegadot
          B2 <- B2 + Q2%*%S2%*%omegadot

        }

        A2 <- A2/length(IDs_not_fold)
        B2 <- B2/length(IDs_not_fold)


        #Compute inverses
        decompB <- svd(B2)
        U <-  decompB$u
        D <-  decompB$d
        V <-  decompB$v


        diagmatinvB <- diag(ifelse(D>=(log(TotT^2)/n)^(1/2), 1/D, 0), length(D), length(D))

        matBinv <-  V%*%diagmatinvB%*%t(U)

        Gammaomegamat[[l]] <- A2%*%matBinv

      }


      for(l in 1:length(I)){
        IDs_not_fold <- seqn[!seqn%in%I[[l]]]

        M2 <- matrix(0,ncol(w), ncol(w))
        L2 <- matrix(0,1,ncol(w))

        for (j1 in IDs_not_fold){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          M2 <- M2 + t(w2)%*%q2%*%w2


          L2 <- L2 + (t(cbind(c(Omega)))%*%H2 -  as.numeric(Gammaomegamat[[l]])%*%Q2)%*%((w2%x%(y2 - w2%*%betalasso[[l]])) + ((y2 - w2%*%betalasso[[l]])%x%w2))

        }

        M2 <- M2/length(IDs_not_fold)


        decompM <- eigen(M2)
        diagmat <- Re(decompM$values)
        U <- decompM$vectors

        diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

        matMinv <- U%*%diagmatinv%*%t(U)  #Compute suitable estimator of the Moore-Penrose inverse of M


        L2 <- -L2/length(IDs_not_fold)

        Gammabetamat[[l]] <- L2%*%matMinv
      }


      gsmlr <- vector("list",n)
      psism0 <- 0


      for(l in 1:length(I)){

        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmlr[[j1]] <- (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]]) - psism0

        }
      }


      gsmestlr <- vector("list",n)

      for(l in 1:length(I)){
        for (j1 in I[[l]]){
          y2 <- y[id == j1]
          w2 <- w[id == j1, , drop = FALSE]
          v2 <- v[id == j1, , drop = FALSE]
          q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
          h2 <- MASS::ginv(v2)

          Q2 <- diag(TotT^2) - (diag(TotT) - q2)%x%(diag(TotT) - q2)
          H2 <- h2%x%h2

          U2 <- (y2 - w2%*%betalasso[[l]])%x%(y2 - w2%*%betalasso[[l]]) - varerrors[[l]]


          gsmestlr[[j1]] <-  (t(cbind(c(Omega)))%*%H2 - as.numeric(Gammaomegamat[[l]])%*%Q2)%*%U2 - as.numeric(Gammabetamat[[l]])%*%t(w2)%*%q2%*%(y2-w2%*%betalasso[[l]])
        }
      }

      gfmestlrv1 <-  gfmestlrv2 <- vector("list",n)
      for(l in 1:length(I)){
        psifmestvect1 <-  psifmestvect2 <- psifmest[[l]]
        for(j1 in I[[l]]){
          gfmestlrv1[[j1]] <- gfmestlr[i,j1] + psifmestvect1[i]
          gfmestlrv2[[j1]] <- gfmestlr[j,j1] + psifmestvect2[j]
        }
      }

      gsmestlr <- unlist(gsmestlr)
      gfmestlrv1 <- unlist(gfmestlrv1)
      gfmestlrv2 <- unlist(gfmestlrv2)


      psivar0 <- 0

      est <- mean(unlist(gsmlr)) - Resultsfm2[i,1]*Resultsfm2[j,1]



      var <- var(gsmestlr) + (mean(gfmestlrv2)^2)*var(gfmestlrv1) + (mean(gfmestlrv1)^2)*var(gfmestlrv2) - 2*mean(gfmestlrv2)*stats::cov(gsmestlr,gfmestlrv1) -2*mean(gfmestlrv1)*stats::cov(gsmestlr,gfmestlrv2) + 2*mean(gfmestlrv1)*mean(gfmestlrv2)*stats::cov(gfmestlrv1,gfmestlrv2)

      seest <- sqrt(var)/sqrt(n)

      Resultscov[[paste0("v", i, "v", j)]] <- list("Est." = est, "S.e." = seest)

    }

  }

}




if(!is.null(C1)){

  if(qr(C1)$rank != ncol(C1)) {
    stop("Error: The rank of C1 is not equal to its number of columns.")
  }

  Resultsbeta1 <- matrix(NA,ncol(C1),2)
  Resultsbeta2 <- matrix(NA,indreg,2)
  psibeta0 <- rep(0,ncol(C1))
  Gammamatbeta <- list()

  # Loop over individuals
  for(l in 1:length(I)){
    WQWmat <- matrix(0,ncol(w), ncol(w))

    IDs_not_fold <- seqn[!seqn%in%I[[l]]]

    for(j1 in IDs_not_fold){
      y2 <- y[id == j1]
      w2 <- w[id == j1, , drop = FALSE]
      v2 <- v[id == j1, , drop = FALSE]
      q2 <- diag(TotT) - v2 %*%MASS::ginv(t(v2)%*%v2)%*%t(v2)
      h2 <- MASS::ginv(v2)

      WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
    }

    Mmat[[l]] <- WQWmat/length(IDs_not_fold)


    #Compute suitable estimator of the Moore-Penrose inverse of M
    decompM <- eigen(Mmat[[l]])
    diagmat <- Re(decompM$values)
    U <- decompM$vectors

    diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

    matMinv <- U%*%diagmatinv%*%t(U)

    Gammamatbeta[[l]] <- t(C1)%*%matMinv

  }




  gbetalr <- matrix(NA,ncol(C1),n)
  QC1 <- diag(ncol(w)) - C1%*%MASS::ginv(t(C1)%*%C1)%*%t(C1)


  for(l in 1:length(I)){

    for(j1 in I[[l]]){

      gbetalr[,j1] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%MASS::ginv(t(C1)%*%C1)%*%psibeta0 - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
    }

  }



  Rbeta <- vector("list",L)
  psibetaest <- vector("list",L)

  for(l in 1:length(I)){
    IDs_not_fold <- seqn[!seqn%in%I[[l]]]
    Rbetainner <- 0

    for(j1 in IDs_not_fold){
      Rbetainner <-  Rbetainner + t(Wmat[[j1]])%*%Qmat[[j1]]%*%ymat[[j1]]
    }
    Rbeta[[l]] <- Rbetainner/length(IDs_not_fold)
    psibetaest[[l]] <- Gammamatbeta[[l]]%*%Rbeta[[l]]
  }



  gbetaestlr <- matrix(NA,ncol(C1),n)


  for(l in 1:length(I)){

    for(j1 in I[[l]]){

      gbetaestlr[,j1] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%MASS::ginv(t(C1)%*%C1)%*%as.numeric(psibetaest[[l]]) - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
    }

  }

  est <- rowMeans(gbetalr)

  meangbetaestlr <- rowMeans(gbetaestlr)

  secondmeangbetaestlr <- gbetaestlr%*%t(gbetaestlr)/n

  varest <- secondmeangbetaestlr -   meangbetaestlr%*%t(meangbetaestlr) #Variance estimator of LR moment for first moment of UH

  seest <- sqrt(diag(varest))/sqrt(n)


  Resultsbeta1[,1] <- round(est,2)
  Resultsbeta1[,2] <- round(seest,2)

  colnames(Resultsbeta1) <- c("Est.", "S.e.")
  rownames(Resultsbeta1) <- paste0("C1.", 1:ncol(C1))


  for(regre in 1:indreg){
    C1 <- rep(0, ncol(w))
    C1[regre] <- 1


    psibeta0 <- 0
    Gammamatbeta <- list()

    # Loop over individuals
    for(l in 1:length(I)){
      WQWmat <- matrix(0,ncol(w), ncol(w))

      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      for(j1 in IDs_not_fold){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
      }

      Mmat[[l]] <- WQWmat/length(IDs_not_fold)


      #Compute suitable estimator of the Moore-Penrose inverse of M
      decompM <- eigen(Mmat[[l]])
      diagmat <- Re(decompM$values)
      U <- decompM$vectors

      diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

      matMinv <- U%*%diagmatinv%*%t(U)

      Gammamatbeta[[l]] <- t(C1)%*%matMinv

    }



    gbetalr <- vector("list", n)

    QC1 <- diag(ncol(w)) - C1%*%MASS::ginv(t(C1)%*%C1)%*%t(C1)

    for(l in 1:length(I)){

      for(j1 in I[[l]]){

        gbetalr[[j1]] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%MASS::ginv(t(C1)%*%C1)*psibeta0 - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
      }

    }



    Rbeta <- vector("list",L)
    psibetaest <- vector("list",L)

    for(l in 1:length(I)){
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]
      Rbetainner <- 0

      for(j1 in IDs_not_fold){
        Rbetainner <-  Rbetainner + t(Wmat[[j1]])%*%Qmat[[j1]]%*%ymat[[j1]]
      }
      Rbeta[[l]] <- Rbetainner/length(IDs_not_fold)
      psibetaest[[l]] <- Gammamatbeta[[l]]%*%Rbeta[[l]]
    }



    gbetaestlr <- vector("list", n)


    for(l in 1:length(I)){

      for(j1 in I[[l]]){

        gbetaestlr[[j1]] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%MASS::ginv(t(C1)%*%C1)*as.numeric(psibetaest[[l]]) - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
      }

    }

    gbetalr <- unlist(gbetalr)
    gbetaestlr <- unlist(gbetaestlr)

    est <- mean(gbetalr)
    varest <- mean((gbetaestlr  - mean(gbetaestlr))^2)

     seest <- sqrt(varest)/sqrt(n)


Resultsbeta2[regre,1] <- round(est,2)
Resultsbeta2[regre,2] <- round(seest,2)
  }

colnames(Resultsbeta2) <- c("Est.", "S.e.")
rownames(Resultsbeta2) <- paste0("W.", 1:indreg)

Resultsbeta <- rbind(Resultsbeta1,Resultsbeta2)

}else{

  Resultsbeta <- matrix(NA,indreg,2)

  for(regre in 1:indreg){
    C1 <- rep(0, ncol(w))
    C1[regre] <- 1


    psibeta0 <- 0
    Gammamatbeta <- list()

    # Loop over individuals
    for(l in 1:length(I)){
      WQWmat <- matrix(0,ncol(w), ncol(w))

      IDs_not_fold <- seqn[!seqn%in%I[[l]]]

      for(j1 in IDs_not_fold){
        y2 <- y[id == j1]
        w2 <- w[id == j1, , drop = FALSE]
        v2 <- v[id == j1, , drop = FALSE]
        q2 <- diag(TotT) - v2 %*% MASS::ginv(t(v2) %*% v2) %*% t(v2)
        h2 <- MASS::ginv(v2)

        WQWmat <- WQWmat + t(w2) %*% q2 %*% w2
      }

      Mmat[[l]] <- WQWmat/length(IDs_not_fold)


      #Compute suitable estimator of the Moore-Penrose inverse of M
      decompM <- eigen(Mmat[[l]])
      diagmat <- Re(decompM$values)
      U <- decompM$vectors

      diagmatinv <- diag(ifelse(diagmat>=(log(ncol(w))/n)^(1/2), 1/diagmat, 0), length(diagmat), length(diagmat))

      matMinv <- U%*%diagmatinv%*%t(U)

      Gammamatbeta[[l]] <- t(C1)%*%matMinv

    }


    gbetalr <- vector("list", n)

    QC1 <- diag(ncol(w)) - C1%*%MASS::ginv(t(C1)%*%C1)%*%t(C1)

    for(l in 1:length(I)){

      for(j1 in I[[l]]){

        gbetalr[[j1]] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%ginv(t(C1)%*%C1)*psibeta0 - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
      }

    }



    Rbeta <- vector("list",L)
    psibetaest <- vector("list",L)

    for(l in 1:length(I)){
      IDs_not_fold <- seqn[!seqn%in%I[[l]]]
      Rbetainner <- 0

      for(j1 in IDs_not_fold){
        Rbetainner <-  Rbetainner + t(Wmat[[j1]])%*%Qmat[[j1]]%*%ymat[[j1]]
      }
      Rbeta[[l]] <- Rbetainner/length(IDs_not_fold)
      psibetaest[[l]] <- Gammamatbeta[[l]]%*%Rbeta[[l]]
    }



    gbetaestlr <- vector("list", n)


    for(l in 1:length(I)){

      for(j1 in I[[l]]){

        gbetaestlr[[j1]] <- Gammamatbeta[[l]]%*%t(Wmat[[j1]])%*%Qmat[[j1]]%*%(ymat[[j1]] - Wmat[[j1]]%*%C1%*%ginv(t(C1)%*%C1)*as.numeric(psibetaest[[l]]) - Wmat[[j1]]%*%QC1%*%betalasso[[l]])
      }

    }

    gbetalr <- unlist(gbetalr)
    gbetaestlr <- unlist(gbetaestlr)

    est <- mean(gbetalr)
    varest <- mean((gbetaestlr  - mean(gbetaestlr))^2)

    seest <- sqrt(varest)/sqrt(n)


    Resultsbeta[regre,1] <- round(est,2)
    Resultsbeta[regre,2] <- round(seest,2)
  }

  colnames(Resultsbeta) <- c("Est.", "S.e.")
  rownames(Resultsbeta) <- paste0("W.", 1:indreg)

}



return(list("FM" = Resultsfm, "SM"=Resultssm, "Cov"= Resultscov, "CP"=Resultsbeta))


}




