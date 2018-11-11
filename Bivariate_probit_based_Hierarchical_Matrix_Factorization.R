#####Bivariate probit based Hierarchical Matrix Factorization#####
options(warn=0)
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(mvtnorm)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####�C�ӂ̕��U�����U�s����쐬������֐�####
##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
s <- 2   #�����ϐ���
k <- 10   #��ꐔ
hh <- 10000   #���[�U�[��
item <- 3000   #�A�C�e����
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #�w���ڐG��
hhpt <- sum(pt)
vec <- rep(1, k)

#ID�̐ݒ�
user_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

##�A�C�e���̊����𐶐�
#�Z�O�����g�����𐶐�
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(1.0, topic))) %*% 1:topic)

#�������z����A�C�e���𐶐�
item_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  item_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, phi[z[user_id[user_list[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}
w <- unlist(lapply(item_list, length))

#�ʂɘa����邽�߂̃X�p�[�X�s��
user_dt <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(item_id, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))

#���������f�[�^������
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="�A�C�e���̍w���p�x", main="�A�C�e���̍w���p�x���z")
gc(); gc()


##�f���x�N�g���𐶐�
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #�f�[�^������
column <- ncol(X)
index_column <- matrix(1:(column*s), nrow=column, ncol=s)

##�K�w���f���̐����ϐ���ݒ�
#���[�U�[�̐����ϐ�
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������
column_u <- ncol(u)
index_k <- matrix(1:(k*s), nrow=k, ncol=s)

#�A�C�e���̐����ϐ�
k1 <- 2; k2 <- 3; k3 <- 4
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #�f�[�^������
column_v <- ncol(v)


####�����ϐ��𐶐�####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)

  ##�K�w���f���̃p�����[�^�𐶐�
  #�K�w���f���̕��U�p�����[�^
  Cov <- Covt <- diag(runif(column*s, 0.01, 0.15))   #�f���x�N�g���̊K�w���f���̕��U
  Cov_u <- Cov_ut <- diag(runif(k*s, 0.01, 0.2))   #���[�U�[�̊K�w���f���̕��U
  Cov_v <- Cov_vt <- diag(runif(k, 0.01, 0.2), k)   #�A�C�e���̊K�w���f���̕��U

  #�K�w���f���̉�A�W����ݒ�
  alpha<- matrix(0, nrow=column_u, ncol=column*s)
  for(j in 1:s){
    alpha[, index_column[, j]] <- cbind(rnorm(column_u, -0.2, 0.5), 
                                        matrix(rnorm((column-1)*column_u, 0, 0.5), nrow=column_u, ncol=column-1))
  }
  alphat <- alpha
  alpha_u <- alpha_ut <- matrix(rnorm(k*column_u*s, 0, 0.5), nrow=column_u, ncol=k*s)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.5), nrow=column_v, ncol=k)
  
  
  ##���f���p�����[�^�𐶐�
  #�f���x�N�g���ƍs�񕪉��̃p�����[�^�𐶐�
  beta <- betat <- u %*% alpha + mvrnorm(hh, rep(0, column*s), Cov)
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k*s), Cov_u)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  #���֍s��𐶐�
  Sigma <- matrix(0.5, nrow=s, ncol=s)
  diag(Sigma) <- 1
  Sigmat <- Sigma
  
  ##���݌��p���牞���ϐ��𐶐�
  #��A���f���̕��ύ\��
  y <- mu <- matrix(0, nrow=hhpt, ncol=s)
  W <- theta_u[user_id, ]
  H <- theta_v[item_id, ]
  for(j in 1:s){
    mu[, j] <- as.numeric((X * beta[user_id, index_column[, j]]) %*% rep(1, column) + (W[, index_k[, j]] * H) %*% vec)
  }
  
  #���ϗʐ��K���z������݌��p�𐶐�
  er <- mvrnorm(hhpt, rep(0, s), Sigma)   #���f���덷
  U <- UT <- mu + er   #���݌��p
  y <- ifelse(U > 0, 1, 0)   #�����ϐ��𐶐�
  
  #�X�g�b�v����
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.15){
    break
  }
}

####�}���R�t�A�������e�J�����@��Bivariate probit based Hierarchical Matrix Factorization�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##���ϗʐ��K���z�̏����t�����Ғl�Ə����t�����U���v�Z����֐�
cdMVN <- function(mean, Cov, dependent, U){
  
  #���U�����U�s��̃u���b�N�s����`
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #�����t�����U�Ə����t�����ς��v�Z
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #�����t�����ς��v�Z
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #�����t�����U���v�Z
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

##�f�[�^�̐ݒ�
#�萔�̐ݒ�
XX_list <- list()
for(i in 1:hh){
  XX_list[[i]] <- t(X[user_list[[i]], ]) %*% X[user_list[[i]], ]
}

#�����̃C���f�b�N�X
xx_index1 <- matrix(1:(s*column), nrow=column, ncol=s, byrow=T)
xx_index2 <- matrix(1:(s*k), nrow=k, ncol=s, byrow=T)

##���O���z�̐ݒ�
#�K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=column_u, ncol=column*s)
Deltabar_u <- matrix(0, nrow=column_u, ncol=k*s)
Deltabar_v <- matrix(0, nrow=column_v, ncol=k)
ADelta <- ADelta_u <- 0.01 * diag(column_u)
ADelta_v <- 0.01 * diag(column_v)
nu <- nu1 <- nu2 <- k + 1
V <- nu * diag(column*s)
V1 <- nu * diag(k*s)
V2 <- nu * diag(k)

#���f���̂̎��O���z
m <- column + 1   #�t�E�B�V���[�g���z�̎��R�x
Rn <- m * diag(s)   #�t�E�B�V���[�g���z�̎��R�x


##�p�����[�^�̐^�l
#�K�w���f���̃p�����[�^
Cov <- Covt; inv_Cov <- solve(Cov)
Cov_u <- Cov_ut; inv_Cov_u <- solve(Cov_u)
Cov_v <- Cov_vt; inv_Cov_v <- solve(Cov_v)
alpha <- alphat; alpha_u <- alpha_ut; alpha_v <- alpha_vt
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u ; v_mu <- v %*% alpha_v

#�f���x�N�g���ƍs�񕪉��̃p�����[�^
beta <- betat 
beta_vec <- beta[user_id, ]
beta_mu <- matrix(0, nrow=hhpt, ncol=s)
for(j in 1:s){
beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
}
theta_u <- theta_ut 
theta_v <- theta_vt 

#���֍s����p�����[�^
Sigma <- Sigmat
inv_Sigma <- solve(Sigma)

#���f���̕��ύ\��
mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
W <- theta_u[user_id, ]; H <- theta_v[item_id, ]
for(j in 1:s){
WH[, j] <- (W[, index_k[, j]] * H) %*% vec
mu[, j] <- as.numeric(beta_mu[, j]) + WH[, j]
}
U <- UT


##�p�����[�^�̏����l
#�K�w���f���̕��U�p�����[�^
Cov <- 0.01 * diag(column*s); inv_Cov <- solve(Cov)
Cov_u <- 0.01 * diag(k*s); inv_Cov_u <- solve(Cov_u)
Cov_v <- 0.01 * diag(k); inv_Cov_v <- solve(Cov_v)

#�K�w���f���̉�A�W����ݒ�
alpha <- matrix(rnorm(column*column_u*s, 0, 0.1), nrow=column_u, ncol=column*s)
alpha_u <- matrix(rnorm(k*column_u*s, 0, 0.5), nrow=column_u, ncol=k*s)
alpha_v <- matrix(rnorm(k*ncol(v), 0, 0.5), nrow=column_v, ncol=k)
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u ; v_mu <- v %*% alpha_v

#�f���x�N�g���ƍs�񕪉��̃p�����[�^���쐬
beta <- u %*% alpha + mvrnorm(hh, rep(0, column*s), Cov)
beta_vec <- beta[user_id, ]
beta_mu <- matrix(0, nrow=hhpt, ncol=s)
for(j in 1:s){
  beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
}
theta_u <- u %*% alpha_u + mvrnorm(hh, rep(0, k*s), Cov_u)
theta_v <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)

#���֍s��𐶐�
Sigma <- diag(s)
inv_Sigma <- solve(Sigma)

#���f���̕��ύ\��
mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
W <- theta_u[user_id, ]; H <- theta_v[item_id, ]
for(j in 1:s){
  WH[, j] <- (W[, index_k[, j]] * H) %*% vec
  mu[, j] <- as.numeric(beta_mu[, j]) + WH[, j]
}
U <- (y)*abs(mu) + (-(1-y)*abs(mu))


##�T���v�����O���ʂ̊i�[�p�z��
COV <- array(0, dim=c(column*s, column*s, R/keep))
COV_U <- array(0, dim=c(k*s, k*s, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
ALPHA <- array(0, dim=c(column_u, column*s, R/keep))
ALPHA_U <- array(0, dim=c(column_u, k*s, R/keep))
ALPHA_V <- array(0, dim=c(column_v, k, R/keep))
BETA <- array(0, dim=c(hh, column*s, R/keep))
THETA_U <- array(0, dim=c(hh, k*s, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
SIGMA <- array(0, dim=c(s, s, R/keep))


##�ΐ��ޓx�̊�l
#1�p�����[�^���f���ł̑ΐ��ޓx
Prob <- matrix(colMeans(y), nrow=hhpt, ncol=s, byrow=T)
print(LLst <- sum(y*log(Prob) + (1-y)*log(1-Prob)))

#�^�l�ł̑ΐ��ޓx
LLi <- rep(0, s)
mut <- matrix(0, nrow=hhpt, ncol=s)
WT <- theta_ut[user_id, ]; HT <- theta_vt[item_id, ]
for(j in 1:s){
  mut[, j] <- as.numeric((X * betat[user_id, index_column[, j]]) %*% rep(1, column) + (WT[, index_k[, j]] * HT) %*% vec)
  Prob <- pnorm(mut[, j], Sigmat[j, j]); Prob[Prob==0] <- 10^-200; Prob[Prob==1] <- 0.99999999   #�����m��
  LLi[j] <- sum(y[, j]*log(Prob) + (1-y[, j])*log(1-Prob))   #�ΐ��ޓx
}
print(LLbest <- sum(LLi))


##�ؒf���K���z�̐ؒf�̈���`
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)
a0 <- ifelse(y==0, -10^-100, 0)
b0 <- ifelse(y==1, 10^-100, 0)


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z������݌��p�𐶐�
  #���f���̕��ύ\��
  mu <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    mu[, j] <- beta_mu[, j] + WH[, j]
  }
  
  #���݌��p�𐶐�
  for(j in 1:s){
    MVR <- cdMVN(mu, Sigma, j, U)   #�����t�����z�ƕ��U�𐶐�
    MVR_U <- as.numeric(MVR$CDmu)
    MVR_S <- as.numeric(sqrt(MVR$CDvar))
    U[, j] <- rtnorm(MVR_U, MVR_S, a[, j], b[, j])   #���p�𐶐�
    U[is.infinite(U[, j])==TRUE, j] <- 0
  }

  ##���[�U�[�ʂɑf���x�N�g�����T���v�����O
  #���f���̉����ϐ���ݒ�
  er <- U - WH
  
  for(i in 1:hh){
    #�f�[�^�̒��o
    er_u <- er[user_list[[i]], ]
    x <- X[user_list[[i]], ]
    
    #��A�W���̎��㕪�z�̃p�����[�^
    XX <- XX_list[[i]]; XU <- as.numeric(t(x) %*% er_u %*% inv_Sigma)
    inv_XVX <- solve(kronecker(inv_Sigma, XX) + inv_Cov)
    beta_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov %*% alpha_mu[i, ]))   #��A�x�N�g���̊��Ғl
    
    #���ϗʐ��K���z����f���x�N�g�����T���v�����O
    beta[i, ] <- mvrnorm(1, beta_mu, inv_XVX)
  }
  
  #�f���x�N�g���̕���
  beta_vec <- beta[user_id, ]
  beta_mu <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
  }
  
  
  ##���[�U�[�̓����s��̃p�����[�^���T���v�����O
  #�����ϐ��̐ݒ�
  er <- U - beta_mu
  
  for(i in 1:hh){
    #�f�[�^�̒��o
    er_w <- er[user_list[[i]], ]
    x <- H[user_list[[i]], ]
    
    #��A�W���̎��㕪�z�̃p�����[�^
    XX <- t(x) %*% x; XU <- as.numeric(t(x) %*% er_w %*% inv_Sigma)
    inv_XVX <- solve(kronecker(inv_Sigma, XX) + inv_Cov_u)
    w_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov_u %*% u_mu[i, ]))   #��A�x�N�g���̊��Ғl
    
    #���ϗʐ��K���z����f���x�N�g�����T���v�����O
    theta_u[i, ] <- mvrnorm(1, w_mu, inv_XVX)
  }
  #���[�U�[�����s���ϊ�
  W <- theta_u[user_id, ]
  
  
  ##�A�C�e���̓����s��̃p�����[�^���T���v�����O
  #�����ϐ��̐ݒ�
  er <- U - beta_mu
  Chol <- chol(inv_Sigma)   #���U�����U�s��̋t�s����R���c�L�[����
  
  for(j in 1:item){
    #�f�[�^�̒��o
    er_h <- as.numeric(t(Chol %*% t(er[item_list[[j]], ])))
    x <- W[item_list[[j]], ]
    x_vec <- cbind(as.numeric(t(x[, index_k[, 1]])), as.numeric(t(x[, index_k[, 2]])))

    #��A�W���̎��㕪�z�̃p�����[�^
    Chol_x <- matrix(x_vec %*% t(Chol), nrow=w[j]*s, ncol=k, byrow=T)
    XX <- t(Chol_x) %*% Chol_x; XU <- t(Chol_x) %*% er_h
    inv_XVX <- solve(XX + inv_Cov_v)
    h_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov_v %*% v_mu[j, ]))   #��A�x�N�g���̊���
    
    #���ϗʐ��K���z����f���x�N�g�����T���v�����O
    theta_v[j, ] <- mvrnorm(1, h_mu, inv_XVX)
  }
  #�A�C�e�������s���ϊ�
  H <- theta_v[item_id, ]
  
  
  ##�t�E�B�V���[�g���z���瑊�֍s����T���v�����O
  #���f���̕��ύ\��
  mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    WH[, j] <- (W[, index_k[, j]] * H) %*% vec
    mu[, j] <- beta_mu[, j] + WH[, j]
  }
  
  #�t�E�B�V���[�g���z�̃p�����[�^
  R_error <- U - mu
  IW_R <- t(R_error) %*% R_error + Rn
  Sn <- hhpt + m
  
  #�t�E�B�V���[�g���z����p�����[�^���T���v�����O
  Sigma_hat <- rwishart(Sn + 1000, solve(IW_R))$IW
  Sigma <- cov2cor(Sigma_hat)   #���֍s��ɕϊ�
  inv_Sigma <- solve(Sigma)
  Sigma
  
  ##�K�w���f���̃p�����[�^�̃T���v�����O
  #���ϗʉ�A���f������f���x�N�g���̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(beta, u, Deltabar, ADelta, nu, V)
  alpha <- out$B
  alpha_mu <- u %*% alpha   #�f���x�N�g���̓����s��̕��ύ\��
  Cov <- diag(diag(out$Sigma))
  inv_Cov <- solve(Cov)
  
  #���ϗʉ�A���f�����烆�[�U�[�����s��K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_u, u, Deltabar_u, ADelta_u, nu1, V1)
  alpha_u <- out$B
  u_mu <- u %*% alpha_u   #���[�U�[�����s��̕��ύ\��
  Cov_u <- diag(diag(out$Sigma))
  inv_Cov_u <- solve(Cov_u)
  
  #���ϗʉ�A���f������A�C�e�������s��K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu2, V2)
  alpha_v <- out$B
  v_mu <- v %*% alpha_v   #�A�C�e�������s��̕��ύ\��
  Cov_v <- diag(diag(out$Sigma))
  inv_Cov_v <- solve(Cov_v)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    COV[, , mkeep] <- Cov
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    ALPHA[, , mkeep] <- alpha
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    BETA[, , mkeep] <- beta
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    SIGMA[, , mkeep] <- Sigma
  }
  
  ##�ΐ��ޓx�֐��̌v�Z�ƃT���v�����O���ʂ̊m�F
  if(rp%%disp==0){
    LLi <- rep(0, s)
    for(j in 1:s){
      Prob <- pnorm(mu[, j], Sigma[j, j]); Prob[Prob==0] <- 10^-200; Prob[Prob==1] <- 0.99999999   #�����m��
      LLi[j] <- sum(y[, j]*log(Prob) + (1-y[, j])*log(1-Prob))   #�ΐ��ޓx
    }
    LL <- sum(LLi)
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLbest, LLst))
    print(Sigma)
    print(round(rbind(diag(Cov_v), diag(Cov_vt)), 3))
  }
}