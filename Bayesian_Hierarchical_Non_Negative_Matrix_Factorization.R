#####Bayesian Hierarchical Non Negative Matrix Factorization#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(NMF)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(78594)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
k <- 10   #��ꐔ
hh <- 5000   #���[�U�[��
item <- 2000  #�A�C�e����
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #�ڐG��
hhpt <- sum(pt)   #���ڐG��
vec_k <- rep(1, k)

#ID��ݒ�
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #�f�[�^�̌���
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}


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
u_col <- ncol(u)

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
v_col <- ncol(v)


##�A�C�e���̊����𐶐�
#�Z�O�����g�����𐶐�
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

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

#���������f�[�^������
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="�A�C�e���̍w���p�x", main="�A�C�e���̍w���p�x���z")
gc(); gc()


####�����ϐ��𐶐�####
rp <- 0
repeat {
  rp <- rp + 1

  ##NMF�̃p�����[�^�𐶐�
  #�K���}���z�̎ړx�p�����[�^��ݒ�
  alpha_u <- alpha_ut <- matrix(rnorm(k*ncol(u), 0, 0.3), nrow=ncol(u), ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.3), nrow=ncol(v), ncol=k)
  lambda_u <- exp(u %*% alpha_u)
  lambda_v <- exp(v %*% alpha_v)
  
  #�K���}���z�̌`��p�����[�^
  beta_u <- beta_ut <- rnorm(1, 1, 0.5)
  beta_v <- beta_vt <- rnorm(1, 1, 0.5)
  
  #�K���}���z����s�񕪉��̃p�����[�^�𐶐�
  theta_u <- theta_ut <- matrix(rgamma(hh*k, as.numeric(lambda_u), beta_u), nrow=hh, ncol=k)
  theta_v <- theta_vt <- matrix(rgamma(item*k, as.numeric(lambda_v), beta_v), nrow=item, ncol=k)
  
  ##�|�A�\�����z���牞���ϐ��𐶐�
  WH <- (theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k   #���Ғl
  y�@<- rpois(hhpt, WH)   
  
  #break����
  print(rp)
  print(max(y))
  if(max(y) < 75 & max(y) > 25){
    break
  }
}

#���������f�[�^�̃q�X�g�O����
hist(y, main="�w���p�x�̕��z", xlab="�w���p�x", col="grey", breaks=50)


####�}���R�t�A�������e�J�����@��Bayesian Hierarchical NMF�𐄒�####
##�ΐ����㕪�z���v�Z����֐�
loglike <- function(beta, alpha, inv_tau, y, y_log, x){
  
  #�K���}��A���f���̑ΐ��ޓx
  lambda <- as.numeric(exp(x %*% beta))   #���Ғl
  Lho <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)   #�ΐ��ޓx�֐�
  
  #���ϗʐ��K���z�̑ΐ����O���z
  log_mvn <- -1/2 * as.numeric(t(beta) %*% inv_tau %*% beta)
  
  #�ΐ����㕪�z
  LL <- Lho + log_mvn
  return(list(LL=LL, Lho=Lho))
}

##HMC�Ŏړx�p�����[�^���T���v�����O���邽�߂̊֐�
#�K���}��A�̑ΐ����㕪�z�̔����֐�
dloglike <- function(beta, alpha, inv_tau, y, y_log, x){ 
  
  #���Ғl�̐ݒ�
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #���Ғl
  
  #�����֐��̐ݒ�
  dlgamma <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #�ړx�p�����[�^�̌��z�x�N�g��
  dmvn <- as.numeric(-inv_tau %*% beta)
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -(dlgamma + dmvn)
  return(LLd)
}

#���[�v�t���b�O�@�������֐�
leapfrog_u <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta_u, inv_tau_u, d, d_log, u) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta_u, inv_tau_u, d, d_log, u) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta_v, inv_tau_v, d, d_log, v) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta_v, inv_tau_v, d, d_log, v) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

#�`��p�����[�^�̑ΐ����㕪�z�̔����֐�
dloglike_alpha <- function(alpha, beta, y, y_log, x, n, k){ 
  #���Ғl�̐ݒ�
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #���Ғl
  
  #���z�x�N�g���̌v�Z
  dlgamma <- (n*k)*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #�`��p�����[�^�̌��z�x�N�g��
  return(dlgamma)
}


##�A���S���Y���̐ݒ�
R <- 5000
keep <- 4
burnin <- 1000/keep
iter <- 0
disp <- 10
e <- 0.001
L <- 3

##���O���z�̐ݒ�
gamma_u <- rep(0, u_col); gamma_v <- rep(0, v_col)
inv_tau_u <- solve(100 * diag(u_col)); inv_tau_v <- solve(100 * diag(v_col))
inv_tau <- solve(100 * diag(k))

##�p�����[�^�̐^�l
#�K���}���z�̎ړx�p�����[�^
alpha_u <- alpha_ut
alpha_v <- alpha_vt
lambda_u <- exp(u %*% alpha_u)
lambda_v <- exp(v %*% alpha_v)

#�K���}���z�̌`��p�����[�^
beta_u <- beta_ut 
beta_v <- beta_vt

#�s�񕪉��̃p�����[�^
theta_u <- theta_ut
theta_v <- theta_vt
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #���Ғl

##�����l�̐ݒ�
#�K���}���z�̎ړx�p�����[�^
alpha_u <- matrix(rnorm(k*ncol(u), 0, 0.1), nrow=ncol(u), ncol=k)
alpha_v <- matrix(rnorm(k*ncol(v), 0, 0.1), nrow=ncol(v), ncol=k)
lambda_u <- exp(u %*% alpha_u)
lambda_v <- exp(v %*% alpha_v)

#�K���}���z�̌`��p�����[�^
beta_u <- 1.0
beta_v <- 1.0

#�s�񕪉��̃p�����[�^
theta_u <- matrix(rgamma(hh*k, as.numeric(lambda_u), beta_u), nrow=hh, ncol=k)
theta_v <- matrix(rgamma(item*k, as.numeric(lambda_v), beta_v), nrow=item, ncol=k)
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #���Ғl


##�T���v�����O���ʂ̕ۑ��p�z��
gamma_u1 <- gamma_v1 <- rep(k)
gamma_u2 <- gamma_v2 <- rep(k)
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
ALPHA_U <- array(0, dim=c(u_col, k, R/keep))
ALPHA_V <- array(0, dim=c(v_col, k, R/keep))
BETA_U <- rep(0, R/keep)
BETA_V <- rep(0, R/keep)

##���[�U�[����уA�C�e���̃C���f�b�N�X���쐬
#�ʂɘa����邽�߂̃X�p�[�X�s��
user_dt <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))

##�ΐ��ޓx�̊�l
LLst <- sum(dpois(y, mean(y), log=TRUE))
LLbest <- sum(dpois(y, as.numeric((theta_ut[user_id, ] * theta_vt[item_id, ]) %*% vec_k), log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##���[�U�[�����s����T���v�����O
  #�⏕�ϐ�lambda���X�V
  theta_vec2 <- theta_v[item_id, ]
  lambda <- (theta_u[user_id, ] * theta_vec2) / WH

  #���[�U�[���Ƃ̃K���}���z�̃p�����[�^��ݒ�
  lambda_y <- lambda * y   #�v�f���Ƃ̊��Ғl
  W1 <- as.matrix(user_dt %*% lambda_y + lambda_u)
  W2 <- as.matrix(user_dt %*% theta_vec2 + beta_u)
  
  #�K���}���z���p�����[�^���T���v�����O
  theta_u <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  #theta_u <- theta_u / matrix(colSums(theta_u), nrow=hh, ncol=k, byrow=T) * hh/(k/2)   #�e��x�N�g���𐳋K��
  
  
  ##�A�C�e�������s����T���v�����O
  #�⏕�ϐ�lambda���X�V
  theta_vec1 <- theta_u[user_id, ]
  WH <- as.numeric((theta_vec1 * theta_vec2) %*% vec_k)
  lambda <- (theta_vec1 * theta_vec2) / WH

  #�A�C�e�����Ƃ̃K���}���z�̃p�����[�^��ݒ�
  lambda_y <- lambda * y   #�v�f���Ƃ̊��Ғl
  H1 <- as.matrix(item_dt %*% lambda_y + lambda_v)
  H2 <- as.matrix(item_dt %*% theta_vec1 + beta_v)
  
  #�K���}���z���p�����[�^���T���v�����O
  theta_v <- matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k)
  WH <- as.numeric((theta_vec1 * theta_v[item_id, ]) %*% vec_k)
  
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  for(j in 1:k){
    
    ##���[�U�[�����s��̎ړx�p�����[�^���T���v�����O
    #HMC�̐V�����p�����[�^�𐶐�
    rold <- as.numeric(mvrnorm(1, rep(0, u_col), diag(u_col)))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
    alphad <- alpha_u[, j]
    
    #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
    d <- theta_u[, j]; d_log <- log(d)
    res <- leapfrog_u(rold, alphad, dloglike, e, L)
    rnew <- res$r
    alphan <- res$z
    
    #�ړ��O�ƈړ���̃n�~���g�j�A��
    Hnew <- -loglike(alphan, beta_u, inv_tau_u, d, d_log, u)$LL + as.numeric(rnew^2 %*% rep(1, u_col))/2
    Hold <- -loglike(alphad, beta_u, inv_tau_u, d, d_log, u)$LL + as.numeric(rold^2 %*% rep(1, u_col))/2
    
    #�p�����[�^�̍̑�������
    rand <- runif(1)   #��l���z���痐���𔭐�
    gamma <- min(c(1, exp(Hold - Hnew)))   #�̑𗦂�����
    gamma_u1[j] <- gamma
    
    #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
    flag <- as.numeric(gamma > rand)
    alpha_u[, j] <- flag*alphan + (1-flag)*alphad
    
    
    ##�A�C�e�������s��̎ړx�p�����[�^���T���v�����O
    #HMC�̐V�����p�����[�^�𐶐�
    rold <- as.numeric(mvrnorm(1, rep(0, v_col), diag(v_col)))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
    alphad <- alpha_v[, j]
    
    #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
    d <- theta_v[, j]; d_log <- log(d)
    res <- leapfrog_v(rold, alphad, dloglike, e, L)
    rnew <- res$r
    alphan <- res$z
  
    #�ړ��O�ƈړ���̃n�~���g�j�A��
    Hnew <- -loglike(alphan, beta_v, inv_tau_v, d, d_log, v)$LL + as.numeric(rnew^2 %*% rep(1, v_col))/2
    Hold <- -loglike(alphad, beta_v, inv_tau_v, d, d_log, v)$LL + as.numeric(rold^2 %*% rep(1, v_col))/2
    
    #�p�����[�^�̍̑�������
    rand <- runif(1)   #��l���z���痐���𔭐�
    gamma <- min(c(1, exp(Hold - Hnew)))   #�̑𗦂�����
    gamma_v1[j] <- gamma
    
    #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
    flag <- as.numeric(gamma > rand)
    alpha_v[, j] <- flag*alphan + (1-flag)*alphad
  }
  
  ##���[�U�[�����s��̌`��p�����[�^���T���v�����O
  #MH�@�̐V�����p�����[�^�𐶐�
  d <- as.numeric(theta_u); d_log <- log(d)
  s <- sign(dloglike_alpha(beta_u, alpha_u, d, d_log, u, hh, k))   #���z�̕���
  betad <- beta_u
  betan <- betad + s*abs(rnorm(1, 0, 0.01))  
  
  #�Ɨ�MH�@�̑ΐ��ޓx
  lognew <- loglike(alpha_u, betan, inv_tau, d, d_log, u)$Lho
  logold <- loglike(alpha_u, betad, inv_tau, d, d_log, u)$Lho 
  
  #�p�����[�^�̍̑�������
  rand <- runif(1)   #��l���z���痐���𔭐�
  gamma <- min(c(1, exp(lognew - logold)))   #�̑𗦂�����
  gamma_u2 <- gamma

  #gamma�̒l�Ɋ�Â��V����alpha���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  beta_u <- flag*betan + (1-flag)*betad
  
  ##�A�C�e�������s��̌`��p�����[�^���T���v�����O
  #MH�@�̐V�����p�����[�^�𐶐�
  d <- as.numeric(theta_v); d_log <- log(d)
  s <- sign(dloglike_alpha(beta_v, alpha_v, d, d_log, v, item, k))   #���z�̕���
  betad <- beta_v
  betan <- betad + s*abs(rnorm(1, 0, 0.01))  
  
  #�Ɨ�MH�@�̑ΐ��ޓx
  lognew <- loglike(alpha_v, betan, inv_tau_v, d, d_log, v)$Lho
  logold <- loglike(alpha_v, betad, inv_tau_v, d, d_log, v)$Lho 
  
  #�p�����[�^�̍̑�������
  rand <- runif(1)   #��l���z���痐���𔭐�
  gamma <- min(c(1, exp(lognew - logold)))   #�̑𗦂�����
  gamma_v2 <- gamma
  
  #gamma�̒l�Ɋ�Â��V����alpha���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  beta_v <- flag*betan + (1-flag)*betad
  
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    BETA_U[mkeep] <- beta_u
    BETA_V[mkeep] <- beta_v
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx�̍X�V
    LL <- sum(dpois(y, WH, log=TRUE))
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(rbind(gamma_u1, gamma_v1), 3))
    print(round(c(gamma_u2, gamma_v2), 3))
    print(round(c(beta_u, beta_v, beta_ut, beta_vt), 3))
  }
}

####�T���v�����O���ʂ̗v��ƓK���x####
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W %*% H))[-index_z1], log=TRUE))
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W0 %*% H0))[-index_z1], log=TRUE))


