#####Probit regression based Tensor Factorization model#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

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
hh <- 10000   #���[�U�[��
item <- 2000   #�A�C�e����
context <- 20   #�R���e�L�X�g��
N0 <- hh*item*context
k <- 10   #��ꐔ
vec <- rep(1, k)

#ID��ݒ�
user_id0 <- rep(1:hh, rep(item*context, hh))
item_id0 <- rep(rep(1:item, context), hh)
context_id0 <- rep(rep(1:context, rep(item, context)), hh)

##�����x�N�g���𐶐�
#�����m���𐶐�
user_prob <- rbeta(hh, 13.0, 47.5)
item_prob <- rbeta(item, 13.5, 52.5)
context_prob <- rbeta(context, 9.0, 45.0)
prob <- user_prob[user_id0]*item_prob[item_id0]*context_prob[context_id0]

#�x���k�[�C���z���猇���x�N�g���𐶐�
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)

#�����x�N�g������id���č\��
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
context_id <- context_id0[z_vec==1]
rm(user_id0); rm(item_id0); rm(context_id0); rm(z_vec); rm(prob)
gc(); gc()

#�w�������J�E���g
freq <- plyr::count(user_id); freq_user <- freq$freq[order(freq$x)]
freq <- plyr::count(item_id); freq_item <- freq$freq[order(freq$x)]
freq <- plyr::count(context_id); freq_context <- freq$freq[order(freq$x)]
hist(freq_user, col="grey", breaks=25, main="���[�U�[���Ƃ̍w����", xlab="�w����")
hist(freq_item, col="grey", breaks=25, main="�A�C�e�����Ƃ̍w����", xlab="�w����")


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


####�����ϐ��𐶐�####
for(rp in 1:1000){
  print(rp)
  
  ##���[�U�[�x�[�X�̊K�w���f���̃p�����[�^
  #�K�w���f���̕��U�p�����[�^
  Cov_ut1 <- Cov_u1 <- runif(1, 0.25, 0.4)
  Cov_ut2 <- Cov_u2 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  Cov_ut3 <- Cov_u3 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  
  #��A�W����ݒ�
  alpha_u1 <- rep(0, ncol(u))
  alpha_u2 <- alpha_u3 <- matrix(0, nrow=ncol(u), ncol=k)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u1[j] <- runif(1, -0.6, -0.2)
      alpha_u2[j, ] <- runif(k, -0.7, -0.1)
      alpha_u3[j, ] <- runif(k, -0.7, -0.1)
      
    } else {
      alpha_u1[j] <- runif(1, -0.6, 0.4)
      alpha_u2[j, ] <- runif(k, -0.6, 0.5)
      alpha_u3[j, ] <- runif(k, -0.6, 0.5)
    }
  }
  alpha_ut1 <- alpha_u1; alpha_ut2 <- alpha_u2; alpha_ut3 <- alpha_u3
  
  #���ϗʉ�A���f�����烆�[�U�[�ʂ̉�A�p�����[�^�𐶐�
  theta_u1 <- theta_ut1 <- as.numeric(u %*% alpha_u1 + rnorm(hh, 0, Cov_u1))   #�ϗʌ��ʂ̃p�����[�^
  theta_u2 <- theta_ut2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)   #�s�񕪉��̃p�����[�^
  theta_u3 <- theta_ut3 <- u %*% alpha_u3 + mvrnorm(hh, rep(0, k), Cov_u3)   #�e���\�������̃p�����[�^
  
  
  ##�A�C�e���x�[�X�̊K�w���f���̃p�����[�^
  #���U�����U�s���ݒ�
  Cov_vt1 <- Cov_v1 <- runif(1, 0.25, 0.4)
  Cov_vt2 <- Cov_v2 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  Cov_vt3 <- Cov_v3 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  
  #��A�W����ݒ�
  alpha_v1 <- rep(0, ncol(v))
  alpha_v2 <- alpha_v3 <- matrix(0, nrow=ncol(v), ncol=k)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v1[j] <- runif(1, -0.5, -0.2)
      alpha_v2[j, ] <- runif(k, -0.7, 0.4)
      alpha_v3[j, ] <- runif(k, -0.7, 0.4)
    } else {
      alpha_v1[j] <- runif(1, -0.7, 0.5)
      alpha_v2[j, ] <- runif(k, -0.7, 0.5)
      alpha_v3[j, ] <- runif(k, -0.7, 0.5)
    }
  }
  alpha_vt1 <- alpha_v1; alpha_vt2 <- alpha_v2; alpha_vt3 <- alpha_v3
  
  #���ϗʉ�A���f�����烆�[�U�[�ʂ̉�A�p�����[�^�𐶐�
  theta_v1 <- theta_vt1 <- as.numeric(v %*% alpha_v1 + rnorm(item, 0, Cov_v1))   #�ϗʌ��ʂ̃p�����[�^
  theta_v2 <- theta_vt2 <- v %*% alpha_v2 + mvrnorm(item, rep(0, k), Cov_v2)   #�s�񕪉��̃p�����[�^
  theta_v3 <- theta_vt3 <- v %*% alpha_v3 + mvrnorm(item, rep(0, k), Cov_v3)   #�e���\�������̃p�����[�^
  
  ##�R���e�L�X�g�x�[�X�̊K�w���f���̃p�����[�^
  alpha_c1 <- alpha_ct1 <- 0
  alpha_c3 <- alpha_ct3 <- alpha_c2 <- alpha_ct2 <- rep(0, k)
  Cov_c1 <- Cov_ct1 <- runif(1, 0.25, 0.4)
  Cov_c2 <- Cov_ct2 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  Cov_c3 <- Cov_ct3 <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  theta_c1 <- theta_ct1 <- rnorm(context, alpha_c1, Cov_c1)
  theta_c2 <- theta_ct2 <- mvrnorm(context, alpha_c2, Cov_c2)
  theta_c3 <- theta_ct3 <- mvrnorm(context, alpha_c3, Cov_c3)


  ##���K���z������p�ƍw���x�N�g���𐶐�
  #�s�񕪉��̃p�����[�^�𐶐�
  uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
  
  #�e���\�������̃p�����[�^�𐶐�
  uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)
  
  #���݌��p�𐶐�
  mu <- theta_u1[user_id] + theta_v1[item_id] + theta_c2[context_id] + uv + uc + vc + uvc   #���Ғl
  U <- mu + rnorm(N, 0, 1)   #�덷�𐶐�
  
  #�w���x�N�g���ɕϊ�
  y <- ifelse(U > 0, 1, 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break����
}


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�e���\�������𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

##�C���f�b�N�X��ݒ�
user_index <- item_index <- context_index <- list()
ui_id <- uc_id <- iu_id <- ic_id <- cu_id <- ci_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  uc_id[[i]] <- context_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  ic_id[[j]] <- context_id[item_index[[j]]]
}
for(j in 1:context){
  context_index[[j]] <- which(context_id==j)
  cu_id[[j]] <- user_id[context_index[[j]]]
  ci_id[[j]] <- item_id[context_index[[j]]]
}
vec <- rep(1, k)
const1 <- hh / 2.0  #���K���萔
const2 <- item / 2.0

##���O���z�̐ݒ�
#�ϗʌ��ʂ̊K�w���f���̎��O���z
gamma_u <- rep(0, ncol(u)); tau_u <- 100 * diag(ncol(u)); inv_tau_u <- solve(tau_u)
gamma_v <- rep(0, ncol(v)); tau_v <- 100 * diag(ncol(v)); inv_tau_v <- solve(tau_v)
gamma_c <- 0; tau_c <- 100; inv_tau_c <- 1/tau_c
v0 <- 1; s0 <- 1

#���[�U�[�̊K�w���f���̎��O���z
Deltabar1 <- matrix(rep(0, ncol(u)*k), nrow=ncol(u), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta1 <- 0.01 * diag(rep(1, ncol(u)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̊K�w���f���̎��O���z
Deltabar2 <- matrix(rep(0, ncol(v)*k), nrow=ncol(v), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta2 <- 0.01 * diag(rep(1, ncol(v)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#�R���e�L�X�g�̊K�w���f���̎��O���z
Deltabar3 <- rep(0, k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta3 <- 0.01 * diag(k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu3 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V3 <- nu3 * diag(rep(1, k))   #�t�E�B�V���[�g���z�̃p�����[�^


##�p�����[�^�̐^�l
alpha_u1 <- alpha_ut1; alpha_u2 <- alpha_ut2; alpha_u3 <- alpha_ut3
user_mu1 <- as.numeric(u %*% alpha_u1); user_mu2 <- u %*% alpha_u2; user_mu3 <- u %*% alpha_u3
alpha_v1 <- alpha_vt1; alpha_v2 <- alpha_vt2; alpha_v3 <- alpha_vt3
item_mu1 <- as.numeric(v %*% alpha_v1); item_mu2 <- v %*% alpha_v2; item_mu3 <- v %*% alpha_v3
alpha_c1 <- alpha_ct1; alpha_c2 <- alpha_ct2; alpha_c3 <- alpha_ct3
Cov_u1 <- Cov_ut1; Cov_u2 <- Cov_ut2; Cov_u3 <- Cov_ut3; inv_Cov_u2 <- solve(Cov_u2); inv_Cov_u3 <- solve(Cov_u3)
Cov_v1 <- Cov_vt1; Cov_v2 <- Cov_vt2; Cov_v3 <- Cov_vt3; inv_Cov_v2 <- solve(Cov_v2); inv_Cov_v3 <- solve(Cov_v3)
Cov_c1 <- Cov_ct1; Cov_c2 <- Cov_ct2; Cov_c3 <- Cov_ct3; inv_Cov_c2 <- solve(Cov_c2); inv_Cov_c3 <- solve(Cov_c3)
theta_u1 <- theta_ut1; theta_u2 <- theta_ut2; theta_u3 <- theta_ut3
theta_v1 <- theta_vt1; theta_v2 <- theta_vt2; theta_v3 <- theta_vt3
theta_c1 <- theta_ct1; theta_c2 <- theta_ct2; theta_c3 <- theta_ct3
sigma <- 1

#�s�񕪉��ƃe���\�������̃p�����[�^
uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)


##�p�����[�^�̏����l
#�K�w���f���̃p�����[�^
alpha_u1 <- runif(ncol(u), -0.3, 0.3)
alpha_u2 <- matrix(runif(ncol(u)*k, -0.3, 0.3), nrow=ncol(u), ncol=k)
alpha_u3 <- matrix(runif(ncol(u)*k, -0.3, 0.3), nrow=ncol(u), ncol=k)
alpha_v1 <- runif(ncol(v), -0.3, 0.3)
alpha_v2 <- matrix(runif(ncol(v)*k, -0.3, 0.3), nrow=ncol(v), ncol=k)
alpha_v3 <- matrix(runif(ncol(v)*k, -0.3, 0.3), nrow=ncol(v), ncol=k)
alpha_c1 <- 0; alpha_c2 <- alpha_c3 <- rep(0, k)
Cov_u1 <- 0.25; Cov_u2 <- 0.05 * diag(k); Cov_u3 <- 0.05 * diag(k)
Cov_v1 <- 0.25; Cov_v2 <- 0.05 * diag(k); Cov_v3 <- 0.05 * diag(k)
Cov_c1 <- 0.25; Cov_c2 <- 0.05 * diag(k); Cov_c3 <- 0.05 * diag(k)

#�ϗʌ��ʂ̃p�����[�^
theta_u1 <- u %*% alpha_u1 + rnorm(hh, 0, Cov_u1)
theta_v1 <- v %*% alpha_v1 + rnorm(item, 0, Cov_v1)
theta_c1 <- rnorm(context, alpha_c1, Cov_c1)

#�s�񕪉��̃p�����[�^
theta_u2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)
theta_v2 <- v %*% alpha_v2 + mvrnorm(item, rep(0, k), Cov_v2)
theta_c2 <- mvrnorm(context, alpha_c2, Cov_c2)
uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)

#�e���\�������̃p�����[�^
theta_u3 <- u %*% alpha_u3 + mvrnorm(hh, rep(0, k), Cov_u3)
theta_v3 <- v %*% alpha_v3 + mvrnorm(item, rep(0, k), Cov_v3)
theta_t3 <- mvrnorm(context, alpha_c3, Cov_c3)
uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)
sigma <- 1


##�T���v�����O���ʂ̃p�����[�^�̕ۑ��p�z��
#���f���p�����[�^�̊i�[�p�z��
THETA_U1 <- matrix(0, nrow=R/keep, ncol=hh)
THETA_U2 <- array(0, dim=c(hh, k, R/keep))
THETA_U3 <- array(0, dim=c(hh, k, R/keep))
THETA_V1 <- matrix(0, nrow=R/keep, ncol=item)
THETA_V2 <- array(0, dim=c(item, k, R/keep))
THETA_V3 <- array(0, dim=c(item, k, R/keep))
THETA_C1 <- matrix(0, nrow=R/keep, ncol=context)
THETA_C2 <- array(0, dim=c(context, k, R/keep))
THETA_C3 <- array(0, dim=c(context, k, R/keep))

#�K�w���f���̊i�[�p�z��
ALPHA_U1 <- matrix(0, nrow=R/keep, ncol=ncol(u))
ALPHA_U2 <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_U3 <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_V1 <- matrix(0, nrow=R/keep, ncol=ncol(v))
ALPHA_V2 <- array(0, dim=c(ncol(v), k, R/keep))
ALPHA_V3 <- array(0, dim=c(ncol(v), k, R/keep))
COV_U1 <- COV_V1 <- COV_C1 <- rep(0, R/keep)
COV_U2 <- COV_V2 <- COV_C2 <- array(0, dim=c(k, k, R/keep))
COV_U3 <- COV_V3 <- COV_C3 <- array(0, dim=c(k, k, R/keep))


##�f�[�^�̒萔��ݒ�
uu <- t(u) %*% u + inv_tau_u; inv_uu <- solve(uu)
vv <- t(v) %*% v + inv_tau_v; inv_vv <- solve(vv)

##�ؒf�̈���`
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##�ΐ��ޓx�̊�l
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #�ΐ��ޓx


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O
for(rp in 1:R){
  
  ##�ؒf���K���z������݌��p�𐶐�
  theta_u_vec1 <- theta_u1[user_id]; theta_v_vec1 <- theta_v1[item_id]; theta_c_vec1 <- theta_c1[context_id]
  mu <- theta_u_vec1 + theta_v_vec1 + theta_c_vec1 + uv + uc + vc + uvc   #���݌��p�̊��Ғl
  U <- extraDistr::rtnorm(N, mu, sigma, a, b)   #���݌��p�𐶐�
  
  
  ##���[�U�[�̕ϗʌ��ʂ��T���v�����O
  #���f���̉����ϐ�
  u_er <- U - theta_v_vec1 - theta_c_vec1 - uv - uc - vc - uvc   
  
  #���[�U�[�̕ϗʌ��ʂ̎��㕪�z�̃p�����[�^
  u_mu <- rep(0, hh)
  for(i in 1:hh){
    u_mu[i] <- mean(u_er[user_index[[i]]])
  }
  weights <- Cov_u1^2 / (sigma^2/freq_user + Cov_u1^2)   #�d�݌W��
  mu_par <- weights*u_mu + (1-weights)*user_mu1   #���㕪�z�̕���
  
  #���K���z��莖�㕪�z���T���v�����O
  theta_u1 <- rnorm(hh, mu_par, sqrt(1 / (1/Cov_u1^2 + freq_user/sigma^2)))
  theta_u_vec1 <- theta_u1[user_id]
  
  ##�A�C�e���̕ϗʌ��ʂ��T���v�����O
  #���f���̉����ϐ�
  v_er <- U - theta_u_vec1 - theta_c_vec1 - uv - uc - vc - uvc
  
  #�A�C�e���̕ϗʌ��ʂ̎��㕪�z�̃p�����[�^
  v_mu <- rep(0, item)
  for(j in 1:item){
    v_mu[j] <- mean(v_er[item_index[[j]]])
  }
  weights <- Cov_v1^2 / (sigma^2/freq_item + Cov_v1^2)   #�d�݌W��
  mu_par <- weights*v_mu + (1-weights)*item_mu1   #���㕪�z�̕���
  
  #���K���z��莖�㕪�z���T���v�����O
  theta_v1 <- rnorm(item, mu_par, sqrt(1 / (1/Cov_v1^2 + freq_item/sigma^2)))
  theta_v_vec1 <- theta_v1[item_id]
  
  ##�R���e�L�X�g�̕ϗʌ��ʂ��T���v�����O
  #���f���̉����ϐ�
  c_er <- U - theta_u_vec1 - theta_v_vec1 - uv - uc - vc - uvc
  
  #�A�C�e���̕ϗʌ��ʂ̎��㕪�z�̃p�����[�^
  c_mu <- rep(0, context)
  for(j in 1:context){
    c_mu[j] <- mean(c_er[context_index[[j]]])
  }
  weights <- Cov_c1^2 / (sigma^2/freq_context + Cov_c1^2)   #�d�݌W��
  mu_par <- weights*c_mu + (1-weights)*alpha_c1   #���㕪�z�̕���
  
  #���K���z��莖�㕪�z���T���v�����O
  theta_c1 <- rnorm(context, mu_par, sqrt(1 / (1/Cov_c1^2 + freq_context/sigma^2)))
  theta_c_vec1 <- theta_c1[context_id]
  
  
  ##���[�U�[�̓����s����T���v�����O
  #���f���̉����ϐ�
  u_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - vc - uvc
  
  for(i in 1:hh){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_v2[ui_id[[i]], ] + theta_c2[uc_id[[i]], ]
    Xy <- t(X) %*% u_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u2)
    mu <- inv_XXV %*% (Xy + inv_Cov_u2 %*% user_mu2[i, ])   #���㕪�z�̕���
    
    #���ϗʐ��K���z���烆�[�U�[�����s����T���v�����O
    theta_u2[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #�s�񕪉��̃p�����[�^���X�V
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  
  ##�A�C�e���̓����s����T���v�����O
  #���f���̉����ϐ�
  v_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uc - uvc
  
  for(j in 1:item){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_u2[iu_id[[j]], ] + theta_c2[ic_id[[j]], ]
    Xy <- t(X) %*% v_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v2)
    mu <- inv_XXV %*% (Xy + inv_Cov_v2 %*% item_mu2[j, ])   #���㕪�z�̕���
    
    #���ϗʐ��K���z���烆�[�U�[�����s����T���v�����O
    theta_v2[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #�s�񕪉��̃p�����[�^���X�V
  uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
  
  ##�R���e�L�X�g�̓����s����T���v�����O
  #���f���̉����ϐ�
  c_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uv - uvc
  
  for(j in 1:context){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_u2[cu_id[[j]], ] + theta_v2[ci_id[[j]], ]
    Xy <- t(X) %*% c_er[context_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_c2)
    mu <- inv_XXV %*% (Xy + inv_Cov_c2 %*% alpha_c2)   #���㕪�z�̕���
    
    #���ϗʐ��K���z���烆�[�U�[�����s����T���v�����O
    theta_c2[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #�s�񕪉��̃p�����[�^���X�V
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
  
  
  ##���[�U�[�̃e���\���̃p�����[�^���T���v�����O
  #���f���̉����ϐ�
  y_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uv - uc - vc
  
  for(i in 1:hh){
    #�����x�N�g���̃p�����[�^
    X <- theta_v3[ui_id[[i]], ] * theta_c3[uc_id[[i]], ]
    Xy <- t(X) %*% y_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_u3 %*% user_mu3[i, ])   #���ϗʐ��K���z�̕��σx�N�g��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_u3[i, ] <- mvrnorm(1, beta_mu, sigma^2*inv_XXV)
  }
  
  ##�A�C�e���̃e���\���̃p�����[�^���T���v�����O
  for(j in 1:item){
    #�����x�N�g���̃p�����[�^
    X <- theta_u3[iu_id[[j]], ] * theta_c3[ic_id[[j]], ]
    Xy <- t(X) %*% y_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_v3 %*% item_mu3[j, ])   #���ϗʐ��K���z�̕��σx�N�g��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_v3[j, ] <- mvrnorm(1, beta_mu, sigma^2*inv_XXV)
  }
  
  ##�R���e�L�X�g�̃e���\���̃p�����[�^���T���v�����O
  for(j in 1:context){
    #�����x�N�g���̃p�����[�^
    X <- theta_u3[cu_id[[j]], ] * theta_v3[ci_id[[j]], ]
    Xy <- t(X) %*% y_er[context_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_c3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_c3 %*% alpha_c3)   #���ϗʐ��K���z�̕��σx�N�g��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_c3[j, ] <- mvrnorm(1, beta_mu, sigma^2*inv_XXV)
  }
  
  
  ##���[�U�[�̕ϗʌ��ʂ̊K�w���f���̃p�����[�^���T���v�����O
  #���[�U�[�̕ϗʌ��ʂ̊K�w���f���̉�A�x�N�g�����X�V
  beta_mu <- inv_uu %*% (t(u) %*% theta_u1 + inv_tau_u %*% gamma_u)
  alpha_u1 <- mvrnorm(1, beta_mu, Cov_u1^2*inv_uu)   #���ϗʐ��K���z�����A�x�N�g�����T���v�����O
  user_mu1 <- as.numeric(u %*% alpha_u1)
  
  #���[�U�[�̊K�w���f���̕W���΍����X�V
  er <- theta_u1 - user_mu1
  s1 <- t(er) %*% er + s0; v1 <- hh + v0   #�K���}���z�̃p�����[�^
  Cov_u1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #�K���}���z����W���΍����T���v�����O
  inv_Cov_u1 <- 1 / Cov_u1
  
  ##�A�C�e���̕ϗʌ��ʂ̊K�w���f���̃p�����[�^���T���v�����O
  #�A�C�e���̕ϗʌ��ʂ̊K�w���f���̉�A�x�N�g�����X�V
  beta_mu <- inv_vv %*% (t(v) %*% theta_v1 + inv_tau_v %*% gamma_v)
  alpha_v1 <- mvrnorm(1, beta_mu, Cov_v1^2*inv_vv)   #���ϗʐ��K���z�����A�x�N�g�����T���v�����O
  item_mu1 <- as.numeric(v %*% alpha_v1)
  
  #���[�U�[�̊K�w���f���̕W���΍����X�V
  er <- theta_v1 - item_mu1
  s1 <- t(er) %*% er + s0; v1 <- item + v0   #�K���}���z�̃p�����[�^
  Cov_v1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #�K���}���z����W���΍����T���v�����O
  inv_Cov_v1 <- 1 / Cov_v1
  
  ##�R���e�L�X�g�̕ϗʌ��ʂ̊K�w���f���̃p�����[�^���T���v�����O
  #�R���e�L�X�g�̕ϗʌ��ʂ̊K�w���f���̕��σx�N�g�����X�V
  mu <- mean(theta_c1)
  weights <- tau_c^2 / (Cov_c1^2/context + tau_c^2)   #�d�݌W��
  mu_par <- weights*mu   #���㕪�z�̕���
  alpha_c1 <- rnorm(1, mu_par, sqrt(1 / (1/tau_c^2 + context/Cov_c1^2)))   #���K���z����ϗʌ��ʂ��T���v�����O
  
  #�R���e�L�X�g�̊K�w���f���̕W���΍����X�V
  er <- theta_c1 - alpha_c1
  s1 <- t(er) %*% er + s0; v1 <- context + v0   #�K���}���z�̃p�����[�^
  Cov_c1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #�K���}���z����W���΍����T���v�����O
  inv_Cov_c1 <- 1 / Cov_c1
  
  
  ##���[�U�[�����s��̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f�����烆�[�U�[�̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_u2, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u2 <- out$B
  user_mu2 <- u %*% alpha_u2   #���[�U�[�����s��̕��ύ\��
  Cov_u2 <- out$Sigma
  inv_cov_u2 <- solve(Cov_u2)
  
  ##�A�C�e�������s��̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������A�C�e���̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_v2, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v2 <- out$B
  item_mu2 <- v %*% alpha_v2   #�A�C�e�������s��̕��ύ\��
  Cov_v2 <- out$Sigma
  inv_cov_v2 <- solve(Cov_v2)
  
  ##�R���e�L�X�g�����s��̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  IW <- t(theta_c2) %*% theta_c2 + solve(V3)
  Sn <- nu3 + context
  Cov_c2 <- rwishart(Sn, solve(IW))$IW   #�t�E�B�V���[�g���z����p�����[�^���T���v�����O
  inv_Cov_c2 <- solve(Cov_c2)
  
  
  ##���[�U�[�e���\���̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f�����烆�[�U�[�̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_u3, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u3 <- out$B
  user_mu3 <- u %*% alpha_u3   #���[�U�[�e���\���̕��ύ\��
  Cov_u3 <- out$Sigma
  inv_cov_u3 <- solve(Cov_u3)
  
  ##�A�C�e���e���\���̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������A�C�e���̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_v3, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v3 <- out$B
  item_mu3 <- v %*% alpha_v3   #�A�C�e�������s��̕��ύ\��
  Cov_v3 <- out$Sigma
  inv_cov_v3 <- solve(Cov_v3)
  
  ##�R���e�L�X�g�e���\���̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  IW <- t(theta_c3) %*% theta_c3 + solve(V3)
  Sn <- nu3 + context
  Cov_c3 <- rwishart(Sn, solve(IW))$IW   #�t�E�B�V���[�g���z����p�����[�^���T���v�����O
  inv_Cov_c3 <- solve(Cov_c3)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep

    #���f���p�����[�^�̊i�[
    THETA_U1[mkeep, ] <- theta_u1 ; THETA_U2[, , mkeep] <- theta_u2; THETA_U3[, , mkeep] <- theta_u3
    THETA_V1[mkeep, ] <- theta_v1; THETA_V2[, , mkeep] <- theta_v2; THETA_V3[, , mkeep] <- theta_v3
    THETA_C1[mkeep, ] <- theta_c1; THETA_C2[, , mkeep] <- theta_c2; THETA_C3[, , mkeep] <- theta_c3
    
    #�K�w���f���̃p�����[�^���i�[
    ALPHA_U1[mkeep, ] <- alpha_u1; ALPHA_U2[, , mkeep] <- alpha_u2; ALPHA_U3[, , mkeep] <- alpha_u3
    ALPHA_V1[mkeep, ] <- alpha_v1; ALPHA_V2[, , mkeep] <- alpha_v2; ALPHA_V3[, , mkeep] <- alpha_v3
    COV_U1[mkeep] <- Cov_u1; COV_U2[, , mkeep] <- Cov_u2; COV_U3[, , mkeep] <- Cov_u3
    COV_V1[mkeep] <- Cov_v1; COV_V2[, , mkeep] <- Cov_v2; COV_V3[, , mkeep] <- Cov_v3
    COV_C1[mkeep] <- Cov_c1; COV_C2[, , mkeep] <- Cov_c2; COV_C3[, , mkeep] <- Cov_c3
  }

  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    mu <- theta_u_vec1 + theta_v_vec1 + theta_c_vec1 + uv + uc + vc + uvc   #���݌��p�̊��Ғl
    prob <- pnorm(mu, 0, sigma)   #�w���m��
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #�ΐ��ޓx
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLst))
  }
}