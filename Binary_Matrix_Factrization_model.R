#####�����l�̂���x�C�W�A���o�C�i���s����q����#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(condMVNorm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)
source("bdiag_m.R")

#set.seed(5897)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k <- 15   #��ꐔ
hh <- 5000   #���[�U�[��
item <- 1500   #�A�C�e����

##�p�����[�^�̐ݒ�
sigma <- 1
A <- A_T <- mvrnorm(hh, rep(-0.5, k), diag(1, k))   #���[�U�[�̓����s��
B <- B_T <- mvrnorm(item, rep(0.5, k), diag(1, k))   #�A�C�e���̓����s��
beta1 <- rbeta(hh, 8.5, 10.0)   #���[�U-�w���m��
beta2 <- rbeta(item, 5.0, 6.0)   #�A�C�e���w���m��


##���f���Ɋ�Â������ϐ��𐶐�
AB <- A %*% t(B)   #���Ғl
U <- matrix(0, nrow=hh, ncol=item)
Y <- matrix(0, nrow=hh, ncol=item)
Z <- matrix(0, nrow=hh, ncol=item)

for(j in 1:item){
  #���݌��p�𐶐�
  u_vec <- rnorm(hh, AB[, j], sigma)   #���K���z������݌��p�𐶐�
  U[, j] <- u_vec
  
  #�����𐶐�
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z[, j] <- deficit   #��������
  
  #�]���x�N�g������
  Y[, j] <- ifelse(u_vec > 0, 1, 0)   #�o�C�i���x�N�g���ɕϊ�
}

##ID�ƕ]���x�N�g����ݒ�
N <- length(Z[Z==1])
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#�]��������v�f�̂ݒ��o
index_z <- index_user <- which(as.numeric(t(Z))==1)
user_id <- user_id0[index_user]
item_id <- item_id0[index_user]
y_vec <- as.numeric(t(Y))
u_vec <- as.numeric(t(U))
mean(y_vec); sum(y_vec)


#�����l�̂���w���x�N�g���ɕϊ�
y <- y_vec[index_user]   #�����̂���]���x�N�g��
u <- u_vec[index_user]
Y <- matrix(y_vec, nrow=hh, ncol=item, byrow=T)   #�w�������̊��S�f�[�^
mean(y); length(y); sum(y)


##�C���f�b�N�X�̍쐬
index_user <- list()
index_item <- list()
for(i in 1:hh){
  index_user[[i]] <- which(user_id==i)
}
for(j in 1:item){
  index_item[[j]] <- which(item_id==j)
}


####�}���R�t�A�������e�J�����@�Ńp�����[�^�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
sigma <- 1
tau0 <- 1/100
Ca <- diag(1, k)
Cb <- diag(1, k)

#�K�w���f���̎��O���z
Bbar1 <- Bbar1 <- rep(0, k)
A1 <- B1 <- 0.01 * diag(1, k)
nu1 <- nu2 <- k/2
V1 <- nu1 * diag(k)
V2 <- nu2 * diag(k)

##�����l�̐ݒ�
Cov_A <- diag(0.5, k); inv_A <- solve(Cov_A)
Cov_B <- diag(0.01, k); inv_B <- solve(Cov_B)
A_mu <- rep(-0.5, k)
B_mu <- rep(0, k)
A <- mvrnorm(hh, rep(0.0, k), Cov_A)
B <- mvrnorm(item, rep(0.0, k), Cov_B)

##�p�����[�^�̊i�[�p�z��
A_SEG <- matrix(0, nrow=hh, ncol=k)
B_SEG <- matrix(0, nrow=item, ncol=k)
MU_A <- matrix(0, nrow=R/keep, ncol=k)
MU_B <- matrix(0, nrow=R/keep, ncol=k)
COV_A <- matrix(0, nrow=R/keep, ncol=k)
COV_B <- matrix(0, nrow=R/keep, ncol=k)

##�ؒf�̈���`
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##�ΐ��ޓx���x�X�g
mu_vec <- as.numeric(t(A_T %*% t(B_T)))[index_z]
LLbest <- sum(dbinom(y, 1, pnorm(mu_vec, 0, 1), log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z������ݕϐ����T���v�����O
  AB_vec <- as.numeric(t(A %*% t(B)))[index_z]   #���݌��p�̕��σx�N�g��
  U <- rtnorm(AB_vec, sigma, a, b)   #�ؒf���K���z������݌��p���T���v�����O
  
  ##���[�U�[�����s��̃p�����[�^���T���v�����O
  for(i in 1:hh){
    index <- item_id[index_user[[i]]]
    
    #�����s��̎��㕪�z�̃p�����[�^
    Xy <- t(B[index, ]) %*% U[index_user[[i]]]
    XXV <- t(B[index, ]) %*% B[index, ] + inv_A
    inv_XXV <- solve(XXV) 
    mu <- inv_XXV %*% (Xy + inv_A %*% A_mu)
    
    #���ϗʐ��K���z���玖�㕽�ς��T���v�����O
    A[i, ] <- mvrnorm(1, mu, inv_XXV)
  }
  A <- A - matrix(colMeans(A) - A_mu, nrow=hh, ncol=k, byrow=T)   #��x�N�g���𐳋K��
  
  ##�A�C�e�������s��̃p�����[�^���T���v�����O
  for(j in 1:item){
    index <- user_id[index_item[[j]]]
    
    #�����s��̎��㕪�z�̃p�����[�^
    Xy <- t(A[index, ]) %*% U[index_item[[j]]]
    
    XXV <- t(A[index, ]) %*% A[index, ] + inv_B
    inv_XXV <- solve(XXV) 
    mu <- inv_XXV %*% (Xy + inv_B %*% B_mu)
    
    #���ϗʐ��K���z���玖�㕽�ς��T���v�����O
    B[j, ] <- mvrnorm(1, mu, inv_XXV)
  }
  
  ##���[�U�[�����s��̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  V_par <- V2 + t(A) %*% A
  Sn <- nu1 + hh
  Cov_A <- bayesm::rwishart(Sn, solve(V_par))$IW  
  inv_A <- solve(Cov_A)
  
  
  ##�A�C�e�������s��̊K�w���f���̃p�����[�^���T���v�����O
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  V_par <- V2 + t(B) %*% B
  Sn <- nu2 + item
  Cov_B <- bayesm::rwishart(Sn, solve(V_par))$IW  
  inv_B <- solve(Cov_B)
  
  #���ϗʐ��K���z���畽�σx�N�g�����T���v�����O
  beta_mu <- item/(item + tau0) * colMeans(B)
  B_mu <- mvrnorm(1, beta_mu, Cov_B/(item + tau0))
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    MU_A[mkeep, ] <- A_mu
    MU_B[mkeep, ] <- B_mu
    COV_A[mkeep, ] <- diag(Cov_A)
    COV_B[mkeep, ] <- diag(Cov_B)
  } 
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    A_SEG <- A_SEG + A
    B_SEG <- B_SEG + B
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    mu_vec <- as.numeric(t(A %*% t(B)))[index_z]
    LL <- sum(dbinom(y, 1, pnorm(mu_vec, 0, 1), log=TRUE))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL, LLbest))
    print(round(B_mu, 3))
    print(round(diag(Cov_B), 3))
  }
}

##���S�f�[�^�ɑ΂�����덷�𐄒�
#�]���x�N�g���ƌ����x�N�g����ݒ�
z_vec <- as.numeric(t(Z))
y_vec <- as.numeric(t(Y))
mu_vec <- as.numeric(t(A %*% t(B)))

#�ϑ��f�[�^�̓��덷
er_obz <- sum((y_vec[z_vec==1] - mu_vec[z_vec==1])^2)
er_obz / sum(z_vec==1)

#�����f�[�^�̓��덷
er_na <- sum((y_vec[z_vec==0] - mu_vec[z_vec==0])^2)
er_na / sum(z_vec==0)
cbind(z_vec, y_vec, round(mu_vec, 2))
