#####�x�C�W�A���񕉒l�e���\�����q����#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

##�f�[�^�̐ݒ�
hh <- 3000   #���[�U�[��
item <- 500   #�A�C�e����
time <- 12   #�ϑ����Ԑ�
k <- 10   #��ꐔ

##�񕉒l�e���\�����q�����̉���ɏ]���f�[�^�𐶐�
#�K���}���z���p�����[�^��ݒ�
alpha01 <- 0.3; beta01 <- 0.9
alpha02 <- 0.3; beta02 <- 0.85
alpha03 <- 0.15; beta03 <- 1.5

#�K���}���z��������s��𐶐�
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
C0 <- matrix(rgamma(time*k, alpha03, beta03), nrow=k, ncol=time)

#�|�A�\�����z�̊��Ғl���v�Z
WHC0 <- array(0, dim=c(hh, item, time))
for(j in 1:k){
  WHC0 <- WHC0 + W0[, j] %o% H0[j, ] %o% C0[j, ]   #�����s��̊O��
}

#�|�A�\�����z����w���ʂ̐���
Data <- array(0, dim=c(hh, item, time))
LLbest <- 0
for(rp in 1:time){
  for(j in 1:item){
    Data[, j, rp] <- rpois(hh, WHC0[, j, rp])
  }
}
sum(Data)

#�ΐ��ޓx�̌v�Z
LLbest <- 0
for(j in 1:k){
  LLbest <- LLbest + sum(dpois(Data[, , j], WHC0[, , j], log=TRUE))
}
LLbest


####�}���R�t�A�������e�J�����@�Ŕ񕉒l�e���\�������𐄒�####
##�A���S���Y���̐ݒ�
R <- 2000
keep <- 2
disp <- 8
iter <- 0

##���O���z�̐ݒ�
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01
alpha3 <- 0.01; beta3 <- 0.01

##�����l�̐ݒ�
W <- matrix(rgamma(hh*k, 0.1, 1.0), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 1.0), nrow=k, ncol=item)
C <- matrix(rgamma(time*k, 0.1, 1.0), nrow=k, ncol=time)

##�T���v�����O���ʂ̕ۑ��p�z��
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
C_array <- array(0, dim=c(k, time, R/keep))
lambda <- array(0, dim=c(hh, item, k))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�K���}���z���W���T���v�����O
  #�|�A�\�����z�̊��Ғl
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #�⏕�ϐ���p����W���T���v�����O
  W_old <- W   #1�O��W�̓����s����i�[
  for(j in 1:k){
    lambda <- W_old[, j] %o% H[j, ] %o% C[j, ] / WHC   #�⏕�ϐ����X�V
    w1 <- alpha1 + rowSums(lambda * Data)
    w2 <- beta1 + sum(H[j, ] %o% C[j, ])
    W[, j] <- rgamma(hh, w1, w2)   #�K���}���z����W���T���v�����O
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #�e��x�N�g���𐳋K��
  
  ##�K���}���z���H���T���v�����O
  #�|�A�\�����z�̊��Ғl
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #�⏕�ϐ���p����H���T���v�����O
  H_old <- H   #1�O��H�̓����s����i�[
  for(j in 1:k){
    lambda <- W[, j] %o% H_old[j, ] %o% C[j, ] / WHC   #�⏕�ϐ����X�V
    h1 <- alpha2 + rowSums(colSums(lambda * Data))
    h2 <- beta2 + sum(W[, j] %o% C[j, ])
    H[j, ] <- rgamma(item, h1, h2)   #�K���}���z����H���T���v�����O
  }
  
  ##�K���}���z���C���T���v�����O
  #�|�A�\�����z�̊��Ғl
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #�⏕�ϐ���p����C���T���v�����O
  C_old <- C    #1�O��C�̓����s����i�[
  for(j in 1:k){
    lambda <- W[, j] %o% H[j, ] %o% C_old[j, ] / WHC   #�⏕�ϐ����X�V
    c1 <- alpha3 + colSums(colSums(lambda * Data))
    c2 <- beta3 + sum(W[, j] %o% H[j, ])
    C[j, ] <- rgamma(time, c1, c2)   #�K���}���z����H���T���v�����O
  }
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep ] <- H
    C_array[, , mkeep] <- C
  }
    
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- 0
    for(j in 1:k){
      LL <- LL + sum(dpois(Data[, , j], WHC[, , j], log=TRUE))
    }
    print(rp)
    print(c(LL, LLbest))
  }
}