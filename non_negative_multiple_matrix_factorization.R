#####�񕉒l���d�s����q����#####
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

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
hh <- 3000   #���[�U�[��
item <- 300   #�A�C�e����
group <- 20   #�O���[�v��
category <- 15   #�J�e�S���[��
k <- 10   #���ݕϐ���

#�Ή��s��̐ݒ�
pr1 <- runif(group, 0.25, 1.0)
pr2 <- runif(category, 0.25, 1.0)
V <- rmnom(hh, 1, pr1)   #���[�U�[�O���[�v�Ή��s��
W <- rmnom(item, 1, pr2)   #���i�J�e�S���Ή��s��


##�񕉒l�s����q�����̉���ɏ]���f�[�^�𐶐�
#�K���}���z���p�����[�^��ݒ�
A0 <- matrix(0, nrow=hh, ncol=k)
B0 <- matrix(0, nrow=k, ncol=item)
alpha01 <- matrix(0, nrow=k, ncol=group)
beta01 <- matrix(0, nrow=k, ncol=group)
alpha02 <- matrix(0, nrow=k, ncol=category)
beta02 <- matrix(0, nrow=k, ncol=category)

for(i in 1:k){
  #���[�U�[�����s��𐶐�
  for(j1 in 1:group){
    alpha01[i, j1] <- runif(1, 0.01, 0.4)
    beta01[i, j1] <- runif(1, 0.6, 1.8)
    A0[V[, j1]==1, i] <- rgamma(sum(V[, j1]), alpha01[i, j1], beta01[i, j1])
  }
  #�A�C�e�������s��𐶐�
  for(j2 in 1:category){
    alpha02[i, j2] <- runif(1, 0.01, 0.3)
    beta02[i, j2] <- runif(1, 0.6, 1.5)
    B0[i, W[, j2]==1] <- rgamma(sum(W[, j2]), alpha02[i, j2], beta02[i, j2])
  }
}

#�|�A�\�����z�̊��Ғl��ݒ�
AB <- A0 %*% B0

#�O���[�v�����s�񂨂�уJ�e�S���[�����s��𐶐�
C0 <- t(t(W) %*% t(B0))
D0 <- t(V) %*% A0 



##�J���o�b�N���C�u���[�_�C�o�[�W�F���X�����Ɋ�Â������ϐ��𐶐�
#���[�U�[���i�w���s��𐶐�
X <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  X[, j] <- rpois(hh, AB[, j])
}
rowSums(X)
sum(X)

#���[�U�[�J�e�S���w���s���ݒ�
Y <- matrix(0, nrow=hh, ncol=category)
for(j in 1:category){
  Y[, j] <- rowSums(X[, W[, j]==1])
}
colSums(Y)

#�O���[�v���i�w���s���ݒ�
Z <- matrix(0, nrow=group, ncol=item)
for(j in 1:group){
  Z[j, ] <- colSums(X[V[, j]==1, ])
}

#�^�̃p�����[�^�ł̑ΐ��ޓx
LL <- sum(dpois(X, A0 %*% B0, log=TRUE))


####�}���R�t�A�������e�J�����@��MMNF�𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2
iter <- 0
disp <- 10

##���O���z�̐ݒ�
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##�����l�̐ݒ�
A <- matrix(rgamma(hh*k, 0.25, 0.5), nrow=hh, ncol=k)
B <- matrix(rgamma(item*k, 0.25, 0.5), nrow=k, ncol=item)
C <- matrix(rgamma(category*k, 0.25, 0.5), nrow=k, ncol=category)
D <- matrix(rgamma(group*k, 0.25, 0.5), nrow=group, ncol=k)


##�T���v�����O���ʂ̕ۑ��p�z��
A_array <- array(0, dim=c(hh, k, R/keep))
B_array <- array(0, dim=c(k, item, R/keep))
C_array <- array(0, dim=c(k, category, R/keep))
D_array <- array(0, dim=c(group, k, R/keep))


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�⏕�ϐ�lambda���X�V
  lambda1 <- array(0, dim=c(hh, item, k))
  lambda2 <- array(0, dim=c(hh, category, k))
  AB <- A %*% B
  AC <- A %*% C
  for(j in 1:k){
    lambda1[, , j] <- A[, j] %*% t(B[j, ]) / AB
    lambda2[, , j] <- A[, j] %*% t(C[j, ]) / AC
  }
  
  ##�K���}���z��胆�[�U�[�����s����T���v�����O
  for(j in 1:k){
    a1 <- alpha1 + (rowSums(lambda1[, , j] * X) + rowSums(lambda2[, , j] * Y))/2
    a2 <- sum(B[j, ])
    A[, j] <- rgamma(hh, a1, a2)   
  }
  A <- A / matrix(colSums(A), nrow=hh, ncol=k, byrow=T) * hh/5   #�e��x�N�g���𐳋K��
  
  ##�⏕�ϐ�lambda���X�V
  lambda1 <- array(0, dim=c(hh, item, k))
  lambda2 <- array(0, dim=c(group, item, k))
  AB <- A %*% B
  BD <- D %*% B
  for(j in 1:k){
    lambda1[, , j] <- A[, j] %*% t(B[j, ]) / AB
    lambda2[, , j] <- D[, j] %*% t(B[j, ]) / BD
  }
  
  ##�K���}���z���A�C�e�������s����T���v�����O
  for(j in 1:k){
    b1 <- alpha2 + (colSums(lambda1[, , j] * X) + colSums(lambda2[, , j] * Z))/2
    b2 <- beta2 + sum(A[, j])
    B[j, ] <- rgamma(item, b1, b2)  
  }
  
  ##���`����������O���[�v�����s�񂨂�уJ�e�S���[�����s����X�V
  C <- t(t(W) %*% t(B))   #�J�e�S���[�����s��̍X�V
  D <- t(V) %*% A   #�O���[�v�����s��̍X�V
  
  
  ##�T���v�����O���ʂ̕ۑ��ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    A_array[, , mkeep] <- A
    B_array[, , mkeep ] <- B
    C_array[, , mkeep] <- C
    D_array[, , mkeep ] <- D
    if(rp%%disp==0){
      print(c(sum(dpois(X, A %*% B, log=T)), LL))
      print(rp)
    }
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2000/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
#���x�N�g��A�̃p�����[�^�̉���
matplot(t(A_array[1:10, 1, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[1:10, 2, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[1:10, 3, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[1:10, 4, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[1:10, 5, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[11:20, 6, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[11:20, 7, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[11:20, 8, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[11:20, 9, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")
matplot(t(A_array[11:20, 10, ]), type="l", xlab="�T���v�����O��", ylab="���x�N�g���̃p�����[�^")

#�����W��B�̃p�����[�^�̉���
matplot(t(B_array[1, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[2, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[3, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[4, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[5, 1:10, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[6, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[7, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[8, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[9, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(B_array[10, 11:20, ]), type="l", xlab="�T���v�����O��", ylab="�����W��")

#�����W��B�̃p�����[�^�̉���
matplot(t(C_array[1, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[2, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[3, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[4, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[5, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[6, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[7, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[8, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[9, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(C_array[10, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")

#�����W��B�̃p�����[�^�̉���
matplot(t(D_array[1, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[2, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[3, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[4, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[5, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[6, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[7, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[8, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[9, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")
matplot(t(D_array[10, , ]), type="l", xlab="�T���v�����O��", ylab="�����W��")