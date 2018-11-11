#####�񕉒l�s����q����#####
library(MASS)
library(NMF)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####�f�[�^�̔���####
seg <- 15   #�Z�O�����g��
n <- 50   #�Z�O�����g���Ƃ̃T���v����
N <- seg*n   #�S�T���v����
category <- 100   #�J�e�S���[��
lambda <- 70   #�w�����̃|�A�\������
y <- rpois(N, lambda)   #�w��������
Y <- data.frame(seg=rep(1:seg, rep(n, seg)), y)

#�w���p�x�s��𔭐�������
X.freq <- matrix(0, 0, ncol=category)
for(i in 1:seg){
  p <- runif(category, 0, 5)
  freq <- t(apply(Y[Y$seg==i, ], 1, function(x) rmultinom(1, x[2], p)))
  X.freq <- rbind(X.freq, freq)
}
w.freq <- colSums(X.freq)

####�񕉒l�s����q�����𐄒�####
#�p�����[�^�̏����l��ݒ�
H <- matrix(runif(N*seg, 0, 3), nrow=N, ncol=seg)
U <- matrix(runif(N*seg, 0, 3), nrow=seg, ncol=category) 
HU <- H %*% U

#�A���S���Y���̐ݒ�
tol <- 0.1
diff <- 100
SE1 <- sum((X.freq - HU)^2)

##�p�����[�^���X�V(���덷�)
while(diff >= tol){
  #�A���S���Y��
  U1 <- U * t(H) %*% X.freq / t(H) %*% H %*% U
  H1  <- H * X.freq %*% t(U1) / H %*% U1 %*% t(U1)
  
  #�p�����[�^���X�V
  HU <- H1 %*% U1
  SE <- sum((X.freq - HU)^2)
  U <- U1
  H <- H1
  diff <- abs(SE1-SE) 
  SE1 <- SE
  print(diff)
}

#���肳�ꂽ�p�����[�^
d <- cbind(as.numeric(X.freq), round(as.numeric(HU), 3))
round(H, 2)
round(U, 2)
round(colSums(HU), 0)
colSums(X.freq)

#�֐��Ő���
res <- nmf(X.freq, seg)
H_fun <- basis(res)
U_fun <- coef(res)
HU_fun <- H_fun %*% U_fun
f <- cbind(as.numeric(X.freq), round(as.numeric(HU), 3), round(as.numeric(HU_fun), 3))