#####Bayesian Multinomial Factorization Machines#####
library(MASS)
library(matrixStats)
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
select <- 10   #�I������
s <- 5   #��ꐔ
N <- 20000   #�T���v����
topic <- 30   #�g�s�b�N��
k <- 3   #�f�[�^�^�C�v��

###�����ϐ����Ó��ɂȂ�܂Ńf�[�^�̐������J��Ԃ�
rp <- 0
repeat { 
  print(rp <- rp + 1)
  
  ##�����ϐ��̐���
  u <- array(0, dim=c(N, topic, k))
  x <- matrix(0, nrow=N, ncol=topic)
  pr <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, k)))
  
  for(j in 1:k){
    #�f�[�^�^�C�v���Ƃ̐����ϐ��𐶐�
    d <- matrix(rgamma(N*topic, 5.5, 15.0), nrow=N, ncol=topic)
    u[, , j] <- d * matrix(rbinom(N*topic, 1, rep(runif(topic, 0.35, 0.7), rep(N, topic))), nrow=N, ncol=topic)
    
    #�������������ϐ��𐶐�
    x <- x + u[, , j] * pr[j]
  }
  
  ##��������ݒ�
  #�������̃C���f�b�N�X���쐬
  v_matrix <- matrix(0, nrow=k-1, ncol=k-1)
  v_matrix[upper.tri(v_matrix)] <- 1
  v_list <- list()
  for(j in 1:(k-1)){
    index <- which(v_matrix[j, ]==1)
    if(length(index) > 0){
      v_list[[j]] <- cbind(j1=j, j2=which(v_matrix[j, ]==1))
    }
  }
  v_index <- do.call(rbind, v_list)
  
  #�������̓��͕ϐ����쐬
  v <- u[, v_index[, 1]+1] * u[, v_index[, 2]+1]
  
  #�����ϐ��̊����C���f�b�N�X
  allocation_index11 <- allocation_index12 <- matrix(0, nrow=k-1, ncol=k-2)
  allocation_index21 <- allocation_index22 <- matrix(0, nrow=k-1, ncol=k-2)
  for(j in 1:(k-1)){
    index <- which(rowSums(v_index==j)==1)
    allocation_index11[j, ] <- allocation_index12[j, ] <- 
      rowSums(matrix(as.numeric(v_index[index, ]!=j), nrow=length(index)) * v_index[index, ])
    allocation_index21[j, ] <- allocation_index22[j, ] <-  index
  }
  allocation_index11[lower.tri(allocation_index11)] <- 0
  allocation_index21[lower.tri(allocation_index21)] <- 0
  vec <- rep(1, s)
  j_data12 <- matrix(1:(k-1), nrow=k-1, ncol=k-2) 
  j_data11 <- j_data12 * (allocation_index11 > 0)
  
  ##�p�����[�^�Ɖ����ϐ��𐶐�
  #�p�����[�^�𐶐�
  sigma <- sigmat <- 0.5
  beta <- betat <- c(5.0, runif(k-1, -1.25, 0.75))
  theta <- thetat <- mvrnorm(k-1, rep(0, s), diag(0.15, s))
  
  #�����ϐ��𐶐�
  u_mu <- u %*% beta
  v_mu <- rep(0, N)
  for(j in 1:(k-1)){
    v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% ((theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec)
  }
  mu <- u_mu + v_mu   #���Ғl���Z�o
  y0 <- rnorm(N, mu, sigma)   #���K���z���牞���ϐ��𐶐�
  
  if(sd(y0) > 1.75 & sd(y0) < 2.25 & min(y0) > -7.5 & max(y0) < 17.5 & mean(y0) > 4.5 & mean(y0) < 6.0){
    break
  }
}

#�����ϐ���1�`10�ɕϊ�����
y <- round(y0)
y[y > 10] <- 10; y[y < 1] <- 1
hist(y0, breaks=25, col="grey", main="�^�̃X�R�A���z", xlab="�X�R�A")
hist(y, breaks=25, col="grey", main="�ؒf���ꂽ�X�R�A���z", xlab="�X�R�A")


####�}���R�t�A�������e�J�����@��FM�𐄒�####
##�A���S���Y���̐ݒ�
R <- 2000
burnin <- 750
keep <- 2
disp <- 10
iter <- 0

##���O���z�̐ݒ�
alpha1 <- rep(0, k)
alpha2 <- rep(0, s)
tau1 <- 100 * diag(k)
tau2 <- 100 * diag(s)
inv_tau1 <- solve(tau1)
inv_tau2 <- solve(tau2)
s0 <- 1.0
v0 <- 1.0

##�^�l�̐ݒ�
beta <- betat
sigma <- sigmat
theta <- thetat

##�����l�̐ݒ�
beta <- as.numeric(solve(t(u) %*% u) %*% t(u) %*% y); beta[1] <- mean(y)
theta <- mvrnorm(k-1, rep(0, s), diag(0.1, s))
sigma <- 1.0


##�p�����[�^�̊i�[�p�z��
m <- 0
BETA <- matrix(0, nrow=R/keep, ncol=k)
SIGMA <- rep(0, R/keep)
THETA <- matrix(0, nrow=k-1, ncol=s)

##�C���f�b�N�X�ƃf�[�^�̒萔��ݒ�
#�C���f�b�N�X��ݒ�
index_list11 <- index_list12 <- list()
index_list21 <- index_list22 <- index_list23 <- list()

for(j in 1:(k-1)){
  #�f�[�^�𒊏o
  index <- (allocation_index11==j) + (j_data11==j)
  
  #����p�����[�^�̃C���f�b�N�X
  j_index <- v_index[rowSums(v_index * cbind(v_index[, 1]==j, v_index[, 2]==j)) > 0, ]
  index_list11[[j]] <- j_index[j_index!=j]
  index_list12[[j]] <- allocation_index21[index==1]
  
  #�Œ�p�����[�^�̃C���f�b�N�X
  index1 <- as.numeric(t(allocation_index11 * (1-index)))
  index2 <- as.numeric(t(allocation_index21 * (1-index)))
  index3 <- as.numeric(t(j_data11 * (1-index)))
  index_list21[[j]] <- index1[index1 > 0]
  index_list22[[j]] <- index2[index2 > 0]
  index_list23[[j]] <- index3[index3 > 0]
}

#�f�[�^�̐ݒ�
uu <- t(u) %*% u
inv_uu <- solve(uu + inv_tau1)
v_array <- array(0, dim=c(N, k-2, k-1))
for(j in 1:(k-1)){
  v_array[, , j] <- v[, index_list12[[j]]]
}
v_mu <- rep(0, N)
for(j in 1:(k-1)){
  v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
}

##�ΐ��ޓx�̊�l
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##��A�W�����T���v�����O
  #�����ϐ��̌덷���Z�o
  er <- as.numeric(y - v_mu)
  
  #���ϗʐ��K���z�̃p�����[�^��ݒ�
  beta_mu <- inv_uu %*% t(u) %*% er   #���ϗʐ��K���z�̕��σx�N�g��
  cov <- sigma^2 * inv_uu   #���ϗʐ��K���z�̕��U
  
  #���ϗʐ��K���z�����A�W�����T���v�����O
  beta <- mvrnorm(1, beta_mu, cov)
  u_mu <- u %*% beta
  
  
  ##���ݍ�p���̓����x�N�g�����T���v�����O
  for(j in 1:(k-1)){
    #�����ϐ���ݒ�
    er <- as.numeric(y - u_mu - v[, index_list22[[j]]] %*% (theta[index_list21[[j]], ] * theta[index_list23[[j]], ]) %*% vec)
    
    #�����x�N�g���̃p�����[�^��ݒ�
    x <- v_array[, , j] %*% theta[index_list11[[j]], ] 
    inv_xxv <- solve(t(x) %*% x + inv_tau2)
    theta_mu <- inv_xxv %*% t(x) %*% er   #���ϗʐ��K���z�̕��ς׃N�g��
    cov <- sigma^2 * inv_xxv   #���ϗʐ��K���z�̕��U
    
    #���ϗʐ��K���z�����A�W�����T���v�����O
    theta[j, ] <- mvrnorm(1, theta_mu, cov)
  }
  
  #���ݍ�p�̕��σx�N�g�����X�V
  v_mu <- rep(0, N)
  for(j in 1:(k-1)){
    v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
  }
  
  ##���f���̕W���΍����T���v�����O
  mu <- u_mu + v_mu
  er <- as.numeric(y - mu)
  
  #�t�K���}���z�̃p�����[�^
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- N + v0
  
  #�t�K���}���z����W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    SIGMA[mkeep] <- sigma
    
    if(rp >= burnin){
      m <- m + 1
      THETA <- THETA + theta
    }
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���Z�o
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLst))
    print(round(rbind(beta, betat), 3))
    print(c(sigma, sigmat))
  }
}

####���茋�ʂ̊m�F�ƓK���x####
##�T���v�����O���ʂ̉���
matplot(BETA, type="l", xlab="�T���v�����O��", main="beta�̃T���v�����O���ʂ̃v���b�g")
plot(1:length(SIGMA), SIGMA, type="l", xlab="�T���v�����O��", main="sigma�̃T���v�����O���ʂ̃v���b�g")

##�p�����[�^�̎��㕽��
#�o�[���C������
RS1 <- burnin / keep
RS2 <- R / keep

#���㕽�ς��v�Z  
beta <- colMeans(BETA[RS1:RS2, ])
theta <- THETA / m
sigma <- mean(SIGMA[RS1:RS2])

#�p�����[�^�̐^�l�Ɣ�r
round(rbind(beta=beta, betat), 3)
round(rbind(theta=theta, thetat), 3)
round(c(sigma, sigmat), 3)