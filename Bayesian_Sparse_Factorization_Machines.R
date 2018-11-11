#####Bayesian Sparse Factorization Machines#####
library(MASS)
library(Matrix)
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
item <- 3000   #�A�C�e����
tag <- 150   #�^�O��
pt <- rtpois(item, rgamma(item, 30.0, 0.225), a=1, b=Inf)   #�w���ڐG��
N <- sum(pt)
n <- rtpois(item, 1.0, a=0, b=5)   #�^�O��
k <- 7   #��ꐔ
vec_k <- rep(1, k)

#ID��ݒ�
item_id <- rep(1:item, pt)
pt_id <- as.numeric(unlist(tapply(1:N, item_id, rank)))
ID <- data.frame(no=1:N, id=item_id, t=pt_id)   #�f�[�^�̌���
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}

##�K�w���f���̐����ϐ���ݒ�
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

##�^�O�𐶐�
#�p�����[�^�̐ݒ�
topic <- 25
omega <- extraDistr::rdirichlet(topic, rep(0.5, tag))
z <- as.numeric(rmnom(item, 1,  extraDistr::rdirichlet(item, rep(1.0, topic))) %*% 1:topic)

#�������z����^�O�𐶐�
max_n <- max(n)
tag_list <- list()
tag_id <- matrix(0, nrow=N, ncol=max_n)
tag_data <- matrix(0, nrow=N, ncol=tag); storage.mode(tag_data) <- "integer"

for(j in 1:item){
  repeat { 
    x <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(x)==1){
      tag_list[[j]] <- (x * 1:tag)[x > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
  tag_data[item_list[[j]], ] <- matrix(x, nrow=length(item_list[[j]]), ncol=tag, byrow=T)
}
tag_id0 <- tag_id
tag_id0[tag_id0==0] <- tag+1

#�g�ݍ��킹���쐬
index_combine <- t(combn(c(1, 2, 3, 4, 5), m=2))
combine_list <- list()
combine_n <- rep(0, max(index_combine[, 1]))
for(j in 1:max(index_combine[, 1])){
  combine_list[[j]] <- index_combine[which(index_combine[, 1]==j), 2]
  combine_n[j] <- length(combine_list[[j]])
}

#�������̃C���f�b�N�X���쐬
index_n <- which(rowSums(tag_id > 0) >= 2)
tag_n <- length(index_n)
tag_dt <- sparseMatrix(rep(1:tag_n, max_n), 1:(tag_n*max_n), x=rep(1, tag_n*max_n), dims=c(tag_n, tag_n*max_n))


####�����ϐ��𐶐�####
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  ##���f���̃p�����[�^
  beta <- betat <- 5.5
  sigma <- sigmat <- 0.75
  
  ##�A�C�e���̊K�w���f���̃p�����[�^
  #�W���΍���ݒ�
  tau_v <- tau_vt <- runif(1, 0.3, 0.75)
  
  #��A�W����ݒ�
  alpha_v <- rep(0, ncol(v))
  for(j in 1:ncol(v)){
    alpha_v[j] <- runif(1, -1.25, 1.25)
  }
  alpha_vt <- alpha_v
  
  #��A���f������A�C�e���ʂ̕ϗʌ��ʂ𐶐�
  theta_v <- theta_vt <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
  theta_vec1 <- theta_v[item_id]
  
  ##�^�O�̃p�����[�^�𐶐�
  #�^�O�ʂ̕ϗʌ��ʂ𐶐�
  tau_r <- tau_rt <- runif(1, 0.25, 0.5)
  theta_r <- theta_rt <- rnorm(tag, 0, tau_r)
  theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
  
  ##���ݍ�p�̃p�����[�^�𐶐�
  #�����x�N�g���̃p�����[�^�𐶐�
  tau_g <- tau_gt <- runif(k, 0.1, 0.4) * diag(k)
  theta_g <- theta_gt <- mvrnorm(tag, rep(0, k), tau_g)
  theta_g0 <- rbind(theta_g, 0)

  #���ݍ�p�̃x�N�g���𐶐�
  WH <- rep(0, N)
  for(j in 1:length(combine_n)){
    W <- theta_g0[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
  }
  
  ##���K���z����]���x�N�g���𐶐�
  mu <- beta + theta_vec1 + theta_vec2 + WH   #���Ғl��ݒ� 
  y0 <- rnorm(N, mu, sigma)   #�����ϐ��𐶐�

  #break����
  if(min(y0) > -3.0 & max(y0) < 14.0 & mean(y0) > 4.5 & mean(y0) < 6.0){
    break
  }
}

#�����ϐ���1�`10�ɕϊ�����
y <- round(y0)
y[y > 10] <- 10; y[y < 1] <- 1
hist(y0, breaks=25, col="grey", main="�^�̃X�R�A���z", xlab="�X�R�A")
hist(y, breaks=25, col="grey", main="�ؒf���ꂽ�X�R�A���z", xlab="�X�R�A")


####�}���R�t�A�������e�J�����@��Bayesian SFM�𐄒�####
##�A���S���Y���̐ݒ�
R <- 2000
burnin <- 500
keep <- 2
disp <- 10
iter <- 0

#�f�[�^�ƃC���f�b�N�X��ݒ�
#�C���f�b�N�X��ݒ�
tag_list1 <- tag_list2 <- dt_list <- list()
tag_n1 <- tag_n2 <- rep(0, tag)
for(i in 1:tag){
  index1 <- index2 <- c()
  for(j in 1:max_n){
    index1 <- c(index1, which(tag_id[, j]==i))
    index2 <- c(index2, index_n[which(tag_id[index_n, j]==i)])
  }
  tag_list1[[i]] <- sort(index1)
  tag_list2[[i]] <- sort(index2)
  tag_n1[i] <- length(tag_list1[[i]])
  tag_n2[i] <- length(tag_list2[[i]])
  dt_list[[i]] <- sparseMatrix(rep(1:tag_n2[i], max_n), 1:(tag_n2[i]*max_n), x=rep(1, tag_n2[i]*max_n), 
                               dims=c(tag_n2[i], tag_n2[i]*max_n))
}


##���O���z�̐ݒ�
#��A�p�����[�^�̎��O���z
alpha1 <- 0
alpha2 <- rep(0, ncol(v))
alpha3 <- rep(0, k)

#���U�̎��O���z
s0 <- 0.1
v0 <- 0.1   
nu <- 1   #�t�E�B�V���[�g���z�̎��R�x
V <- 0.01 * diag(k)    #�t�E�B�V���[�g���z�̃p�����[�^
tau1 <- 100
tau2 <- 100 * diag(ncol(v))
inv_tau2 <- solve(tau2)
tau3 <- 100 * diag(k)
inv_tau3 <- solve(tau3)

##�^�l�̐ݒ�
#���f���p�����[�^
beta <- betat
sigma <- sigmat

#�K�w���f���̃p�����[�^
alpha_v <- alpha_vt
tau_v <- tau_vt
tau_r <- tau_rt
tau_g <- tau_gt
inv_tau_g <- solve(tau_g)

#�ϗʌ��ʂ̃p�����[�^
theta_v <- theta_vt
theta_r <- theta_rt
theta_g <- theta_gt
theta_g0 <- rbind(theta_g, 0)

#�]���x�N�g���̊��Ғl
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
y_mu <- beta + theta_vec1 + theta_vec2 + WH   #���Ғl��ݒ� 

##�����l�̐ݒ�
#���f���p�����[�^
beta <- mean(y)
sigma <- 1.0

#�K�w���f���̃p�����[�^
alpha_v <- rep(0, ncol(v))
tau_v <- 0.1
tau_r <- 0.1
tau_g <- 0.2 * diag(k)

#�ϗʌ��ʂ̃p�����[�^
theta_v <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
theta_r <- rnorm(tag, 0, tau_r)
theta_g <- mvrnorm(tag, rep(0, k), tau_g)
theta_g0 <- rbind(theta_g, 0)

#�]���x�N�g���̊��Ғl
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
mu <- beta + theta_vec1 + theta_vec2 + WH   #���Ғl��ݒ� 


##�p�����[�^�̊i�[�p�z��
#���f���p�����[�^
BETA <- rep(0, R/keep)
SIGMA <- rep(0, R/keep)

#�K�w���f���̃p�����[�^
ALPHA_V <- matrix(0, nrow=R/keep, ncol=ncol(v))
TAU_V <- rep(0, R/keep)
TAU_R <- rep(0, R/keep)
TAU_G <- array(0, dim=c(k, k, R/keep))

#�ϗʌ��ʂ̃p�����[�^
THETA_V <- matrix(0, nrow=R/keep, ncol=item)
THETA_R <- matrix(0, nrow=R/keep, ncol=tag)
THETA_G <- array(0, dim=c(tag, k, R/keep))


##�ΐ��ޓx�̊�l
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))
LLbest <- sum(dnorm(y, y_mu, sigmat, log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���f�����ς��T���v�����O
  #���f���덷��ݒ�
  er <- y - theta_vec1 - theta_vec2 - WH 
  
  #���K���z�̃p�����[�^
  omega <- (N/sigma^2) / (1/tau1 + N/sigma^2)   #�d�݌W��
  beta_mu <- omega * mean(er)   #���K���z�̕���
  cov <- 1 / (1/tau1 + N/sigma^2)   #���K���z�̕W���΍�
  
  #���K���z����p�����[�^���T���v�����O
  beta <- rnorm(1, beta_mu, cov)   
  
  ##���f���̕W���΍����T���v�����O
  #���f���덷��ݒ�
  er <- y - beta - theta_vec1 - theta_vec2 - WH
  
  #�t�K���}���z�̃p�����[�^
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- N + v0
  
  #�t�K���}���z����W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##�A�C�e���̕ϗʌ��ʂ��T���v�����O
  #���f���덷��ݒ�
  er <- y - beta - theta_vec2 - WH
  alpha_mu <- as.numeric(v %*% alpha_v)
  
  #���K���z�̃p�����[�^
  omega <- (pt/sigma^2) / (1/tau_v^2 + pt/sigma^2)   #�d�݌W��
  cov <- 1 / (1/tau_v^2 + pt/sigma^2)   #���K���z�̕W���΍�
  theta_mu <- rep(0, item)
  for(j in 1:item){
    theta_mu[j] <- (1-omega[j])*alpha_mu[j] + omega[j]*mean(er[item_list[[j]]])   #���K���z�̕���
  }
  
  #���K���z����p�����[�^���T���v�����O
  theta_v <- rnorm(item, theta_mu, cov)
  theta_vec1 <- theta_v[item_id]

  
  ##�A�C�e���^�O�̕ϗʌ��ʂ��T���v�����O
  #���f���덷��ݒ�
  er <- y - beta - theta_vec1 - WH

  #���K���z�̃p�����[�^
  omega <- (tag_n1/sigma^2) / (1/tau_r^2 + tag_n1/sigma^2)   #�d�݌W��
  cov <- 1 / (1/tau_r^2 + tag_n1/sigma^2)   #���K���z�̕W���΍�
  theta_mu <- rep(0, tag)
  for(j in 1:tag){
    #�T���v�����O�ΏۊO�̃A�C�e���^�O�̕ϗʌ���
    theta_r0 <- c(theta_r, 0); theta_r0[j] <- 0
    r <- as.numeric(matrix(theta_r0[tag_id0[tag_list1[[j]], ]], nrow=tag_n1[j], ncol=max_n) %*% rep(1, max_n))
    
    #���K���z�̕���
    er_y <- er[tag_list1[[j]]] - r
    theta_mu[j] <- omega[j]*mean(er_y)   
  }

  #���K���z����p�����[�^���T���v�����O
  theta_r <- rnorm(tag, theta_mu, cov)
  theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
  
  
  ##���ݍ�p�̓����x�N�g�����T���v�����O
  #���f���덷��ݒ�
  er <- y - beta - theta_vec1 - theta_vec2
  for(i in 1:tag){
    #�f�[�^�̐ݒ�
    theta_g0 <- rbind(theta_g, 0); theta_g0[i, ] <- 0
    index <- tag_list2[[i]]
    tag_vec <- tag_id0[index, ]
    
    #�T���v�����O�ΏۊO�̌��ݍ�p�x�N�g���̃p�����[�^
    wh <- rep(0, length(index))
    for(j in 1:length(combine_n)){
      w <- theta_g0[tag_vec[, j], ]
      h <- as.matrix(dt_list[[i]][, 1:(length(index)*combine_n[j])] %*% theta_g0[tag_vec[, combine_list[[j]]], ])
      wh <- wh + as.numeric((w * h) %*% vec_k)
    }
    
    #���ϗʐ��K���z�̃p�����[�^
    er_y <- er[index] - wh
    x <- as.matrix(dt_list[[i]] %*% theta_g0[tag_vec, ])
    inv_xxv <- solve(t(x) %*% x + inv_tau_g)
    theta_mu <- as.numeric(inv_xxv %*% t(x) %*% er_y)   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- sigma^2 * inv_xxv   #���ϗʐ��K���z�̕��U
  
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_g[i, ] <- mvrnorm(1, theta_mu, cov)
  }
    
  #���ݍ�p�x�N�g���̊��Ғl
  theta_g0 <- rbind(theta_g, 0)
  WH <- rep(0, N)
  for(j in 1:length(combine_n)){
    W <- theta_g0[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
  }
  
  ##�A�C�e���̕ϗʌ��ʂ̊K�w���f�����T���v�����O
  #���ϗʉ�A���f������p�����[�^���T���v�����O
  inv_xxv <- solve(t(v) %*% v + inv_tau2) 
  mu <- as.numeric(inv_xxv %*% t(v) %*% theta_v)
  alpha_v <- mvrnorm(1, mu, tau_v*inv_xxv)   #���ϗʐ��K���z����p�����[�^���T���v�����O
  alpha_mu <- as.numeric(v %*% alpha_v)
  
  #�t�K���}���z���番�U���T���v�����O
  er <- theta_v - alpha_mu
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- item + v0
  tau_v <- sqrt(1/rgamma(1, v1/2, s1/2))
  inv_tau_v <- solve(tau_v)
  
  ##�^�O�̕ϗʌ��ʂ̊K�w���f�����T���v�����O
  #�t�K���}���z���番�U���T���v�����O
  s1 <- as.numeric(t(theta_r) %*% theta_r) + s0
  v1 <- tag + v0
  tau_r <- sqrt(1/rgamma(1, v1/2, s1/2))
  inv_tau_r <- solve(tau_r)
  
  ##�^�O�̓����x�N�g���̊K�w���f�����T���v�����O
  #�t�E�B�V���[�g���z���番�U���T���v�����O
  IW_R <- t(theta_g) %*% theta_g + V
  Sn <- tag + nu
  tau_g <- diag(diag(rwishart(Sn, solve(IW_R))$IW))
  inv_tau_g <- solve(tau_g)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�p�����[�^���i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    
    #���f���p�����[�^���i�[
    BETA[mkeep] <- beta
    SIGMA[mkeep] <- sigma
    
    #�K�w���f���̃p�����[�^���i�[
    ALPHA_V[mkeep, ] <- alpha_v
    TAU_V[mkeep] <- tau_v
    TAU_R[mkeep] <- tau_r
    TAU_G[, , mkeep] <- tau_g
    
    #�ϗʌ��ʂ̃p�����[�^���i�[
    THETA_V[mkeep, ] <- theta_v
    THETA_R[mkeep, ] <- theta_r
    THETA_G[, , mkeep] <- theta_g
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���Z�o
    mu <- beta + theta_vec1 + theta_vec2 + WH   #���Ғl��ݒ� 
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(c(beta, betat), 3))
    print(c(sigma, sigmat))
  }
}


####���茋�ʂ̊m�F�ƓK���x####
##�T���v�����O���ʂ̉���
plot(1:(R/keep), BETA, type="l", xlab="�T���v�����O��", main="beta�̃T���v�����O���ʂ̃v���b�g")
matplot(THETA_V[, 1:10], type="l")
matplot(t(THETA_G[100, , ]), type="l")
t(THETA_G[100, , ])

plot(1:length(SIGMA), SIGMA, type="l", xlab="�T���v�����O��", main="sigma�̃T���v�����O���ʂ̃v���b�g")
THETA_V[, 1]

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

