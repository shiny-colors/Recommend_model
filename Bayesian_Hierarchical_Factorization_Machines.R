#####Bayesian Hierarchical Factorization Machines#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

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
s1 <- 10; vec1 <- rep(1, s1)   #�s�񕪉��̊�ꐔ
s2 <- 5; vec2 <- rep(1, s2)   #���ݍ�p�̊�ꐔ
hh <- 5000   #���[�U�[��
item <- 3000   #�A�C�e����
pt <- rtpois(hh, rgamma(hh, 35.0, 0.2), a=1, b=Inf)   #�w���ڐG��
f <- sum(pt)
vec_s1 <- rep(1, s1)
vec_s2 <- rep(1, s2)

#ID��ݒ�
user_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:f, user_id, rank)))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

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

#�X�p�[�X�s����쐬
user_data <- sparseMatrix(1:f, user_id, x=rep(1, f), dims=c(f, hh))
user_data_T <- t(user_data)
item_data <- sparseMatrix(1:f, item_id, x=rep(1, f), dims=c(f, item))
item_data_T <- t(item_data)

#���������f�[�^������
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="�A�C�e���̍w���p�x", main="�A�C�e���̍w���p�x���z")
gc(); gc()


##�����ϐ��̐���
#���f���̐����ϐ��𐶐�
k1 <- 3; k2 <- 4; k3 <- 4
k <- k1 + k2 + k3
x1 <- matrix(0, nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  x1[, j] <- rbeta(f, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������
z <- cbind(x1, x2 ,x3)

#���[�U�[�̐����ϐ��𐶐�
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(0, nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  u1[, j] <- rbeta(hh, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������


#�A�C�e���̐����ϐ��𐶐�
k1 <- 3; k2 <- 2; k3 <- 4
v1 <- matrix(0, nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  v1[, j] <- rbeta(item, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #�f�[�^������


##�������̐ݒ�
#�g�ݍ��킹���쐬
index_combine <- t(combn(c(1:ncol(z)), m=2))
combine_list <- list()
combine_n <- rep(0, max(index_combine[, 1]))
for(j in 1:max(index_combine[, 1])){
  combine_list[[j]] <- index_combine[which(index_combine[, 1]==j), 2]
  combine_n[j] <- length(combine_list[[j]])
}

#�������̔z���ݒ�
z_list <- list()
for(j in 1:length(combine_list)){
  z_list[[j]] <- z[, j] * z[, combine_list[[j]]]
}
zz <- do.call(cbind, z_list)


####�����ϐ��𐶐�####
rp <- 0
repeat { 
  print(rp <- rp + 1)
  
  ##�p�����[�^�Ɖ����ϐ��𐶐�
  #���f���̕W���΍�
  sigma <- sigmat <- 1.0
  
  #�K�w���f���̕��U�p�����[�^
  Cov_x <- Cov_xt <- runif(ncol(x), 0.05, 0.25) * diag(ncol(x))
  Cov_u <- Cov_ut <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_v <- Cov_vt <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_z <- Cov_zt <- runif(s2, 0.05, 0.15) * diag(s2)
  
  #�K�w���f���̉�A�W����ݒ�
  alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x))
  alpha_u <- matrix(0, nrow=ncol(u), ncol=s1)
  alpha_v <- matrix(0, nrow=ncol(v), ncol=s1)
  alpha_z <- array(0, dim=c(ncol(u), s2, k-1))
  
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.4, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.225), nrow=s2, ncol=k-1)
    } else {
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.35, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.225), nrow=s2, ncol=k-1)
    }
  }
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(s2, -0.5, 0.4)
    } else {
      alpha_v[j, ] <- runif(s2, -0.4, 0.4)
    }
  }
  alpha_xt <- alpha_x; alpha_ut <- alpha_u; alpha_vt <- alpha_v; alpha_zt <- alpha_z   #�^�l���i�[
  
  #���ϗʉ�A���f�����烆�[�U�[�ʂ̉�A�p�����[�^�𐶐�
  theta_x <- theta_xt <- u %*% alpha_x + mvrnorm(hh, rep(0, ncol(x)), Cov_x)   #�ϗʌ��ʂ̃p�����[�^
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, s1), Cov_u)   #���[�U�[�̍s�񕪉��̃p�����[�^
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, s1), Cov_v)   #�A�C�e���̍s�񕪉��̃p�����[�^
  
  theta_z <- array(0, c(hh, s2, k-1))
  for(j in 1:(k-1)){
    theta_z[, , j] <- u %*% alpha_z[, , j] + mvrnorm(hh, rep(0, s2), Cov_z)   #���ݍ�p�̃p�����[�^
  }
  theta_xt <- theta_x; theta_ut <- theta_u; theta_vt <- theta_v; theta_zt <- theta_z
  
  ##���K���z������p�ƍw���x�N�g���𐶐�
  #�ϗʌ��ʂ̃p�����[�^
  x_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  #�s�񕪉��̃p�����[�^
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

  #���ݍ�p�̃p�����[�^
  theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
  for(i in 1:length(combine_n)){
    vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
    for(j in 1:combine_n[[i]]){
      vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
    }
    z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
  }
  
  #���݌��p�Ɖ����ϐ��𐶐�
  mu <- mut <- x_mu + z_vec + uv   #���Ғl
  U <- mu + rnorm(f, 0, sigma)   #���݌��p�𐶐�
  
  #�w���x�N�g���ɕϊ�
  y <- ifelse(U > 0, 1, 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break����
}

#�������������ϐ����m�F
mean(y)   #�w���m��
prob <- pnorm(U, 0, sigma)   #�����m��
mean(prob[y==1]); mean(prob[y==0])   #�w���L���ʂ̉����m��
hist(U, col="grey", main="���݌��p�̕��z", xlab="���݌��p", breaks=25)


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�YFactorization Machines�𐄒�####
##�ؒf���K���z�̗����𔭐�������֐�
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 1000
keep <- 2  
iter <- 0
burnin <- 100
disp <- 10

##�C���f�b�N�X�ƃf�[�^�̒萔��ݒ�
#�s�񕪉��̃C���f�b�N�X���쐬
user_index <- item_index <- list()
ui_id <- iu_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
}

#���ݍ�p���̃C���f�b�N�X���쐬
index_allocation1 <- rep(1:(ncol(z)-1), combine_n)
index_allocation2 <- unlist(combine_list)
index_allocation <- matrix(1:ncol(z), nrow=ncol(z), ncol=ncol(z), byrow=T); diag(index_allocation) <- 0

index_z <- index_list1 <- index_list2 <- list()
for(j in 1:(k-1)){
  index_z[[j]] <- sort(c(which(index_allocation1==j), which(index_allocation2==j)))
  index_list1[[j]] <- cbind(index_allocation1[index_z[[j]]], index_allocation2[index_z[[j]]])
  index_list2[[j]] <- cbind(index_allocation1[-index_z[[j]]], index_allocation2[-index_z[[j]]])
}

#���ݍ�p���̃f�[�^�̐ݒ�
zz_array1 <- array(0, dim=c(f, k-2, k-1))
zz_array2 <- array(0, dim=c(f, ncol(zz)-(k-2), k-1))
for(j in 1:(k-1)){
  zz_array1[, , j] <- zz[, index_z[[j]]]
  zz_array2[, , j] <- zz[, -index_z[[j]]]
}

#���͕ϐ��̒萔��ݒ�
xx_list <- list()
for(i in 1:hh){
  xx_list[[i]] <- t(x[user_index[[i]], ]) %*% x[user_index[[i]], ]
}


##���O���z�̐ݒ�
#�ϗʌ��ʂ̊K�w���f���̎��O���z
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕���
ADelta1 <- 0.01 * diag(1, ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#���[�U�[�̍s�񕪉��̊K�w���f���̎��O���z
Deltabar2 <- matrix(0, nrow=ncol(u), ncol=s1)   #�K�w���f���̉�A�W���̎��O���z�̕���
ADelta2 <- 0.01 * diag(1, ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- s1 + 1   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(rep(1, s1)) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̍s�񕪉��̊K�w���f���̎��O���z
Deltabar3 <- matrix(0, nrow=ncol(v), ncol=s1)   #�K�w���f���̉�A�W���̎��O���z�̕���
ADelta3 <- 0.01 * diag(1, ncol(v))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu3 <- s1 + 1   #�t�E�B�V���[�g���z�̎��R�x
V3 <- nu3 * diag(rep(1, s1)) #�t�E�B�V���[�g���z�̃p�����[�^

#���ݍ�p���̊K�w���f���̎��O���z
Deltabar4 <- matrix(0, nrow=ncol(u), ncol=s2)   #�K�w���f���̉�A�W���̎��O���z�̕���
ADelta4 <- 0.01 * diag(1, ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu4 <- s2 + 1   #�t�E�B�V���[�g���z�̎��R�x
V4 <- nu4 * diag(rep(1, s2)) #�t�E�B�V���[�g���z�̃p�����[�^

#�ϗʌ��ʂ̊K�w���f���̎��O���z
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕���
ADelta1 <- 0.01 * diag(1, ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^


##�p�����[�^�̐^�l
#�K�w���f���̃p�����[�^�̐^�l
Cov_x <- Cov_xt; Cov_u <- Cov_ut; Cov_v <- Cov_vt; Cov_z <- array(Cov_zt, dim=c(s2, s2, k-1))
inv_Cov_x <- solve(Cov_x); inv_Cov_u <- solve(Cov_u); inv_Cov_v <- solve(Cov_v)
inv_Cov_z <- array(0, dim=c(s2, s2, k-1))
for(j in 1:(k-1)){
  inv_Cov_z[, , j] <- solve(Cov_z[, , j])
}
alpha_x <- alpha_xt; alpha_u <- alpha_ut; alpha_v <- alpha_vt; alpha_z <- alpha_zt
x_mu <- u %*% alpha_x; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
z_mu <- array(0, c(hh, s2, k-1))
for(j in 1:(k-1)){
  z_mu[, , j] <- u %*% alpha_z[, , j]
}

#���f���p�����[�^�̐^�l
sigma <- sigmat
theta_x <- theta_xt; theta_u <- theta_ut; theta_v <- theta_vt; theta_z <- theta_zt
user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))

#�s�񕪉��̃p�����[�^
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

#���ݍ�p���̃p�����[�^
theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
for(i in 1:length(combine_n)){
  vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
  for(j in 1:combine_n[[i]]){
    vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
  }
  z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
}


##�p�����[�^�̏����l
#�K�w���f���̏����l
Cov_x <- diag(0.01, ncol(x)); Cov_u <- Cov_v <- diag(0.01, s1); Cov_z <- array(diag(0.01, s2), dim=c(s2, s2, k-1))
alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x)); x_mu <- u %*% alpha_x
alpha_u <- matrix(0, nrow=ncol(u), ncol=s1); u_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=s1); v_mu <- v %*% alpha_v
alpha_z <- array(0, dim=c(ncol(u), s2, k-1)); z_mu <- array(0, dim=c(hh, s2, k-1))

#���f���p�����[�^�̏����l
sigma <- 1
theta_x <- mvrnorm(hh, as.numeric(solve(t(x) %*% x) %*% t(x) %*% y), Cov_x) 
theta_u <- mvrnorm(hh, rep(0, s1), Cov_u)
theta_v <- mvrnorm(item, rep(0, s1), Cov_v)
theta_z <- array(0, dim=c(hh, s2, k-1))
for(j in 1:(k-1)){
  theta_z[, , j] <- mvrnorm(hh, rep(0, s2), Cov_z[, , j])
}

#�s�񕪉��̃p�����[�^
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

#���ݍ�p���̃p�����[�^
theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
for(i in 1:length(combine_n)){
  vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
  for(j in 1:combine_n[[i]]){
    vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
  }
  z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
}


##�p�����[�^�̊i�[�p�z��
#���f���p�����[�^�̊i�[�p�z��
d <- 0
THETA_X <- array(0, dim=c(hh, k, R/keep))
THETA_U <- array(0, dim=c(hh, s1, R/keep))
THETA_V <- array(0, dim=c(item, s1, R/keep))
THETA_Z <- array(0, dim=c(hh, s2, k-1))

#�K�w���f���̊i�[�p�z��
ALPHA_X <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), s1, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), s1, R/keep))
ALPHA_Z <- array(0, dim=c(ncol(u), s2, k-1, R/keep))
COV_X <- array(0, dim=c(k, k, R/keep))
COV_U <- COV_V <- array(0, dim=c(s1, s1, R/keep))
COV_Z <- array(0, dim=c(s2, s2, k-1, R/keep))


##�ؒf�̈���`
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##�ΐ��ޓx�̊�l
#1�p�����[�^���f���̑ΐ��ޓx
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #�ΐ��ޓx

#�x�X�g���f���̑ΐ��ޓx
prob <- pnorm(mut, 0, sigmat)   #�w���m��
prob[prob==1] <- 0.9999999; prob[prob==0] <- 0.0000001
LLbest <- sum(y*log(prob) + (1-y)*log(1-prob))   #�ΐ��ޓx


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
    
  ##�ؒf���K���z�����݌��p�𐶐�
  mu <- user_mu + uv + z_vec   #���݌��p�̊��Ғl
  U <- extraDistr::rtnorm(f, mu, sigma, a, b)   #���݌��p�𐶐�
  
  ##���[�U�[�̉�A�x�N�g�����T���v�����O
  #���f���̉����ϐ�
  y_er <- U - uv - z_vec
  
  for(i in 1:hh){
    #��A�x�N�g���̎��㕪�z�̃p�����[�^
    XX <- xx_list[[i]]
    Xy <- t(x[user_index[[i]], ]) %*% y_er[user_index[[i]]]
    inv_XXV <- solve(XX + inv_Cov_x)
    mu <- inv_XXV %*% (Xy + inv_Cov_x %*% x_mu[i, ])   #���㕪�z�̕���
    
    #���ϗʐ��K���z�����A�x�N�g�����T���v�����O
    theta_x[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  
  ##���[�U�[�̓����s����T���v�����O
  #���f���̉����ϐ�
  u_er <- U - user_mu - z_vec
  
  for(i in 1:hh){
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_v[ui_id[[i]], ]
    Xy <- t(X) %*% u_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u)
    mu <- inv_XXV %*% (Xy + inv_Cov_u %*% u_mu[i, ])
    
    #���ϗʐ��K���z��������x�N�g�����T���v�����O
    theta_u[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  
  ##�A�C�e���̓����s����T���v�����O
  for(j in 1:item){
    
    #�����x�N�g���̎��㕪�z�̃p�����[�^
    X <- theta_u[iu_id[[j]], ]
    Xy <- t(X) %*% u_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v)
    mu <- inv_XXV %*% (Xy + inv_Cov_v %*% v_mu[j, ])
    
    #���ϗʐ��K���z��������x�N�g�����T���v�����O
    theta_v[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #�s�񕪉��̃p�����[�^���X�V
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
  
  
  ##���ݍ�p���̓����x�N�g�����T���v�����O
  #���f���̉����ϐ�
  z_er <- U - user_mu - uv

  for(i in 1:hh){
    #�f�[�^�̒��o
    zz1 <- zz_array1[user_index[[i]], , ]; zz2 <- zz_array2[user_index[[i]], , ]
    theta <- t(theta_z[i, , ])
    er_vec <- z_er[user_index[[i]]]
  
    for(j in 1:(k-1))
      #�����ϐ��̐ݒ�
      er <- er_vec - as.numeric(zz2[, , j] %*% (theta[index_list2[[j]][, 1], ] * theta[index_list2[[j]][, 2], ]) %*% vec_s2)
      
      #���ݍ�p���̎��㕪�z�̃p�����[�^
      X <- zz1[, , j] %*% theta[index_list1[[j]][, 2], ]
      Xy <- t(X) %*% er
      inv_XXV <- solve(t(X) %*% X + inv_Cov_z[, , j])
      mu <- as.numeric(inv_XXV %*% (Xy + inv_Cov_z[, , j] %*% z_mu[i, , j]))   #���㕪�z�̕���
      
      #���ϗʐ��K���z������ݍ�p�����T���v�����O
      theta_z[i, , j] <- mvrnorm(1, mu, sigma^2*inv_XXV)
      theta[j, ] <- theta_z[i, , j]
  }
  
  #���ݍ�p���̃p�����[�^���X�V
  theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
  for(i in 1:length(combine_n)){
    vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
    for(j in 1:combine_n[[i]]){
      vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
    }
    z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
  }
  
  
  ##���[�U�[�̉�A�x�N�g���̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������p�����[�^���T���v�����O
  out <- rmultireg(theta_x, u, Deltabar1, ADelta1, nu1, V1)
  alpha_x <- out$B; x_mu <- u %*% alpha_x   
  Cov_x <- out$Sigma; inv_Cov_x <- solve(Cov_x)
  
  ##���[�U�[�����s��̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������p�����[�^���T���v�����O
  out <- rmultireg(theta_u, u, Deltabar2, ADelta2, nu2, V2)
  alpha_u <- out$B; u_mu <- u %*% alpha_u   
  Cov_u <- out$Sigma; inv_Cov_u <- solve(Cov_u)
  
  ##�A�C�e���̓����s��̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������p�����[�^���T���v�����O
  out <- rmultireg(theta_v, v, Deltabar3, ADelta3, nu3, V3)
  alpha_v <- out$B; v_mu <- v %*% alpha_v   
  Cov_v <- out$Sigma; inv_Cov_v <- solve(Cov_v)
  
  ##���ݍ�p���̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f���p�����[�^���T���v�����O
  for(j in 1:(k-1)){
    out <- rmultireg(theta_z[, , j], u, Deltabar4, ADelta4, nu4, V4)
    alpha_z[, , j] <- out$B; z_mu[, , j] <- u %*% alpha_z[, , j]
    Cov_z[, , j] <- out$Sigma; inv_Cov_z[, , j] <- solve(Cov_z[, , j])
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̊i�[
  if(rp%%keep==0){
    mkeep <- rp/keep

    #���f���p�����[�^�̊i�[
    THETA_X[, , mkeep] <- theta_x
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    
    #�K�w���f���̊i�[�p�z��
    ALPHA_X[, , mkeep] <- alpha_x
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    ALPHA_Z[, , , mkeep] <- alpha_z
    COV_X[, , mkeep] <- Cov_x
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    COV_Z[, , , mkeep] <- Cov_z
    
    if(rp >= burnin){
      d <- d + 1
      THETA_Z <- THETA_Z + theta_z
    }
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    mu <- user_mu + uv + z_vec   #���݌��p�̊��Ғl
    prob <- pnorm(mu, 0, sigma)   #�w���m��
    prob[prob==1] <- 0.9999999; prob[prob==0] <- 0.0000001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))   #�ΐ��ޓx
    
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(c(LL, LLbest, LLst))
  }
}

