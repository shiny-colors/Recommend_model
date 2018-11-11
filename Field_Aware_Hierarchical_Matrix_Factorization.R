#####Field Aware Hierarchical Matrix Factorization#####
options(warn=0)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(stringr)
library(extraDistr)
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
k <- 10   #��ꐔ
hh <- 10000   #���[�U�[��
item <- 3000   #�A�C�e����
company <- 1000   #���[�J�[��
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #�w���ڐG��
hhpt <- sum(pt)
vec <- rep(1, k)

#ID��ݒ�
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #�f�[�^�̌���
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}


##�f���x�N�g���𐶐�
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #�f�[�^������


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

#���[�J�[�̐����ϐ�
k1 <- 2; k2 <- 2; k3 <- 4
d1 <- matrix(runif(company*k1, 0, 1), nrow=company, ncol=k1)
d2 <- matrix(0, nrow=company, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.2, 0.6)
  d2[, j] <- rbinom(company, 1, pr)
}
d3 <- rmnom(company, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(d3))]
d0 <- cbind(1, d1, d2, d3)   #�f�[�^������


##�A�C�e���ƃ��[�J�[�̊����𐶐�
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

#�A�C�e���������烁�[�J�[�𐶐�
z <- rmnom(item, 1, as.numeric(extraDistr::rdirichlet(1, rep(2.5, company))))
index_z <- as.numeric(z[, colSums(z)!=0] %*% 1:sum(colSums(z)!=0))  
company_id <- left_join(data.frame(id=item_id), data.frame(id=1:item, no=index_z), by="id")$no   #id��ݒ�
company <- length(unique(company_id)); d <- d0[1:company, ]
company_list <- list()
for(j in 1:company){
  company_list[[j]] <- which(company_id==j)
}

#�ʂɘa����邽�߂̃X�p�[�X�s��
user_vec <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_vec <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))
company_vec <- sparseMatrix(sort(company_id), unlist(company_list), x=rep(1, hhpt), dims=c(company, hhpt))


#���������f�[�^������
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
freq_company <- plyr::count(item_id); freq_company$x <- as.character(freq_company$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="�A�C�e���̍w���p�x", main="�A�C�e���̍w���p�x���z")
hist(freq_company$freq, breaks=25, col="grey", xlab="��ƒP�ʂ̍w���p�x", main="��ƒP�ʂ̍w���p�x���z")
gc(); gc()


####�����ϐ��𐶐�####
for(rp in 1:1000){
  print(rp)
  
  ##�f���x�N�g���̃p�����[�^
  beta <- betat <- c(-0.6, runif(ncol(x)-1, -1.0, 1.0))
  
  ##�K�w���f���̃p�����[�^�𐶐�
  #�K�w���f���̕��U�p�����[�^
  Cov_u1 <- Cov_ut1 <- diag(runif(k, 0.01, 0.25), k)   #���[�U�[-�A�C�e���̊K�w���f���̕��U
  Cov_u2 <- Cov_ut2 <- diag(runif(k, 0.01, 0.25), k)   #���[�U�[-���[�J�[�̊K�w���f���̕��U
  Cov_v <- Cov_vt <- diag(runif(k, 0.01, 0.25), k)   #�A�C�e���̊K�w���f���̕��U
  Cov_d <- Cov_dt <- diag(runif(k, 0.01, 0.25), k)   #���[�J�[�̊K�w���f���̕��U
  
  #�K�w���f���̉�A�W����ݒ�
  alpha_u1 <- alpha_ut1 <- matrix(rnorm(k*ncol(u), 0, 0.5), nrow=ncol(u), ncol=k)
  alpha_u2 <- alpha_ut2 <- matrix(rnorm(k*ncol(u), 0, 0.5), nrow=ncol(u), ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.5), nrow=ncol(v), ncol=k)
  alpha_d <- alpha_dt <- matrix(rnorm(k*ncol(d), 0, 0.5), nrow=ncol(d), ncol=k)
  
  ##�s�񕪉��̃p�����[�^�𐶐�
  theta_u1 <- theta_ut1 <- u %*% alpha_u1 + mvrnorm(hh, rep(0, k), Cov_u1)
  theta_u2 <- theta_ut2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  theta_d <- theta_dt <- d %*% alpha_d + mvrnorm(company, rep(0, k), Cov_d)
  
  
  ##���W�b�g���f������w���x�N�g���𐶐�
  #���W�b�g�ƍw���m����ݒ�
  uv <- uvt <- rowSums(theta_u1[user_id, ] * theta_v[item_id, ])
  ud <- udt <- rowSums(theta_u2[user_id, ] * theta_d[company_id, ])
  mu <- as.numeric(x %*% beta)
  logit <- mu + uv + ud
  prob <- exp(logit) / (1 + exp(logit))
  
  #�x���k�[�C���z����w���x�N�g���𐶐�
  y <- rbinom(hhpt, 1, prob)
  print(mean(y))
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break����
}

#�w�������m�F
sum(y); mean(y)
mean(prob[y==0]); mean(prob[y==1])
hist(prob, main="�w���m���̐^�l�̕��z", xlab="�w���m��", col="grey", breaks=25)

  
  
####�}���R�t�A�������e�J�����@��FFMF���f���𐄒�####
##�ΐ����㕪�z���v�Z����֐�
#�ΐ��ޓx���v�Z����֐�
LLho <- function(mu, uv, ud, y, dt){
  
  #���W�b�g���f���̑ΐ��ޓx
  logit_exp <- exp(mu + uv + ud)   #���W�b�g�̊��Ғl�̎w��
  prob <- logit_exp / (1 + logit_exp)
  LLi_logit <- y*log(prob) + (1-y)*log(1-prob)

  #���[�U�[���Ƃ̑ΐ��ޓx
  LLi <- as.numeric(dt %*% LLi_logit)
  return(LLi)
}

#�ΐ��ޓx���v�Z����֐�
loglike <- function(mu, uv, ud, y){
  
  #���W�b�g���f���̑ΐ��ޓx
  logit_exp <- exp(mu + uv + ud)   #���W�b�g�̊��Ғl�̎w��
  prob <- logit_exp / (1 + logit_exp)
  LLi_logit <- y*log(prob) + (1-y)*log(1-prob)
  LL <- sum(LLi_logit)   #�ΐ��ޓx�̘a
  return(LL)
}

#���ϗʐ��K���z�̑ΐ����O���z
dmv <- function(er, inv_Cov, k){
  Li <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k))
  return(Li)
}

##�f���x�N�g���̃T���v�����O�ɕK�v�Ȋ֐�
#�f���x�N�g���̑ΐ����㕪�z�̔����֐�
dloglike <- function(beta, uv, ud, y, x, inv_Cov){
  #�����m���̐ݒ�
  logit_exp <- exp(as.numeric(x %*% beta) + uv + ud)   #���W�b�g�̎w���֐�
  prob <- logit_exp / (1 + logit_exp)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  dlogit <- y*x - x*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% beta)   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�

  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -as.numeric((colSums(dlogit) + dmvn))
  return(LLd)
}

#�f���x�N�g���̃��[�v�t���b�O�@�������֐�
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, uv, ud, y, x, inv_tau) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, uv, ud, y, x, inv_tau) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##���[�U�[�̃A�C�e���ɑ΂���s�񕪉��̃p�����[�^�̃T���v�����O�ɕK�v�Ȋ֐�
#�s�񕪉��̃p�����[�^�̑ΐ����㕪�z�̔����֐�
dloglike_u1 <- function(theta_u1, mu, theta_vec, ud, y, mu_u1, inv_Cov, user_vec, user_id, item_id){
  
  #�����m���̐ݒ�
  uv <- as.numeric((theta_u1[user_id, ] * theta_vec) %*% vec)
  logit_exp <- exp(mu + uv + ud)   #���W�b�g�̎w���֐�
  prob <- logit_exp / (1 + logit_exp)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  er <- theta_u1 - mu_u1
  dlogit <- y*theta_vec - theta_vec*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -as.matrix(user_vec %*% dlogit + dmvn)
  return(LLd)
}


#�s�񕪉��̃p�����[�^�̃��[�v�t���b�O�@�������֐�
leapfrog_u1 <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, mu, theta_vec, ud, y, mu_u1, inv_Cov_u1, user_vec, user_id, item_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, mu, theta_vec, ud, y, mu_u1, inv_Cov_u1, user_vec, user_id, item_id) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##���[�U�[�̃��[�J�[�ɑ΂���s�񕪉��̃p�����[�^�̃T���v�����O�ɕK�v�Ȋ֐�
#�s�񕪉��̃p�����[�^�̑ΐ����㕪�z�̔����֐�
dloglike_u2 <- function(theta_u2, mu, theta_vec, uv, y, mu_u2, inv_Cov, user_vec, user_id, company_id){
  
  #�����m���̐ݒ�
  ud <- as.numeric((theta_u2[user_id, ] * theta_vec) %*% vec)
  logit_exp <- exp(mu + ud + uv)   #���W�b�g�̎w���֐�
  prob <- logit_exp / (1 + logit_exp)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  er <- theta_u2 - mu_u2
  dlogit <- y*theta_vec - theta_vec*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -as.matrix(user_vec %*% dlogit + dmvn)
  return(LLd)
}

#�s�񕪉��̃p�����[�^�̃��[�v�t���b�O�@�������֐�
leapfrog_u2 <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, mu, theta_vec, uv, y, mu_u2, inv_Cov_u2, user_vec, user_id, company_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, mu, theta_vec, uv, y, mu_u2, inv_Cov_u2, user_vec, user_id, company_id) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##�A�C�e���̍s�񕪉��̃p�����[�^�̃T���v�����O�ɕK�v�Ȋ֐�
#�s�񕪉��̃p�����[�^�̑ΐ����㕪�z�̔����֐�
dloglike_v <- function(theta_v, mu, theta_vec, ud, y, mu_v, inv_Cov, item_vec, user_id, item_id){
  
  #�����m���̐ݒ�
  uv <- as.numeric((theta_vec * theta_v[item_id, ]) %*% vec)
  logit_exp <- exp(mu + uv + ud)   #���W�b�g�̎w���֐�
  prob <- logit_exp / (1 + logit_exp)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  er <- theta_v - mu_v
  dlogit <- y*theta_vec - theta_vec*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -as.matrix(item_vec %*% dlogit + dmvn)
  return(LLd)
}

#�s�񕪉��̃p�����[�^�̃��[�v�t���b�O�@�������֐�
leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, mu, theta_vec, ud, y, mu_v, inv_Cov_v, item_vec, user_id, item_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, mu, theta_vec, ud, y, mu_v, inv_Cov_v, item_vec, user_id, item_id) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##���[�J�[�̍s�񕪉��̃p�����[�^�̃T���v�����O�ɕK�v�Ȋ֐�
#�s�񕪉��̃p�����[�^�̑ΐ����㕪�z�̔����֐�
dloglike_d <- function(theta_d, mu, theta_vec, uv, y, mu_d, inv_Cov, company_vec, user_id, company_id){
  
  #�����m���̐ݒ�
  ud <- as.numeric((theta_vec * theta_d[company_id, ]) %*% vec)
  logit_exp <- exp(mu + uv + ud)   #���W�b�g�̎w���֐�
  prob <- logit_exp / (1 + logit_exp)   #�m���̌v�Z
  
  #�����֐��̐ݒ�
  er <- theta_d - mu_d
  dlogit <- y*theta_vec - theta_vec*prob   #���W�X�e�B�b�N��A�̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -as.matrix(company_vec %*% dlogit + dmvn)
  return(LLd)
}

#�s�񕪉��̃p�����[�^�̃��[�v�t���b�O�@�������֐�
leapfrog_d <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, mu, theta_vec, uv, y, mu_d, inv_Cov_d, company_vec, user_id, company_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, mu, theta_vec, uv, y, mu_d, inv_Cov_d, company_vec, user_id, company_id) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##�A���S���Y���̐ݒ�
LL1 <- -100000000   #�ΐ��ޓx�̏����l
R <- 2000
keep <- 2  
iter <- 0
burnin <- 100/keep
disp <- 10
e1 <- 0.001
e2 <- 0.1
e3 <- 0.05
L <- 3

##���O���z�̐ݒ�
#�f���x�N�g���̎��O���z
delta <- rep(0, ncol(x))
inv_tau <- solve(100 * diag(ncol(x)))

#���[�U�[�̊K�w���f���̎��O���z
Deltabar1 <- matrix(rep(0, ncol(u)*k), nrow=ncol(u), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta1 <- 0.01 * diag(rep(1, ncol(u)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- 1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̊K�w���f���̎��O���z
Deltabar2 <- matrix(rep(0, ncol(v)*k), nrow=ncol(v), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta2 <- 0.01 * diag(rep(1, ncol(v)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- 1   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#���[�J�[�̊K�w���f���̎��O���z
Deltabar3 <- matrix(rep(0, ncol(d)*k), nrow=ncol(d), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta3 <- 0.01 * diag(rep(1, ncol(d)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu3 <- 1   #�t�E�B�V���[�g���z�̎��R�x
V3 <- nu2 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^


##�p�����[�^�̐^�l
#�f���x�N�g���̃p�����[�^
beta <- betat

#�K�w���f���̃p�����[�^
alpha_u1 <- alpha_ut1; Cov_u1 <- Cov_ut1; inv_Cov_u1 <- solve(Cov_u1)
alpha_u2 <- alpha_ut2; Cov_u2 <- Cov_ut2; inv_Cov_u2 <- solve(Cov_u2)
alpha_v <- alpha_vt; Cov_v1 <- Cov_vt; inv_Cov_v <- solve(Cov_v)
alpha_d <- alpha_dt; Cov_d1 <- Cov_dt; inv_Cov_d <- solve(Cov_d)

#�s�񕪉��̃p�����[�^
theta_u1 <- theta_ut1
theta_u2 <- theta_ut2
theta_v <- theta_vt
theta_d <- theta_dt
uv <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec)
ud <- as.numeric((theta_u2[user_id, ] * theta_d[company_id, ]) %*% vec)


##�����l�̐ݒ�
#�f���x�N�g���̃p�����[�^
res <- glm(y ~ x[, -1], family=binomial)
beta <- as.numeric(res$coefficients)

##�K�w���f���̃p�����[�^�𐶐�
#�K�w���f���̕��U�p�����[�^
Cov_u1 <- 0.1 * diag(k); inv_Cov_u1 <- solve(Cov_u1)   #���[�U�[-�A�C�e���̊K�w���f���̕��U
Cov_u2 <- 0.1 * diag(k); inv_Cov_u2 <- solve(Cov_u2)   #���[�U�[-���[�J�[�̊K�w���f���̕��U
Cov_v <- 0.1 * diag(k); inv_Cov_v <- solve(Cov_v)   #�A�C�e���̊K�w���f���̕��U
Cov_d <- 0.1 * diag(k); inv_Cov_d <- solve(Cov_d)   #���[�J�[�̊K�w���f���̕��U

#�K�w���f���̉�A�W����ݒ�
alpha_u1 <- matrix(0, nrow=ncol(u), ncol=k); mu_u1 <- u %*% alpha_u1
alpha_u2 <- matrix(0, nrow=ncol(u), ncol=k); mu_u2 <- u %*% alpha_u2
alpha_v <- matrix(0, nrow=ncol(v), ncol=k); mu_v <- v %*% alpha_v
alpha_d <- matrix(0, nrow=ncol(d), ncol=k); mu_d <- d %*% alpha_d

#�s�񕪉��̃p�����[�^�𐶐�
theta_u1 <- u %*% alpha_u1 + mvrnorm(hh, rep(0, k), Cov_u1)
theta_u2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u1)
theta_v <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
theta_d <- d %*% alpha_d + mvrnorm(company, rep(0, k), Cov_d)
uv <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec)
ud <- as.numeric((theta_u2[user_id, ] * theta_v[item_id, ]) %*% vec)

##�f�[�^�̐ݒ�
x_col <- ncol(x)
u_col <- ncol(u)
v_col <- ncol(v)
d_col <- ncol(d)

##�T���v�����O���ʂ̃p�����[�^�̕ۑ��p�z��
#���f���p�����[�^�̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=x_col)
THETA_U1 <- array(0, dim=c(hh, k, R/keep))
THETA_U2 <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
THETA_D <- array(0, dim=c(company, k, R/keep))

#�K�w���f���̊i�[�p�z��
ALPHA_U1 <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_U2 <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k, R/keep))
ALPHA_D <- array(0, dim=c(ncol(d), k, R/keep))
COV_U1 <- COV_U2 <- COV_V <- COV_D <- array(0, dim=c(k, k, R/keep))

##�ΐ��ޓx�̊�l
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #�ΐ��ޓx
LLglm <- as.numeric(logLik(res))
exp_logit <- exp(as.numeric(x %*% betat) + uvt + udt)
prob <- exp_logit / (1 + exp_logit)
LLt <- sum(y*log(prob)) + sum((1-y)*log(1-prob))


####HMC�Ńp�����[�^���T���v�����O
for(rp in 1:R){
  
  ##�f���x�N�g���̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(1, rep(0, x_col), diag(x_col))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  betad <- beta
  mu_old <- mu
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog(rold, betad, dloglike, e1, L)
  rnew <- res$r
  betan <- res$z
  mu_new <- x %*% betan
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -(loglike(mu_new, uv, ud, y) + dmv(betan, inv_tau, x_col))# + as.numeric(rnew^2 %*% rep(1, x_col))/2
  Hold <- -(loglike(mu_old, uv, ud, y) + dmv(betad, inv_tau, x_col))# + as.numeric(rold^2 %*% rep(1, x_col))/2

  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(1) #��l���z���痐���𔭐�
  gamma <- min(1, exp(Hold - Hnew))   #�̑𗦂�����
  gamma_beta <- gamma
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- gamma > rand
  beta <- flag*betan + (1-flag)*betad
  mu <- as.numeric(x %*% beta)

  
  ##���[�U�[�̃A�C�e���ɑ΂���s�񕪉��̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(hh, rep(0, k), diag(k))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad <- theta_u1
  uv_old <- uv

  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  theta_vec <- theta_v[item_id, ]
  res <- leapfrog_u1(rold, thetad, dloglike_u1, e2, L)
  rnew <- res$r
  thetan <- res$z
  uv_new <- as.numeric((thetan[user_id, ] * theta_vec) %*% vec)
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  er_new <- thetan - mu_u1
  er_old <- thetad - mu_u1
  Hnew <- -(LLho(mu, uv_new, ud, y, user_vec) + dmv(er_new, inv_Cov_u1, k)) + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -(LLho(mu, uv_old, ud, y, user_vec) + dmv(er_old, inv_Cov_u1, k)) + as.numeric(rold^2 %*% rep(1, k))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma_u1 <- mean(gamma)
  

  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_u1 <- flag*thetan + (1-flag)*thetad
  uv <- as.numeric((theta_u1[user_id, ] * theta_vec) %*% vec)
  
  
  ##�A�C�e���̍s�񕪉��̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  theta_vec <- theta_u1[user_id, ]
  rold <- mvrnorm(item, rep(0, k), diag(k))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad <- theta_v
  uv_old <- uv
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog_v(rold, thetad, dloglike_v, e2, L)
  rnew <- res$r
  thetan <- res$z
  uv_new <- as.numeric((theta_vec * thetan[item_id, ]) %*% vec)
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  er_new <- thetan - mu_v
  er_old <- thetad - mu_v
  Hnew <- -(LLho(mu, uv_new, ud, y, item_vec) + dmv(er_new, inv_Cov_v, k)) + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -(LLho(mu, uv_old, ud, y, item_vec) + dmv(er_old, inv_Cov_v, k)) + as.numeric(rold^2 %*% rep(1, k))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(item)   #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma_v <- mean(gamma)
  

  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_v <- flag*thetan + (1-flag)*thetad
  uv <- as.numeric((theta_vec * theta_v[item_id, ]) %*% vec)
  
  
  ##���[�U�[�̃��[�J�[�ɑ΂���s�񕪉��̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(hh, rep(0, k), diag(k))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad <- theta_u2
  ud_old <- ud
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  theta_vec <- theta_d[company_id, ]
  res <- leapfrog_u2(rold, thetad, dloglike_u2, e2, L)
  rnew <- res$r
  thetan <- res$z
  ud_new <- as.numeric((thetan[user_id, ] * theta_vec) %*% vec)
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  er_new <- thetan - mu_u2
  er_old <- thetad - mu_u2
  Hnew <- -(LLho(mu, uv, ud_new, y, user_vec) + dmv(er_new, inv_Cov_u2, k)) + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -(LLho(mu, uv, ud_old, y, user_vec) + dmv(er_old, inv_Cov_u2, k)) + as.numeric(rold^2 %*% rep(1, k))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh)   #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma_u2 <- mean(gamma)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_u2 <- flag*thetan + (1-flag)*thetad
  ud <- as.numeric((theta_u2[user_id, ] * theta_vec) %*% vec)
  
  
  ##���[�J�[�̍s�񕪉��̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(company, rep(0, k), diag(k))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad <- theta_d
  ud_old <- ud
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  theta_vec <- theta_u2[user_id, ]
  res <- leapfrog_d(rold, thetad, dloglike_d, e3, L)
  rnew <- res$r
  thetan <- res$z
  ud_new <- as.numeric((theta_vec * thetan[company_id, ]) %*% vec)
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  er_new <- thetan - mu_d
  er_old <- thetad - mu_d
  Hnew <- -(LLho(mu, uv, ud_new, y, company_vec) + dmv(er_new, inv_Cov_d, k)) + as.numeric(rnew^2 %*% rep(1, k))/2
  Hold <- -(LLho(mu, uv, ud_old, y, company_vec) + dmv(er_old, inv_Cov_d, k)) + as.numeric(rold^2 %*% rep(1, k))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(company)   #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma_d <- mean(gamma)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_d <- flag*thetan + (1-flag)*thetad; theta_d / 2
  ud <- as.numeric((theta_vec * theta_d[company_id, ]) %*% vec)
  
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #���[�U�[�̍s�񕪉��̃p�����[�^���T���v�����O
  out_u1 <- rmultireg(theta_u1, u, Deltabar1, ADelta1, nu1, V1)
  out_u2 <- rmultireg(theta_u2, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u1 <- out_u1$B; alpha_u2 <- out_u2$B
  Cov_u1 <- out_u1$Sigma; Cov_u2 <- out_u2$Sigma
  mu_u1 <- u %*% alpha_u1; mu_u2 <- u %*% alpha_u2
  inv_Cov_u1 <- solve(Cov_u1); inv_Cov_u2 <- solve(Cov_u2)
  
  #�A�C�e���̍s�񕪉��̃p�����[�^���T���v�����O
  out_v <- rmultireg(theta_v, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v <- out_v$B
  Cov_v <- out_v$Sigma
  mu_v <- v %*% alpha_v
  inv_Cov_v <- solve(Cov_v)
  
  #���[�J�[�̍s�񕪉��̃p�����[�^���T���v�����O
  out_d <- rmultireg(theta_d, d, Deltabar3, ADelta3, nu3, V3)
  alpha_d <- out_d$B
  Cov_d <- out_d$Sigma
  mu_d <- d %*% alpha_d
  inv_Cov_d <- solve(Cov_d)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    #���f���̃p�����[�^���i�[
    BETA[mkeep, ] <- beta
    THETA_U1[, , mkeep] <- theta_u1
    THETA_U2[, , mkeep] <- theta_u2
    THETA_V[, , mkeep] <- theta_v
    THETA_D[, , mkeep] <- theta_d
    
    #�K�w���f���̃p�����[�^
    ALPHA_U1[, , mkeep] <- alpha_u1 
    ALPHA_U2[, , mkeep] <- alpha_u2
    ALPHA_D[, , mkeep] <- alpha_v 
    ALPHA_D[, , mkeep] <- alpha_d 
    COV_U1[, , mkeep] <- Cov_u1
    COV_U2[, , mkeep] <- Cov_u2
    COV_V[, , mkeep] <- Cov_v
    COV_D[, , mkeep] <- Cov_d
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    logit_exp <- exp(mu + uv + ud)   #���W�b�g�̎w���֐�
    prob <- logit_exp / (1 + logit_exp)   #�w���m��
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))   #�ΐ��ޓx
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLt, LLglm, LLst))
    print(round(c(gamma_beta, gamma_u1, gamma_u2, gamma_v, gamma_d), 3))
    print(rbind(beta, betat))
  }
}
