#####Multinomial logit based Latent Factor Model#####
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
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
#set.seed(34027)
##�f�[�^�̐ݒ�
r <- 7   #��ꐔ
select <- 7   #�I������
hh <- 5000   #����Ґ�
item <- 1500   #�A�C�e����
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #�w���ڐG��
hhpt <- sum(pt)   #�S�T���v����
allocation_index <- matrix(1:(r*(select-1)), nrow=select-1, ncol=r, byrow=T)   #���̊���

#ID��ݒ�
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #�f�[�^�̌���


##�����ϐ��̔���
#�����t�������ϐ�
k11 <- 3; k12 <- 3; k13 <- 4
k1 <- k11 + k12 + k13-1
X1 <- array(0, dim=c(hhpt, k1, select))
Data1 <- matrix(0, nrow=hhpt*select, ncol=k1)

for(j in 1:select){
  x1 <- matrix(runif(hhpt*k11, 0, 1), nrow=hhpt, ncol=k11)
  x2 <- matrix(rbinom(hhpt*k12, 1, rep(runif(k12, 0.15, 0.6), rep(hhpt, k12))), nrow=hhpt, ncol=k12)
  x3 <- rmnom(hhpt, 1, extraDistr::rdirichlet(1, rep(1.5, k13))); x3 <- x3[, -which.min(colSums(x3))]
  X1[, , j] <- cbind(x1, x2, x3)
} 
for(j in 1:k1){
  Data1[, j] <- as.numeric(t(X1[, j, ]))   #�����ϐ����x�N�g���ϊ�
}
user_id_vec <- rep(user_id, rep(select, hhpt))


#�����������ϐ�
k21 <- 2; k22 <- 3; k23 <- 4
x1 <- matrix(runif(hhpt*k21, 0, 1), nrow=hhpt, ncol=k21)
x2 <- matrix(0, nrow=hhpt, ncol=k22)
for(j in 1:k22){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k23, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
data <- cbind(1, x1, x2, x3)   #�f�[�^�̌���

#�����������ϐ����x�N�g���ϊ�
Data2 <- c()
for(j in 1:ncol(data)){
  print(j)
  data_list <- list()
  for(i in 1:hhpt){
    data_list[[i]] <- rbind(diag(data[i, j], select-1), 0)
  }
  Data2 <- cbind(Data2, do.call(rbind, data_list))
}
Data <- cbind(Data2, Data1)
sparse_data <- as(Data, "CsparseMatrix")   #�X�p�[�X�s��ɕϊ�


##���[�U�[�̊K�w���f���̐����ϐ��̔���
#���[�U�[�̐����ϐ�
k1 <- 1; k2 <- 3; k3 <- 5
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #�f�[�^������

#�A�C�e���̐����ϐ�
k1 <- 2; k2 <- 2; k3 <- 4
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #�f�[�^������


##�A�C�e���̊����𐶐�
#�Z�O�����g�����𐶐�
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

#�������z����A�C�e���𐶐�
item_id <- as.numeric(rmnom(hhpt, 1, phi[z[user_id], ]) %*% 1:item)
freq <- plyr::count(item_id); freq$x <- as.character(freq$x)
hist(freq$freq, breaks=25, col="grey", xlab="�A�C�e���̍w���p�x", main="�A�C�e���̍w���p�x���z")
gc(); gc()


####�����ϐ��̔���####
rp <- 0
repeat { 
    
  ##�p�����[�^�𐶐�
  #���W�b�g�̃p�����[�^�𐶐�
  k1 <- ncol(Data)   #�p�����[�^��
  alpha <- alphat <- matrix(rnorm(k1*ncol(u), 0, 0.5), nrow=ncol(u), ncol=k1)   #�K�w���f���̉�A�W��
  Cov <- Covt <- runif(k1, 0.05, 0.2) * diag(k1)   #�K�w���f���̕��U
  theta <- thetat <- u %*% alpha + mvrnorm(hh, rep(0, k1), Cov)   #���[�U�[�ʉ�A�W��
  
  #���[�U�[�̍s�񕪉��̃p�����[�^
  k2 <- r*(select-1)
  alpha_u <- alpha_ut <- matrix(rnorm(ncol(u)*k2, 0, 0.5), nrow=ncol(u), ncol=k2)
  Cov_u <- Cov_ut <- runif(k2, 0.01, 0.15) * diag(k2)
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k2), Cov_u)
  
  #�A�C�e���̍s�񕪉��̃p�����[�^
  alpha_v <- alpha_vt <- matrix(rnorm(ncol(v)*k2, 0, 0.5), nrow=ncol(v), ncol=k2)
  Cov_v <- Cov_vt <- runif(k2, 0.01, 0.15) * diag(k2)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k2), Cov_v)
  
  ##�����ϐ��𐶐�
  #���W�b�g��ݒ�
  mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #�f���x�N�g���̊��Ғl
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #�s�񕪉��̊��Ғl
  } 
  logit_exp <- exp(mu + uv)   #���W�b�g�̊��Ғl�̎w��
  
  #�I�����ʂ𐶐�
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))   #�I���m��
  y <- rmnom(hhpt, 1, prob)   #�������z���牞���ϐ��𐶐�
  
  #��~����
  print(round(c(mean(rowMaxs(prob)), min(colSums(y))), 3))
  if(mean(rowMaxs(prob)) > 0.725 & min(colSums(y)) > hhpt/(select*5)) break
}
colSums(y)
hist(prob[y==1], main="�I�����ʂ̑I���m��", xlab="�I���m��", breaks=25, col="grey")
hist(rowMaxs(prob), main="�����Ƃ������I���m��", xlab="�I���m��", breaks=25, col="grey")

#�I�u�W�F�N�g�̏���
rm(data_list); rm(Data1); rm(Data2)
gc(); gc()


####�n�~���g�j�A�������e�J�����@�Ńp�����[�^���T���v�����O####
##�ΐ����㕪�z���v�Z����֐�
#�ΐ��ޓx���v�Z����֐�
loglike <- function(mu, uv, y, hh, select, dt){
  
  #���W�b�g���f���̑ΐ��ޓx
  logit_exp <- exp(mu + uv)   #���W�b�g�̊��Ғl�̎w��
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))   #�I���m��
  LLi_logit <- log(as.numeric((y * prob) %*% rep(1, select)))
  
  #���[�U�[���Ƃ̑ΐ��ޓx
  LLi <- as.numeric(dt %*% LLi_logit)
  return(LLi)
}

#���ϗʐ��K���z�̑ΐ����O���z
dmv <- function(er, inv_Cov, k){
  Li <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k))
  return(Li)
}

##�f���x�N�g���̃T���v�����O�ɕK�v�Ȋ֐�
#�f���x�N�g���̑ΐ����㕪�z�̔����֐�
dloglike <- function(theta, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt){
  
  #�����m���̐ݒ�
  logit_exp <- exp(matrix((Data * theta[user_id_vec, ]) %*% rep(1, k1), nrow=hhpt, ncol=select, byrow=T) + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  prob_vec <- as.numeric(t(prob))
  
  #�����֐��̐ݒ�
  er <- theta - alpha_mu
  dlogit <- (y_vec - prob_vec) * Data   #���W�b�g���f���̑ΐ��ޓx�̔����֐�
  dmvn <- -t(inv_Cov %*% t(er))   #���ϗʐ��K���z�̑ΐ����O���z�̔����֐�

  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -(user_dt %*% dlogit + dmvn)
  return(LLd)
}

#�f���x�N�g���̃��[�v�t���b�O�@�������֐�
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##���[�U�[�̓����s��̃T���v�����O�ɕK�v�Ȋ֐�
#���[�U�[�̍s�񕪉��̑ΐ����㕪�z�̔����֐�
dloglike_u <- function(theta_u, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, r, k2,
                       user_id, user_vec, allocation_index){
  
  #�����m���̐ݒ�
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_u[user_id, ] * theta_dv
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  logit_exp <- exp(mu + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  
  #�����֐��̐ݒ�
  er <- theta_u - u_mu
  dlogit <- matrix(0, nrow=hhpt, ncol=k2)
  for(j in 1:(select-1)){
    dlogit[, allocation_index[j, ]] <- (y[, j] - prob[, j]) * theta_dv[, allocation_index[j, ]]   #���W�b�g�̔����֐�
  }
  dmvn <- -t(inv_Cov_u %*% t(er))   #���ϗʐ��K���z�̔����֐�

  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -(as.matrix(user_vec %*% dlogit) + dmvn)
  return(LLd)
}

#���[�U�[�̍s�񕪉��̃��[�v�t���b�O�@�������֐�
leapfrog_u <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, g, k2, user_id, user_vec, allocation_index) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, g, k2, user_id, user_vec, allocation_index) / 2
    list(r=r2, z=z2) # 1��̈ړ���̉^���ʂƍ��W
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##�A�C�e���̓����s��̃T���v�����O�ɕK�v�Ȋ֐�
#�A�C�e���̍s�񕪉��̑ΐ����㕪�z�̔����֐�
dloglike_v <- function(theta_v, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, r, k2,
                       item_id, item_vec, allocation_index){
  
  #�����m���̐ݒ�
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_du * theta_v[item_id, ]
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  logit_exp <- exp(mu + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  
  #�����֐��̐ݒ�
  er <- theta_v - v_mu
  dlogit <- matrix(0, nrow=hhpt, ncol=k2)
  for(j in 1:(select-1)){
    dlogit[, allocation_index[j, ]] <- (y[, j] - prob[, j]) * theta_du[, allocation_index[j, ]]
  }
  dmvn <- -t(inv_Cov_v %*% t(er))
  
  #�ΐ����㕪�z�̔����֐��̘a
  LLd <- -(as.matrix(item_vec %*% dlogit) + dmvn)
  return(LLd)
}

#�A�C�e���̍s�񕪉��̃��[�v�t���b�O�@�������֐�
leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, g, k2, item_id, item_vec, allocation_index) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, g, k2, item_id, item_vec, allocation_index) / 2
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
R <- 1000
keep <- 2  
iter <- 0
burnin <- 100/keep
disp <- 5
e1 <- 0.1
e2 <- 0.05
L1 <- 3
L2 <- 5

##�C���f�b�N�X�ƃf�[�^�̐ݒ�
#�C���f�b�N�X�̐ݒ�
user_vec <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
user_dt <- sparseMatrix(user_id_vec, 1:length(user_id_vec), x=rep(1, length(user_id_vec)), dims=c(hh, length(user_id_vec)))
item_vec <- sparseMatrix(item_id, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))

#�f�[�^�̐ݒ�
g <- r
k1 <- ncol(Data)
k2 <- r*(select-1)
y_vec <- as.numeric(t(y))

##���O���z�̐ݒ�
#�f���x�N�g���̊K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=ncol(u), ncol=k1)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta <- 0.01 * diag(ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu <- select-1   #�t�E�B�V���[�g���z�̎��R�x
V <- nu * diag(k1) #�t�E�B�V���[�g���z�̃p�����[�^

#���[�U�[�̍s�񕪉��̃p�����[�^�̊K�w���f���̎��O���z
Deltabar_u <- matrix(0, nrow=ncol(u), ncol=k2)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta_u <- 0.01 * diag(ncol(u))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- select-1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(k2) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̍s�񕪉��̃p�����[�^�̊K�w���f���̎��O���z
Deltabar_v <- matrix(0, nrow=ncol(v), ncol=k2)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta_v <- 0.01 * diag(ncol(v))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- select-1   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(k2) #�t�E�B�V���[�g���z�̃p�����[�^


##�p�����[�^�̐^�l
#�K�w���f���̃p�����[�^
alpha <- alphat; alpha_u <- alpha_ut; alpha_v <- alpha_vt
Cov <- Covt; Cov_u <- Cov_ut; Cov_v <- Cov_vt
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
inv_Cov <- solve(Cov); inv_Cov_u <- solve(Cov_u); inv_Cov_v <- solve(Cov_v)

#���f���̃p�����[�^
theta <- thetat
theta_u <- theta_ut
theta_v <- theta_vt

#�p�����[�^�̐ݒ�
mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #�f���x�N�g���̊��Ғl
uv <- matrix(0, nrow=hhpt, ncol=select)
uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
for(j in 1:(select-1)){
  uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #�s�񕪉��̊��Ғl
} 

#�x�X�g�ȑΐ��ޓx
logit_exp <- exp(mu + uv)
prob <- logit_exp / rowSums(logit_exp)
LL_best <- sum(y * log(prob))

##�����l�̐ݒ�
#�K�w���f���̃p�����[�^�̏����l
alpha <- matrix(0, nrow=ncol(u), ncol=k1); alpha_mu <- u %*% alpha
alpha_u <- matrix(0, nrow=ncol(u), ncol=k2); alpha_u <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=k2); alpha_v <- v %*% alpha_v
Cov <- diag(0.01, k1); inv_Cov <- solve(Cov)
Cov_u <- diag(0.01, k2); inv_Cov_u <- solve(Cov_u)
Cov_v <- diag(0.01, k2); inv_Cov_v <- solve(Cov_v)

#���f���̃p�����[�^�̏����l
theta <- alpha_mu + mvrnorm(hh, rep(0, k1), Cov)
theta_u <- u_mu + mvrnorm(hh, rep(0, k2), Cov_u)
theta_v <- v_mu + mvrnorm(item, rep(0, k2), Cov_v)

#�p�����[�^�̐ݒ�
mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #�f���x�N�g���̊��Ғl
uv <- matrix(0, nrow=hhpt, ncol=select)
uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
for(j in 1:(select-1)){
  uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #�s�񕪉��̊��Ғl
} 

##�p�����[�^�̊i�[�p�z��
ALPHA <- array(0, dim=c(ncol(u), k1, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), k2, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k2, R/keep))
COV <- array(0, dim=c(k1, k1, R/keep))
COV_U <- array(0, dim=c(k2, k2, R/keep))
COV_V <- array(0, dim=c(k2, k2, R/keep))

#���f���̃p�����[�^�̊i�[�p�z��
d <- 0
THETA <- matrix(0, nrow=hh, ncol=k1)
THETA_U <- matrix(0, nrow=hh, ncol=k2)
THETA_V <- matrix(0, nrow=item, ncol=k2)


##�ΐ��ޓx�̊�l
#�����l�̑ΐ��ޓx
logit_exp <- exp(mu + uv)
prob <- logit_exp / rowSums(logit_exp)
LL1 <- sum(y * log(prob))

#���σ��f���̑ΐ��ޓx
LLst <- sum(log(colMeans(y)[as.numeric(y %*% 1:select)]))


####HMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##�f���x�N�g���̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(hh, rep(0, k1), diag(k1))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad <- theta
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  res <- leapfrog(rold, thetad, dloglike, e1, L1)
  rnew <- res$r
  thetan <- as.matrix(res$z)
  
  #���W�b�g�̃p�����[�^�̐ݒ�
  mu_old <- matrix(as.numeric((Data * thetad[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  mu_new <- matrix(as.numeric((Data * thetan[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  er_old <- thetad - alpha_mu
  er_new <- thetan - alpha_mu
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -(loglike(mu_new, uv, y, hh, select, user_vec) + dmv(er_new, inv_Cov, k1)) + as.numeric(rnew^2 %*% rep(1, k1))/2
  Hold <- -(loglike(mu_old, uv, y, hh, select, user_vec) + dmv(er_old, inv_Cov, k1)) + as.numeric(rold^2 %*% rep(1, k1))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh) #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma1 <- mean(gamma)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta <- flag*thetan + (1-flag)*thetad
  
  
  ##���[�U�[�̍s�񕪉��̃p�����[�^���T���v�����O
  #�f���x�N�g���̃��W�b�g
  mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(hh, rep(0, k2), diag(k2))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad_u <- theta_u
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  theta_dv <- theta_v[item_id, ]
  res <- leapfrog_u(rold, thetad_u, dloglike_u, e1, L1)
  rnew <- res$r
  thetan_u <- res$z
  
  #���W�b�g�̃p�����[�^�̐ݒ�
  uv_old <- uv
  uv_new <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- thetan_u[user_id, ] * theta_dv
  for(j in 1:(select-1)){
    uv_new[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  er_old <- thetad_u - u_mu
  er_new <- thetan_u - u_mu
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -(loglike(mu, uv_new, y, hh, select, user_vec) + dmv(er_new, inv_Cov_u, k2)) + as.numeric(rnew^2 %*% rep(1, k2))/2
  Hold <- -(loglike(mu, uv_old, y, hh, select, user_vec) + dmv(er_old, inv_Cov_u, k2)) + as.numeric(rold^2 %*% rep(1, k2))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(hh) #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma2 <- mean(gamma)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_u <- flag*thetan_u + (1-flag)*thetad_u
  
  #�s�񕪉��̃p�����[�^���X�V
  theta_du <- theta_u[user_id, ]
  uv_dt <- theta_du * theta_dv
  uv <- matrix(0, nrow=hhpt, ncol=select)
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  
  
  ##�A�C�e���̍s�񕪉��̃p�����[�^���T���v�����O
  #HMC�̐V�����p�����[�^�𐶐�
  rold <- mvrnorm(item, rep(0, k2), diag(k2))   #�W�����ϗʐ��K���z����p�����[�^�𐶐�
  thetad_v <- theta_v
  
  #���[�v�t���b�O�@�ɂ��1�X�e�b�v�ړ�
  theta_du <- theta_u[user_id, ]
  res <- leapfrog_v(rold, thetad_v, dloglike_v, e2, L2)
  rnew <- res$r
  thetan_v <- res$z
  
  #���W�b�g�̃p�����[�^�̐ݒ�
  uv_old <- uv
  uv_new <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_du * thetan_v[item_id, ]
  for(j in 1:(select-1)){
    uv_new[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  er_old <- thetad_v - v_mu
  er_new <- thetan_v - v_mu
  
  #�ړ��O�ƈړ���̃n�~���g�j�A��
  Hnew <- -(loglike(mu, uv_new, y, item, select, item_vec) + dmv(er_new, inv_Cov_v, k2)) + as.numeric(rnew^2 %*% rep(1, k2))/2
  Hold <- -(loglike(mu, uv_old, y, item, select, item_vec) + dmv(er_old, inv_Cov_v, k2)) + as.numeric(rold^2 %*% rep(1, k2))/2
  
  #HMC�@�ɂ��p�����[�^�̍̑�������
  rand <- runif(item)   #��l���z���痐���𔭐�
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #�̑𗦂�����
  gamma3 <- mean(gamma)
  
  #alpha�̒l�Ɋ�Â��V����beta���̑����邩�ǂ���������
  flag <- as.numeric(gamma > rand)
  theta_v <- flag*thetan_v + (1-flag)*thetad_v
  
  #�s�񕪉��̃p�����[�^���X�V
  uv_dt <- theta_du * theta_v[item_id, ]
  uv <- matrix(0, nrow=hhpt, ncol=select)
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������f���x�N�g���̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta, u, Deltabar, ADelta, nu, V)
  alpha <- out$B
  alpha_mu <- u %*% alpha   #�f���x�N�g���̊��Ғl
  Cov <- out$Sigma
  inv_Cov <- solve(Cov)
  
  #���ϗʉ�A���f�����烆�[�U�[�����s��̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_u, u, Deltabar_u, ADelta_u, nu1, V1)
  alpha_u <- out$B
  u_mu <- u %*% alpha_u   #���[�U�[�����s��̊��Ғl
  Cov_u <- out$Sigma
  inv_Cov_u <- solve(Cov_u)
  
  #���ϗʉ�A���f������A�C�e�������s��̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu2, V2)
  alpha_v <- out$B
  v_mu <- v %*% alpha_v   #�A�C�e�������s��̊��Ғl
  Cov_v <- out$Sigma
  inv_Cov_v <- solve(Cov_v)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�p�����[�^�̊i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    ALPHA[, , mkeep] <- alpha
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    COV[, , mkeep] <- Cov
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
  }
  if(rp >= burnin & rp%%keep==0){
    d <- d + 1
    THETA <- THETA + theta
    THETA_U <- THETA_U + theta_u
    THETA_V <- THETA_V + theta_v
  }
  
  #�T���v�����O���ʂ̕\��
  if(rp%%disp==0){
    #�ΐ��ޓx�̌v�Z�ƕ\��
    logit_exp <- exp(mu + uv)
    prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
    LL <- sum(log((y * prob) %*% rep(1, select)))
  
    #���ʂ̊m�F
    print(rp)
    print(c(LL, LL_best, LL1, LLst))
    print(round(cbind(alpha[, 1:(select-1)], alphat[, 1:(select-1)]), 3))
    print(round(rbind(diag(Cov)[1:((select-1)*3)], diag(Covt)[1:((select-1)*3)]), 3))
    print(round(c(gamma1, gamma2, gamma3), 3))
  }
}

