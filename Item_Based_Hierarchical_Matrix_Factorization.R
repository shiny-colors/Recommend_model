#####Item Based Hierarchical Matrix Factorization#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(monomvn)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

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
item <- 40000   #���A�C�e����
trigar_item <- 4000   #�g���K�[�̃A�C�e����
display <- 10   #�\������
w <- extraDistr::rtpois(trigar_item, rgamma(trigar_item, 1.85, 0.075), a=0, b=Inf)   #�A�C�e�����Ƃ̉{����

##�A�C�e���̓����s��𐶐�
#�g�s�b�N������ݒ�
topic <- 250
theta <- extraDistr::rdirichlet(1, rep(2.5, topic))
Zi <- as.numeric(rmnom(item, 1, theta) %*% 1:topic)

#�p�����[�^�𐶐�
const <- 10
topic1 <- 60; topic2 <- 80
k1 <- 9; k2 <- 8
alpha1 <- extraDistr::rdirichlet(topic, rep(0.1, topic1))*const
alpha2 <- extraDistr::rdirichlet(topic, rep(0.1, topic2))*const
alpha3 <- matrix(rgamma(topic*k1, 4.5, 11.0), nrow=topic, ncol=k1)
alpha4 <- matrix(rbeta(topic*k2, 2.0, 3.0), nrow=topic, ncol=k2)

#�����x�N�g���𐶐�
topic_data1 <- matrix(0, nrow=item, ncol=topic1)
topic_data2 <- matrix(0, nrow=item, ncol=topic2)
data1 <- matrix(0, nrow=item, ncol=k1); data2 <- matrix(0, nrow=item, ncol=k2)

#�g�s�b�N���ƂɃf�[�^�𐶐�
for(j in 1:topic){
  index_z <- which(Zi==j); n <- length(index_z)
  topic_data1[index_z, ] <- extraDistr::rdirichlet(n, alpha1[j, ])
  topic_data2[index_z, ] <- extraDistr::rdirichlet(n, alpha2[j, ])
  data1[index_z, ] <- abs(mvrnorm(n, alpha3[j, ], runif(k1, 0.01, 0.2) * diag(k1)))
  data2[index_z, ] <- matrix(rbinom(n*k2, 1, alpha4[j, ]), nrow=n, ncol=k2, byrow=T)
}
u <- cbind(1, topic_data1, topic_data2, data1, data2)   #�f�[�^�̌���
rm(topic_data1); rm(topic_data2); rm(data1); rm(data2)
gc(); gc()


##�^�[�Q�b�g�A�C�e���̐ݒ�
#�A�C�e���Ԃ̃R�T�C���ގ��x
ab <- u[1:trigar_item, ] %*% t(u)   #�x�N�g���Ԃ̓���
b <- sqrt(rowSums(u^2))   #���K���W��

#�R�T�C���ގ��x����^�[�Q�b�g�A�C�e�����擾
get <- 50; type <- 6
trigar_id_list <-  target_id_list <- no_id_list <- list()

for(i in 1:trigar_item){
  if(i%%100==0){
    print(i)
  }
  a <- sqrt(sum(u[i, ]^2))
  score <- ab[i, ] / (a * b)
  sortlist <- order(score, decreasing = TRUE)
  get_id <- sortlist[1:(get+1)][-1]
  
  allocation <- matrix(0, nrow=type, ncol=display)
  for(j in 1:type){
    repeat{
      target_vec <- rmnom(1, display, seq(get, 15, length=get))
      if(max(target_vec)==1){
        break
      }
    }
    allocation[j, ] <- get_id[target_vec==1]
  }
  z <- as.numeric(rmnom(w[i], 1, rep(1, type)) %*% 1:type)   
  target_id_list[[i]] <- as.numeric(t(allocation[z, ]))
  trigar_id_list[[i]] <- rep(i, length(target_id_list[[i]]))
  no_id_list[[i]] <- rep(1:display, length(z))
}
#���X�g��ϊ�
trigar_id0 <- unlist(trigar_id_list)
target_id0 <- unlist(target_id_list)
trigar_item <- length(unique(trigar_id0))
target_item <- length(unique(target_id0))
u_trigar <- u[unique(trigar_id0), ]
u_target <- u[unique(target_id0), ]


#ID���Đݒ�
trigar_id <- trigar_id0
target_id <- left_join(data.frame(id=target_id0),
                       data.frame(id=unique(target_id0), no=1:length(unique(target_id0))), by="id")$no
joint_id <- rep(1:sum(w), rep(display, sum(w)))
no_id <- unlist(no_id_list)
N <- length(trigar_id)

##�f���x�N�g���𐶐�
#�g�ݍ��킹��id���쐬
comb_id0 <- paste(trigar_id, target_id)
comb_id <- left_join(data.frame(id=comb_id0, stringsAsFactors = FALSE),
                     data.frame(id=unique(comb_id0), no=1:length(unique(comb_id0)), stringsAsFactors = FALSE), by="id")$no

#�A�C�e���ގ��x�𐶐�
x1 <- rbeta(length(unique(comb_id0)), 2.0, 5.0)[comb_id]

#�A�C�e���^�C�v�̗ގ��x�𐶐�
v <- as.numeric(rmnom(item, 1, extraDistr::rdirichlet(item, rep(2.0, 150))) %*% 1:150)
v_id0 <- paste(v[trigar_id0], v[target_id0])
v_id <- left_join(data.frame(id=v_id0, stringsAsFactors = FALSE),
                  data.frame(id=unique(v_id0), no=1:length(unique(v_id0)), stringsAsFactors = FALSE), by="id")$no
x2 <- rbeta(length(unique(v_id)), 1.5, 3.5)[v_id]

#�f�[�^������
x <- cbind(1, x1, x2)
colnames(x) <- c("intercept", "item", "type")


##�����ϐ��𐶐�
repeat{
  
  #�f���x�N�g���̃p�����[�^�𐶐�
  beta <- betat <- c(-3.2, 0.5, 0.6)
  sigma <- sigamt <- 1.0
  
  #�K�w���f���̃p�����[�^�𐶐�
  alpha_u1 <- alpha_ut1 <- matrix(rnorm(ncol(u)*k, 0, 0.25), nrow=ncol(u_trigar), ncol=k)
  alpha_u2 <- alpha_ut2 <- matrix(rnorm(ncol(u)*k, 0, 0.25), nrow=ncol(u_target), ncol=k)
  Cov_u1 <- Cov_ut1 <- runif(k, 0.01, 0.1) * diag(k)
  Cov_u2 <- Cov_ut2 <- runif(k, 0.01, 0.1) * diag(k)
  
  #�����s��̃p�����[�^�𐶐�
  W <- WT<- u_trigar %*% alpha_u1 + mvrnorm(trigar_item, rep(0, k), Cov_u1)
  H <- HT <- u_target %*% alpha_u2 + mvrnorm(target_item, rep(0, k), Cov_u2)
  
  #���݌��p����R���o�[�W�����𐶐�
  mu <- as.numeric(x %*% beta + (W[trigar_id, ] * H[target_id, ]) %*% rep(1, k))
  U <- mu + rnorm(N, 0, sigma)   #���݌��p
  y <- ifelse(U > 0, 1, 0)   #�R���o�[�W�����x�N�g���ɕϊ�
  if(mean(y) > 0.05 & mean(y) < 0.1){
    break
  }
}

####�}���R�t�A�������e�J�����@��Item Based Hierarchical Matrix Factorization�𐄒�####
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
trigar_index <- target_index <- list()
t_id1 <- t_id2 <- list()

for(i in 1:trigar_item){
  trigar_index[[i]] <- which(trigar_id==i)
  t_id1[[i]] <- target_id[trigar_index[[i]]]
}
for(i in 1:target_item){
  target_index[[i]] <- which(target_id==i)
  t_id2[[i]] <- trigar_id[target_index[[i]]]
}
vec <- rep(1, k)

##���O���z�̐ݒ�
#�K�w���f���̎��O���z
Deltabar1 <- Deltabar2 <- matrix(0, nrow=ncol(u), ncol=k)
ADelta1 <- 0.1 * diag(rep(1, ncol(u))); ADelta2 <- 0.01 * diag(rep(1, ncol(u)))
nu1 <- nu2 <- k + 1 
V1 <- V2 <- nu1 * diag(rep(1, k))

#�f���x�N�g���̎��O���z
gamma <- rep(0, ncol(x))
tau <- 100 * diag(ncol(x))
inv_tau <- solve(tau)

##�p�����[�^�̐^�l
beta <- betat; names(betat) <- colnames(x)
alpha_u1 <- alpha_ut1; alpha_u2 <- alpha_ut2
trigar_mu <- u_trigar %*% alpha_u1; target_mu <- u_target %*% alpha_u2
Cov_u1 <- Cov_ut1; Cov_u2 <- Cov_ut2
inv_Cov_u1 <- solve(Cov_u1); inv_Cov_u2 <- solve(Cov_u2)
W <- WT; H <- HT

#�s�񕪉��̃p�����[�^
WH <- as.numeric((W[trigar_id, ] * H[target_id, ]) %*% vec)


##�p�����[�^�̏����l
beta <- c(-3.0, 0.0, 0.0)
alpha_u1 <- matrix(rnorm(ncol(u)*k, 0, 0.1), nrow=ncol(u), ncol=k)
alpha_u2 <- matrix(rnorm(ncol(u)*k, 0, 0.1), nrow=ncol(u), ncol=k)
Cov_u1 <- Cov_u2 <- 0.01 * diag(k)
W <- u_trigar %*% alpha_u1 + mvrnorm(trigar_item, rep(0, k), Cov_u1)
H <- u_target %*% alpha_u2 + mvrnorm(target_item, rep(0, k), Cov_u2)

#�s�񕪉��̃p�����[�^
WH <- as.numeric((W[trigar_id, ] * H[target_id, ]) %*% vec)


##�T���v�����O���ʂ̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(x))
W_array <- array(0, dim=c(trigar_item, k, R/keep))
H_array <- array(0, dim=c(target_item, k, R/keep))
ALPHA_U1 <- ALPHA_U2 <- array(0, dim=c(ncol(u), k, R/keep))
COV_U1 <- COV_U2 <- array(0, dim=c(k, k, R/keep))

##�f�[�^�̒萔��ݒ�
xx <- t(x) %*% x + inv_tau; inv_xx <- solve(xx)

##�ؒf�̈���`
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##�ΐ��ޓx�̊�l
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #�ΐ��ޓx


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�ؒf���K���z������݌��p���T���v�����O
  mu <- as.numeric(x %*% beta + WH)   #���݌��p�̊��Ғl
  U <- extraDistr::rtnorm(N, mu, sigma, a, b)   #���݌��p�̐���
  
  ##�f���x�N�g���̃p�����[�^���T���v�����O
  #�����ϐ���ݒ�
  u_er <- U - WH   #���f���덷
  
  #���ϗʐ��K���z�̃p�����[�^��ݒ�
  xy <- t(x) %*% u_er
  mu <- as.numeric(inv_xx %*% xy)   #���ϗʐ��K���z�̕��σx�N�g��
  
  #���ϗʐ��K���z����f���x�N�g���̃p�����[�^���T���v�����O
  beta <- mvrnorm(1, mu, sigma^2*inv_xx)
  
  ##�g���K�[�̓����s����T���v�����O
  #�����ϐ���ݒ�
  trigar_er <- U - as.numeric(x %*% beta)   #���f���덷
  
  for(j in 1:trigar_item){
    #�����s��̎��㕪�z�̃p�����[�^��ݒ�
    X <- H[t_id1[[j]], ]
    Xy <- t(X) %*% trigar_er[trigar_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u1)
    mu <- inv_XXV %*% (Xy + inv_Cov_u1 %*% trigar_mu[j, ])
    
    #���ϗʐ��K���z��������x�N�g�����T���v�����O
    W[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  
  ##�^�[�Q�b�g�̓����s����T���v�����O
  #�����ϐ���ݒ�
  target_er <- U - as.numeric(x %*% beta)   #���f���덷
  
  for(j in 1:target_item){
    #�����s��̎��㕪�z�̃p�����[�^��ݒ�
    X <- W[t_id2[[j]], , drop=FALSE]
    Xy <- t(X) %*% target_er[target_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u2)
    mu <- inv_XXV %*% (Xy + inv_Cov_u2 %*% target_mu[j, ])
    
    #���ϗʐ��K���z��������x�N�g�����T���v�����O
    H[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  
  #�s�񕪉��̃p�����[�^���X�V
  WH <- as.numeric((W[trigar_id, ] * H[target_id, ]) %*% vec)
 
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������g���K�[�̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(W, u_trigar, Deltabar1, ADelta1, nu1, V1)
  alpha_u1 <- out$B
  trigar_mu <- u_trigar %*% alpha_u1   #�g���K�[�A�C�e���̓����s��̕��ύ\��
  Cov_u1 <- out$Sigma
  inv_cov_u1 <- solve(Cov_u1)
  
  ##�A�C�e�������s��̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʉ�A���f������^�[�Q�b�g�̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(H, u_target, Deltabar2, ADelta2, nu2, V2)
  alpha_u2 <- out$B
  target_mu <- u_target %*% alpha_u2   #�^�[�Q�b�g�A�C�e���̓����s��̕��ύ\��
  Cov_u2 <- out$Sigma
  inv_cov_u2 <- solve(Cov_u2)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    
    #���f���p�����[�^�̊i�[
    BETA[mkeep, ] <- beta
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    
    #�K�w���f���̃p�����[�^���i�[
    ALPHA_U1[, , mkeep] <- alpha_u1
    ALPHA_U2[, , mkeep] <- alpha_u2
    COV_U1[, , mkeep] <- Cov_u1
    COV_U2[, , mkeep] <- Cov_u2
  } 
  
  if(rp%%disp==0){
    #�ΐ��ޓx���Z�o
    mu <- as.numeric(x %*% beta) + WH   #���݌��p�̊��Ғl
    prob <- pnorm(mu, 0, sigma)   #�w���m��
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #�ΐ��ޓx
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLst))
    print(round(c(beta, betat), 3))
  }
}

round(prob, 3)
mean(y)
cbind(y, prob)
