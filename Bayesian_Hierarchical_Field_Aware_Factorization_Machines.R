#####Bayesian Hierarchical Field Aware Factorization Machines#####
library(MASS)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
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
hh <- 5000   #���[�U�[��
item <- 2000   #�A�C�e����
tag <- 150   #�^�O��
pt <- rtpois(hh, rgamma(hh, 30.0, 0.25), a=1, b=Inf)   #�w���ڐG��
hhpt <- sum(pt)
n <- rtpois(item, 0.75, a=0, b=Inf)   #�^�O��
k <- 10   #��ꐔ
vec_k <- rep(1, k)

#ID��ݒ�
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #�f�[�^�̌���
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

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
item_data <- colSums(sparseMatrix(1:hhpt, item_id, x=rep(1, hhpt), dims=c(hhpt, item)))

##�^�O�𐶐�
#�p�����[�^�̐ݒ�
omega <- extraDistr::rdirichlet(topic, rep(0.5, tag))
z <- as.numeric(rmnom(item, 1, t(phi) / rowSums(t(phi))) %*% 1:topic)

#�������z����^�O�𐶐�
max_n <- max(n)
tag_list <- list()
tag_id <- matrix(0, nrow=hhpt, ncol=max_n)
for(j in 1:item){
  repeat { 
    x <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(x)==1){
      tag_list[[j]] <- (x * 1:tag)[x > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
  
}
tag_data <- sparseMatrix(matrix(rep(1:hhpt, max_n), nrow=hhpt, ncol=max_n)[tag_id!=0], tag_id[tag_id!=0], 
                         x=rep(1, sum(tag_id!=0)), dims=c(hhpt, tag))



####�����ϐ��𐶐�####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  ##�K�w���f���̃p�����[�^�𐶐�
  ##���[�U�[�x�[�X�̊K�w���f���̃p�����[�^
  #���U�����U�s���ݒ�
  Cov_ut <- Cov_u <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.0025, 0.25)$covariance
  
  #��A�W����ݒ�
  alpha_u <- matrix(0, nrow=ncol(u), ncol=k)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u[j, ] <- runif(k, -1.3, -0.5)
    } else {
      alpha_u[j, ] <- runif(k, -0.6, 0.7)
    }
  }
  alpha_ut <- alpha_u
  
  #���ϗʉ�A���f�����烆�[�U�[�ʂ̉�A�p�����[�^�𐶐�
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k), Cov_u)   #�e���\�������̃p�����[�^
  
  
  ##�A�C�e���x�[�X�̊K�w���f���̃p�����[�^
  #���U�����U�s���ݒ�
  Cov_vt <- Cov_v <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.025, 0.25)$covariance
  
  #��A�W����ݒ�
  alpha_v <- matrix(0, nrow=ncol(v), ncol=k)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(k, -1.2, -0.4)
    } else {
      alpha_v[j, ] <- runif(k, -0.7, 0.8)
    }
  }
  alpha_vt <- alpha_v
  
  #���ϗʉ�A���f������A�C�e���ʂ̉�A�p�����[�^�𐶐�
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  
  ##���ԃx�[�X�̊K�w���f���̃p�����[�^
  alpha_t <- alpha_tt <- rep(-0.25, k)
  Cov_t <- Cov_tt <- diag(0.25, k)
  theta_t <- theta_tt <- mvrnorm(time, alpha_t, Cov_t)
  
  
  ##���K���z������p�ƍw���x�N�g���𐶐�
  #���݌��p�𐶐�
  WHC0 <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC0 <- WHC0 + theta_u[, j] %o% t(theta_v)[j, ] %o% t(theta_t)[j, ]
  }
  #whc0 <- as.numeric((W0[user_id0, ] * t(H0)[item_id0, ] * t(C0)[time_id0, ]) %*% vec)   #������ł�ok
  whc0 <- as.numeric(WHC0)   #�e���\�����x�N�g���ɕϊ�
  u_vec0 <- whc0 + rnorm(hh*item*time, 0, 1)   #�덷�𐶐�
  
  #�w���x�N�g���ɕϊ�
  y0 <- ifelse(u_vec0 > 0, 1, 0)
  if(mean(y0) > 0.3 & mean(y0) < 0.4) break   #break����
}

##�����x�N�g���𐶐�
#�����m���𐶐�
user_prob <- rbeta(hh, 10, 50)
item_prob <- rbeta(item, 15, 55)
time_prob <- rbeta(time, 60, 140)
prob <- user_prob[user_id0]*item_id0[item_id0]*time_prob[time_id0]

#�x���k�[�C���z���猇���x�N�g���𐶐�
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)
y <- y0[z_vec==1]; u_vec <- u_vec0[z_vec==1]; whc <- whc0[z_vec==1]
hist(u_vec, breaks=25, col="grey", main="���ݓI�ȃX�R�A���z", xlab="�X�R�A")

#�����x�N�g������id���č\��
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
time_id <- time_id0[z_vec==1]


####�}���R�t�A�������e�J�����@�ŊK�w�x�C�Y�e���\�������𐄒�####
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
user_index <- item_index <- time_index <- list()
ui_id <- ut_id <- iu_id <- it_id <- tu_id <- ti_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  ut_id[[i]] <- time_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  it_id[[j]] <- time_id[item_index[[j]]]
}
for(j in 1:time){
  time_index[[j]] <- which(time_id==j)
  tu_id[[j]] <- user_id[time_index[[j]]]
  ti_id[[j]] <- item_id[time_index[[j]]]
}
vec <- rep(1, k)
const1 <- hh / 2.0  #���K���萔
const2 <- item / 2.0

##���O���z�̐ݒ�
#���[�U�[�̊K�w���f���̎��O���z
Deltabar1 <- matrix(rep(0, ncol(u)*k), nrow=ncol(u), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta1 <- 0.01 * diag(rep(1, ncol(u)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu1 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V1 <- nu1 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#�A�C�e���̊K�w���f���̎��O���z
Deltabar2 <- matrix(rep(0, ncol(v)*k), nrow=ncol(v), ncol=k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta2 <- 0.01 * diag(rep(1, ncol(v)))   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu2 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V2 <- nu2 * diag(rep(1, k)) #�t�E�B�V���[�g���z�̃p�����[�^

#���Ԃ̊K�w���f���̎��O���z
Aelta <- rep(0, k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
ADelta3 <- 0.01 * diag(k)   #�K�w���f���̉�A�W���̎��O���z�̕��U
nu3 <- k + 1   #�t�E�B�V���[�g���z�̎��R�x
V3 <- nu3 * diag(rep(1, k))   #�t�E�B�V���[�g���z�̃p�����[�^

##�p�����[�^�̐^�l
alpha_u <- alpha_ut; user_mu <- u %*% alpha_u
alpha_v <- alpha_vt; item_mu <- v %*% alpha_v
alpha_t <- alpha_tt
Cov_u <- Cov_ut; inv_Cov1 <- solve(Cov_u)
Cov_v <- Cov_vt; inv_Cov2 <- solve(Cov_v)
Cov_t <- Cov_tt; inv_Cov3 <- solve(Cov_t)
theta_u <- theta_ut
theta_v <- theta_vt
theta_t <- theta_tt
sigma <- 1

##�p�����[�^�̏����l
#�e���\�������̃p�����[�^
theta_u <- mvrnorm(hh, rep(0, k), 0.25 * diag(k))
theta_v <- mvrnorm(item, rep(0, k), 0.25 * diag(k))
theta_t <- mvrnorm(time, rep(0, k), 0.25 * diag(k))
sigma <- 1

#�K�w���f���̃p�����[�^
alpha_u <- matrix(0, nrow=ncol(u), ncol=k); user_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=k); item_mu <- v %*% alpha_v
alpha_t <- rep(0, k)
Cov_u <- diag(0.2, k); inv_Cov1 <- solve(Cov_u)
Cov_v <- diag(0.2, k); inv_Cov2 <- solve(Cov_v)
Cov_t <- diag(0.2, k); inv_Cov3 <- solve(Cov_t)

##�T���v�����O���ʂ̃p�����[�^�̕ۑ��p�z��
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
THETA_T <- array(0, dim=c(time, k, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k, R/keep))
ALPHA_T <- matrix(0, nrow=R/keep, ncol=k)
COV_U <- array(0, dim=c(k, k, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
COV_T <- array(0, dim=c(k, k, R/keep))

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
  
  ##�ؒf���K���z������݌��p�𐶐�
  mu <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ] * theta_t[time_id, ]) %*% vec)   #���݌��p�̊��Ғl
  U <- rtnorm(mu, sigma, a, b)   #���݌��p�𐶐�
  
  ##�e���\�������̓����s��̃p�����[�^���T���v�����O
  ##���[�U�[�����s����T���v�����O
  for(i in 1:hh){
    #�����x�N�g���̃p�����[�^
    X <- theta_v[ui_id[[i]], ] * theta_t[ut_id[[i]], ]
    Xy <- t(X) %*% U[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov1)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov1 %*% user_mu[i, ])   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- inv_XXV   #���ϗʐ��K���z�̕��U�����U�s��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_u[i, ] <- mvrnorm(1, beta_mu, cov)
  }
  #theta_u <- theta_u / matrix(colSums(theta_u), nrow=hh, ncol=k, byrow=T) * const1
  
  ##�A�C�e�������s����T���v�����O
  for(j in 1:item){
    #�����x�N�g���̃p�����[�^
    X <- theta_u[iu_id[[j]], ] * theta_t[it_id[[j]], ]
    Xy <- t(X) %*% U[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X +inv_Cov2)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov2 %*% item_mu[j, ])   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- inv_XXV   #���ϗʐ��K���z�̕��U�����U�s��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_v[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  #theta_v <- theta_v / matrix(colSums(theta_v), nrow=item, ncol=k, byrow=T) * const2
  
  ##���Ԃ̓����s����T���v�����O
  for(j in 1:time){
    #�����x�N�g���̃p�����[�^
    X <- theta_u[tu_id[[j]], ] * theta_v[ti_id[[j]], ]
    Xy <- t(X) %*% U[time_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov3 %*% alpha_t)   #���ϗʐ��K���z�̕��σx�N�g��
    cov <- inv_XXV   #���ϗʐ��K���z�̕��U�����U�s��
    
    #���ϗʐ��K���z����p�����[�^���T���v�����O
    theta_t[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  
  ##�K�w���f���̃p�����[�^���T���v�����O
  ##���ϗʉ�A���f�����烆�[�U�[�̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_u, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u <- out$B
  user_mu <- u %*% alpha_u   #���[�U�[�����s��̕��ύ\��
  Cov_u <- out$Sigma
  inv_cov1 <- solve(Cov_u)
  
  ##���ϗʉ�A���f������A�C�e���̊K�w���f���̃p�����[�^���T���v�����O
  out <- rmultireg(theta_v, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v <- out$B
  item_mu <- v %*% alpha_v   #�A�C�e�������s��̕��ύ\��
  Cov_v <- out$Sigma
  inv_cov2 <- solve(Cov_v)
  
  ##���Ԃ̊K�w���f���̃p�����[�^���T���v�����O
  #���ϗʐ��K���z���畽�σx�N�g�����T���v�����O
  alpha_mu <- solve(ADelta3 + time*inv_Cov3) %*% (time*inv_Cov3 %*% alpha_t)
  alpha_t <- mvrnorm(1, alpha_mu, Cov_t/time)
  
  #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
  er <- theta_t - matrix(alpha_t, nrow=time, ncol=k, byrow=T)   #�덷���Z�o
  IW <- solve(V3) + t(er) %*% er  
  Sn <- nu3 + time 
  Cov_t <- rwishart(Sn, solve(IW))$IW   #�t�E�B�V���[�ƕ��z����cov���T���v�����O
  inv_Cov3 <- solve(Cov_t)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    THETA_T[, , mkeep] <- theta_t
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    ALPHA_T[mkeep, ] <- alpha_t
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    COV_T[, , mkeep] <- Cov_t
  }
  
  #�ΐ��ޓx���v�Z
  if(rp%%disp==0){
    mu <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ] * theta_t[time_id, ]) %*% vec)   #���S�f�[�^�̊��Ғl
    prob <- pnorm(mu, 0, sigma)   #�w���m��
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #�ΐ��ޓx
    
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(c(LL, LLst))
    print(round(alpha_u, 3))
  }
}
mean(prob[y==0])

matplot(t(THETA_U[1, , ]), type="l")
matplot(t(THETA_V[1, , ]), type="l")
matplot(t(THETA_T[1, , ]), type="l")

round(prob, 3)

