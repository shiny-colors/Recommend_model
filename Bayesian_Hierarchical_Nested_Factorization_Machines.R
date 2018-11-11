#####Bayesian Hierarchical Nested Factorization Machines#####
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

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
r <- 5   #��ꐔ
hh <- 5000   #���[�U�[��
w <- rpois(hh, (rgamma(hh, 12.5, 0.1)))   #���[�U�[���Ƃ̃T���v����
f <- sum(w)   #���T���v����

##ID��ݒ�
u_id <- rep(1:hh, w)
t_id <- as.numeric(unlist(tapply(1:f, u_id, rank)))

id_list <- list()
for(i in 1:hh){
  id_list[[i]] <- which(u_id==i)
}

##�����ϐ����Ó��ɂȂ�܂Ńp�����[�^�̐������J��Ԃ�
rp <- 0
repeat {
  print(rp <- rp + 1)
  
  ##�f���x�N�g���𐶐�
  m1 <- 2; m2 <- 3; m3 <- 4
  z1 <- matrix(runif(f*m1, 0, 1), nrow=f, ncol=m1)
  z2 <- matrix(0, nrow=f, ncol=m2)
  for(j in 1:m2){
    pr <- runif(1, 0.25, 0.55)
    z2[, j] <- rbinom(f, 1, pr)
  }
  z3 <- rmnom(f, 1, runif(m3, 0.2, 1.25)); z3 <- z3[, -which.min(colSums(z3))]
  z <- cbind(1, z1, z2, z3)   #�f�[�^������
  m <- ncol(z)
  
  ##�K�w���f���̐����ϐ�
  m1 <- 1; m2 <- 3; m3 <- 5
  u1 <- matrix(runif(hh*m1, 0, 1), nrow=hh, ncol=m1)
  u2 <- matrix(0, nrow=hh, ncol=m2)
  for(j in 1:m2){
    pr <- runif(1, 0.25, 0.55)
    u2[, j] <- rbinom(hh, 1, pr)
  }
  u3 <- rmnom(hh, 1, runif(m3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
  u <- cbind(1, u1, u2, u3)   #�f�[�^������
  
  
  ##�e���\���𐶐�
  k1 <- 15; k2 <- 5
  X_list <- X1_list <- X2_list <- list()
  x1_vec <- x2_vec <- list()
  allocation_vec <- list()
  
  for(i in 1:hh){
    #�p�����[�^�𐶐�
    lambda <- runif(k1, 0.1, 2.0)
    par <- extraDistr::rdirichlet(k1, rep(1.5, k2))
    
    #�f�[�^�𐶐�
    n <- w[i]*k1
    x1 <- matrix(rtpois(n, rep(lambda, w[i]), a=-Inf, b=8), nrow=w[i], ncol=k1, byrow=T)
    x2 <- matrix(rmnom(n, 1, par) %*% 1:k2, nrow=w[i], ncol=k1, byrow=T)
  
    #�f�[�^���i�[
    storage.mode(x1) <- "integer"
    storage.mode(x2) <- "integer"
    X1_list[[i]] <- x1; x1_vec[[i]] <- as.numeric(t(x1))
    X2_list[[i]] <- x2; x2_vec[[i]] <- as.numeric(t(x2))
    allocation_vec[[i]] <- rep(1:k1, w[i])
  }
  
  #���X�g��ϊ�
  Data1 <- do.call(rbind, X1_list)
  Data2 <- do.call(rbind, X2_list)
  
  ##�p�����[�^�𐶐�
  #�f���x�N�g���̃p�����[�^
  beta <- betat  <- c(5.0, rnorm(m-1, 0.25, 0.4))
  sigma <- sigmat <- 0.2
  
  ##�K�w���f���ƃe���\���̃p�����[�^�𐶐�
  #1�K�w�ڂ̃��f���̃p�����[�^
  alpha1 <- array(0, dim=c(ncol(u), r, k1))
  Cov1 <- array(0, dim=c(r, r, k1))
  theta1 <- array(0, dim=c(k1, r, hh))
  for(j in 1:k1){
    alpha1[, , j] <- rbind(mvrnorm(1, rep(0, r), diag(runif(r, 0.005, 0.075))), 
                           mvrnorm(ncol(u)-1, rep(0.0, r), diag(runif(r, 0.005, 0.075))))
    Cov1[, , j] <- diag(runif(r, 0.005, 0.05))
    theta1[j, , ] <- t(u %*% alpha1[, , j] + mvrnorm(hh, rep(0, r), Cov1[, , j]))
  }
  alphat1 <- alpha1; Covt1 <- Cov1; thetat1 <- theta1
  
  #2�K�w�ڂ̃��f���̃p�����[�^
  alpha2 <- array(0, dim=c(ncol(u), r, k2))
  Cov2 <- array(0, dim=c(r, r, k2))
  theta2 <- array(0, dim=c(k2, r, hh))
  for(j in 1:k2){
    alpha2[, , j] <- rbind(mvrnorm(1, rep(0, r), diag(runif(r, 0.005, 0.075))), 
                           mvrnorm(ncol(u)-1, rep(0.0, r), diag(runif(r, 0.005, 0.075))))
    Cov2[, , j] <- diag(runif(r, 0.005, 0.05))
    theta2[j, , ] <- t(u %*% alpha2[, , j] + mvrnorm(hh, rep(0, r), Cov2[, , j]))
  }
  alphat2 <- alpha2; Covt2 <- Cov2; thetat2 <- theta2

  ##�����ϐ��𐶐�
  y0 <- rep(0, f)
  for(i in 1:hh){
    #�f�[�^�𒊏o
    x1 <- x1_vec[[i]]; x2 <- x2_vec[[i]]
    index <- allocation_vec[[i]]
    
    #�e���\���̕��σx�N�g��
    tensor_dt <- matrix((x1 * theta1[index, , i] * theta2[x2, , i]) %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)
    tensor_mu <- as.numeric(tensor_dt %*% rep(1, k1))
    
    #���f���̕��σx�N�g��
    mu <- as.numeric(z[id_list[[i]], ] %*% beta) + tensor_mu
    y0[id_list[[i]]] <- rnorm(w[i], mu, sigma)
  }
  
  #�����ϐ���break����
  print(c(mean(y0), sd(y0)))
  if(sd(y0) > 1.50 & sd(y0) < 2.25 & mean(y0) > 5.0 & mean(y0) < 6.0){
    break
  }
}

#���������X�R�A��]���f�[�^�ɕϊ�
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y <- round(y0_censor, 0)   #�X�R�A���ۂ߂�
hist(y0, breaks=25, col="grey", xlab="�X�R�A�̐^�l", main="�X�R�A���z�̐^�l")
hist(y, breaks=25, col="grey", xlab="�]���X�R�A", main="�]���X�R�A���z")


####�M�u�X�T���v�����O��Bayesian Hierarchical Nested Factorization Machines�𐄒�####
##�A���S���Y���̐ݒ�
R <- 2000
burnin <- 500
keep <- 2
disp <- 4
iter <- 0
LL <- -1000000000

##�l�X�g�ϐ��̃C���f�b�N�X���쐬
DT1 <- matrix(1:k1, nrow=f, ncol=k1, byrow=T)
index_data2 <- freq_data2 <- index_dt1 <- index_data1 <- list()

for(i in 1:hh){
  #�f�[�^�̐ݒ�
  dt1 <- DT1[id_list[[i]], ]
  data1 <- Data1[id_list[[i]], ]
  data2 <- Data2[id_list[[i]], ]
  
  #�K�v�ȃC���f�b�N�X���쐬
  index1 <- index2 <- index3 <- index4 <- index5 <- list()
  for(j in 1:k2){
    index1[[j]] <- index <- which(as.numeric(t(data2==j))==TRUE)
    index2[[j]] <- rowSums(data2==j)
    index3[[j]] <- as.numeric(t(dt1))[index]
    index4[[j]] <- as.numeric(t(data1))[index]
  }
  #�C���f�b�N�X���i�[
  index_data2[[i]] <- index1
  freq_data2[[i]] <- index2
  index_dt1[[i]] <- index3
  index_data1[[i]] <- index4
}

#�O���[�v���Ƃɘa���Ƃ邽�߂̍s����쐬
dt_freq <- list()

for(i in 1:hh){
  freq <- freq_data2[[i]]
  dt_list <- list()
  
  for(j1 in 1:k2){
    #�C���f�b�N�X���쐬
    vec <- freq[[j1]]
    vec_cumsum2 <- cumsum(vec)
    vec_cumsum1 <- c(1, vec_cumsum2[-w[i]]+1)
    dt <- matrix(0, nrow=sum(freq[[j1]]), ncol=w[i])
    
    for(j2 in 1:w[i]){
      #�s����쐬
      if(vec[j2]==0){
        next
      }
      dt[vec_cumsum1[j2]:vec_cumsum2[j2], j2] <- 1
    }
    dt_list[[j1]] <- as(t(dt), "CsparseMatrix")
  }
  dt_freq[[i]] <- dt_list
}

##���O���z�̐ݒ�
#�K�w���f���̎��O���z
Deltabar <- matrix(0, nrow=ncol(u), ncol=r)
ADelta <- 0.01 * diag(rep(1, ncol(u)))
nu <- r + 1 
V <- nu * diag(rep(1, r))

#�f���x�N�g���̎��O���z
gamma <- rep(0, ncol(z))
tau <- 100 * diag(ncol(z))
inv_tau <- solve(tau)

#�W���΍��̎��O���z
s0 <- 1.0
v0 <- 1.0

##�p�����[�^�̐^�l
#�K�w���f���̃p�����[�^
alpha1 <- alphat1
alpha2 <- alphat2
Cov1 <- Covt1
Cov2 <- Covt2

#�K�w���f���̊��Ғl��ݒ�
alpha_mu1 <- array(0, dim=c(hh, r, k1))
inv_Cov1 <- array(0, dim=c(r, r, k1))
for(j in 1:k1){
  alpha_mu1[, , j] <- u %*% alpha1[, , j]
  inv_Cov1[, , j] <- solve(Cov1[, , j])
}
alpha_mu2 <- array(0, dim=c(hh, r, k2))
inv_Cov2 <- array(0, dim=c(r, r, k2))
for(j in 1:k2){
  alpha_mu2[, , j] <- u %*% alpha2[, , j]
  inv_Cov2[, , j] <- solve(Cov2[, , j])
}

#�e���\���̃p�����[�^
theta1 <- thetat1
theta2 <- thetat2

#�f���x�N�g���̃p�����[�^
beta <- betat
sigma <- sigmat
beta_mu <- as.numeric(z %*% beta)


##�����l�̐ݒ�
#1�K�w�ڂ̃��f���̃p�����[�^
alpha1 <- array(0, dim=c(ncol(u), r, k1))
Cov1 <- inv_Cov1 <- array(0, dim=c(r, r, k1))
theta1 <- array(0, dim=c(k1, r, hh))
for(j in 1:k1){
  alpha1[, , j] <- matrix(0, nrow=ncol(u), ncol=r)
  Cov1[, , j] <- 0.005 * diag(r)
  inv_Cov1[, , j] <- solve(Cov1[, , j])
  theta1[j, , ] <- t(u %*% alpha1[, , j] + mvrnorm(hh, rep(0, r), Cov1[, , j]))
}

#2�K�w�ڂ̃��f���̃p�����[�^
alpha2 <- array(0, dim=c(ncol(u), r, k2))
Cov2 <- inv_Cov2 <- array(0, dim=c(r, r, k2))
theta2 <- array(0, dim=c(k2, r, hh))
for(j in 1:k2){
  alpha2[, , j] <- matrix(0, nrow=ncol(u), ncol=r)
  Cov2[, , j] <- 0.005 * diag(r)
  inv_Cov2[, , j] <- solve(Cov2[, , j])
  theta2[j, , ] <- t(u %*% alpha2[, , j] + mvrnorm(hh, rep(0, r), Cov2[, , j]))
}

#�f���x�N�g���̏����l  
beta <- solve(t(z) %*% z) %*% t(z) %*% y
sigma <- sd(y - z %*% beta)
beta_mu <- as.numeric(z %*% beta)


##�p�����[�^�̊i�[�p�z��
#�K�w���f���ƃe���\���̊i�[�p�z��
ALPHA1 <- array(0, dim=c(ncol(u), r, k1, R/keep))
ALPHA2 <- array(0, dim=c(ncol(u), r, k2, R/keep))
COV1 <- array(0, dim=c(r, r, k1, R/keep))
COV2 <- array(0, dim=c(r, r, k2, R/keep))
THETA1 <- array(0, dim=c(k1, r, hh, R/keep))
THETA2 <- array(0, dim=c(k2, r, hh, R/keep))

#�f���x�N�g���̊i�[�p�z��
BETA <- matrix(0, nrow=R/keep, ncol=ncol(z))
SIGMA <- rep(0, R/keep)
gc(); gc()

##�f�[�^�ƃC���f�b�N�X�̐ݒ�
#�f�[�^�̐ݒ�
id_vec <- rep(1:f, rep(k1, f))
y_vec <- y[id_vec]
inv_ZZV <- solve(t(z) %*% z + inv_tau)


#�C���f�b�N�X�̐ݒ�
index_vec <- matrix(1:(k1*max(w)), nrow=max(w), ncol=k1, byrow=T)

##�ΐ��ޓx�̊�l
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))
LLsq <- sum(dnorm(y, beta_mu, sd(y - beta_mu), log=TRUE))



####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�e���\���̃p�����[�^���T���v�����O
  #���f���덷���牞���ϐ���ݒ�
  er <- y - beta_mu
  tensor_mu_list <- list()
  
  for(i in 1:hh){
    #�f�[�^�𒊏o
    er_vec <- er[id_list[[i]]]
    x1 <- x1_vec[[i]]; x2 <- x2_vec[[i]]
    gamma2 <- theta2[x2, , i]
    index <- allocation_vec[[i]]
    index_w <- index_vec[1:w[i], ]
    
    ##1�K�w�ڂ̃e���\���̃p�����[�^���T���v�����O
    for(j in 1:k1){
      
      #�e���\���̐ݒ�
      gamma1 <- theta1[, , i]; gamma1[j, ] <- 1
      tensor_dt <- x1 * gamma1[index, ] * gamma2
      tensor_mu <- as.numeric((matrix(tensor_dt %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)[, -j]) %*% rep(1, k1-1))
      tensor_er <- er_vec - tensor_mu   #�����ϐ���ݒ�
      
      #���ϗʐ��K���z�̃p�����[�^��ݒ�
      X <- tensor_dt[index_w[, j], ]   #���͕ϐ���ݒ�
      Xy <- t(X) %*% tensor_er
      inv_XXV <- solve(t(X) %*% X + inv_Cov1[, , j])
      theta_mu <- inv_XXV %*% (Xy + inv_Cov1[, , j] %*% alpha_mu1[i, , j])   #���ϗʐ��K���z�̕��σx�N�g��
      
      #���ϗʐ��K���z����p�����[�^���T���v�����O
      theta1[j, , i] <- mvrnorm(1, theta_mu, sigma^2*inv_XXV)
    }
    
    ##2�K�w�ڂ̃e���\���̃p�����[�^���T���v�����O
    gamma2_dt <- theta2[x2_vec[[i]], , i]
    gamma1 <- theta1[index, , i]
    tensor_dt <- x1 * gamma1 
    
    for(j in 1:k2){
      
      #�C���f�b�N�X�ƒ萔�𒊏o
      index_x1 <- index_dt1[[i]][[j]]
      index_x2 <- index_data2[[i]][[j]]
      x_vec <- index_data1[[i]][[j]]
      dt <- dt_freq[[i]][[j]]
      
      #�����ϐ���ݒ�
      gamma2 <- gamma2_dt 
      gamma2[index_x2, ] <- 0
      tensor_mu <- as.numeric(matrix((tensor_dt * gamma2) %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T) %*% rep(1, k1))
      tensor_er <- er_vec - tensor_mu   #���f���덷
      
      #���͕ϐ���ݒ�
      X <- dt %*% (x_vec * theta1[index_x1, , i])
      
      #���ϗʐ��K���z�̃p�����[�^��ݒ�
      Xy <- t(X) %*% tensor_er
      inv_XXV <- solve(t(X) %*% X + inv_Cov2[, , j])
      theta_mu <- inv_XXV %*% (Xy + inv_Cov2[, , j] %*% alpha_mu2[i, , j])   #���ϗʐ��K���z�̕��σx�N�g��
      
      #���ϗʐ��K���z����p�����[�^���T���v�����O
      theta2[j, , i] <- as.numeric(mvrnorm(1, theta_mu, sigma^2*inv_XXV))
    }
    #�e���\���̊��Ғl
    mu <- matrix((x1 * theta1[index, , i] * theta2[x2, , i]) %*% rep(1, r), nrow=w[i], ncol=k1, byrow=T)
    tensor_mu_list[[i]] <- as.numeric(mu %*% rep(1, k1))
  }

  ##�f���x�N�g���̃p�����[�^���T���v�����O
  #�����ϐ���ݒ�
  tensor_mu <- unlist(tensor_mu_list)
  y_er <- y - tensor_mu
  
  #���ϗʐ��K���z�̃p�����[�^
  Zy <- t(z) %*% y_er
  par_mu <- as.numeric(inv_ZZV %*% Zy)
  
  #���ϗʐ��K���z����f���x�N�g�����T���v�����O
  beta <- mvrnorm(1, par_mu, sigma^2*inv_ZZV)
  beta_mu <- as.numeric(z %*% beta)
  
  ##���f���̕W���΍����T���v�����O
  #���f���̊��Ғl��ݒ�
  mu <- beta_mu + tensor_mu
  er <- y - mu
  
  #�K���}���z�̃p�����[�^
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- f + v0

  #�t�K���}���z����W���΍����T���v�����O
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##���ϗʉ�A���f���ŊK�w���f���̃p�����[�^���T���v�����O
  #1�K�w�ڂ̃e���\���̊K�w���f���̃p�����[�^���T���v�����O
  for(j in 1:k1){
    out <- rmultireg(t(theta1[j, , ]), u, Deltabar, ADelta, nu, V)
    alpha1[, , j] <- out$B
    alpha1[, , j] <- alphat1[, , j]
    Cov1[, , j] <- out$Sigma  
    Cov1[, , j] <- Covt1[, , j]
    alpha_mu1[, , j] <- u %*% alpha1[, , j]
    inv_Cov1[, , j] <- solve(Cov1[, , j])
  }
  
  #2�K�w�ڂ̃e���\���̊K�w���f���̃p�����[�^���T���v�����O
  for(j in 1:k2){
    out <- rmultireg(t(theta2[j, , ]), u, Deltabar, ADelta, nu, V)
    alpha2[, , j] <- out$B
    alpha2[, , j] <- alphat2[, , j]
    Cov2[, , j] <- out$Sigma  
    Cov2[, , j] <- Covt2[, , j]
    alpha_mu2[, , j] <- u %*% alpha2[, , j]
    inv_Cov2[, , j] <- solve(Cov2[, , j])
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�p�����[�^�̊i�[
  if(rp%%keep==0){
    mkeep <- rp/keep
    mkeep <- 1
    BETA[mkeep, ] <- beta
    SIGMA[mkeep] <- sigma
    THETA1[, , , mkeep] <- theta1
    THETA2[, , , mkeep] <- theta2
    ALPHA1[, , , mkeep] <- alpha1
    ALPHA2[, , , mkeep] <- alpha2
    COV1[, , , mkeep] <- Cov1
    COV2[, , , mkeep] <- Cov2
  }
  
  #�T���v�����O���ʂ̕\��
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))   #��ă��f���̑ΐ��ޓx
    
    #�T���v�����O���ʂ�\��
    print(rp)
    print(c(LL, LLsq, LLst))
    print(round(rbind(beta=as.numeric(beta), betat), 3))
    print(sigma)
  }
}

