#####Truncated Poisson based Hierarchical Matrix Factorization#####
options(warn=0)
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##データの設定
k <- 7   #基底数
hh <- 10000   #ユーザー数
item <- 3000   #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)
vec <- rep(1, k)

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}


##素性ベクトルを生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合


##階層モデルの説明変数を設定
#ユーザーの説明変数
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合

#アイテムの説明変数
k1 <- 2; k2 <- 3; k3 <- 4
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #データを結合

##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

#多項分布からアイテムを生成
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

#個別に和を取るためのスパース行列
user_vec <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_vec <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


####応答変数を生成####
rp <- 0
repeat {
  
  ##素性ベクトルのパラメータ
  #階層モデルのパラメータ
  Cov <- Covt <- diag(c(runif(1, 0.01, 0.1), runif(ncol(x)-1, 0.005, 0.05)), ncol(x))
  alpha <- alphat <- matrix(rnorm(ncol(x)*ncol(u), 0, 0.4), nrow=ncol(u), ncol=ncol(x))
  alpha[1, 1] <- 0.3
  
  #多変量正規分布からパラメータを生成
  beta <- betat <- u %*% alpha + mvrnorm(hh, rep(0, ncol(x)), Cov)
  
  ##行列分解のパラメータを生成
  #階層モデルの分散パラメータ
  Cov_u <- Cov_ut <- diag(runif(k, 0.005, 0.05), k)   #ユーザー-アイテムの階層モデルの分散
  Cov_v <- Cov_vt <- diag(runif(k, 0.05, 0.05), k)   #アイテムの階層モデルの分散
  
  #階層モデルの回帰係数を設定
  alpha_u <- alpha_ut <- matrix(rnorm(k*ncol(u), 0, 0.2), nrow=ncol(u), ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.2), nrow=ncol(v), ncol=k)
  
  #多変量正規分布からパラメータを生成
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k), Cov_u)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  ##切断ポアソン分布から購買ベクトルを生成
  #期待値の設定
  uv <- uvt <- rowSums(theta_u[user_id, ] * theta_v[item_id, ])
  mu <- as.numeric((x * beta[user_id, ]) %*% rep(1, ncol(x)))
  lambda <- exp(mu + uv)
  
  #購買ベクトルを生成
  y <- extraDistr::rtpois(hhpt, lambda, a=0, b=Inf)
  
  print(c(mean(y), max(y)))
  if(mean(y) > 2.5 & max(y) < 200) break   #break条件
}

#購買数を確認
sum(y); max(y)
summary(y)
hist(y, main="購買頻度の真値の分布", xlab="購買頻度", col="grey", breaks=25)


####ハミルトニアンモンテカルロ法でTruncated Poisson based HMFを推定####
##対数事後分布を計算する関数
#対数尤度を計算する関数
loglike <- function(mu, uv, y, x, y_lfactorial, dt){
  
  #切断ポアソン回帰モデルの対数尤度
  lambda <- exp(mu + uv)   #期待値
  Lho <- y*log(lambda) - lambda - log(1-exp(-lambda)) - y_lfactorial   #対数尤度関数
  
  #対数事後分布
  LLi <- as.numeric(dt %*% Lho)
  return(LLi)
}

#多変量正規分布の対数事前分布
dmv <- function(er, inv_Cov, k){
  Li <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k))
  return(Li)
}

##素性ベクトルのサンプリングに必要な関数
#素性ベクトルの対数事後分布の微分関数
dloglike <- function(beta, uv, y, x, alpha_mu, inv_Cov, const, dt){ 
  
  #期待値の設定
  lambda <- as.numeric(exp(as.numeric((x * beta[user_id, ]) %*% vec_x) + uv))
  lambda_exp <- exp(-lambda)
  lambda_x <- x*lambda
  
  #微分関数の設定
  er <- beta - alpha_mu
  dltpois <- const - lambda_x - lambda_exp * (lambda_x) / (1-lambda_exp)   #切断ポアソン微分関数の設定
  dmvn <- -t(inv_Cov %*% t(er))   #多変量正規分布の対数事前分布の微分関数
  
  #対数事後分布の微分関数の和
  LLd <- -((dt %*% dltpois) + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, uv, y, x, alpha_mu, inv_Cov, const, user_vec) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, uv, y, x, alpha_mu, inv_Cov, const, user_vec) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##ユーザー特徴行列のサンプリングに必要な関数
#ユーザー特徴行列の対数事後分布の微分関数
dloglike_u <- function(theta_u, theta_vec, mu, y, x, mu_u, inv_Cov, const, dt, user_id){ 
  
  #期待値の設定
  uv <- as.numeric((theta_u[user_id, ] * theta_vec) %*% vec)
  lambda <- as.numeric(exp(mu + uv))
  lambda_exp <- exp(-lambda)
  lambda_x <- theta_vec*lambda
  
  #微分関数の設定
  er <- theta_u - mu_u
  dltpois <- const - lambda_x - lambda_exp * (lambda_x) / (1-lambda_exp)   #切断ポアソン微分関数の設定
  dmvn <- -t(inv_Cov %*% t(er))   #多変量正規分布の対数事前分布の微分関数
  
  #対数事後分布の微分関数の和
  LLd <- -((dt %*% dltpois) + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog_u <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_vec, mu, y, x, mu_u, inv_Cov_u, const_u, user_vec, user_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_vec, mu, y, x, mu_u, inv_Cov_u, const_u, user_vec, user_id) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##アイテム特徴行列のサンプリングに必要な関数
#アイテム特徴行列の対数事後分布の微分関数
dloglike_v <- function(theta_v, theta_vec, mu, y, x, mu_v, inv_Cov, const, dt, item_id){ 
  
  #期待値の設定
  uv <- as.numeric((theta_v[item_id, ] * theta_vec) %*% vec)
  lambda <- as.numeric(exp(mu + uv))
  lambda_exp <- exp(-lambda)
  lambda_x <- theta_vec*lambda
  
  #微分関数の設定
  er <- theta_v - mu_v
  dltpois <- const - lambda_x - lambda_exp * (lambda_x) / (1-lambda_exp)   #切断ポアソン微分関数の設定
  dmvn <- -t(inv_Cov %*% t(er))   #多変量正規分布の対数事前分布の微分関数
  
  #対数事後分布の微分関数の和
  LLd <- -((dt %*% dltpois) + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_vec, mu, y, x, mu_v, inv_Cov_v, const_v, item_vec, item_id) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_vec, mu, y, x, mu_v, inv_Cov_v, const_v, item_vec, item_id) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}


##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 100/keep
disp <- 10
e <- 0.025
L <- 3

##事前分布の設定
#素性ベクトルの階層モデルの事前分布
Deltabar <- matrix(rep(0, ncol(u)*ncol(x)), nrow=ncol(u), ncol=ncol(x))   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, ncol(u)))   #階層モデルの回帰係数の事前分布の分散
nu <- 1   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, ncol(x))) #逆ウィシャート分布のパラメータ

#ユーザーの階層モデルの事前分布
Deltabar_u <- matrix(rep(0, ncol(u)*k), nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta_u <- 0.01 * diag(rep(1, ncol(u)))   #階層モデルの回帰係数の事前分布の分散
nu_u <- 1   #逆ウィシャート分布の自由度
V_u <- nu_u * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#アイテムの階層モデルの事前分布
Deltabar_v <- matrix(rep(0, ncol(v)*k), nrow=ncol(v), ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta_v <- 0.01 * diag(rep(1, ncol(v)))   #階層モデルの回帰係数の事前分布の分散
nu_v <- 1   #逆ウィシャート分布の自由度
V_v <- nu_v * diag(rep(1, k)) #逆ウィシャート分布のパラメータ


##パラメータの真値
#素性ベクトルのパラメータ
beta <- betat
mu <- rowSums(x * beta[user_id, ])

#階層モデルのパラメータ
alpha <- alphat; Cov <- Covt; inv_Cov <- solve(Cov); alpha_mu <- u %*% alpha
alpha_u <- alpha_ut; Cov_u <- Cov_ut; inv_Cov_u <- solve(Cov_u); mu_u <- u %*% alpha_u
alpha_v <- alpha_vt; Cov_v <- Cov_vt; inv_Cov_v <- solve(Cov_v); mu_v <- v %*% alpha_v

#行列分解のパラメータ
theta_u <- theta_ut
theta_v <- theta_vt
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec)

##初期値の設定
#素性ベクトルの階層モデルのパラメータ
Cov <- 0.01 * diag(ncol(x)); inv_Cov <- solve(Cov)   
alpha <- matrix(0, nrow=ncol(u), ncol=ncol(x)); alpha_mu <- u %*% alpha

#素性ベクトルのパラメータ
res <- glm(y ~ x[, -1], family=poisson)
beta <- mvrnorm(hh, rep(0, ncol(x)), Cov)

##階層モデルのパラメータを生成
#階層モデルの分散パラメータ
Cov_u <- 0.01 * diag(k); inv_Cov_u <- solve(Cov_u)   #ユーザー-アイテムの階層モデルの分散
Cov_v <- 0.01 * diag(k); inv_Cov_v <- solve(Cov_v)   #アイテムの階層モデルの分散

#階層モデルの回帰係数を設定
alpha_u <- matrix(0, nrow=ncol(u), ncol=k); mu_u <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=k); mu_v <- v %*% alpha_v

#行列分解のパラメータを生成
theta_u <- mvrnorm(hh, rep(0, k), Cov_u)
theta_v <- mvrnorm(item, rep(0, k), Cov_v)
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec)

##定数とデータの設定
#定数の設定
y_lfactorial <- lfactorial(y)
const <- y * x

#データの設定
x_col <- ncol(x)
u_col <- ncol(u)
vec_x <- rep(1, x_col)
vec <- rep(1, k)


##サンプリング結果のパラメータの保存用配列
#モデルパラメータの格納用配列
BETA <- array(0, dim=c(hh, x_col, R/keep))
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))

#階層モデルの格納用配列
ALPHA <- array(0, dim=c(ncol(u), x_col, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k, R/keep))
COV <- array(0, dim=c(x_col, x_col, R/keep))
COV_U <- COV_V <- array(0, dim=c(k, k, R/keep))

##対数尤度の基準値
#1パラメータモデルの対数尤度
LLst <- sum(extraDistr::dtpois(y, mean(y), a=0, b=Inf, log=TRUE))

#真値での対数尤度
beta_mu <- as.numeric((x * betat[user_id, ]) %*% vec_x)
uv <- as.numeric((theta_ut[user_id, ] * theta_vt[item_id, ]) %*% vec)
lambda <- exp(beta_mu + uv)
LLbest <- sum(extraDistr::dtpois(y, lambda, a=0, b=Inf, log=TRUE))


####HMCでパラメータをサンプリング
for(rp in 1:R){
  
  ##素性ベクトルのパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(hh, rep(0, x_col), diag(x_col))   #標準多変量正規分布からパラメータを生成
  betad <- beta
  mu_old <- mu
  
  #リープフロッグ法による1ステップ移動
  res <- leapfrog(rold, betad, dloglike, e, L)
  rnew <- res$r
  betan <- res$z
  mu_new <- as.numeric((x * betan[user_id, ]) %*% vec_x)
  
  #移動前と移動後のハミルトニアン
  er_new <- betan - alpha_mu
  er_old <- betad - alpha_mu
  Hnew <- -(loglike(mu_new, uv, y, x, y_lfactorial, user_vec) + dmv(er_new, inv_Cov, x_col)) + as.numeric(rnew^2 %*% vec_x)/2
  Hold <- -(loglike(mu_old, uv, y, x, y_lfactorial, user_vec) + dmv(er_old, inv_Cov, x_col)) + as.numeric(rold^2 %*% vec_x)/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(hh) #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma_beta <- mean(gamma, na.rm=TRUE)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- gamma > rand
  beta <- as.matrix(flag*betan + (1-flag)*betad)
  na_beta <- is.na(beta); if(sum(na_beta) > 0){beta[na_beta] <- betad}
  mu <- as.numeric((x * beta[user_id, ]) %*% vec_x)
  
  ##ユーザーのアイテムに対する行列分解のパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(hh, rep(0, k), diag(k))   #標準多変量正規分布からパラメータを生成
  thetad <- theta_u
  uv_old <- uv
  
  #リープフロッグ法による1ステップ移動
  theta_vec <- theta_v[item_id, ]
  const_u <- y * theta_vec
  res <- leapfrog_u(rold, thetad, dloglike_u, e, L)
  rnew <- res$r
  thetan <- res$z
  uv_new <- as.numeric((thetan[user_id, ] * theta_vec) %*% vec)
  
  #移動前と移動後のハミルトニアン
  er_new <- thetan - mu_u
  er_old <- thetad - mu_u
  Hnew <- -(loglike(mu, uv_new, y, x, y_lfactorial, user_vec) + dmv(er_new, inv_Cov_u, k)) + as.numeric(rnew^2 %*% vec)/2
  Hold <- -(loglike(mu, uv_old, y, x, y_lfactorial, user_vec) + dmv(er_old, inv_Cov_u, k)) + as.numeric(rold^2 %*% vec)/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(hh)   #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma_u <- mean(gamma)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  theta_u <- as.matrix(flag*thetan + (1-flag)*thetad)
  na_theta_u <- is.na(theta_u); if(sum(na_theta_u) > 0){theta_u[na_theta_u] <- thetad}
  uv <- as.numeric((theta_u[user_id, ] * theta_vec) %*% vec)  
  
  
  ##アイテムの行列分解のパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(item, rep(0, k), diag(k))   #標準多変量正規分布からパラメータを生成
  thetad <- theta_v
  uv_old <- uv
  
  #リープフロッグ法による1ステップ移動
  theta_vec <- theta_u[user_id, ]
  const_v <- y * theta_vec
  res <- leapfrog_v(rold, thetad, dloglike_v, e, L)
  rnew <- res$r
  thetan <- res$z
  uv_new <- as.numeric((theta_vec * thetan[item_id, ]) %*% vec)
  
  #移動前と移動後のハミルトニアン
  er_new <- thetan - mu_v
  er_old <- thetad - mu_v
  Hnew <- -(loglike(mu, uv_new, y, x, y_lfactorial, item_vec) + dmv(er_new, inv_Cov_v, k)) + as.numeric(rnew^2 %*% vec)/2
  Hold <- -(loglike(mu, uv_old, y, x, y_lfactorial, item_vec) + dmv(er_old, inv_Cov_v, k)) + as.numeric(rold^2 %*% vec)/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(item)   #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma_v <- mean(gamma)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  theta_v <- as.matrix(flag*thetan + (1-flag)*thetad)
  na_theta_v <- is.na(theta_v); if(sum(na_theta_v) > 0){theta_v[na_theta_v] <- thetad}
  uv <- as.numeric((theta_vec * theta_v[item_id, ]) %*% vec)
  
  ##階層モデルのパラメータをサンプリング
  #素性ベクトルのパラメータをサンプリング
  out <- rmultireg(beta, u, Deltabar, ADelta, nu, V)
  alpha <- out$B
  Cov <- out$Sigma
  alpha_mu <- u %*% alpha
  inv_Cov <- solve(Cov)
  
  #ユーザーの行列分解のパラメータをサンプリング
  out_u <- rmultireg(theta_u, u, Deltabar_u, ADelta_u, nu_u, V_u)
  alpha_u <- out_u$B
  Cov_u <- out_u$Sigma
  mu_u <- u %*% alpha_u
  inv_Cov_u <- solve(Cov_u)
  
  #アイテムの行列分解のパラメータをサンプリング
  out_v <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu_v, V_v)
  alpha_v <- out_v$B
  Cov_v <- out_v$Sigma
  mu_v <- v %*% alpha_v
  inv_Cov_v <- solve(Cov_v)
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    #モデルのパラメータを格納
    BETA[, , mkeep] <- beta
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    
    #階層モデルのパラメータ
    ALPHA[, , mkeep] <- alpha
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v 
    COV[, , mkeep] <- Cov
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    LL <- sum(loglike(mu, uv, y, x, y_lfactorial, user_vec))   #対数尤度
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(c(gamma_beta, gamma_u, gamma_v), 3))
    print(c(sum(rowSums(na_beta) > 0), sum(rowSums(na_theta_u) > 0), sum(rowSums(na_theta_v) > 0)))
    print(round(rbind(diag(Cov), diag(Covt)), 3))
  }
}


