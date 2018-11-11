#####Bayesian Hierarchical Non Negative Matrix Factorization#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(NMF)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(78594)

####データの発生####
#データの設定
k <- 10   #基底数
hh <- 5000   #ユーザー数
item <- 2000  #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #接触数
hhpt <- sum(pt)   #総接触数
vec_k <- rep(1, k)

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}


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
u_col <- ncol(u)

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
v_col <- ncol(v)


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

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1

  ##NMFのパラメータを生成
  #ガンマ分布の尺度パラメータを設定
  alpha_u <- alpha_ut <- matrix(rnorm(k*ncol(u), 0, 0.3), nrow=ncol(u), ncol=k)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.3), nrow=ncol(v), ncol=k)
  lambda_u <- exp(u %*% alpha_u)
  lambda_v <- exp(v %*% alpha_v)
  
  #ガンマ分布の形状パラメータ
  beta_u <- beta_ut <- rnorm(1, 1, 0.5)
  beta_v <- beta_vt <- rnorm(1, 1, 0.5)
  
  #ガンマ分布から行列分解のパラメータを生成
  theta_u <- theta_ut <- matrix(rgamma(hh*k, as.numeric(lambda_u), beta_u), nrow=hh, ncol=k)
  theta_v <- theta_vt <- matrix(rgamma(item*k, as.numeric(lambda_v), beta_v), nrow=item, ncol=k)
  
  ##ポアソン分布から応答変数を生成
  WH <- (theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k   #期待値
  y　<- rpois(hhpt, WH)   
  
  #break条件
  print(rp)
  print(max(y))
  if(max(y) < 75 & max(y) > 25){
    break
  }
}

#生成したデータのヒストグラム
hist(y, main="購買頻度の分布", xlab="購買頻度", col="grey", breaks=50)


####マルコフ連鎖モンテカルロ法でBayesian Hierarchical NMFを推定####
##対数事後分布を計算する関数
loglike <- function(beta, alpha, inv_tau, y, y_log, x){
  
  #ガンマ回帰モデルの対数尤度
  lambda <- as.numeric(exp(x %*% beta))   #期待値
  Lho <- sum(alpha * as.numeric(-y/lambda - x %*% beta) + alpha*log(alpha) - lgamma(alpha) + (alpha-1)*y_log)   #対数尤度関数
  
  #多変量正規分布の対数事前分布
  log_mvn <- -1/2 * as.numeric(t(beta) %*% inv_tau %*% beta)
  
  #対数事後分布
  LL <- Lho + log_mvn
  return(list(LL=LL, Lho=Lho))
}

##HMCで尺度パラメータをサンプリングするための関数
#ガンマ回帰の対数事後分布の微分関数
dloglike <- function(beta, alpha, inv_tau, y, y_log, x){ 
  
  #期待値の設定
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #期待値
  
  #微分関数の設定
  dlgamma <- colSums((y-lambda) / (lambda^2/alpha) * lambda * x)   #尺度パラメータの勾配ベクトル
  dmvn <- as.numeric(-inv_tau %*% beta)
  
  #対数事後分布の微分関数の和
  LLd <- -(dlgamma + dmvn)
  return(LLd)
}

#リープフロッグ法を解く関数
leapfrog_u <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta_u, inv_tau_u, d, d_log, u) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta_u, inv_tau_u, d, d_log, u) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, beta_v, inv_tau_v, d, d_log, v) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, beta_v, inv_tau_v, d, d_log, v) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

#形状パラメータの対数事後分布の微分関数
dloglike_alpha <- function(alpha, beta, y, y_log, x, n, k){ 
  #期待値の設定
  mu <- as.numeric(x %*% beta)
  lambda <- exp(mu)   #期待値
  
  #勾配ベクトルの計算
  dlgamma <- (n*k)*(log(alpha) - digamma(alpha)) + sum(1 - y/lambda + log(y/lambda))   #形状パラメータの勾配ベクトル
  return(dlgamma)
}


##アルゴリズムの設定
R <- 5000
keep <- 4
burnin <- 1000/keep
iter <- 0
disp <- 10
e <- 0.001
L <- 3

##事前分布の設定
gamma_u <- rep(0, u_col); gamma_v <- rep(0, v_col)
inv_tau_u <- solve(100 * diag(u_col)); inv_tau_v <- solve(100 * diag(v_col))
inv_tau <- solve(100 * diag(k))

##パラメータの真値
#ガンマ分布の尺度パラメータ
alpha_u <- alpha_ut
alpha_v <- alpha_vt
lambda_u <- exp(u %*% alpha_u)
lambda_v <- exp(v %*% alpha_v)

#ガンマ分布の形状パラメータ
beta_u <- beta_ut 
beta_v <- beta_vt

#行列分解のパラメータ
theta_u <- theta_ut
theta_v <- theta_vt
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #期待値

##初期値の設定
#ガンマ分布の尺度パラメータ
alpha_u <- matrix(rnorm(k*ncol(u), 0, 0.1), nrow=ncol(u), ncol=k)
alpha_v <- matrix(rnorm(k*ncol(v), 0, 0.1), nrow=ncol(v), ncol=k)
lambda_u <- exp(u %*% alpha_u)
lambda_v <- exp(v %*% alpha_v)

#ガンマ分布の形状パラメータ
beta_u <- 1.0
beta_v <- 1.0

#行列分解のパラメータ
theta_u <- matrix(rgamma(hh*k, as.numeric(lambda_u), beta_u), nrow=hh, ncol=k)
theta_v <- matrix(rgamma(item*k, as.numeric(lambda_v), beta_v), nrow=item, ncol=k)
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #期待値


##サンプリング結果の保存用配列
gamma_u1 <- gamma_v1 <- rep(k)
gamma_u2 <- gamma_v2 <- rep(k)
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
ALPHA_U <- array(0, dim=c(u_col, k, R/keep))
ALPHA_V <- array(0, dim=c(v_col, k, R/keep))
BETA_U <- rep(0, R/keep)
BETA_V <- rep(0, R/keep)

##ユーザーおよびアイテムのインデックスを作成
#個別に和を取るためのスパース行列
user_dt <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))

##対数尤度の基準値
LLst <- sum(dpois(y, mean(y), log=TRUE))
LLbest <- sum(dpois(y, as.numeric((theta_ut[user_id, ] * theta_vt[item_id, ]) %*% vec_k), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){

  ##ユーザー特徴行列をサンプリング
  #補助変数lambdaを更新
  theta_vec2 <- theta_v[item_id, ]
  lambda <- (theta_u[user_id, ] * theta_vec2) / WH

  #ユーザーごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y   #要素ごとの期待値
  W1 <- as.matrix(user_dt %*% lambda_y + lambda_u)
  W2 <- as.matrix(user_dt %*% theta_vec2 + beta_u)
  
  #ガンマ分布よりパラメータをサンプリング
  theta_u <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  #theta_u <- theta_u / matrix(colSums(theta_u), nrow=hh, ncol=k, byrow=T) * hh/(k/2)   #各列ベクトルを正規化
  
  
  ##アイテム特徴行列をサンプリング
  #補助変数lambdaを更新
  theta_vec1 <- theta_u[user_id, ]
  WH <- as.numeric((theta_vec1 * theta_vec2) %*% vec_k)
  lambda <- (theta_vec1 * theta_vec2) / WH

  #アイテムごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y   #要素ごとの期待値
  H1 <- as.matrix(item_dt %*% lambda_y + lambda_v)
  H2 <- as.matrix(item_dt %*% theta_vec1 + beta_v)
  
  #ガンマ分布よりパラメータをサンプリング
  theta_v <- matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k)
  WH <- as.numeric((theta_vec1 * theta_v[item_id, ]) %*% vec_k)
  
  
  ##階層モデルのパラメータをサンプリング
  for(j in 1:k){
    
    ##ユーザー特徴行列の尺度パラメータをサンプリング
    #HMCの新しいパラメータを生成
    rold <- as.numeric(mvrnorm(1, rep(0, u_col), diag(u_col)))   #標準多変量正規分布からパラメータを生成
    alphad <- alpha_u[, j]
    
    #リープフロッグ法による1ステップ移動
    d <- theta_u[, j]; d_log <- log(d)
    res <- leapfrog_u(rold, alphad, dloglike, e, L)
    rnew <- res$r
    alphan <- res$z
    
    #移動前と移動後のハミルトニアン
    Hnew <- -loglike(alphan, beta_u, inv_tau_u, d, d_log, u)$LL + as.numeric(rnew^2 %*% rep(1, u_col))/2
    Hold <- -loglike(alphad, beta_u, inv_tau_u, d, d_log, u)$LL + as.numeric(rold^2 %*% rep(1, u_col))/2
    
    #パラメータの採択を決定
    rand <- runif(1)   #一様分布から乱数を発生
    gamma <- min(c(1, exp(Hold - Hnew)))   #採択率を決定
    gamma_u1[j] <- gamma
    
    #alphaの値に基づき新しいbetaを採択するかどうかを決定
    flag <- as.numeric(gamma > rand)
    alpha_u[, j] <- flag*alphan + (1-flag)*alphad
    
    
    ##アイテム特徴行列の尺度パラメータをサンプリング
    #HMCの新しいパラメータを生成
    rold <- as.numeric(mvrnorm(1, rep(0, v_col), diag(v_col)))   #標準多変量正規分布からパラメータを生成
    alphad <- alpha_v[, j]
    
    #リープフロッグ法による1ステップ移動
    d <- theta_v[, j]; d_log <- log(d)
    res <- leapfrog_v(rold, alphad, dloglike, e, L)
    rnew <- res$r
    alphan <- res$z
  
    #移動前と移動後のハミルトニアン
    Hnew <- -loglike(alphan, beta_v, inv_tau_v, d, d_log, v)$LL + as.numeric(rnew^2 %*% rep(1, v_col))/2
    Hold <- -loglike(alphad, beta_v, inv_tau_v, d, d_log, v)$LL + as.numeric(rold^2 %*% rep(1, v_col))/2
    
    #パラメータの採択を決定
    rand <- runif(1)   #一様分布から乱数を発生
    gamma <- min(c(1, exp(Hold - Hnew)))   #採択率を決定
    gamma_v1[j] <- gamma
    
    #alphaの値に基づき新しいbetaを採択するかどうかを決定
    flag <- as.numeric(gamma > rand)
    alpha_v[, j] <- flag*alphan + (1-flag)*alphad
  }
  
  ##ユーザー特徴行列の形状パラメータをサンプリング
  #MH法の新しいパラメータを生成
  d <- as.numeric(theta_u); d_log <- log(d)
  s <- sign(dloglike_alpha(beta_u, alpha_u, d, d_log, u, hh, k))   #勾配の方向
  betad <- beta_u
  betan <- betad + s*abs(rnorm(1, 0, 0.01))  
  
  #独立MH法の対数尤度
  lognew <- loglike(alpha_u, betan, inv_tau, d, d_log, u)$Lho
  logold <- loglike(alpha_u, betad, inv_tau, d, d_log, u)$Lho 
  
  #パラメータの採択を決定
  rand <- runif(1)   #一様分布から乱数を発生
  gamma <- min(c(1, exp(lognew - logold)))   #採択率を決定
  gamma_u2 <- gamma

  #gammaの値に基づき新しいalphaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  beta_u <- flag*betan + (1-flag)*betad
  
  ##アイテム特徴行列の形状パラメータをサンプリング
  #MH法の新しいパラメータを生成
  d <- as.numeric(theta_v); d_log <- log(d)
  s <- sign(dloglike_alpha(beta_v, alpha_v, d, d_log, v, item, k))   #勾配の方向
  betad <- beta_v
  betan <- betad + s*abs(rnorm(1, 0, 0.01))  
  
  #独立MH法の対数尤度
  lognew <- loglike(alpha_v, betan, inv_tau_v, d, d_log, v)$Lho
  logold <- loglike(alpha_v, betad, inv_tau_v, d, d_log, v)$Lho 
  
  #パラメータの採択を決定
  rand <- runif(1)   #一様分布から乱数を発生
  gamma <- min(c(1, exp(lognew - logold)))   #採択率を決定
  gamma_v2 <- gamma
  
  #gammaの値に基づき新しいalphaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  beta_v <- flag*betan + (1-flag)*betad
  
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    BETA_U[mkeep] <- beta_u
    BETA_V[mkeep] <- beta_v
  }
  
  if(rp%%disp==0){
    #対数尤度の更新
    LL <- sum(dpois(y, WH, log=TRUE))
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(rbind(gamma_u1, gamma_v1), 3))
    print(round(c(gamma_u2, gamma_v2), 3))
    print(round(c(beta_u, beta_v, beta_ut, beta_vt), 3))
  }
}

####サンプリング結果の要約と適合度####
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W %*% H))[-index_z1], log=TRUE))
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W0 %*% H0))[-index_z1], log=TRUE))



