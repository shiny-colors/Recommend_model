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

####データの発生####
#set.seed(34027)
##データの設定
r <- 7   #基底数
select <- 7   #選択肢数
hh <- 5000   #消費者数
item <- 1500   #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)   #全サンプル数
allocation_index <- matrix(1:(r*(select-1)), nrow=select-1, ncol=r, byrow=T)   #基底の割当

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合


##説明変数の発生
#条件付き説明変数
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
  Data1[, j] <- as.numeric(t(X1[, j, ]))   #説明変数をベクトル変換
}
user_id_vec <- rep(user_id, rep(select, hhpt))


#無条件説明変数
k21 <- 2; k22 <- 3; k23 <- 4
x1 <- matrix(runif(hhpt*k21, 0, 1), nrow=hhpt, ncol=k21)
x2 <- matrix(0, nrow=hhpt, ncol=k22)
for(j in 1:k22){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k23, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
data <- cbind(1, x1, x2, x3)   #データの結合

#無条件説明変数をベクトル変換
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
sparse_data <- as(Data, "CsparseMatrix")   #スパース行列に変換


##ユーザーの階層モデルの説明変数の発生
#ユーザーの説明変数
k1 <- 1; k2 <- 3; k3 <- 5
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合

#アイテムの説明変数
k1 <- 2; k2 <- 2; k3 <- 4
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
item_id <- as.numeric(rmnom(hhpt, 1, phi[z[user_id], ]) %*% 1:item)
freq <- plyr::count(item_id); freq$x <- as.character(freq$x)
hist(freq$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


####応答変数の発生####
rp <- 0
repeat { 
    
  ##パラメータを生成
  #ロジットのパラメータを生成
  k1 <- ncol(Data)   #パラメータ数
  alpha <- alphat <- matrix(rnorm(k1*ncol(u), 0, 0.5), nrow=ncol(u), ncol=k1)   #階層モデルの回帰係数
  Cov <- Covt <- runif(k1, 0.05, 0.2) * diag(k1)   #階層モデルの分散
  theta <- thetat <- u %*% alpha + mvrnorm(hh, rep(0, k1), Cov)   #ユーザー別回帰係数
  
  #ユーザーの行列分解のパラメータ
  k2 <- r*(select-1)
  alpha_u <- alpha_ut <- matrix(rnorm(ncol(u)*k2, 0, 0.5), nrow=ncol(u), ncol=k2)
  Cov_u <- Cov_ut <- runif(k2, 0.01, 0.15) * diag(k2)
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k2), Cov_u)
  
  #アイテムの行列分解のパラメータ
  alpha_v <- alpha_vt <- matrix(rnorm(ncol(v)*k2, 0, 0.5), nrow=ncol(v), ncol=k2)
  Cov_v <- Cov_vt <- runif(k2, 0.01, 0.15) * diag(k2)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k2), Cov_v)
  
  ##応答変数を生成
  #ロジットを設定
  mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #素性ベクトルの期待値
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #行列分解の期待値
  } 
  logit_exp <- exp(mu + uv)   #ロジットの期待値の指数
  
  #選択結果を生成
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))   #選択確率
  y <- rmnom(hhpt, 1, prob)   #多項分布から応答変数を生成
  
  #停止条件
  print(round(c(mean(rowMaxs(prob)), min(colSums(y))), 3))
  if(mean(rowMaxs(prob)) > 0.725 & min(colSums(y)) > hhpt/(select*5)) break
}
colSums(y)
hist(prob[y==1], main="選択結果の選択確率", xlab="選択確率", breaks=25, col="grey")
hist(rowMaxs(prob), main="もっとも高い選択確率", xlab="選択確率", breaks=25, col="grey")

#オブジェクトの消去
rm(data_list); rm(Data1); rm(Data2)
gc(); gc()


####ハミルトニアンモンテカルロ法でパラメータをサンプリング####
##対数事後分布を計算する関数
#対数尤度を計算する関数
loglike <- function(mu, uv, y, hh, select, dt){
  
  #ロジットモデルの対数尤度
  logit_exp <- exp(mu + uv)   #ロジットの期待値の指数
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))   #選択確率
  LLi_logit <- log(as.numeric((y * prob) %*% rep(1, select)))
  
  #ユーザーごとの対数尤度
  LLi <- as.numeric(dt %*% LLi_logit)
  return(LLi)
}

#多変量正規分布の対数事前分布
dmv <- function(er, inv_Cov, k){
  Li <- -1/2 * as.numeric((er %*% inv_Cov * er) %*% rep(1, k))
  return(Li)
}

##素性ベクトルのサンプリングに必要な関数
#素性ベクトルの対数事後分布の微分関数
dloglike <- function(theta, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt){
  
  #応答確率の設定
  logit_exp <- exp(matrix((Data * theta[user_id_vec, ]) %*% rep(1, k1), nrow=hhpt, ncol=select, byrow=T) + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  prob_vec <- as.numeric(t(prob))
  
  #微分関数の設定
  er <- theta - alpha_mu
  dlogit <- (y_vec - prob_vec) * Data   #ロジットモデルの対数尤度の微分関数
  dmvn <- -t(inv_Cov %*% t(er))   #多変量正規分布の対数事前分布の微分関数

  #対数事後分布の微分関数の和
  LLd <- -(user_dt %*% dlogit + dmvn)
  return(LLd)
}

#素性ベクトルのリープフロッグ法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, uv, Data, y_vec, alpha_mu, inv_Cov, hh, hhpt, select, k1, user_id_vec, user_dt) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##ユーザーの特徴行列のサンプリングに必要な関数
#ユーザーの行列分解の対数事後分布の微分関数
dloglike_u <- function(theta_u, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, r, k2,
                       user_id, user_vec, allocation_index){
  
  #応答確率の設定
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_u[user_id, ] * theta_dv
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  logit_exp <- exp(mu + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  
  #微分関数の設定
  er <- theta_u - u_mu
  dlogit <- matrix(0, nrow=hhpt, ncol=k2)
  for(j in 1:(select-1)){
    dlogit[, allocation_index[j, ]] <- (y[, j] - prob[, j]) * theta_dv[, allocation_index[j, ]]   #ロジットの微分関数
  }
  dmvn <- -t(inv_Cov_u %*% t(er))   #多変量正規分布の微分関数

  #対数事後分布の微分関数の和
  LLd <- -(as.matrix(user_vec %*% dlogit) + dmvn)
  return(LLd)
}

#ユーザーの行列分解のリープフロッグ法を解く関数
leapfrog_u <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, g, k2, user_id, user_vec, allocation_index) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_dv, mu, y, u_mu, inv_Cov_u, hh, hhpt, select, g, k2, user_id, user_vec, allocation_index) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##アイテムの特徴行列のサンプリングに必要な関数
#アイテムの行列分解の対数事後分布の微分関数
dloglike_v <- function(theta_v, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, r, k2,
                       item_id, item_vec, allocation_index){
  
  #応答確率の設定
  uv <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_du * theta_v[item_id, ]
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  logit_exp <- exp(mu + uv)
  prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
  
  #微分関数の設定
  er <- theta_v - v_mu
  dlogit <- matrix(0, nrow=hhpt, ncol=k2)
  for(j in 1:(select-1)){
    dlogit[, allocation_index[j, ]] <- (y[, j] - prob[, j]) * theta_du[, allocation_index[j, ]]
  }
  dmvn <- -t(inv_Cov_v %*% t(er))
  
  #対数事後分布の微分関数の和
  LLd <- -(as.matrix(item_vec %*% dlogit) + dmvn)
  return(LLd)
}

#アイテムの行列分解のリープフロッグ法を解く関数
leapfrog_v <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, g, k2, item_id, item_vec, allocation_index) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, theta_du, mu, y, v_mu, inv_Cov_v, item, hhpt, select, g, k2, item_id, item_vec, allocation_index) / 2
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
R <- 1000
keep <- 2  
iter <- 0
burnin <- 100/keep
disp <- 5
e1 <- 0.1
e2 <- 0.05
L1 <- 3
L2 <- 5

##インデックスとデータの設定
#インデックスの設定
user_vec <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
user_dt <- sparseMatrix(user_id_vec, 1:length(user_id_vec), x=rep(1, length(user_id_vec)), dims=c(hh, length(user_id_vec)))
item_vec <- sparseMatrix(item_id, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))

#データの設定
g <- r
k1 <- ncol(Data)
k2 <- r*(select-1)
y_vec <- as.numeric(t(y))

##事前分布の設定
#素性ベクトルの階層モデルの事前分布
Deltabar <- matrix(0, nrow=ncol(u), ncol=k1)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu <- select-1   #逆ウィシャート分布の自由度
V <- nu * diag(k1) #逆ウィシャート分布のパラメータ

#ユーザーの行列分解のパラメータの階層モデルの事前分布
Deltabar_u <- matrix(0, nrow=ncol(u), ncol=k2)   #階層モデルの回帰係数の事前分布の分散
ADelta_u <- 0.01 * diag(ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- select-1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(k2) #逆ウィシャート分布のパラメータ

#アイテムの行列分解のパラメータの階層モデルの事前分布
Deltabar_v <- matrix(0, nrow=ncol(v), ncol=k2)   #階層モデルの回帰係数の事前分布の分散
ADelta_v <- 0.01 * diag(ncol(v))   #階層モデルの回帰係数の事前分布の分散
nu2 <- select-1   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(k2) #逆ウィシャート分布のパラメータ


##パラメータの真値
#階層モデルのパラメータ
alpha <- alphat; alpha_u <- alpha_ut; alpha_v <- alpha_vt
Cov <- Covt; Cov_u <- Cov_ut; Cov_v <- Cov_vt
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
inv_Cov <- solve(Cov); inv_Cov_u <- solve(Cov_u); inv_Cov_v <- solve(Cov_v)

#モデルのパラメータ
theta <- thetat
theta_u <- theta_ut
theta_v <- theta_vt

#パラメータの設定
mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #素性ベクトルの期待値
uv <- matrix(0, nrow=hhpt, ncol=select)
uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
for(j in 1:(select-1)){
  uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #行列分解の期待値
} 

#ベストな対数尤度
logit_exp <- exp(mu + uv)
prob <- logit_exp / rowSums(logit_exp)
LL_best <- sum(y * log(prob))

##初期値の設定
#階層モデルのパラメータの初期値
alpha <- matrix(0, nrow=ncol(u), ncol=k1); alpha_mu <- u %*% alpha
alpha_u <- matrix(0, nrow=ncol(u), ncol=k2); alpha_u <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=k2); alpha_v <- v %*% alpha_v
Cov <- diag(0.01, k1); inv_Cov <- solve(Cov)
Cov_u <- diag(0.01, k2); inv_Cov_u <- solve(Cov_u)
Cov_v <- diag(0.01, k2); inv_Cov_v <- solve(Cov_v)

#モデルのパラメータの初期値
theta <- alpha_mu + mvrnorm(hh, rep(0, k1), Cov)
theta_u <- u_mu + mvrnorm(hh, rep(0, k2), Cov_u)
theta_v <- v_mu + mvrnorm(item, rep(0, k2), Cov_v)

#パラメータの設定
mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   #素性ベクトルの期待値
uv <- matrix(0, nrow=hhpt, ncol=select)
uv_dt <- theta_u[user_id, ] * theta_v[item_id, ]
for(j in 1:(select-1)){
  uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)   #行列分解の期待値
} 

##パラメータの格納用配列
ALPHA <- array(0, dim=c(ncol(u), k1, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), k2, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k2, R/keep))
COV <- array(0, dim=c(k1, k1, R/keep))
COV_U <- array(0, dim=c(k2, k2, R/keep))
COV_V <- array(0, dim=c(k2, k2, R/keep))

#モデルのパラメータの格納用配列
d <- 0
THETA <- matrix(0, nrow=hh, ncol=k1)
THETA_U <- matrix(0, nrow=hh, ncol=k2)
THETA_V <- matrix(0, nrow=item, ncol=k2)


##対数尤度の基準値
#初期値の対数尤度
logit_exp <- exp(mu + uv)
prob <- logit_exp / rowSums(logit_exp)
LL1 <- sum(y * log(prob))

#平均モデルの対数尤度
LLst <- sum(log(colMeans(y)[as.numeric(y %*% 1:select)]))


####HMCでパラメータをサンプリング####
for(rp in 1:R){

  ##素性ベクトルのパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(hh, rep(0, k1), diag(k1))   #標準多変量正規分布からパラメータを生成
  thetad <- theta
  
  #リープフロッグ法による1ステップ移動
  res <- leapfrog(rold, thetad, dloglike, e1, L1)
  rnew <- res$r
  thetan <- as.matrix(res$z)
  
  #ロジットのパラメータの設定
  mu_old <- matrix(as.numeric((Data * thetad[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  mu_new <- matrix(as.numeric((Data * thetan[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  er_old <- thetad - alpha_mu
  er_new <- thetan - alpha_mu
  
  #移動前と移動後のハミルトニアン
  Hnew <- -(loglike(mu_new, uv, y, hh, select, user_vec) + dmv(er_new, inv_Cov, k1)) + as.numeric(rnew^2 %*% rep(1, k1))/2
  Hold <- -(loglike(mu_old, uv, y, hh, select, user_vec) + dmv(er_old, inv_Cov, k1)) + as.numeric(rold^2 %*% rep(1, k1))/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(hh) #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma1 <- mean(gamma)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  theta <- flag*thetan + (1-flag)*thetad
  
  
  ##ユーザーの行列分解のパラメータをサンプリング
  #素性ベクトルのロジット
  mu <- matrix(as.numeric((Data * theta[user_id_vec, ]) %*% rep(1, k1)), nrow=hhpt, ncol=select, byrow=T)   
  
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(hh, rep(0, k2), diag(k2))   #標準多変量正規分布からパラメータを生成
  thetad_u <- theta_u
  
  #リープフロッグ法による1ステップ移動
  theta_dv <- theta_v[item_id, ]
  res <- leapfrog_u(rold, thetad_u, dloglike_u, e1, L1)
  rnew <- res$r
  thetan_u <- res$z
  
  #ロジットのパラメータの設定
  uv_old <- uv
  uv_new <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- thetan_u[user_id, ] * theta_dv
  for(j in 1:(select-1)){
    uv_new[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  er_old <- thetad_u - u_mu
  er_new <- thetan_u - u_mu
  
  #移動前と移動後のハミルトニアン
  Hnew <- -(loglike(mu, uv_new, y, hh, select, user_vec) + dmv(er_new, inv_Cov_u, k2)) + as.numeric(rnew^2 %*% rep(1, k2))/2
  Hold <- -(loglike(mu, uv_old, y, hh, select, user_vec) + dmv(er_old, inv_Cov_u, k2)) + as.numeric(rold^2 %*% rep(1, k2))/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(hh) #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma2 <- mean(gamma)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  theta_u <- flag*thetan_u + (1-flag)*thetad_u
  
  #行列分解のパラメータを更新
  theta_du <- theta_u[user_id, ]
  uv_dt <- theta_du * theta_dv
  uv <- matrix(0, nrow=hhpt, ncol=select)
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  
  
  ##アイテムの行列分解のパラメータをサンプリング
  #HMCの新しいパラメータを生成
  rold <- mvrnorm(item, rep(0, k2), diag(k2))   #標準多変量正規分布からパラメータを生成
  thetad_v <- theta_v
  
  #リープフロッグ法による1ステップ移動
  theta_du <- theta_u[user_id, ]
  res <- leapfrog_v(rold, thetad_v, dloglike_v, e2, L2)
  rnew <- res$r
  thetan_v <- res$z
  
  #ロジットのパラメータの設定
  uv_old <- uv
  uv_new <- matrix(0, nrow=hhpt, ncol=select)
  uv_dt <- theta_du * thetan_v[item_id, ]
  for(j in 1:(select-1)){
    uv_new[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  er_old <- thetad_v - v_mu
  er_new <- thetan_v - v_mu
  
  #移動前と移動後のハミルトニアン
  Hnew <- -(loglike(mu, uv_new, y, item, select, item_vec) + dmv(er_new, inv_Cov_v, k2)) + as.numeric(rnew^2 %*% rep(1, k2))/2
  Hold <- -(loglike(mu, uv_old, y, item, select, item_vec) + dmv(er_old, inv_Cov_v, k2)) + as.numeric(rold^2 %*% rep(1, k2))/2
  
  #HMC法によりパラメータの採択を決定
  rand <- runif(item)   #一様分布から乱数を発生
  gamma <- rowMins(cbind(1, exp(Hold - Hnew)))   #採択率を決定
  gamma3 <- mean(gamma)
  
  #alphaの値に基づき新しいbetaを採択するかどうかを決定
  flag <- as.numeric(gamma > rand)
  theta_v <- flag*thetan_v + (1-flag)*thetad_v
  
  #行列分解のパラメータを更新
  uv_dt <- theta_du * theta_v[item_id, ]
  uv <- matrix(0, nrow=hhpt, ncol=select)
  for(j in 1:(select-1)){
    uv[, j] <- uv_dt[, allocation_index[j, ]] %*% rep(1, r)
  }
  
  ##階層モデルのパラメータをサンプリング
  #多変量回帰モデルから素性ベクトルの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta, u, Deltabar, ADelta, nu, V)
  alpha <- out$B
  alpha_mu <- u %*% alpha   #素性ベクトルの期待値
  Cov <- out$Sigma
  inv_Cov <- solve(Cov)
  
  #多変量回帰モデルからユーザー特徴行列の階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_u, u, Deltabar_u, ADelta_u, nu1, V1)
  alpha_u <- out$B
  u_mu <- u %*% alpha_u   #ユーザー特徴行列の期待値
  Cov_u <- out$Sigma
  inv_Cov_u <- solve(Cov_u)
  
  #多変量回帰モデルからアイテム特徴行列の階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu2, V2)
  alpha_v <- out$B
  v_mu <- v %*% alpha_v   #アイテム特徴行列の期待値
  Cov_v <- out$Sigma
  inv_Cov_v <- solve(Cov_v)
  
  ##パラメータの格納とサンプリング結果の表示
  #パラメータの格納
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
  
  #サンプリング結果の表示
  if(rp%%disp==0){
    #対数尤度の計算と表示
    logit_exp <- exp(mu + uv)
    prob <- logit_exp / as.numeric(logit_exp %*% rep(1, select))
    LL <- sum(log((y * prob) %*% rep(1, select)))
  
    #結果の確認
    print(rp)
    print(c(LL, LL_best, LL1, LLst))
    print(round(cbind(alpha[, 1:(select-1)], alphat[, 1:(select-1)]), 3))
    print(round(rbind(diag(Cov)[1:((select-1)*3)], diag(Covt)[1:((select-1)*3)]), 3))
    print(round(c(gamma1, gamma2, gamma3), 3))
  }
}


