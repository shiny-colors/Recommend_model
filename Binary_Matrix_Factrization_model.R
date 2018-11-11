#####欠損値のあるベイジアンバイナリ行列因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(condMVNorm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)
source("bdiag_m.R")

#set.seed(5897)

####データの発生####
##データの設定
k <- 15   #基底数
hh <- 5000   #ユーザー数
item <- 1500   #アイテム数

##パラメータの設定
sigma <- 1
A <- A_T <- mvrnorm(hh, rep(-0.5, k), diag(1, k))   #ユーザーの特徴行列
B <- B_T <- mvrnorm(item, rep(0.5, k), diag(1, k))   #アイテムの特徴行列
beta1 <- rbeta(hh, 8.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 5.0, 6.0)   #アイテム購買確率


##モデルに基づき応答変数を生成
AB <- A %*% t(B)   #期待値
U <- matrix(0, nrow=hh, ncol=item)
Y <- matrix(0, nrow=hh, ncol=item)
Z <- matrix(0, nrow=hh, ncol=item)

for(j in 1:item){
  #潜在効用を生成
  u_vec <- rnorm(hh, AB[, j], sigma)   #正規分布から潜在効用を生成
  U[, j] <- u_vec
  
  #欠損を生成
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z[, j] <- deficit   #欠損を代入
  
  #評価ベクトルを代入
  Y[, j] <- ifelse(u_vec > 0, 1, 0)   #バイナリベクトルに変換
}

##IDと評価ベクトルを設定
N <- length(Z[Z==1])
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#評価がある要素のみ抽出
index_z <- index_user <- which(as.numeric(t(Z))==1)
user_id <- user_id0[index_user]
item_id <- item_id0[index_user]
y_vec <- as.numeric(t(Y))
u_vec <- as.numeric(t(U))
mean(y_vec); sum(y_vec)


#欠損値のある購買ベクトルに変換
y <- y_vec[index_user]   #欠損のある評価ベクトル
u <- u_vec[index_user]
Y <- matrix(y_vec, nrow=hh, ncol=item, byrow=T)   #購買履歴の完全データ
mean(y); length(y); sum(y)


##インデックスの作成
index_user <- list()
index_item <- list()
for(i in 1:hh){
  index_user[[i]] <- which(user_id==i)
}
for(j in 1:item){
  index_item[[j]] <- which(item_id==j)
}


####マルコフ連鎖モンテカルロ法でパラメータを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
sigma <- 1
tau0 <- 1/100
Ca <- diag(1, k)
Cb <- diag(1, k)

#階層モデルの事前分布
Bbar1 <- Bbar1 <- rep(0, k)
A1 <- B1 <- 0.01 * diag(1, k)
nu1 <- nu2 <- k/2
V1 <- nu1 * diag(k)
V2 <- nu2 * diag(k)

##初期値の設定
Cov_A <- diag(0.5, k); inv_A <- solve(Cov_A)
Cov_B <- diag(0.01, k); inv_B <- solve(Cov_B)
A_mu <- rep(-0.5, k)
B_mu <- rep(0, k)
A <- mvrnorm(hh, rep(0.0, k), Cov_A)
B <- mvrnorm(item, rep(0.0, k), Cov_B)

##パラメータの格納用配列
A_SEG <- matrix(0, nrow=hh, ncol=k)
B_SEG <- matrix(0, nrow=item, ncol=k)
MU_A <- matrix(0, nrow=R/keep, ncol=k)
MU_B <- matrix(0, nrow=R/keep, ncol=k)
COV_A <- matrix(0, nrow=R/keep, ncol=k)
COV_B <- matrix(0, nrow=R/keep, ncol=k)

##切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##対数尤度をベスト
mu_vec <- as.numeric(t(A_T %*% t(B_T)))[index_z]
LLbest <- sum(dbinom(y, 1, pnorm(mu_vec, 0, 1), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在変数をサンプリング
  AB_vec <- as.numeric(t(A %*% t(B)))[index_z]   #潜在効用の平均ベクトル
  U <- rtnorm(AB_vec, sigma, a, b)   #切断正規分布から潜在効用をサンプリング
  
  ##ユーザー特徴行列のパラメータをサンプリング
  for(i in 1:hh){
    index <- item_id[index_user[[i]]]
    
    #特徴行列の事後分布のパラメータ
    Xy <- t(B[index, ]) %*% U[index_user[[i]]]
    XXV <- t(B[index, ]) %*% B[index, ] + inv_A
    inv_XXV <- solve(XXV) 
    mu <- inv_XXV %*% (Xy + inv_A %*% A_mu)
    
    #多変量正規分布から事後平均をサンプリング
    A[i, ] <- mvrnorm(1, mu, inv_XXV)
  }
  A <- A - matrix(colMeans(A) - A_mu, nrow=hh, ncol=k, byrow=T)   #列ベクトルを正規化
  
  ##アイテム特徴行列のパラメータをサンプリング
  for(j in 1:item){
    index <- user_id[index_item[[j]]]
    
    #特徴行列の事後分布のパラメータ
    Xy <- t(A[index, ]) %*% U[index_item[[j]]]
    
    XXV <- t(A[index, ]) %*% A[index, ] + inv_B
    inv_XXV <- solve(XXV) 
    mu <- inv_XXV %*% (Xy + inv_B %*% B_mu)
    
    #多変量正規分布から事後平均をサンプリング
    B[j, ] <- mvrnorm(1, mu, inv_XXV)
  }
  
  ##ユーザー特徴行列の階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  V_par <- V2 + t(A) %*% A
  Sn <- nu1 + hh
  Cov_A <- bayesm::rwishart(Sn, solve(V_par))$IW  
  inv_A <- solve(Cov_A)
  
  
  ##アイテム特徴行列の階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  V_par <- V2 + t(B) %*% B
  Sn <- nu2 + item
  Cov_B <- bayesm::rwishart(Sn, solve(V_par))$IW  
  inv_B <- solve(Cov_B)
  
  #多変量正規分布から平均ベクトルをサンプリング
  beta_mu <- item/(item + tau0) * colMeans(B)
  B_mu <- mvrnorm(1, beta_mu, Cov_B/(item + tau0))
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    MU_A[mkeep, ] <- A_mu
    MU_B[mkeep, ] <- B_mu
    COV_A[mkeep, ] <- diag(Cov_A)
    COV_B[mkeep, ] <- diag(Cov_B)
  } 
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    A_SEG <- A_SEG + A
    B_SEG <- B_SEG + B
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    mu_vec <- as.numeric(t(A %*% t(B)))[index_z]
    LL <- sum(dbinom(y, 1, pnorm(mu_vec, 0, 1), log=TRUE))
    
    #サンプリング結果を確認
    print(rp)
    print(c(LL, LLbest))
    print(round(B_mu, 3))
    print(round(diag(Cov_B), 3))
  }
}

##完全データに対する二乗誤差を推定
#評価ベクトルと欠損ベクトルを設定
z_vec <- as.numeric(t(Z))
y_vec <- as.numeric(t(Y))
mu_vec <- as.numeric(t(A %*% t(B)))

#観測データの二乗誤差
er_obz <- sum((y_vec[z_vec==1] - mu_vec[z_vec==1])^2)
er_obz / sum(z_vec==1)

#欠損データの二乗誤差
er_na <- sum((y_vec[z_vec==0] - mu_vec[z_vec==0])^2)
er_na / sum(z_vec==0)
cbind(z_vec, y_vec, round(mu_vec, 2))

