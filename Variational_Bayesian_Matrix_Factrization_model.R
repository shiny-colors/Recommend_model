#####欠損値のある変分ベイズ行列因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(5897)

####データの発生####
##データの設定
k <- 15   #基底数
hh <- 5000   #ユーザー数
item <- 1500   #アイテム数

##パラメータの設定
sigma <- 1
A <- A_T <- mvrnorm(hh, rep(0, k), diag(1, k))   #ユーザーの特徴行列
B <- B_T <- mvrnorm(item, rep(0, k), diag(1, k))   #アイテムの特徴行列
beta1 <- rbeta(hh, 8.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 5.0, 6.0)   #アイテム購買確率


##モデルに基づき応答変数を生成
AB <- A %*% t(B)   #期待値
Y <- matrix(0, nrow=hh, ncol=item)
Z <- matrix(0, nrow=hh, ncol=item)

for(j in 1:item){
  #評価ベクトルを生成
  y_vec <- rnorm(hh, AB[, j], 1)   #正規分布から評価ベクトルを生成
  
  #欠損を生成
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z[, j] <- deficit   #欠損を代入
  
  #評価ベクトルを代入
  Y[, j] <- y_vec
}

##IDと評価ベクトルを設定
N <- length(Z[Z==1])
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

#評価がある要素のみ抽出
index_user <- which(as.numeric(t(Z))==1)
user_id <- user_id0[index_user]
item_id <- item_id0[index_user]
y_vec <- as.numeric(t(Y))

#評価ベクトルを1～5の間の離散値に収める
score_mu <- 3   #平均スコア
y0 <- as.numeric(round(scale(y_vec) + score_mu))   #平均3の整数値評価ベクトル
y0[y0 < 1] <- 1
y0[y0 > 5] <- 5
y <- y0[index_user]   #欠損のある評価ベクトル
Y <- matrix(y0, nrow=hh, ncol=item, byrow=T)   #評価ベクトルの完全データ


##インデックスの作成
index_user <- list()
index_item <- list()
for(i in 1:hh){
  index_user[[i]] <- which(user_id==i)
}
for(j in 1:item){
  index_item[[j]] <- which(item_id==j)
}

####変分ベイズ法でパラメータを推定####
##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
tol <- 1
iter <- 1
dl <- 100

##事前分布の設定
sigma <- 1
Ca <- diag(1, k)
Cb <- diag(1, k)

##初期値の設定
Cov_A <- array(diag(1, k), dim=c(k, k, hh))
Cov_B <- array(diag(1, k), dim=c(k, k, item))
A <- mvrnorm(hh, rep(0.0, k), diag(0.5, k))
B <- mvrnorm(item, rep(0.0, k), diag(0.5, k))


####変分ベイズ法でパラメータを更新####
while(abs(dl) > tol){   #dlがtol以上なら繰り返す

  ##ユーザー特徴行列のパラメータを更新
  A <- matrix(0, nrow=hh, ncol=k)
  Cov_A <- array(0, dim=c(k, k, hh))
  
  #分散成分を更新
  for(i in 1:hh){
    index <- item_id[index_user[[i]]]
    Cov_sum <- matrix(0, k, k)
    for(j in 1:length(index)){
      Cov_sum <- Cov_sum + Cov_B[, , index[j]]
    }
    Cov_A[, , i] <- sigma^2 * (solve((t(B[index, ]) %*% B[index, ] + Cov_sum) + sigma^2 * solve(Ca)))
  }
  
  #ユーザーごとに特徴ベクトルを更新
  for(i in 1:hh){
    A[i, ] <- sigma^-2 * Cov_A[, , i] %*% colSums(y[index_user[[i]]] * B[item_id[index_user[[i]]], ])
  }
  
  ##アイテム特徴行列のパラメータを更新
  #分散成分を更新
  B <- matrix(0, nrow=item, ncol=k)
  Cov_B <- array(0, dim=c(k, k, item))
  
  #分散成分を更新
  for(j in 1:item){
    index <- user_id[index_item[[j]]]
    Cov_sum <- matrix(0, k, k)
    for(l in 1:length(index)){
      Cov_sum <- Cov_sum + Cov_A[, , index[l]]
    }
    Cov_B[, , j] <- sigma^2 * solve((t(A[index, ]) %*% A[index, ] + Cov_sum) + sigma^2 * solve(Cb))
  }
  
  #アイテムごとに特徴ベクトルを更新
  for(j in 1:item){
    B[j, ] <- sigma^-2 * Cov_B[, , j] %*% colSums(y[index_item[[j]]] * A[user_id[index_item[[j]]], ])
  }
  
  ##ハイパーパラメータを更新
  #パラメータの事前分布を更新
  Ca <- diag((colSums(A^2) + diag(apply(Cov_A, c(1, 2), sum))) / hh)
  Cb <- diag((colSums(B^2) + diag(apply(Cov_B, c(1, 2), sum))) / item)
  
  #標準偏差を更新
  score <- rowSums(A[user_id, ] * B[item_id, ])
  sigma <- sqrt(sum((y - score)^2) / length(y))
  
  
  ##アルゴリズムの収束判定
  LL <- sum(dnorm(y, score, sigma, log=TRUE))   #対数尤度を更新
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
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
