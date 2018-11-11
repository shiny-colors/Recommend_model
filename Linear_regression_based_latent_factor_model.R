#####線形回帰ベース潜在因子モデル#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)


####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

##相関行列から分散共分散行列を作成する関数を定義
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
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
##データの設定
k <- 10   #基底数
hh <- 5000   #ユーザー数
item <- 2000   #アイテム数

##IDの設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)


##素性ベクトルを生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hh*item*k1, 0, 1), nrow=hh*item, ncol=k1)
x2 <- matrix(0, nrow=hh*item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hh*item, 1, pr)
}
x3 <- rmnom(hh*item, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x0 <- cbind(x1, x2, x3)   #データを結合


##階層モデルの説明変数を生成
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


##応答変数が妥当になるまでパラメータの生成を繰り返す
for(rp in 1:1000){
  print(rp)
  
  ##素性ベクトルの回帰係数を生成
  beta <- rep(0, ncol(x0))
  for(j in 1:ncol(x0)){
    beta[j] <- runif(1, -0.6, 1.6)
  }
  betat <- beta
  
  ##階層モデルのパラメータを生成
  ##ユーザーベースの階層モデルのパラメータ
  #分散共分散行列を設定
  sigma_ut <- sigma_u <- 0.4
  Cov_ut <- Cov_u <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  
  #回帰係数を設定
  alpha_u <- matrix(0, nrow=ncol(u), ncol=k+1)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u[j, ] <- runif(k+1, -0.55, 1.3)
    } else {
      alpha_u[j, ] <- runif(k+1, -0.4, 0.5)
    }
  }
  alpha_ut <- alpha_u
  
  #多変量回帰モデルからユーザー個別の回帰パラメータを生成
  theta_u <- u %*% alpha_u + cbind(rnorm(hh, 0, sigma_u), mvrnorm(hh, rep(0, k), Cov_u))
  theta_ut1 <- theta_u1 <- theta_u[, 1]   #ランダム効果のパラメータ
  theta_ut2 <- theta_u2 <- theta_u[, -1]   #行列分解のパラメータ
  
  
  ##アイテムベースの階層モデルのパラメータ
  #分散共分散行列を設定
  sigma_vt <- sigma_v <- 0.4
  Cov_vt <- Cov_v <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.01, 0.25)$covariance
  
  #回帰係数を設定
  alpha_v <- matrix(0, nrow=ncol(v), ncol=k+1)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(k+1, -0.55, 1.3)
    } else {
      alpha_v[j, ] <- runif(k+1, -0.4, 0.5)
    }
  }
  alpha_vt <- alpha_v
  
  #多変量回帰モデルからアイテム個別の回帰パラメータを生成
  theta_v <- v %*% alpha_v + cbind(rnorm(item, 0, sigma_v), mvrnorm(item, rep(0, k), Cov_v))
  theta_vt1 <- theta_v1 <- theta_v[, 1]   #ランダム効果のパラメータ
  theta_vt2 <- theta_v2 <- theta_v[, -1]   #行列分解のパラメータ
  
  
  ##LRBLFモデルの応答変数を生成
  #モデルの観測誤差
  sigmat <- sigma <- 0.75   
  
  #モデルの平均構造から応答変数を生成
  mu <- x0 %*% beta + theta_u1[user_id0] + theta_v1[item_id0] + as.numeric(t(theta_u2 %*% t(theta_v2)))   #平均構造
  y0 <- rnorm(hh*item, mu, sigma)
  
  #応答変数のbreak条件
  if(mean(y0) > 4.5 & mean(y0) < 6.5 & min(y0) > -7.5 & max(y0) < 17.5) break
}

#生成したスコアを評価データに変換
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y_full <- round(y0_censor, 0)   #スコアを丸める


##欠損ベクトルを生成
#欠損有無のベータ分布のパラメータを設定
beta1 <- rbeta(hh, 8.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 6.5, 8.0)   #アイテム購買確率

#欠損がある購買データを生成
Z <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  deficit <- rbinom(hh, 1, beta1 * beta2[j])
  Z[, j] <- deficit   #欠損を代入
}

#欠損インデックス
z_vec <- as.numeric(t(Z))
index_z1 <- which(z_vec==1)
index_z0 <- which(z_vec==0)
N <- length(index_z1)

#欠損ベクトルに応じてデータを抽出
user_id <- user_id0[index_z1]
item_id <- item_id0[index_z1]
x <- x0[index_z1, ]
y <- y_full[index_z1]
n1 <- plyr::count(user_id)$freq
n2 <- plyr::count(item_id)$freq


#生成した応答変数のヒストグラム
hist(y0, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #元データ
hist(y_full, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #完全データのスコア分布
hist(y, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #購買データのスコア分布


####モンテカルロEMアルゴリズムでLRBLFモデルを推定####
##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
tol <- 0.25
iter <- 1
dl <- 100
L <- 500   #モンテカルロサンプリング数

##初期値の設定
#モデルパラメータの初期値
beta <- rep(0, ncol(x))
beta_mu <- x %*% beta   #固定効果の期待値   
sigma <- 0.5

#ユーザーベースのパラメータの初期値
sigma_u <- 0.2 
Cov_u <- 0.1 * diag(k)
inv_Cov_u <- solve(Cov_u)
alpha_u <- matrix(0, nrow=ncol(u), ncol=k+1)
theta_mu1 <- u %*% alpha_u
theta_u1 <- theta_mu11 <- theta_mu1[, 1]
theta_mu12 <- theta_mu1[, -1]
theta_u2 <- theta_mu12 + mvrnorm(hh, rep(0, k), Cov_u)

#アイテムベースのパラメータの初期値
sigma_v <- 0.2
Cov_v <- 0.1 * diag(k)
inv_Cov_v <- solve(Cov_v)
alpha_v <- matrix(0, nrow=ncol(v), ncol=k+1)
theta_mu2 <- v %*% alpha_v
theta_v1 <- theta_mu21 <- theta_mu2[, 1]
theta_mu22 <- t(theta_mu2[, -1])
theta_v2 <- theta_mu22 + t(mvrnorm(item, rep(0, k), Cov_v))

#行列分解の初期値
uv <- as.numeric(t(theta_u2 %*% theta_v2))[index_z1]


##データの設定
#インデックスの作成
user_list <- item_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}
#定数の設定
inv_xx <- solve(t(x) %*% x)
inv_uu <- solve(t(u) %*% u)
inv_vv <- solve(t(v) %*% v)

#対数尤度の初期値
Mu <- beta_mu + theta_u1[user_id] + theta_v1[item_id] + uv   #完全データの平均構造
LL <- sum(dnorm(y, Mu, sigma, log=TRUE))   #完全データの対数尤度を更新
print(LL)


####モンテカルロEMアルゴリズムをパラメータを推定####
while(dl > 0){   #dlがtol以上なら繰り返す
  
  ###モンテカルロEステップで潜在変数をサンプリング
  ##ユーザーのランダム効果をサンプリング
  #データの設定
  u_er <- as.numeric(y - beta_mu - theta_v1[item_id] - uv)   #ユーザーのランダム効果の誤差
  
  #事後分布のパラメータを設定
  u_mu <- rep(0, hh)
  for(i in 1:hh){
    u_mu[i] <- mean(u_er[user_list[[i]]])
  }
  weights <- sigma_u^2 / (sigma^2/n1 + sigma_u^2)   #重み係数
  mu_par <- weights*u_mu + (1-weights)*theta_mu11   #事後分布の平均
  
  
  #正規分布より事後分布をサンプリング
  theta_u_data <- matrix(rnorm(hh*L, mu_par, sqrt(1 / (1/sigma_u^2 + n1/sigma^2))), nrow=hh, ncol=L)
  theta_u1 <- rowMeans(theta_u_data)
  u1_vars <- rowVars(theta_u_data)
  
  
  ##アイテムのランダム効果をサンプリング
  #データの設定
  i_er <- as.numeric(y - beta_mu - theta_u1[user_id] - uv)
  
  #事後分布のパラメータ
  i_mu <- rep(0, item)
  for(j in 1:item){
    i_mu[j] <- mean(i_er[item_list[[j]]])
  }
  weights <- sigma_v^2 / (sigma^2/n2 + sigma_v^2)   #重み係数
  mu_par <- weights*i_mu + (1-weights)*theta_mu21   #事後分布の平均
  
  #正規分布より事後分布をサンプリング
  theta_v1_data <- matrix(rnorm(item*L, mu_par, sqrt(1 / (1/sigma_v^2 + n2/sigma^2))), nrow=item, ncol=L)
  theta_v1 <- rowMeans(theta_v1_data)
  v1_vars <- rowVars(theta_v1_data)
  
  
  ##ユーザー特徴行列のパラメータをサンプリング
  #データの設定
  theta_u_vec <- theta_u1[user_id]; theta_v_vec <- theta_v1[item_id]
  uv_er <- as.numeric(y - beta_mu - theta_u_vec - theta_v_vec)
  theta_v2_T <- t(theta_v2)
  theta_u2 <- matrix(0, nrow=hh, ncol=k)
  u2_vars <- matrix(0, nrow=k, ncol=k)
  
  #ユーザーごとに特徴ベクトルをサンプリング
  for(i in 1:hh){
    
    #特徴ベクトルの事後分布のパラメータ
    index <- item_id[user_list[[i]]]   #アイテムインデックス
    Xy <- t(theta_v2_T[index, , drop=FALSE]) %*% uv_er[user_list[[i]]]
    XXV <- (t(theta_v2_T[index, , drop=FALSE]) %*% theta_v2_T[index, , drop=FALSE]) + inv_Cov_u
    inv_XXV <- solve(XXV)
    mu <- inv_XXV %*% (Xy + inv_Cov_u %*% theta_mu12[i, ])   #事後分布の平均
    
    #多変量正規分布からユーザー特徴ベクトルをサンプリング
    theta_u2_data <- mvrnorm(L, mu, sigma^2*inv_XXV)
    theta_u2[i, ] <- colMeans(theta_u2_data)   #モンテカルロ平均
    u2_vars <- u2_vars + var(theta_u2_data)
  }
  
  
  ##アイテム特徴行列のパラメータをサンプリング
  #アイテムごとに特徴ベクトルをサンプリング
  theta_v2 <- matrix(0, nrow=k, ncol=item)
  v2_vars <- matrix(0, nrow=k, ncol=k)
  for(j in 1:item){
    
    #特徴ベクトルの事後分布のパラメータ
    index <- user_id[item_list[[j]]]   #アイテムインデックス
    Xy <- t(theta_u2[index, , drop=FALSE]) %*% uv_er[item_list[[j]]]
    XXV <- (t(theta_u2[index, , drop=FALSE]) %*% theta_u2[index, , drop=FALSE]) + inv_Cov_v
    inv_XXV <- solve(XXV)
    mu <- inv_XXV %*% (Xy + inv_Cov_v %*% theta_mu22[, j])   #事後分布の平均
    
    #多変量正規分布からアイテム特徴ベクトルをサンプリング
    theta_v2_data <- mvrnorm(L, mu, sigma^2*inv_XXV)
    theta_v2[, j] <- colMeans(theta_v2_data)   #モンテカルロ平均
    v2_vars <- v2_vars + var(theta_v2_data)
  }
  
  
  ###Mステップで完全データの尤度を最大化
  #行列分解のパラメータ
  uv <- as.numeric(t(theta_u2 %*% theta_v2))[index_z1]
  
  ##素性ベクトルのパラメータを更新
  y_er <- y - theta_u_vec - theta_v_vec - uv   #応答変数を設定
  beta <- inv_xx %*% t(x) %*% y_er   #最小二乗法で素性ベクトルを更新
  beta_mu <- x %*% beta   #素性ベクトルの平均構造
  
  ##観測モデルの誤差パラメータを更新
  er <- y - beta_mu - theta_u_vec - theta_v_vec - uv
  sigma <- sd(er)
  
  
  ##階層モデルのパラメータを更新
  #ユーザーのランダム効果の階層モデルのパラメータを更新
  alpha_u[, 1] <- inv_uu %*% t(u) %*% theta_u1
  theta_mu11 <- as.numeric(u %*% alpha_u[, 1])
  sigma_u <- sqrt((sum(u1_vars) + sum((theta_u1 - theta_mu11)^2)) / hh)
  
  
  #アイテムのランダム効果の階層モデルのパラメータを更新
  alpha_v[, 1] <- inv_vv %*% t(v) %*% theta_v1
  theta_mu21 <- as.numeric(v %*% alpha_v[, 1])
  sigma_v <- sqrt((sum(v1_vars) + sum((theta_v1 - theta_mu21)^2)) / item)
  
  #ユーザー特徴行列の階層モデルのパラメータを更新
  alpha_u[, -1] <- inv_uu %*% t(u) %*% theta_u2
  theta_mu12 <- u %*% alpha_u[, -1]
  Cov_u <- (u2_vars + t(theta_u2 - theta_mu12) %*% (theta_u2 - theta_mu12)) / hh
  inv_Cov_u
  
  #アイテム特徴行列の階層モデルのパラメータを更新
  alpha_v[, -1] <- inv_vv %*% t(v) %*% t(theta_v2)
  theta_mu22 <- t(v %*% alpha_v[, -1])
  Cov_v <- (v2_vars + (theta_v2 - theta_mu22) %*% t(theta_v2 - theta_mu22)) / item
  inv_Cov_v <- solve(Cov_v)
  
  ##アルゴリズムの収束判定
  Mu <- beta_mu + theta_u_vec + theta_v_vec + uv   #完全データの平均構造
  LL <- sum(dnorm(y, Mu, sigma, log=TRUE))   #完全データの対数尤度を更新
  iter <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}

####推定結果の確認と適合度####
##比較対象モデルの対数尤度
#最小二乗法での対数尤度
X <- cbind(1, x, u[user_id, -1], v[item_id, -1])
b <- solve(t(X) %*% X) %*% t(X) %*% y
tau <- sd(y - X %*% b)
LLsq <- sum(dnorm(y, X %*% b, tau, log=TRUE))

#真値での対数尤度
uvt <- as.numeric(t(theta_ut2 %*% t(theta_vt2)))[index_z1]
beta_mut <- x %*% betat
Mut <- beta_mut + theta_ut1[user_id] + theta_vt1[item_id] + uvt   #完全データの平均構造
LLt <- sum(dnorm(y, Mut, sd(y-Mut), log=TRUE))

#学習データに対する対数尤度の比較
c(LL, LLt, LLsq)


##欠損データ(テストデータ)に対する対数尤度と二乗誤差を計算
#推定結果での平均構造
beta_mu <- x0 %*% beta
uv <- as.numeric(t(theta_u2 %*% theta_v2))
Mu <- as.numeric(beta_mu + theta_u1[user_id0] + theta_v1[item_id0] + uv)
sigmaf <- sd(y_full[index_z0] - Mu[index_z0])

#真値での平均構造
beta_mut <- x0 %*% betat
uvt <- as.numeric(t(theta_ut2 %*% t(theta_vt2)))
Mut <- as.numeric(beta_mut + theta_ut1[user_id0] + theta_vt1[item_id0] + uvt)   #完全データの平均構造

#最小二乗法での平均構造
X <- cbind(1, x0, u[user_id0, -1], v[item_id0, -1])
Mu_sq <- X %*% b

#テストデータに対する対数尤度
sum(dnorm(y_full[index_z0], Mu[index_z0], sigmaf, log=TRUE))   #推定結果の対数尤度
sum(dnorm(y_full[index_z0], Mut[index_z0], sigmat, log=TRUE))   #真値の対数尤度
sum(dnorm(y_full[index_z0], Mu_sq[index_z0], tau, log=TRUE))   #最小二乗法の対数尤度

#テストデータに対する二乗誤差
sum((y_full[index_z0] - Mu[index_z0])^2)   #推定結果の二乗誤差
sum((y_full[index_z0] - Mut[index_z0])^2)   #真値の二乗誤差
sum((y_full[index_z0] - Mu_sq[index_z0])^2)   #最小二乗法の二乗誤差

#結果を比較
result <- round(data.table(z_vec, y_full, Mu, Mut, Mu_sq), 3)
colnames(result) <- c("z_vec", "y_full", "Mu", "Mut", "Mu_sq")
