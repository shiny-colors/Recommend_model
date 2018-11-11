#####Constrained Multi Matrix Factorization#####
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
context <- 200   #コンテキスト数

#コンテキストを生成
context_list <- list()
prob <- as.numeric(extraDistr::rdirichlet(1, rep(20.0, context)))
for(j in 1:item){
  if(j%%100==0){
    print(j)
  }
  context_list[[j]] <- as.numeric(rmnom(hh, 1, prob) %*% 1:context)
}

##IDの設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)
context_id0 <- unlist(context_list)


##応答変数が妥当になるまでパラメータの生成を繰り返す
for(rp in 1:1000){
  print(rp)
  
  ##CMMFモデルの応答変数を生成
  sigmat <- sigma <- 0.5   #モデルの観測誤差 
  
  #制約付きの特徴行列
  W <- WT <- mvrnorm(hh, rep(0.675, k), diag(0.25, k))
  H <- HT <- mvrnorm(item, rep(0.5, k), diag(0.2, k))
  V <- VT <- mvrnorm(context, rep(0.25, k), diag(0.15, k))
  
  #モデルの平均構造から応答変数を生成
  mu <- rowSums(W[user_id0, ] * H[item_id0, ]) + rowSums(W[user_id0, ] * V[context_id0, ])
  y0 <- rtnorm(hh*item, mu, sigma, a=1, b=10)
  
  #応答変数のbreak条件
  if(mean(y0) < 6.5 & mean(y0) > 4.5) break
}

#生成したスコアを評価データに変換
y0_censor <- ifelse(y0 < 1, 1, ifelse(y0 > 10, 10, y0)) 
y_vec <- round(y0_censor, 0)   #スコアを丸める

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
context_id <- context_id0[index_z1]
y <- y_vec[index_z1]
n1 <- plyr::count(user_id)$freq
n2 <- plyr::count(item_id)$freq
n3 <- plyr::count(context_id)$freq

#生成した応答変数のヒストグラム
hist(y0, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #元データ
hist(y_vec, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #完全データのスコア分布
hist(y, col="grey", xlab="スコア", main="ユーザー×アイテムのスコア分布")   #購買データのスコア分布


####マルコフ連鎖モンテカルロ法でCMMFを推定####
##アルゴリズムの設定
R <- 2000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
theta <- rep(0, k)
tau <- 100 * diag(k)
inv_tau <- solve(tau)
s0 <- 1.0
v0 <- 1.0

##真値の設定
W <- WT
H <- HT
V <- VT
sigma <- sigmat

##初期値の設定
sigma <- sd(y)
W <- mvrnorm(hh, rep(0, k), 0.2 * diag(k))
H <- mvrnorm(item, rep(0, k), 0.2 * diag(k))
V <- mvrnorm(context, rep(0, k), 0.2 * diag(k))

##パラメータの格納用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(item, k, R/keep))
V_array <- array(0, dim=c(context, k, R/keep))
SIGMA <- rep(0, R/keep)


##インデックスを設定
user_index <- item_index <- context_index <- list()
ui_id <- uc_id <- iu_id <- ic_id <- cu_id <- ci_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  uc_id[[i]] <- context_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  ic_id[[j]] <- context_id[item_index[[j]]]
}
for(j in 1:context){
  context_index[[j]] <- which(context_id==j)
  cu_id[[j]] <- user_id[context_index[[j]]]
  ci_id[[j]] <- item_id[context_index[[j]]]
}
vec <- rep(1, k)
const <- hh / 100


##対数尤度の基準値
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザーの特徴行列をサンプリング
  for(i in 1:hh){
    #特徴ベクトルのパラメータ
    X <- H[ui_id[[i]], ] + V[uc_id[[i]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% y[user_index[[i]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    W[i, ] <- mvrnorm(1, beta_mu, cov)  
  }

  
  ##アイテムの特徴行列をサンプリング
  W_vec <- W[user_id, ]
  er <- as.numeric(y - (W_vec * V[context_id, ]) %*% vec)   #誤差を設定
  
  for(j in 1:item){
    #特徴ベクトルのパラメータ
    X <- W[iu_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% er[item_index[[j]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    H[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  ##コンテキストの特徴行列をサンプリング
  er <- as.numeric(y - (W_vec * H[item_id, ]) %*% vec)   #誤差を設定
  for(j in 1:context){
    #特徴ベクトルのパラメータ
    X <- W[cu_id[[j]], ]
    inv_XXV <- solve(t(X) %*% X + inv_tau)
    beta_mu <- inv_XXV %*% t(X) %*% er[context_index[[j]]]   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_XXV
    
    #多変量正規分布からパラメータをサンプリング
    V[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  ##モデルパラメータをサンプリング
  #モデルの標準偏差をサンプリング
  WHV <- as.numeric((W_vec * H[item_id, ]) %*% vec + (W_vec * V[context_id, ]) %*% vec)
  er <- y - WHV
  s <- s0 + t(er) %*% er
  v <- v0 + N
  sigma <- sqrt(1/(rgamma(1, v/2, s/2)))   #逆ガンマ分布からsigmaをサンプリング
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    #パラメータの格納
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    V_array[, , mkeep] <- V
    SIGMA[mkeep] <- sigma
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    LL <- sum(dnorm(y, WHV, sigma, log=TRUE))
    print(rp)
    print(c(LL, LLst))
    print(round(c(sigma, sigmat), 3))
  }
}

matplot(t(W_array[1, , ]), type="l")

