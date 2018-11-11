#####Bayesian Multinomial Factorization Machines#####
library(MASS)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

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
select <- 10   #選択肢数
s <- 5   #基底数
N <- 20000   #サンプル数
topic <- 30   #トピック数
k <- 3   #データタイプ数

###応答変数が妥当になるまでデータの生成を繰り返す
rp <- 0
repeat { 
  print(rp <- rp + 1)
  
  ##説明変数の生成
  u <- array(0, dim=c(N, topic, k))
  x <- matrix(0, nrow=N, ncol=topic)
  pr <- as.numeric(extraDistr::rdirichlet(1, rep(3.0, k)))
  
  for(j in 1:k){
    #データタイプごとの説明変数を生成
    d <- matrix(rgamma(N*topic, 5.5, 15.0), nrow=N, ncol=topic)
    u[, , j] <- d * matrix(rbinom(N*topic, 1, rep(runif(topic, 0.35, 0.7), rep(N, topic))), nrow=N, ncol=topic)
    
    #統合した説明変数を生成
    x <- x + u[, , j] * pr[j]
  }
  
  ##交差項を設定
  #交差項のインデックスを作成
  v_matrix <- matrix(0, nrow=k-1, ncol=k-1)
  v_matrix[upper.tri(v_matrix)] <- 1
  v_list <- list()
  for(j in 1:(k-1)){
    index <- which(v_matrix[j, ]==1)
    if(length(index) > 0){
      v_list[[j]] <- cbind(j1=j, j2=which(v_matrix[j, ]==1))
    }
  }
  v_index <- do.call(rbind, v_list)
  
  #交差項の入力変数を作成
  v <- u[, v_index[, 1]+1] * u[, v_index[, 2]+1]
  
  #説明変数の割当インデックス
  allocation_index11 <- allocation_index12 <- matrix(0, nrow=k-1, ncol=k-2)
  allocation_index21 <- allocation_index22 <- matrix(0, nrow=k-1, ncol=k-2)
  for(j in 1:(k-1)){
    index <- which(rowSums(v_index==j)==1)
    allocation_index11[j, ] <- allocation_index12[j, ] <- 
      rowSums(matrix(as.numeric(v_index[index, ]!=j), nrow=length(index)) * v_index[index, ])
    allocation_index21[j, ] <- allocation_index22[j, ] <-  index
  }
  allocation_index11[lower.tri(allocation_index11)] <- 0
  allocation_index21[lower.tri(allocation_index21)] <- 0
  vec <- rep(1, s)
  j_data12 <- matrix(1:(k-1), nrow=k-1, ncol=k-2) 
  j_data11 <- j_data12 * (allocation_index11 > 0)
  
  ##パラメータと応答変数を生成
  #パラメータを生成
  sigma <- sigmat <- 0.5
  beta <- betat <- c(5.0, runif(k-1, -1.25, 0.75))
  theta <- thetat <- mvrnorm(k-1, rep(0, s), diag(0.15, s))
  
  #応答変数を生成
  u_mu <- u %*% beta
  v_mu <- rep(0, N)
  for(j in 1:(k-1)){
    v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% ((theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec)
  }
  mu <- u_mu + v_mu   #期待値を算出
  y0 <- rnorm(N, mu, sigma)   #正規分布から応答変数を生成
  
  if(sd(y0) > 1.75 & sd(y0) < 2.25 & min(y0) > -7.5 & max(y0) < 17.5 & mean(y0) > 4.5 & mean(y0) < 6.0){
    break
  }
}

#応答変数を1～10に変換する
y <- round(y0)
y[y > 10] <- 10; y[y < 1] <- 1
hist(y0, breaks=25, col="grey", main="真のスコア分布", xlab="スコア")
hist(y, breaks=25, col="grey", main="切断されたスコア分布", xlab="スコア")


####マルコフ連鎖モンテカルロ法でFMを推定####
##アルゴリズムの設定
R <- 2000
burnin <- 750
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
alpha1 <- rep(0, k)
alpha2 <- rep(0, s)
tau1 <- 100 * diag(k)
tau2 <- 100 * diag(s)
inv_tau1 <- solve(tau1)
inv_tau2 <- solve(tau2)
s0 <- 1.0
v0 <- 1.0

##真値の設定
beta <- betat
sigma <- sigmat
theta <- thetat

##初期値の設定
beta <- as.numeric(solve(t(u) %*% u) %*% t(u) %*% y); beta[1] <- mean(y)
theta <- mvrnorm(k-1, rep(0, s), diag(0.1, s))
sigma <- 1.0


##パラメータの格納用配列
m <- 0
BETA <- matrix(0, nrow=R/keep, ncol=k)
SIGMA <- rep(0, R/keep)
THETA <- matrix(0, nrow=k-1, ncol=s)

##インデックスとデータの定数を設定
#インデックスを設定
index_list11 <- index_list12 <- list()
index_list21 <- index_list22 <- index_list23 <- list()

for(j in 1:(k-1)){
  #データを抽出
  index <- (allocation_index11==j) + (j_data11==j)
  
  #推定パラメータのインデックス
  j_index <- v_index[rowSums(v_index * cbind(v_index[, 1]==j, v_index[, 2]==j)) > 0, ]
  index_list11[[j]] <- j_index[j_index!=j]
  index_list12[[j]] <- allocation_index21[index==1]
  
  #固定パラメータのインデックス
  index1 <- as.numeric(t(allocation_index11 * (1-index)))
  index2 <- as.numeric(t(allocation_index21 * (1-index)))
  index3 <- as.numeric(t(j_data11 * (1-index)))
  index_list21[[j]] <- index1[index1 > 0]
  index_list22[[j]] <- index2[index2 > 0]
  index_list23[[j]] <- index3[index3 > 0]
}

#データの設定
uu <- t(u) %*% u
inv_uu <- solve(uu + inv_tau1)
v_array <- array(0, dim=c(N, k-2, k-1))
for(j in 1:(k-1)){
  v_array[, , j] <- v[, index_list12[[j]]]
}
v_mu <- rep(0, N)
for(j in 1:(k-1)){
  v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
}

##対数尤度の基準値
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##回帰係数をサンプリング
  #応答変数の誤差を算出
  er <- as.numeric(y - v_mu)
  
  #多変量正規分布のパラメータを設定
  beta_mu <- inv_uu %*% t(u) %*% er   #多変量正規分布の平均ベクトル
  cov <- sigma^2 * inv_uu   #多変量正規分布の分散
  
  #多変量正規分布から回帰係数をサンプリング
  beta <- mvrnorm(1, beta_mu, cov)
  u_mu <- u %*% beta
  
  
  ##交互作用項の特徴ベクトルをサンプリング
  for(j in 1:(k-1)){
    #応答変数を設定
    er <- as.numeric(y - u_mu - v[, index_list22[[j]]] %*% (theta[index_list21[[j]], ] * theta[index_list23[[j]], ]) %*% vec)
    
    #特徴ベクトルのパラメータを設定
    x <- v_array[, , j] %*% theta[index_list11[[j]], ] 
    inv_xxv <- solve(t(x) %*% x + inv_tau2)
    theta_mu <- inv_xxv %*% t(x) %*% er   #多変量正規分布の平均べクトル
    cov <- sigma^2 * inv_xxv   #多変量正規分布の分散
    
    #多変量正規分布から回帰係数をサンプリング
    theta[j, ] <- mvrnorm(1, theta_mu, cov)
  }
  
  #交互作用の平均ベクトルを更新
  v_mu <- rep(0, N)
  for(j in 1:(k-1)){
    v_mu <- v_mu + v[, allocation_index21[j, ], drop=FALSE] %*% (theta[j_data11[j, ], ] * theta[allocation_index11[j, ], ]) %*% vec
  }
  
  ##モデルの標準偏差をサンプリング
  mu <- u_mu + v_mu
  er <- as.numeric(y - mu)
  
  #逆ガンマ分布のパラメータ
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- N + v0
  
  #逆ガンマ分布から標準偏差をサンプリング
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##パラメータの格納とサンプリング結果の表示
  #パラメータを格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    SIGMA[mkeep] <- sigma
    
    if(rp >= burnin){
      m <- m + 1
      THETA <- THETA + theta
    }
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLst))
    print(round(rbind(beta, betat), 3))
    print(c(sigma, sigmat))
  }
}

####推定結果の確認と適合度####
##サンプリング結果の可視化
matplot(BETA, type="l", xlab="サンプリング回数", main="betaのサンプリング結果のプロット")
plot(1:length(SIGMA), SIGMA, type="l", xlab="サンプリング回数", main="sigmaのサンプリング結果のプロット")

##パラメータの事後平均
#バーンイン期間
RS1 <- burnin / keep
RS2 <- R / keep

#事後平均を計算  
beta <- colMeans(BETA[RS1:RS2, ])
theta <- THETA / m
sigma <- mean(SIGMA[RS1:RS2])

#パラメータの真値と比較
round(rbind(beta=beta, betat), 3)
round(rbind(theta=theta, thetat), 3)
round(c(sigma, sigmat), 3)
