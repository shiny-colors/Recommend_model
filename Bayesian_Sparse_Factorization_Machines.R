#####Bayesian Sparse Factorization Machines#####
library(MASS)
library(Matrix)
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
item <- 3000   #アイテム数
tag <- 150   #タグ数
pt <- rtpois(item, rgamma(item, 30.0, 0.225), a=1, b=Inf)   #購買接触数
N <- sum(pt)
n <- rtpois(item, 1.0, a=0, b=5)   #タグ数
k <- 7   #基底数
vec_k <- rep(1, k)

#IDを設定
item_id <- rep(1:item, pt)
pt_id <- as.numeric(unlist(tapply(1:N, item_id, rank)))
ID <- data.frame(no=1:N, id=item_id, t=pt_id)   #データの結合
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}

##階層モデルの説明変数を設定
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

##タグを生成
#パラメータの設定
topic <- 25
omega <- extraDistr::rdirichlet(topic, rep(0.5, tag))
z <- as.numeric(rmnom(item, 1,  extraDistr::rdirichlet(item, rep(1.0, topic))) %*% 1:topic)

#多項分布からタグを生成
max_n <- max(n)
tag_list <- list()
tag_id <- matrix(0, nrow=N, ncol=max_n)
tag_data <- matrix(0, nrow=N, ncol=tag); storage.mode(tag_data) <- "integer"

for(j in 1:item){
  repeat { 
    x <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(x)==1){
      tag_list[[j]] <- (x * 1:tag)[x > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
  tag_data[item_list[[j]], ] <- matrix(x, nrow=length(item_list[[j]]), ncol=tag, byrow=T)
}
tag_id0 <- tag_id
tag_id0[tag_id0==0] <- tag+1

#組み合わせを作成
index_combine <- t(combn(c(1, 2, 3, 4, 5), m=2))
combine_list <- list()
combine_n <- rep(0, max(index_combine[, 1]))
for(j in 1:max(index_combine[, 1])){
  combine_list[[j]] <- index_combine[which(index_combine[, 1]==j), 2]
  combine_n[j] <- length(combine_list[[j]])
}

#交差項のインデックスを作成
index_n <- which(rowSums(tag_id > 0) >= 2)
tag_n <- length(index_n)
tag_dt <- sparseMatrix(rep(1:tag_n, max_n), 1:(tag_n*max_n), x=rep(1, tag_n*max_n), dims=c(tag_n, tag_n*max_n))


####応答変数を生成####
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  ##モデルのパラメータ
  beta <- betat <- 5.5
  sigma <- sigmat <- 0.75
  
  ##アイテムの階層モデルのパラメータ
  #標準偏差を設定
  tau_v <- tau_vt <- runif(1, 0.3, 0.75)
  
  #回帰係数を設定
  alpha_v <- rep(0, ncol(v))
  for(j in 1:ncol(v)){
    alpha_v[j] <- runif(1, -1.25, 1.25)
  }
  alpha_vt <- alpha_v
  
  #回帰モデルからアイテム個別の変量効果を生成
  theta_v <- theta_vt <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
  theta_vec1 <- theta_v[item_id]
  
  ##タグのパラメータを生成
  #タグ個別の変量効果を生成
  tau_r <- tau_rt <- runif(1, 0.25, 0.5)
  theta_r <- theta_rt <- rnorm(tag, 0, tau_r)
  theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
  
  ##交互作用のパラメータを生成
  #特徴ベクトルのパラメータを生成
  tau_g <- tau_gt <- runif(k, 0.1, 0.4) * diag(k)
  theta_g <- theta_gt <- mvrnorm(tag, rep(0, k), tau_g)
  theta_g0 <- rbind(theta_g, 0)

  #交互作用のベクトルを生成
  WH <- rep(0, N)
  for(j in 1:length(combine_n)){
    W <- theta_g0[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
  }
  
  ##正規分布から評価ベクトルを生成
  mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 
  y0 <- rnorm(N, mu, sigma)   #応答変数を生成

  #break条件
  if(min(y0) > -3.0 & max(y0) < 14.0 & mean(y0) > 4.5 & mean(y0) < 6.0){
    break
  }
}

#応答変数を1～10に変換する
y <- round(y0)
y[y > 10] <- 10; y[y < 1] <- 1
hist(y0, breaks=25, col="grey", main="真のスコア分布", xlab="スコア")
hist(y, breaks=25, col="grey", main="切断されたスコア分布", xlab="スコア")


####マルコフ連鎖モンテカルロ法でBayesian SFMを推定####
##アルゴリズムの設定
R <- 2000
burnin <- 500
keep <- 2
disp <- 10
iter <- 0

#データとインデックスを設定
#インデックスを設定
tag_list1 <- tag_list2 <- dt_list <- list()
tag_n1 <- tag_n2 <- rep(0, tag)
for(i in 1:tag){
  index1 <- index2 <- c()
  for(j in 1:max_n){
    index1 <- c(index1, which(tag_id[, j]==i))
    index2 <- c(index2, index_n[which(tag_id[index_n, j]==i)])
  }
  tag_list1[[i]] <- sort(index1)
  tag_list2[[i]] <- sort(index2)
  tag_n1[i] <- length(tag_list1[[i]])
  tag_n2[i] <- length(tag_list2[[i]])
  dt_list[[i]] <- sparseMatrix(rep(1:tag_n2[i], max_n), 1:(tag_n2[i]*max_n), x=rep(1, tag_n2[i]*max_n), 
                               dims=c(tag_n2[i], tag_n2[i]*max_n))
}


##事前分布の設定
#回帰パラメータの事前分布
alpha1 <- 0
alpha2 <- rep(0, ncol(v))
alpha3 <- rep(0, k)

#分散の事前分布
s0 <- 0.1
v0 <- 0.1   
nu <- 1   #逆ウィシャート分布の自由度
V <- 0.01 * diag(k)    #逆ウィシャート分布のパラメータ
tau1 <- 100
tau2 <- 100 * diag(ncol(v))
inv_tau2 <- solve(tau2)
tau3 <- 100 * diag(k)
inv_tau3 <- solve(tau3)

##真値の設定
#モデルパラメータ
beta <- betat
sigma <- sigmat

#階層モデルのパラメータ
alpha_v <- alpha_vt
tau_v <- tau_vt
tau_r <- tau_rt
tau_g <- tau_gt
inv_tau_g <- solve(tau_g)

#変量効果のパラメータ
theta_v <- theta_vt
theta_r <- theta_rt
theta_g <- theta_gt
theta_g0 <- rbind(theta_g, 0)

#評価ベクトルの期待値
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
y_mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 

##初期値の設定
#モデルパラメータ
beta <- mean(y)
sigma <- 1.0

#階層モデルのパラメータ
alpha_v <- rep(0, ncol(v))
tau_v <- 0.1
tau_r <- 0.1
tau_g <- 0.2 * diag(k)

#変量効果のパラメータ
theta_v <- as.numeric(v %*% alpha_v) + rnorm(item, 0, tau_v)
theta_r <- rnorm(tag, 0, tau_r)
theta_g <- mvrnorm(tag, rep(0, k), tau_g)
theta_g0 <- rbind(theta_g, 0)

#評価ベクトルの期待値
theta_vec1 <- theta_v[item_id]
theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
WH <- rep(0, N)
for(j in 1:length(combine_n)){
  W <- theta_g0[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
}
mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 


##パラメータの格納用配列
#モデルパラメータ
BETA <- rep(0, R/keep)
SIGMA <- rep(0, R/keep)

#階層モデルのパラメータ
ALPHA_V <- matrix(0, nrow=R/keep, ncol=ncol(v))
TAU_V <- rep(0, R/keep)
TAU_R <- rep(0, R/keep)
TAU_G <- array(0, dim=c(k, k, R/keep))

#変量効果のパラメータ
THETA_V <- matrix(0, nrow=R/keep, ncol=item)
THETA_R <- matrix(0, nrow=R/keep, ncol=tag)
THETA_G <- array(0, dim=c(tag, k, R/keep))


##対数尤度の基準値
LLst <- sum(dnorm(y, mean(y), sd(y), log=TRUE))
LLbest <- sum(dnorm(y, y_mu, sigmat, log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##モデル平均をサンプリング
  #モデル誤差を設定
  er <- y - theta_vec1 - theta_vec2 - WH 
  
  #正規分布のパラメータ
  omega <- (N/sigma^2) / (1/tau1 + N/sigma^2)   #重み係数
  beta_mu <- omega * mean(er)   #正規分布の平均
  cov <- 1 / (1/tau1 + N/sigma^2)   #正規分布の標準偏差
  
  #正規分布からパラメータをサンプリング
  beta <- rnorm(1, beta_mu, cov)   
  
  ##モデルの標準偏差をサンプリング
  #モデル誤差を設定
  er <- y - beta - theta_vec1 - theta_vec2 - WH
  
  #逆ガンマ分布のパラメータ
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- N + v0
  
  #逆ガンマ分布から標準偏差をサンプリング
  sigma <- sqrt(1/rgamma(1, v1/2, s1/2))
  
  
  ##アイテムの変量効果をサンプリング
  #モデル誤差を設定
  er <- y - beta - theta_vec2 - WH
  alpha_mu <- as.numeric(v %*% alpha_v)
  
  #正規分布のパラメータ
  omega <- (pt/sigma^2) / (1/tau_v^2 + pt/sigma^2)   #重み係数
  cov <- 1 / (1/tau_v^2 + pt/sigma^2)   #正規分布の標準偏差
  theta_mu <- rep(0, item)
  for(j in 1:item){
    theta_mu[j] <- (1-omega[j])*alpha_mu[j] + omega[j]*mean(er[item_list[[j]]])   #正規分布の平均
  }
  
  #正規分布からパラメータをサンプリング
  theta_v <- rnorm(item, theta_mu, cov)
  theta_vec1 <- theta_v[item_id]

  
  ##アイテムタグの変量効果をサンプリング
  #モデル誤差を設定
  er <- y - beta - theta_vec1 - WH

  #正規分布のパラメータ
  omega <- (tag_n1/sigma^2) / (1/tau_r^2 + tag_n1/sigma^2)   #重み係数
  cov <- 1 / (1/tau_r^2 + tag_n1/sigma^2)   #正規分布の標準偏差
  theta_mu <- rep(0, tag)
  for(j in 1:tag){
    #サンプリング対象外のアイテムタグの変量効果
    theta_r0 <- c(theta_r, 0); theta_r0[j] <- 0
    r <- as.numeric(matrix(theta_r0[tag_id0[tag_list1[[j]], ]], nrow=tag_n1[j], ncol=max_n) %*% rep(1, max_n))
    
    #正規分布の平均
    er_y <- er[tag_list1[[j]]] - r
    theta_mu[j] <- omega[j]*mean(er_y)   
  }

  #正規分布からパラメータをサンプリング
  theta_r <- rnorm(tag, theta_mu, cov)
  theta_vec2 <- as.numeric(matrix(c(theta_r, 0)[tag_id0], nrow=N, ncol=max_n) %*% rep(1, max_n))
  
  
  ##交互作用の特徴ベクトルをサンプリング
  #モデル誤差を設定
  er <- y - beta - theta_vec1 - theta_vec2
  for(i in 1:tag){
    #データの設定
    theta_g0 <- rbind(theta_g, 0); theta_g0[i, ] <- 0
    index <- tag_list2[[i]]
    tag_vec <- tag_id0[index, ]
    
    #サンプリング対象外の交互作用ベクトルのパラメータ
    wh <- rep(0, length(index))
    for(j in 1:length(combine_n)){
      w <- theta_g0[tag_vec[, j], ]
      h <- as.matrix(dt_list[[i]][, 1:(length(index)*combine_n[j])] %*% theta_g0[tag_vec[, combine_list[[j]]], ])
      wh <- wh + as.numeric((w * h) %*% vec_k)
    }
    
    #多変量正規分布のパラメータ
    er_y <- er[index] - wh
    x <- as.matrix(dt_list[[i]] %*% theta_g0[tag_vec, ])
    inv_xxv <- solve(t(x) %*% x + inv_tau_g)
    theta_mu <- as.numeric(inv_xxv %*% t(x) %*% er_y)   #多変量正規分布の平均ベクトル
    cov <- sigma^2 * inv_xxv   #多変量正規分布の分散
  
    #多変量正規分布からパラメータをサンプリング
    theta_g[i, ] <- mvrnorm(1, theta_mu, cov)
  }
    
  #交互作用ベクトルの期待値
  theta_g0 <- rbind(theta_g, 0)
  WH <- rep(0, N)
  for(j in 1:length(combine_n)){
    W <- theta_g0[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% theta_g0[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((H * W) %*% vec_k)
  }
  
  ##アイテムの変量効果の階層モデルをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  inv_xxv <- solve(t(v) %*% v + inv_tau2) 
  mu <- as.numeric(inv_xxv %*% t(v) %*% theta_v)
  alpha_v <- mvrnorm(1, mu, tau_v*inv_xxv)   #多変量正規分布からパラメータをサンプリング
  alpha_mu <- as.numeric(v %*% alpha_v)
  
  #逆ガンマ分布から分散をサンプリング
  er <- theta_v - alpha_mu
  s1 <- as.numeric(t(er) %*% er) + s0
  v1 <- item + v0
  tau_v <- sqrt(1/rgamma(1, v1/2, s1/2))
  inv_tau_v <- solve(tau_v)
  
  ##タグの変量効果の階層モデルをサンプリング
  #逆ガンマ分布から分散をサンプリング
  s1 <- as.numeric(t(theta_r) %*% theta_r) + s0
  v1 <- tag + v0
  tau_r <- sqrt(1/rgamma(1, v1/2, s1/2))
  inv_tau_r <- solve(tau_r)
  
  ##タグの特徴ベクトルの階層モデルをサンプリング
  #逆ウィシャート分布から分散をサンプリング
  IW_R <- t(theta_g) %*% theta_g + V
  Sn <- tag + nu
  tau_g <- diag(diag(rwishart(Sn, solve(IW_R))$IW))
  inv_tau_g <- solve(tau_g)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #パラメータを格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    
    #モデルパラメータを格納
    BETA[mkeep] <- beta
    SIGMA[mkeep] <- sigma
    
    #階層モデルのパラメータを格納
    ALPHA_V[mkeep, ] <- alpha_v
    TAU_V[mkeep] <- tau_v
    TAU_R[mkeep] <- tau_r
    TAU_G[, , mkeep] <- tau_g
    
    #変量効果のパラメータを格納
    THETA_V[mkeep, ] <- theta_v
    THETA_R[mkeep, ] <- theta_r
    THETA_G[, , mkeep] <- theta_g
  }
  
  if(rp%%disp==0){
    #対数尤度を算出
    mu <- beta + theta_vec1 + theta_vec2 + WH   #期待値を設定 
    LL <- sum(dnorm(y, mu, sigma, log=TRUE))
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(c(beta, betat), 3))
    print(c(sigma, sigmat))
  }
}


####推定結果の確認と適合度####
##サンプリング結果の可視化
plot(1:(R/keep), BETA, type="l", xlab="サンプリング回数", main="betaのサンプリング結果のプロット")
matplot(THETA_V[, 1:10], type="l")
matplot(t(THETA_G[100, , ]), type="l")
t(THETA_G[100, , ])

plot(1:length(SIGMA), SIGMA, type="l", xlab="サンプリング回数", main="sigmaのサンプリング結果のプロット")
THETA_V[, 1]

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


