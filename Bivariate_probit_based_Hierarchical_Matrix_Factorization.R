#####Bivariate probit based Hierarchical Matrix Factorization#####
options(warn=0)
library(MASS)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(stringr)
library(extraDistr)
library(mvtnorm)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

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
s <- 2   #応答変数数
k <- 10   #基底数
hh <- 10000   #ユーザー数
item <- 3000   #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)
vec <- rep(1, k)

#IDの設定
user_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(1.0, topic))) %*% 1:topic)

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
w <- unlist(lapply(item_list, length))

#個別に和を取るためのスパース行列
user_dt <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(item_id, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


##素性ベクトルを生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
X <- cbind(1, x1, x2, x3)   #データを結合
column <- ncol(X)
index_column <- matrix(1:(column*s), nrow=column, ncol=s)

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
column_u <- ncol(u)
index_k <- matrix(1:(k*s), nrow=k, ncol=s)

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
column_v <- ncol(v)


####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)

  ##階層モデルのパラメータを生成
  #階層モデルの分散パラメータ
  Cov <- Covt <- diag(runif(column*s, 0.01, 0.15))   #素性ベクトルの階層モデルの分散
  Cov_u <- Cov_ut <- diag(runif(k*s, 0.01, 0.2))   #ユーザーの階層モデルの分散
  Cov_v <- Cov_vt <- diag(runif(k, 0.01, 0.2), k)   #アイテムの階層モデルの分散

  #階層モデルの回帰係数を設定
  alpha<- matrix(0, nrow=column_u, ncol=column*s)
  for(j in 1:s){
    alpha[, index_column[, j]] <- cbind(rnorm(column_u, -0.2, 0.5), 
                                        matrix(rnorm((column-1)*column_u, 0, 0.5), nrow=column_u, ncol=column-1))
  }
  alphat <- alpha
  alpha_u <- alpha_ut <- matrix(rnorm(k*column_u*s, 0, 0.5), nrow=column_u, ncol=k*s)
  alpha_v <- alpha_vt <- matrix(rnorm(k*ncol(v), 0, 0.5), nrow=column_v, ncol=k)
  
  
  ##モデルパラメータを生成
  #素性ベクトルと行列分解のパラメータを生成
  beta <- betat <- u %*% alpha + mvrnorm(hh, rep(0, column*s), Cov)
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k*s), Cov_u)
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  #相関行列を生成
  Sigma <- matrix(0.5, nrow=s, ncol=s)
  diag(Sigma) <- 1
  Sigmat <- Sigma
  
  ##潜在効用から応答変数を生成
  #回帰モデルの平均構造
  y <- mu <- matrix(0, nrow=hhpt, ncol=s)
  W <- theta_u[user_id, ]
  H <- theta_v[item_id, ]
  for(j in 1:s){
    mu[, j] <- as.numeric((X * beta[user_id, index_column[, j]]) %*% rep(1, column) + (W[, index_k[, j]] * H) %*% vec)
  }
  
  #多変量正規分布から潜在効用を生成
  er <- mvrnorm(hhpt, rep(0, s), Sigma)   #モデル誤差
  U <- UT <- mu + er   #潜在効用
  y <- ifelse(U > 0, 1, 0)   #応答変数を生成
  
  #ストップ判定
  if(max(colMeans(y)) < 0.4 & min(colMeans(y)) > 0.15){
    break
  }
}

####マルコフ連鎖モンテカルロ法でBivariate probit based Hierarchical Matrix Factorizationを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

##データの設定
#定数の設定
XX_list <- list()
for(i in 1:hh){
  XX_list[[i]] <- t(X[user_list[[i]], ]) %*% X[user_list[[i]], ]
}

#割当のインデックス
xx_index1 <- matrix(1:(s*column), nrow=column, ncol=s, byrow=T)
xx_index2 <- matrix(1:(s*k), nrow=k, ncol=s, byrow=T)

##事前分布の設定
#階層モデルの事前分布
Deltabar <- matrix(0, nrow=column_u, ncol=column*s)
Deltabar_u <- matrix(0, nrow=column_u, ncol=k*s)
Deltabar_v <- matrix(0, nrow=column_v, ncol=k)
ADelta <- ADelta_u <- 0.01 * diag(column_u)
ADelta_v <- 0.01 * diag(column_v)
nu <- nu1 <- nu2 <- k + 1
V <- nu * diag(column*s)
V1 <- nu * diag(k*s)
V2 <- nu * diag(k)

#モデルのの事前分布
m <- column + 1   #逆ウィシャート分布の自由度
Rn <- m * diag(s)   #逆ウィシャート分布の自由度


##パラメータの真値
#階層モデルのパラメータ
Cov <- Covt; inv_Cov <- solve(Cov)
Cov_u <- Cov_ut; inv_Cov_u <- solve(Cov_u)
Cov_v <- Cov_vt; inv_Cov_v <- solve(Cov_v)
alpha <- alphat; alpha_u <- alpha_ut; alpha_v <- alpha_vt
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u ; v_mu <- v %*% alpha_v

#素性ベクトルと行列分解のパラメータ
beta <- betat 
beta_vec <- beta[user_id, ]
beta_mu <- matrix(0, nrow=hhpt, ncol=s)
for(j in 1:s){
beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
}
theta_u <- theta_ut 
theta_v <- theta_vt 

#相関行列をパラメータ
Sigma <- Sigmat
inv_Sigma <- solve(Sigma)

#モデルの平均構造
mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
W <- theta_u[user_id, ]; H <- theta_v[item_id, ]
for(j in 1:s){
WH[, j] <- (W[, index_k[, j]] * H) %*% vec
mu[, j] <- as.numeric(beta_mu[, j]) + WH[, j]
}
U <- UT


##パラメータの初期値
#階層モデルの分散パラメータ
Cov <- 0.01 * diag(column*s); inv_Cov <- solve(Cov)
Cov_u <- 0.01 * diag(k*s); inv_Cov_u <- solve(Cov_u)
Cov_v <- 0.01 * diag(k); inv_Cov_v <- solve(Cov_v)

#階層モデルの回帰係数を設定
alpha <- matrix(rnorm(column*column_u*s, 0, 0.1), nrow=column_u, ncol=column*s)
alpha_u <- matrix(rnorm(k*column_u*s, 0, 0.5), nrow=column_u, ncol=k*s)
alpha_v <- matrix(rnorm(k*ncol(v), 0, 0.5), nrow=column_v, ncol=k)
alpha_mu <- u %*% alpha; u_mu <- u %*% alpha_u ; v_mu <- v %*% alpha_v

#素性ベクトルと行列分解のパラメータを作成
beta <- u %*% alpha + mvrnorm(hh, rep(0, column*s), Cov)
beta_vec <- beta[user_id, ]
beta_mu <- matrix(0, nrow=hhpt, ncol=s)
for(j in 1:s){
  beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
}
theta_u <- u %*% alpha_u + mvrnorm(hh, rep(0, k*s), Cov_u)
theta_v <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)

#相関行列を生成
Sigma <- diag(s)
inv_Sigma <- solve(Sigma)

#モデルの平均構造
mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
W <- theta_u[user_id, ]; H <- theta_v[item_id, ]
for(j in 1:s){
  WH[, j] <- (W[, index_k[, j]] * H) %*% vec
  mu[, j] <- as.numeric(beta_mu[, j]) + WH[, j]
}
U <- (y)*abs(mu) + (-(1-y)*abs(mu))


##サンプリング結果の格納用配列
COV <- array(0, dim=c(column*s, column*s, R/keep))
COV_U <- array(0, dim=c(k*s, k*s, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
ALPHA <- array(0, dim=c(column_u, column*s, R/keep))
ALPHA_U <- array(0, dim=c(column_u, k*s, R/keep))
ALPHA_V <- array(0, dim=c(column_v, k, R/keep))
BETA <- array(0, dim=c(hh, column*s, R/keep))
THETA_U <- array(0, dim=c(hh, k*s, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
SIGMA <- array(0, dim=c(s, s, R/keep))


##対数尤度の基準値
#1パラメータモデルでの対数尤度
Prob <- matrix(colMeans(y), nrow=hhpt, ncol=s, byrow=T)
print(LLst <- sum(y*log(Prob) + (1-y)*log(1-Prob)))

#真値での対数尤度
LLi <- rep(0, s)
mut <- matrix(0, nrow=hhpt, ncol=s)
WT <- theta_ut[user_id, ]; HT <- theta_vt[item_id, ]
for(j in 1:s){
  mut[, j] <- as.numeric((X * betat[user_id, index_column[, j]]) %*% rep(1, column) + (WT[, index_k[, j]] * HT) %*% vec)
  Prob <- pnorm(mut[, j], Sigmat[j, j]); Prob[Prob==0] <- 10^-200; Prob[Prob==1] <- 0.99999999   #応答確率
  LLi[j] <- sum(y[, j]*log(Prob) + (1-y[, j])*log(1-Prob))   #対数尤度
}
print(LLbest <- sum(LLi))


##切断正規分布の切断領域を定義
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)
a0 <- ifelse(y==0, -10^-100, 0)
b0 <- ifelse(y==1, 10^-100, 0)


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  #モデルの平均構造
  mu <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    mu[, j] <- beta_mu[, j] + WH[, j]
  }
  
  #潜在効用を生成
  for(j in 1:s){
    MVR <- cdMVN(mu, Sigma, j, U)   #条件付き分布と分散を生成
    MVR_U <- as.numeric(MVR$CDmu)
    MVR_S <- as.numeric(sqrt(MVR$CDvar))
    U[, j] <- rtnorm(MVR_U, MVR_S, a[, j], b[, j])   #効用を生成
    U[is.infinite(U[, j])==TRUE, j] <- 0
  }

  ##ユーザー別に素性ベクトルをサンプリング
  #モデルの応答変数を設定
  er <- U - WH
  
  for(i in 1:hh){
    #データの抽出
    er_u <- er[user_list[[i]], ]
    x <- X[user_list[[i]], ]
    
    #回帰係数の事後分布のパラメータ
    XX <- XX_list[[i]]; XU <- as.numeric(t(x) %*% er_u %*% inv_Sigma)
    inv_XVX <- solve(kronecker(inv_Sigma, XX) + inv_Cov)
    beta_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov %*% alpha_mu[i, ]))   #回帰ベクトルの期待値
    
    #多変量正規分布から素性ベクトルをサンプリング
    beta[i, ] <- mvrnorm(1, beta_mu, inv_XVX)
  }
  
  #素性ベクトルの平均
  beta_vec <- beta[user_id, ]
  beta_mu <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    beta_mu[, j] <- (X * beta_vec[, index_column[, j]]) %*% rep(1, column)
  }
  
  
  ##ユーザーの特徴行列のパラメータをサンプリング
  #応答変数の設定
  er <- U - beta_mu
  
  for(i in 1:hh){
    #データの抽出
    er_w <- er[user_list[[i]], ]
    x <- H[user_list[[i]], ]
    
    #回帰係数の事後分布のパラメータ
    XX <- t(x) %*% x; XU <- as.numeric(t(x) %*% er_w %*% inv_Sigma)
    inv_XVX <- solve(kronecker(inv_Sigma, XX) + inv_Cov_u)
    w_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov_u %*% u_mu[i, ]))   #回帰ベクトルの期待値
    
    #多変量正規分布から素性ベクトルをサンプリング
    theta_u[i, ] <- mvrnorm(1, w_mu, inv_XVX)
  }
  #ユーザー特徴行列を変換
  W <- theta_u[user_id, ]
  
  
  ##アイテムの特徴行列のパラメータをサンプリング
  #応答変数の設定
  er <- U - beta_mu
  Chol <- chol(inv_Sigma)   #分散共分散行列の逆行列をコレツキー分解
  
  for(j in 1:item){
    #データの抽出
    er_h <- as.numeric(t(Chol %*% t(er[item_list[[j]], ])))
    x <- W[item_list[[j]], ]
    x_vec <- cbind(as.numeric(t(x[, index_k[, 1]])), as.numeric(t(x[, index_k[, 2]])))

    #回帰係数の事後分布のパラメータ
    Chol_x <- matrix(x_vec %*% t(Chol), nrow=w[j]*s, ncol=k, byrow=T)
    XX <- t(Chol_x) %*% Chol_x; XU <- t(Chol_x) %*% er_h
    inv_XVX <- solve(XX + inv_Cov_v)
    h_mu <- as.numeric(inv_XVX %*% (XU + inv_Cov_v %*% v_mu[j, ]))   #回帰ベクトルの期待
    
    #多変量正規分布から素性ベクトルをサンプリング
    theta_v[j, ] <- mvrnorm(1, h_mu, inv_XVX)
  }
  #アイテム特徴行列を変換
  H <- theta_v[item_id, ]
  
  
  ##逆ウィシャート分布から相関行列をサンプリング
  #モデルの平均構造
  mu <- WH <- matrix(0, nrow=hhpt, ncol=s)
  for(j in 1:s){
    WH[, j] <- (W[, index_k[, j]] * H) %*% vec
    mu[, j] <- beta_mu[, j] + WH[, j]
  }
  
  #逆ウィシャート分布のパラメータ
  R_error <- U - mu
  IW_R <- t(R_error) %*% R_error + Rn
  Sn <- hhpt + m
  
  #逆ウィシャート分布からパラメータをサンプリング
  Sigma_hat <- rwishart(Sn + 1000, solve(IW_R))$IW
  Sigma <- cov2cor(Sigma_hat)   #相関行列に変換
  inv_Sigma <- solve(Sigma)
  Sigma
  
  ##階層モデルのパラメータのサンプリング
  #多変量回帰モデルから素性ベクトルの階層モデルのパラメータをサンプリング
  out <- rmultireg(beta, u, Deltabar, ADelta, nu, V)
  alpha <- out$B
  alpha_mu <- u %*% alpha   #素性ベクトルの特徴行列の平均構造
  Cov <- diag(diag(out$Sigma))
  inv_Cov <- solve(Cov)
  
  #多変量回帰モデルからユーザー特徴行列階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_u, u, Deltabar_u, ADelta_u, nu1, V1)
  alpha_u <- out$B
  u_mu <- u %*% alpha_u   #ユーザー特徴行列の平均構造
  Cov_u <- diag(diag(out$Sigma))
  inv_Cov_u <- solve(Cov_u)
  
  #多変量回帰モデルからアイテム特徴行列階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu2, V2)
  alpha_v <- out$B
  v_mu <- v %*% alpha_v   #アイテム特徴行列の平均構造
  Cov_v <- diag(diag(out$Sigma))
  inv_Cov_v <- solve(Cov_v)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    COV[, , mkeep] <- Cov
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    ALPHA[, , mkeep] <- alpha
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    BETA[, , mkeep] <- beta
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    SIGMA[, , mkeep] <- Sigma
  }
  
  ##対数尤度関数の計算とサンプリング結果の確認
  if(rp%%disp==0){
    LLi <- rep(0, s)
    for(j in 1:s){
      Prob <- pnorm(mu[, j], Sigma[j, j]); Prob[Prob==0] <- 10^-200; Prob[Prob==1] <- 0.99999999   #応答確率
      LLi[j] <- sum(y[, j]*log(Prob) + (1-y[, j])*log(1-Prob))   #対数尤度
    }
    LL <- sum(LLi)
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLbest, LLst))
    print(Sigma)
    print(round(rbind(diag(Cov_v), diag(Cov_vt)), 3))
  }
}
