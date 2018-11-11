#####Bias Smoothed Hierarchical Probit Regression model#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(bayesm)
library(MCMCpack)
library(condMVNorm)
library(extraDistr)
library(reshape2)
library(actuar)
library(extraDistr)
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
hh <- 5000   #ユーザー数
item <- 2000   #アイテム数
context <- 200   #コンテキスト数
N0 <- hh*item


##欠損ベクトルとIDを設定
#IDを仮設定
item_id0 <- rep(1:item, rep(context, item))
context_id0 <- rep(1:context, item)
n <- length(item_id0)

#要素ごとの出現確率
beta1 <- rbeta(hh, 3.0, 36.0)
par2 <- rbeta(item, 3.0, 15.0)
par3 <- rbeta(context, 0.6, 8.0)


#ユーザーごとにIDを作成
user_id_list <- item_id_list <- context_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  #確率を生成
  beta2 <- rbeta(item, par2*2.0, (1-par2)*2.0)
  beta3 <- rbeta(context, par3*1.5, (1-par3)*1.5)
  beta_vec2 <- beta2[item_id0]
  beta_vec3 <- beta3[context_id0]
  
  #欠損ベクトルを生成
  prob <- beta1[i] * beta_vec2 * beta_vec3
  deficit <- rbinom(n, 1, prob)
  index_z <- which(deficit==1)
  
  #IDを設定
  user_id_list[[i]] <- rep(i, n)[index_z]
  item_id_list[[i]] <- item_id0[index_z]
  context_id_list[[i]] <- context_id0[index_z]
}
#リストを変換
user_id <- unlist(user_id_list)
item_id <- unlist(item_id_list)
context_id <- unlist(context_id_list)
N <- length(user_id)


#ユーザー×コンテキストのID
uw_index <- paste(user_id, context_id, sep="-")
uw_id <- left_join(data.frame(id=uw_index, no_vec=1:length(uw_index)),
                   data.frame(id=unique(uw_index), no=1:length(unique(uw_index))), by="id")$no

#アイテム×コンテキストのID
vw_index <- paste(item_id, context_id, sep="-")
vw_id <- left_join(data.frame(id=vw_index, no_vec=1:length(vw_index)),
                   data.frame(id=unique(vw_index), no=1:length(unique(vw_index))), by="id")$no


##応答変数が妥当になるまでパラメータの生成を繰り返す
for(rp in 1:1000){
  print(rp)
  
  ##素性ベクトルを生成
  k1 <- 2; k2 <- 3; k3 <- 4
  x1 <- matrix(runif(N*k1, 0, 1), nrow=N, ncol=k1)
  x2 <- matrix(0, nrow=N, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    x2[, j] <- rbinom(N, 1, pr)
  }
  x3 <- rmnom(N, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
  x <- cbind(x1, x2, x3)   #データを結合
  
  
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
  
  #コンテキストの説明変数
  k1 <- 2; k2 <- 2; k2 <- 3
  w1 <- matrix(runif(item*k1, 0, 1), nrow=context, ncol=k1)
  w2 <- matrix(0, nrow=context, ncol=k2)
  for(j in 1:k2){
    pr <- runif(1, 0.25, 0.55)
    w2[, j] <- rbinom(context, 1, pr)
  }
  w3 <- rmnom(context, 1, runif(k2, 0.2, 1.25)); w3 <- w3[, -which.min(colSums(w2))]
  w <- cbind(1, w1, w2, w3)   #データを結合
  
  
  #素性ベクトルの回帰係数を生成
  beta <- rep(0, ncol(x))
  for(j in 1:ncol(x)){
    beta[j] <- runif(1, -0.8, 1.2)
  }
  betat <- beta
  
  
  ##階層モデルのパラメータを生成
  ##ユーザーベースの階層モデルのパラメータ
  tau_u <- tau_ut <- 0.5   #標準偏差
  
  #回帰係数を設定
  alpha_u <- rep(0, ncol(u))
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u[j] <- runif(1, -1.3, -0.5)
    } else {
      alpha_u[j] <- runif(1, -0.4, 0.6)
    }
  }
  alpha_ut <- alpha_u
  
  #回帰モデルからユーザー個別の回帰パラメータを生成
  theta_ut <- theta_u <- as.numeric(u %*% alpha_u + rnorm(hh, 0, tau_u))
  
  
  ##アイテムベースの階層モデルのパラメータ
  tau_v <- tau_vt <- 0.7   #標準偏差
  
  #回帰係数を設定
  alpha_v <- rep(0, ncol(v))
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j] <- runif(1, -1.2, -0.3)
    } else {
      alpha_v[j] <- runif(1, -0.6, 0.7)
    }
  }
  alpha_vt <- alpha_v
  
  #回帰モデルからアイテム個別の回帰パラメータを生成
  theta_vt <- theta_v <- as.numeric(v %*% alpha_v + rnorm(item, 0, tau_v))
  
  
  ##コンテキストベースの階層モデルのパラメータ
  tau_w <- tau_wt <- 0.4   #標準偏差
  
  #回帰係数を設定
  kw <- 2
  alpha_w <- matrix(0, nrow=ncol(w), ncol=kw)
  for(j in 1:ncol(w)){
    if(j==1){
      alpha_w[j, ] <- runif(kw, 0.2, 0.5)
    } else {
      alpha_w[j, ] <- runif(kw, 0.2, 0.7)
    }
  }
  alpha_wt <- alpha_w
  
  #回帰モデルからコンテキスト個別の回帰パラメータを生成
  theta_wt <- theta_w <- w %*% alpha_w + mvrnorm(context, rep(0, kw), tau_w^2 * diag(kw))
  theta_wt1 <- theta_w1 <- theta_w[, 1]
  theta_wt2 <- theta_w2 <- theta_w[, 2]
  
  ##コンテキスト依存のユーザーおよびアイテムバイアスのパラメータを生成
  #コンテキスト依存の変量効果
  tau_uwt <- tau_uw <- 0.4; tau_vwt <- tau_vw <- 0.4
  alpha_uwt <- alpha_uw <- rnorm(unique(uw_id), 0, tau_uw)
  alpha_vwt <- alpha_vw <- rnorm(unique(vw_id), 0, tau_vw)
  
  #コンテキスト依存バイアス平滑化パラメータ
  theta_uwt <- theta_uw <- alpha_uw[uw_id] + theta_u[user_id] * theta_w1[context_id] 
  theta_vwt <- theta_vw <- alpha_vw[vw_id] + theta_v[item_id] * theta_w2[context_id] 
  
  
  ##プロビットモデルから応答変数を生成
  #潜在効用の生成
  mu <- x %*% beta + theta_uw + theta_vw
  U <- rnorm(N, mu, 1)
  
  #応答変数を生成
  y <- as.numeric(U > 0)
  if(mean(y) > 0.3 & mean(y) < 0.4) break   #break条件
}

#####マルコフ連鎖モンテカルロ法でBSHPモデルを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b, L){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- matrix(0, nrow=length(mu), ncol=L)
  for(j in 1:L){
    par[, j] <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  }
  return(par)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 5

##インデックスを作成
user_index <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
}
item_index <- list()
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
}
context_index <- list()
for(j in 1:context){
  context_index[[j]] <- which(context_id==j)
}

#コンテキスト依存ユーザーインデックス
n_uw <- rep(0, length(unique(uw_id)))
uw_index <- list()
context_u <- rep(0, length(n_uw))
context_w1 <- rep(0, length(n_uw))

for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  id <- uw_id[user_index[[i]]]
  min_id <- min(id); max_id <- max(id)
  for(j in min_id:max_id){
    uw_index[[j]] <- user_index[[i]][id==j]
    context_u[j] <- user_id[uw_index[[j]]][1]
    context_w1[j] <- context_id[uw_index[[j]]][1]
    n_uw[j] <- length(uw_index[[j]])
  }
}
N_uw <- length(n_uw)

#コンテキスト依存アイテムインデックス
n_vw <- rep(0, length(unique(vw_id)))
vw_index <- list()
context_v <- rep(0, length(n_vw))
context_w2 <- rep(0, length(n_vw))

for(i in 1:item){
  if(i%%100==0){
    print(i)
  }
  id <- vw_id[item_index[[i]]]
  unique_id <- unique(id)
  for(j in 1:length(unique_id)){
    index <- unique_id[j]
    vw_index[[index]] <- item_index[[i]][id==index]
    context_v[index] <- item_id[vw_index[[index]]][1]
    context_w2[index] <- context_id[vw_index[[index]]][1]
    n_vw[index] <- length(vw_index[[index]])
  }
}
N_vw <- length(n_vw)


##階層モデルのインデックス
#ユーザーインデックス
user_index <- list()
user_n <- rep(0, hh)
for(i in 1:hh){
  user_index[[i]] <- which(context_u==i)
  user_n[i] <- length(user_index[[i]])
}
#アイテムインデックス
item_index <- list()
item_n <- list()
for(j in 1:item){
  item_index[[j]] <- which(context_v==j)
  item_n[j] <- length(item_index[[j]])
}
#コンテキストインデックス
context_index1 <- context_index2 <- list()
context_n1 <- context_n2 <- rep(0, context)
for(j in 1:context){
  context_index1[[j]] <- which(context_w1==j)
  context_index2[[j]] <- which(context_w2==j)
  context_n1[j] <- length(context_index1[[j]])
  context_n2[j] <- length(context_index2[[j]])
}

##データの設定
#素性ベクトルの説明変数の設定
xx <- t(x) %*% x
inv_xx <- solve(xx)

#ユーザーの階層モデルの説明変数の設定
uu <- t(u) %*% u
inv_uu <- solve(uu)

#アイテムの階層モデルの説明変数の設定
vv <- t(v) %*% v
inv_vv <- solve(vv)

#コンテキストの階層モデルの説明変数の設定
ww <- t(w) %*% w
inv_ww <- solve(ww)


##事前分布を設定
#逆ガンマ分布の事前分布
s0 <- 1
v0 <- 1

#素性回帰ベクトルの事前分布
tau1 <- diag(100, ncol(x))
tau_inv1 <- solve(tau1)
mu1 <- rep(0, ncol(x))

#ユーザーの階層モデルの事前分布
tau2 <- diag(100, ncol(u))
tau_inv2 <- solve(tau2)
mu2 <- rep(0, ncol(u))
s01 <- 100; v01 <- 1

#アイテムの階層モデルの事前分布
tau3 <- diag(100, ncol(v))
tau_inv3 <- solve(tau3)
mu3 <- rep(0, ncol(v))
s02 <- 1; v02 <- 1
 
#コンテキストの階層モデルの事前分布
Deltabar <- matrix(rep(0, 2*ncol(w)), nrow=ncol(w), ncol=2)   #階層モデルの回帰係数の事前分布の分散
ADelta <- 0.01 * diag(rep(1, ncol(w)))   #階層モデルの回帰係数の事前分布の分散
nu <- 1   #逆ウィシャート分布の自由度
V <- nu * diag(rep(1, 2)) #逆ウィシャート分布のパラメータ
s03 <- 1; v03 <- 1


##切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)


##パラメータの真値
#素性ベクトルの回帰係数
sigma <- 1.0
beta <- solve(t(x) %*% x) %*% t(x) %*% y

#変量効果のパラメータ
theta_u <- theta_ut   #ユーザーの変量効果
theta_v <- theta_vt   #アイテムの変量効果
theta_w1 <- theta_wt[, 1]   #コンテキストの変量効果
theta_w2 <- theta_wt[, 2]
theta_uw <- theta_uwt   #コンテキスト依存のユーザーの変量効果
theta_vw <- theta_vwt   #コンテキスト依存のアイテムの変量効果

#変量効果のパラメータをベクトル化
theta_u_vec <- theta_u[context_u]
theta_w1_vec <- theta_w1[context_w1]
theta_v_vec <- theta_v[context_v]
theta_w2_vec <- theta_w2[context_w2]

#階層モデルのパラメータを生成
tau_u <- tau_ut   #ユーザーの階層モデルの標準偏差
alpha_u <- alpha_ut   #ユーザーの階層モデルの回帰係数
tau_v <- tau_vt   #アイテムの階層モデルの標準偏差
alpha_v <- alpha_vt   #アイテムの階層モデルの回帰係数
tau_w <- tau_wt   #コンテキストの階層モデルの標準偏差
alpha_w <- alpha_wt   #コンテキストの階層モデルの回帰係数
tau_uw <- tau_uwt   #コンテキスト依存のユーザーバイアスの標準偏差
tau_vw <- tau_vwt   #コンテキスト依存のユーザーバイアスの標準偏差

#階層モデルの回帰モデルの平均構造
u_mu <- as.numeric(u %*% alpha_u)
v_mu <- as.numeric(v %*% alpha_v)
w_mu <- w %*% alpha_w

#モデルの変量効果
theta_uw <- rep(0, N_uw)
theta_vw <- rep(0, N_vw)
for(i in 1:N_uw){
  theta_uw[i] <- theta_uwt[uw_index[[i]]][1]
}
for(i in 1:N_vw){
  theta_vw[i] <- theta_vwt[vw_index[[i]]][1]
}
theta_vwt0 <- theta_vw; theta_uwt0 <- theta_uw
theta_uw_vec <- theta_uw[uw_id]; theta_vw_vec <- theta_vw[vw_id]


##初期値の設定
#素性ベクトルの回帰係数
sigma <- 1.0
beta <- solve(t(x) %*% x) %*% t(x) %*% y 

#階層モデルのパラメータを生成
tau_u <- tau_v <- tau_w <- tau_uw <- tau_vw <- 0.75
tau_w <- 0.4
alpha_u <- runif(ncol(u), -0.1, 0.1)   #ユーザーの階層モデルの回帰係数
alpha_v <- runif(ncol(v), -0.1, 0.1)   #アイテムの階層モデルの回帰係数差
alpha_w <- matrix(runif(ncol(w)*2, -0.1, 0.1), nrow=ncol(w), ncol=2)   #コンテキストの階層モデルの回帰係数

#階層モデルの回帰モデルの平均構造
u_mu <- as.numeric(u %*% alpha_u)
v_mu <- as.numeric(v %*% alpha_v)
w_mu <- w %*% alpha_w

#変量効果のパラメータ
theta_u <- rnorm(hh, 0, tau_u)
theta_v <- rnorm(item, 0, tau_v)
theta_w1 <- rnorm(context, 0, tau_w); theta_w2 <- rnorm(context, 0, tau_w)
theta_uw <- rnorm(N_uw, 0, tau_uw)
theta_vw <- rnorm(N_vw, 0, tau_vw)
theta_uw_vec <- theta_uw[uw_id]
theta_vw_vec <- theta_vw[vw_id]

#変量効果のパラメータをベクトル化
theta_u_vec <- theta_u[context_u]
theta_w1_vec <- theta_w1[context_w1]
theta_v_vec <- theta_v[context_v]
theta_w2_vec <- theta_w2[context_w2]


##パラメータの格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=ncol(x))
ALPHA_U <- matrix(0, nrow=R/keep, ncol=ncol(u))
ALPHA_V <- matrix(0, nrow=R/keep, ncol=ncol(v))
ALPHA_W <- matrix(0, nrow=R/keep, ncol=ncol(w)*2)
THETA_UW <- rep(0, N_uw)
THETA_VW <- rep(0, N_vw)
THETA_U <- matrix(0, nrow=R/keep, ncol=hh)
THETA_V <- matrix(0, nrow=R/keep, ncol=item)
THETA_W <- matrix(0, nrow=R/keep, ncol=context*2)
COV <- matrix(0, nrow=R/keep, ncol=length(c(tau_uw, tau_vw, tau_u, tau_v, tau_w)))
rkeep <- c()


##対数尤度の基準値
beta_mu <- as.numeric(x %*% beta)
mu <- beta_mu + theta_uw_vec + theta_vw_vec   #潜在効用の期待値
prob <- pnorm(mu, 0, sigma)   #購買確率
LL1 <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #対数尤度
LLst <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
print(c(LL1, LLst))



####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  beta_mu <- as.numeric(x %*% beta)
  mu <- beta_mu + theta_uw_vec + theta_vw_vec   #潜在効用の期待値
  U <- as.numeric(rtnorm(mu, sigma, a, b, 1))   #潜在効用を生成
  
  ###変量効果のパラメータをサンプリング
  ##コンテキスト依存のユーザー変量効果の期待値をサンプリング
  #データの設定
  uw_er <- U - beta_mu - theta_vw_vec   #誤差を設定
  
  #階層モデルのパラメータ
  theta_vec <- theta_u_vec * theta_w1_vec 
  
  #事後分布のパラメータを設定
  uw_mu <- rep(0, N_uw)
  for(i in 1:N_uw){
    uw_mu[i] <- mean(uw_er[uw_index[[i]]])
  }
  weights <- tau_uw^2 / (sigma^2/n_uw + tau_uw^2)    #重み係数
  mu_par <- weights*uw_mu + (1-weights)*theta_vec   #事後分布の平均
  tau <- sqrt(1 / (1/tau_uw^2 + n_uw/sigma^2))
  
  #正規分布より事後分布をサンプリング
  theta_uw <- rnorm(N_uw, mu_par, tau)
  theta_uw_vec <- theta_uw[uw_id]
  
  ##コンテキスト依存のユーザー変量効果の分散をサンプリング
  #逆ガンマ分布より分散をサンプリング
  s1 <- s0 + sum((theta_uw - theta_vec)^2)
  v1 <- v0 + N_uw
  tau_uw <- sqrt(1/(rgamma(1, v1/2, s1/2)))
  

  ##コンテキスト依存のアイテム変量効果をサンプリング
  #データの設定
  vw_er <- U - beta_mu - theta_uw_vec   #誤差を設定
  
  #階層モデルのパラメータ
  theta_vec <- theta_v_vec * theta_w2_vec 
  
  #事後分布のパラメータを設定
  vw_mu <- rep(0, N_vw)
  for(i in 1:N_vw){
    vw_mu[i] <- mean(vw_er[vw_index[[i]]])
  }
  weights <- tau_vw^2 / (sigma/n_vw + tau_vw^2)   #重み係数
  mu_par <- weights*vw_mu + (1-weights)*theta_vec   #事後分布の平均
  tau <- sqrt(1 / (1/tau_vw^2 + n_vw/sigma^2))
  
  #正規分布より事後分布をサンプリング
  theta_vw <- rnorm(N_vw, mu_par, tau)
  theta_vw_vec <- theta_vw[vw_id]
  
  ##コンテキスト依存のユーザー変量効果の分散をサンプリング
  #逆ガンマ分布より分散をサンプリング
  s1 <- s0 + sum((theta_vw - theta_vec)^2)
  v1 <- v0 + N_vw
  tau_vw <- sqrt(1/(rgamma(1, v1/2, s1/2)))

  
  ##ユーザー変量効果をサンプリング
  #ユーザーごとに事後分布のパラメータを設定
  inv_tau_u <- 1/tau_u^2
  mu_par <- rep(0, hh)
  sigma_par <- rep(0, hh)
  
  for(i in 1:hh){
    index <- context_w1[user_index[[i]]]
    X <- theta_w1[index]
    Xy <- (t(X) %*% theta_uw[user_index[[i]]])
    XXV <- (t(X) %*% X) + inv_tau_u
    inv_XXV <- 1/XXV
    sigma_par[i] <- tau_uw^2 * inv_XXV
    mu_par[i] <- inv_XXV %*% (Xy + inv_tau_u %*% u_mu[i])
  }
  
  #正規分布から事後分布をサンプリング
  theta_u <- rnorm(hh, mu_par, sqrt(sigma_par))
  
  
  ##アイテム変量効果を推定
  #ユーザーごとに事後分布のパラメータを設定
  inv_tau_v <- 1/tau_v^2
  mu_par <- rep(0, item)
  sigma_par <- rep(0, item)
  
  for(i in 1:item){
    index <- context_w2[item_index[[i]]]
    X <- theta_w2[index]
    Xy <- (t(X) %*% theta_vw[item_index[[i]]])
    XXV <- (t(X) %*% X) + inv_tau_v
    inv_XXV <- 1/XXV
    sigma_par[i] <- tau_vw^2 * inv_XXV
    mu_par[i] <- inv_XXV %*% (Xy + inv_tau_v %*% v_mu[i])
  }
  #正規分布から事後分布をサンプリング
  theta_v <- rnorm(item, mu_par, sqrt(sigma_par))
  
  
  ##コンテキスト変量効果をサンプリング
  #コンテキストごとに事後分布のパラメータを設定
  inv_tau_w <- 1/tau_w^2
  mu_par <- matrix(0, nrow=context, ncol=2)
  sigma_par <- matrix(0, nrow=context, ncol=2)
  
  for(i in 1:context){
    index1 <- context_u[context_index1[[i]]]; index2 <- context_v[context_index2[[i]]]
    X1 <- theta_u[index1]; X2 <- theta_v[index2]
    Xy1 <- (t(X1) %*% theta_uw[context_index1[[i]]])
    Xy2 <- (t(X2) %*% theta_vw[context_index2[[i]]])
    XXV1 <- (t(X1) %*% X1) + inv_tau_w
    XXV2 <- (t(X2) %*% X2) + inv_tau_w
    inv_XXV1 <- 1/XXV1; inv_XXV2 <- 1/XXV2
    sigma_par[i, 1] <- tau_uw^2 * inv_XXV1
    sigma_par[i, 2] <- tau_vw^2 * inv_XXV2
    mu_par[i, 1] <- inv_XXV1 %*% (Xy1 + inv_tau_w %*% w_mu[i, 1])
    mu_par[i, 2] <- inv_XXV2 %*% (Xy2 + inv_tau_w %*% w_mu[i, 2])
  }
  
  #正規分布から事後分布をサンプリング
  theta_w1 <- rnorm(context, mu_par[, 1], sqrt(sigma_par[, 1]))
  theta_w2 <- rnorm(context, mu_par[, 2], sqrt(sigma_par[, 2]))
  theta_w <- cbind(theta_w1, theta_w2)
  
  
  ##パラメータを観測ベクトルのレコード数に変換
  theta_u_vec <- theta_u[context_u]
  theta_w1_vec <- theta_w1[context_w1]
  theta_v_vec <- theta_v[context_v]
  theta_w2_vec <- theta_w2[context_w2]
  
  
  ###素性回帰モデルと階層モデルのパラメータをサンプリング
  ##回帰パラメータをサンプリング
  #回帰ベクトルのパラメータ
  u_er <- U - theta_uw_vec - theta_vw_vec   #応答変数の設定
  inv_XXV <- solve(xx + tau_inv1)
  Xy <- t(x) %*% u_er
  mu_par <- inv_XXV %*% (Xy + tau_inv1 %*% mu1)
  
  #正規分布から回帰ベクトルをサンプリング
  beta <- mvrnorm(1, mu_par, sigma^2*inv_XXV)
  beta_mu <- as.numeric(x %*% beta)   #素性ベクトルの平均構造
  
  
  ##ユーザーの階層モデルのパラメータをサンプリング
  ##回帰パラメータをサンプリング
  #回帰ベクトルのパラメータ
  inv_XXV <- solve(uu + tau_inv2)
  Xy <- t(u) %*% theta_u
  mu_par <- inv_XXV %*% (Xy + tau_inv2 %*% mu2)
  
  #正規分布から回帰ベクトルをサンプリング
  alpha_u <- mvrnorm(1, mu_par, tau_u^2*inv_XXV)
  u_mu <- as.numeric(u %*% alpha_u)   #素性ベクトルの平均構造
  
  ##ユーザー変量効果の分散をサンプリング
  #逆ガンマ分布より分散をサンプリング
  s1 <- s01 + sum((theta_u - u_mu)^2)
  v1 <- v01 + hh
  tau_u <- sqrt(1/(rgamma(1, v1/2, s1/2)))
  

  ##アイテムの階層モデルのパラメータをサンプリング
  ##回帰パラメータをサンプリング
  #回帰ベクトルのパラメータ
  inv_XXV <- solve(vv + tau_inv3)
  Xy <- t(v) %*% theta_v
  mu_par <- inv_XXV %*% (Xy + tau_inv3 %*% mu3)
  
  #正規分布から回帰ベクトルをサンプリング
  alpha_v <- mvrnorm(1, mu_par, tau_v^2*inv_XXV)
  v_mu <- as.numeric(v %*% alpha_v)   #素性ベクトルの平均構造
  
  ##アイテム変量効果の分散をサンプリング
  #逆ガンマ分布より分散をサンプリング
  s1 <- s02 + sum((theta_v - v_mu)^2)
  v1 <- v02 + item
  tau_v <- sqrt(1/(rgamma(1, v1/2, s1/2)))
  
  
  ##コンテキストの階層モデルのパラメータをサンプリング
  #多変量回帰モデルから回帰ベクトルをサンプリング
  out <- rmultireg(theta_w, w, Deltabar, ADelta, nu, V)
  alpha_w <- out$B
  w_mu <- w %*% alpha_w
  
  ##サンプリング結果を保存と結果の表示
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    ALPHA_U[mkeep, ] <- alpha_u
    ALPHA_V[mkeep, ] <- alpha_v
    ALPHA_W[mkeep, ] <- as.numeric(alpha_w)
    THETA_U[mkeep, ] <- theta_u
    THETA_V[mkeep, ] <- theta_v
    THETA_W[mkeep, ] <- as.numeric(theta_w)
    COV[mkeep, ] <- c(tau_uw, tau_vw, tau_u, tau_v, tau_w)
  }
  
  #コンテキスト依存変量効果はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    rkeep <- c(rkeep, rp)
    THETA_UW <- THETA_UW + theta_uw
    THETA_VW <- THETA_VW + theta_vw
  }
    
  if(rp%%disp==0){
    #対数尤度を推定
    Mu <- beta_mu + theta_uw_vec + theta_vw_vec   #完全データの平均構造
    prob <- pnorm(Mu, 0, sigma)   #購買確率
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #対数尤度
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LL1, LLst))
    print(round(c(tau_u, tau_v, tau_w, tau_uw, tau_vw), 3))
    print(round(rbind(beta, betat), 3))
  }
}


matplot(COV, type="l")

tau_u
tau_v
tau_uw
tau_vw
mean(prob[y==1])
mean((1-prob)[y==0])
563474.1 

round(cbind(theta_uw, theta_uwt0, n_uw), 3)
round(cbind(theta_vw, theta_vwt0, n_vw), 3)
round(cbind(theta_u, theta_ut), 3)
round(cbind(theta_v, theta_vt), 3)
round(cbind(theta_w, theta_wt), 3)


