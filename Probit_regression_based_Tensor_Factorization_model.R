#####Probit regression based Tensor Factorization model#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(gtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####データの発生####
##データの設定
hh <- 10000   #ユーザー数
item <- 3000   #アイテム数
context <- 24   #コンテキスト数
pt <- rpois(hh, rgamma(hh, 25, 0.20))   #ユーザーあたりの接触数
N <- sum(pt)   #総レコード数
k <- 10   #基底数
vec <- rep(1, k)

##IDとインデックスを設定
#IDの設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:N, user_id, rank)))

#インデックスの設定
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}
user_data <- sparseMatrix(1:N, user_id, x=rep(1, N), dims=c(N, hh))


##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(0.75, topic))) %*% 1:topic)

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
item_data <- sparseMatrix(1:N, item_id, x=rep(1, N), dims=c(N, item))
item_n <- colSums(item_data)


##コンテキストの割当を生成
#セグメント割当を生成
topic <- 15
phi <- extraDistr::rdirichlet(topic, rep(0.2, context))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(0.5, topic))) %*% 1:topic)

#多項分布からアイテムを生成
context_id_list <- list()
for(i in 1:hh){
  if(i%%1000==0){
    print(i)
  }
  context_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, phi[z[user_id[user_list[[i]]]], ]) %*% 1:context)
}
context_id <- unlist(context_id_list)
context_list <- list()
for(j in 1:context){
  context_list[[j]] <- which(context_id==j)
}
context_data <- sparseMatrix(1:N, context_id, x=rep(1, N), dims=c(N, context))
context_n <- colSums(context_data)


#購買数をカウント
freq <- plyr::count(user_id); freq_user <- freq$freq[order(freq$x)]
freq <- plyr::count(item_id); freq_item <- freq$freq[order(freq$x)]
freq <- plyr::count(context_id); freq_context <- freq$freq[order(freq$x)]
hist(freq_user, col="grey", breaks=25, main="ユーザーごとの購買数", xlab="購買数")
hist(freq_item, col="grey", breaks=25, main="アイテムごとの購買数", xlab="購買数")


##階層モデルの説明変数を設定
#ユーザーの説明変数
k1 <- 4; k2 <- 4; k3 <- 4
u1 <- matrix(runif(hh*k1, 0, 1), nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合
col_u <- ncol(u)

#アイテムの説明変数
k1 <- 3; k2 <- 5; k3 <- 5
v1 <- matrix(runif(item*k1, 0, 1), nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #データを結合
col_v <- ncol(v)


####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##ユーザーベースの階層モデルのパラメータ
  #階層モデルの分散パラメータ
  Cov_ut1 <- Cov_u1 <- runif(1, 0.25, 0.4)
  Cov_ut2 <- Cov_u2 <- diag(runif(k, 0.01, 0.25))
  Cov_ut3 <- Cov_u3 <- diag(runif(k, 0.01, 0.25)) 
  
  #回帰係数を設定
  alpha_u1 <- rep(0, ncol(u))
  alpha_u2 <- alpha_u3 <- matrix(0, nrow=ncol(u), ncol=k)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u1[j] <- runif(1, -0.6, -0.2)
      alpha_u2[j, ] <- runif(k, -0.7, -0.1)
      alpha_u3[j, ] <- runif(k, -0.7, -0.1)
      
    } else {
      alpha_u1[j] <- runif(1, -0.6, 0.4)
      alpha_u2[j, ] <- runif(k, -0.6, 0.5)
      alpha_u3[j, ] <- runif(k, -0.6, 0.5)
    }
  }
  alpha_ut1 <- alpha_u1; alpha_ut2 <- alpha_u2; alpha_ut3 <- alpha_u3
  
  #ユーザーの変量効果と特徴行列を生成
  theta_u1 <- theta_ut1 <- as.numeric(u %*% alpha_u1 + rnorm(hh, 0, Cov_u1))   #変量効果のパラメータ
  theta_u2 <- theta_ut2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)   #行列分解のパラメータ
  theta_u3 <- theta_ut3 <- u %*% alpha_u3 + mvrnorm(hh, rep(0, k), Cov_u3)   #テンソル分解のパラメータ
  
  
  ##アイテムベースの階層モデルのパラメータ
  #分散共分散行列を設定
  Cov_vt1 <- Cov_v1 <- runif(1, 0.25, 0.4)
  Cov_vt2 <- Cov_v2 <- diag(runif(k, 0.01, 0.25))
  Cov_vt3 <- Cov_v3 <- diag(runif(k, 0.01, 0.25))
  
  #回帰係数を設定
  alpha_v1 <- rep(0, ncol(v))
  alpha_v2 <- alpha_v3 <- matrix(0, nrow=ncol(v), ncol=k)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v1[j] <- runif(1, -0.5, -0.2)
      alpha_v2[j, ] <- runif(k, -0.7, 0.4)
      alpha_v3[j, ] <- runif(k, -0.7, 0.4)
    } else {
      alpha_v1[j] <- runif(1, -0.7, 0.5)
      alpha_v2[j, ] <- runif(k, -0.7, 0.5)
      alpha_v3[j, ] <- runif(k, -0.7, 0.5)
    }
  }
  alpha_vt1 <- alpha_v1; alpha_vt2 <- alpha_v2; alpha_vt3 <- alpha_v3
  
  #アイテムの変量効果と特徴行列を生成
  theta_v1 <- theta_vt1 <- as.numeric(v %*% alpha_v1 + rnorm(item, 0, Cov_v1))   #変量効果のパラメータ
  theta_v2 <- theta_vt2 <- v %*% alpha_v2 + mvrnorm(item, rep(0, k), Cov_v2)   #行列分解のパラメータ
  theta_v3 <- theta_vt3 <- v %*% alpha_v3 + mvrnorm(item, rep(0, k), Cov_v3)   #テンソル分解のパラメータ
  
  
  ##コンテキストベースの階層モデルのパラメータ
  #階層モデルのパラメータを生成
  alpha_c1 <- alpha_ct1 <- 0
  alpha_c3 <- alpha_ct3 <- rep(0, k)
  alpha_c2 <- alpha_ct2 <- rep(0, k)
  Cov_c1 <- Cov_ct1 <- runif(1, 0.25, 0.4)
  Cov_c2 <- Cov_ct2 <- diag(runif(k, 0.01, 0.25))
  Cov_c3 <- Cov_ct3 <- diag(runif(k, 0.01, 0.25))
  
  #コンテキストの変量効果と特徴行列を生成
  theta_c1 <- theta_ct1 <- rnorm(context, alpha_c1, Cov_c1)
  theta_c2 <- theta_ct2 <- mvrnorm(context, alpha_c2, Cov_c2)
  theta_c3 <- theta_ct3 <- mvrnorm(context, alpha_c3, Cov_c3)


  ##正規分布から効用と購買ベクトルを生成
  #モデル誤差
  Sigma <- Sigmat <- 1
  
  #行列分解のパラメータを生成
  uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
  
  #テンソル分解のパラメータを生成
  uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)
  
  #潜在効用を生成
  mu <- mut <-  theta_u1[user_id] + theta_v1[item_id] + theta_c2[context_id] + uv + uc + vc + uvc   #期待値
  U <- mu + rnorm(N, 0, Sigma)   #誤差を生成
  
  #購買ベクトルに変換
  y <- ifelse(U > 0, 1, 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break条件
}


####マルコフ連鎖モンテカルロ法で階層ベイズテンソル分解を推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 2000
keep <- 2  
iter <- 0
burnin <- 500/keep
disp <- 10

##事前分布の設定
#変量効果の階層モデルの事前分布
gamma_u <- rep(0, col_u)
tau_u <- 100 * diag(col_u); inv_tau_u <- solve(tau_u)
gamma_v <- rep(0, col_v)
tau_v <- 100 * diag(col_v); inv_tau_v <- solve(tau_v)
gamma_c <- 0
tau_c <- 100; inv_tau_c <- 1/tau_c
v0 <- 1; s0 <- 1

#ユーザーの階層モデルの事前分布
Deltabar1 <- matrix(rep(0, col_u*k), nrow=col_u, ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta1 <- 0.01 * diag(rep(1, col_u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#アイテムの階層モデルの事前分布
Deltabar2 <- matrix(rep(0, col_v*k), nrow=col_v, ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta2 <- 0.01 * diag(rep(1, col_v))   #階層モデルの回帰係数の事前分布の分散
nu2 <- k + 1   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#コンテキストの階層モデルの事前分布
Deltabar3 <- rep(0, k)   #階層モデルの回帰係数の事前分布の分散
ADelta3 <- 0.01 * diag(k)   #階層モデルの回帰係数の事前分布の分散
nu3 <- k + 1   #逆ウィシャート分布の自由度
V3 <- nu3 * diag(rep(1, k))   #逆ウィシャート分布のパラメータ


##パラメータの真値
#階層モデルの真値
alpha_u1 <- alpha_ut1; alpha_u2 <- alpha_ut2; alpha_u3 <- alpha_ut3
user_mu1 <- as.numeric(u %*% alpha_u1); user_mu2 <- u %*% alpha_u2; user_mu3 <- u %*% alpha_u3
alpha_v1 <- alpha_vt1; alpha_v2 <- alpha_vt2; alpha_v3 <- alpha_vt3
item_mu1 <- as.numeric(v %*% alpha_v1); item_mu2 <- v %*% alpha_v2; item_mu3 <- v %*% alpha_v3
alpha_c1 <- alpha_ct1; alpha_c2 <- alpha_ct2; alpha_c3 <- alpha_ct3
Cov_u1 <- Cov_ut1; Cov_u2 <- Cov_ut2; Cov_u3 <- Cov_ut3
inv_Cov_u2 <- solve(Cov_u2); inv_Cov_u3 <- solve(Cov_u3)
Cov_v1 <- Cov_vt1; Cov_v2 <- Cov_vt2; Cov_v3 <- Cov_vt3
inv_Cov_v2 <- solve(Cov_v2); inv_Cov_v3 <- solve(Cov_v3)
Cov_c1 <- Cov_ct1; Cov_c2 <- Cov_ct2; Cov_c3 <- Cov_ct3
inv_Cov_c2 <- solve(Cov_c2); inv_Cov_c3 <- solve(Cov_c3)

#変量効果と特徴行列の真値
Sigma <- Sigmat
theta_u1 <- theta_ut1; theta_u2 <- theta_ut2; theta_u3 <- theta_ut3
theta_v1 <- theta_vt1; theta_v2 <- theta_vt2; theta_v3 <- theta_vt3
theta_c1 <- theta_ct1; theta_c2 <- theta_ct2; theta_c3 <- theta_ct3

#行列分解とテンソル分解のパラメータ
theta_u_vec1 <- theta_u1[user_id]
theta_v_vec1 <- theta_v1[item_id]
theta_c_vec1 <- theta_c1[context_id]
uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)


##パラメータの初期値
#階層モデルのパラメータ
alpha_u1 <- rep(0, col_u)
alpha_u2 <- matrix(0, nrow=col_u, ncol=k)
alpha_u3 <- matrix(0, nrow=col_u, ncol=k)
alpha_v1 <- rep(0, col_v)
alpha_v2 <- matrix(0, nrow=col_v, ncol=k)
alpha_v3 <- matrix(0, nrow=col_v, ncol=k)
alpha_c1 <- 0
alpha_c2 <- alpha_c3 <- rep(0, k)
Cov_u1 <- 0.2; Cov_u2 <- 0.01 * diag(k); Cov_u3 <- 0.01 * diag(k)
Cov_v1 <- 0.2; Cov_v2 <- 0.01 * diag(k); Cov_v3 <- 0.01 * diag(k)
Cov_c1 <- 0.2; Cov_c2 <- 0.01 * diag(k); Cov_c3 <- 0.01 * diag(k)

#変量効果のパラメータ
Sigma <- 1
theta_u1 <- u %*% alpha_u1 + rnorm(hh, 0, Cov_u1)
theta_v1 <- v %*% alpha_v1 + rnorm(item, 0, Cov_v1)
theta_c1 <- rnorm(context, alpha_c1, Cov_c1)
theta_u_vec1 <- theta_u1[user_id]
theta_v_vec1 <- theta_v1[item_id]
theta_c_vec1 <- theta_c1[context_id]

#行列分解のパラメータ
theta_u2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)
theta_v2 <- v %*% alpha_v2 + mvrnorm(item, rep(0, k), Cov_v2)
theta_c2 <- mvrnorm(context, alpha_c2, Cov_c2)
uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)

#テンソル分解のパラメータ
theta_u3 <- u %*% alpha_u3 + mvrnorm(hh, rep(0, k), Cov_u3)
theta_v3 <- v %*% alpha_v3 + mvrnorm(item, rep(0, k), Cov_v3)
theta_t3 <- mvrnorm(context, alpha_c3, Cov_c3)
uvc <- as.numeric((theta_u3[user_id, ] * theta_v3[item_id, ] * theta_c3[context_id, ]) %*% vec)


##データとインデックスを設定
#切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

#データの設定
uu <- t(u) %*% u + inv_tau_u; inv_uu <- solve(uu)
vv <- t(v) %*% v + inv_tau_v; inv_vv <- solve(vv)

#インデックスの設定
ui_id <- uc_id <- iu_id <- ic_id <- cu_id <- ci_id <- list()
for(i in 1:hh){
  ui_id[[i]] <- item_id[user_list[[i]]]
  uc_id[[i]] <- context_id[user_list[[i]]]
}
for(j in 1:item){
  iu_id[[j]] <- user_id[item_list[[j]]]
  ic_id[[j]] <- context_id[item_list[[j]]]
}
for(j in 1:context){
  cu_id[[j]] <- user_id[context_list[[j]]]
  ci_id[[j]] <- item_id[context_list[[j]]]
}


##サンプリング結果のパラメータの保存用配列
#モデルパラメータの格納用配列
THETA_U1 <- matrix(0, nrow=R/keep, ncol=hh)
THETA_U2 <- array(0, dim=c(hh, k, R/keep))
THETA_U3 <- array(0, dim=c(hh, k, R/keep))
THETA_V1 <- matrix(0, nrow=R/keep, ncol=item)
THETA_V2 <- array(0, dim=c(item, k, R/keep))
THETA_V3 <- array(0, dim=c(item, k, R/keep))
THETA_C1 <- matrix(0, nrow=R/keep, ncol=context)
THETA_C2 <- array(0, dim=c(context, k, R/keep))
THETA_C3 <- array(0, dim=c(context, k, R/keep))

#階層モデルの格納用配列
ALPHA_U1 <- matrix(0, nrow=R/keep, ncol=col_u)
ALPHA_U2 <- array(0, dim=c(col_u, k, R/keep))
ALPHA_U3 <- array(0, dim=c(col_u, k, R/keep))
ALPHA_V1 <- matrix(0, nrow=R/keep, ncol=col_v)
ALPHA_V2 <- array(0, dim=c(col_v, k, R/keep))
ALPHA_V3 <- array(0, dim=c(col_v, k, R/keep))
COV_U1 <- COV_V1 <- COV_C1 <- rep(0, R/keep)
COV_U2 <- COV_V2 <- COV_C2 <- array(0, dim=c(k, k, R/keep))
COV_U3 <- COV_V3 <- COV_C3 <- array(0, dim=c(k, k, R/keep))


##対数尤度の基準値
#1パラメータモデルの対数尤度
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #対数尤度

#真値での対数尤度
prob <- pnorm(mut, 0, Sigma)   #購買確率
LLbest <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #対数尤度


####ギブスサンプリングでパラメータをサンプリング
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  mu <- theta_u_vec1 + theta_v_vec1 + theta_c_vec1 + uv + uc + vc + uvc   #潜在効用の期待値
  U <- extraDistr::rtnorm(N, mu, Sigma, a, b)   #潜在効用を生成
  
  
  ##ユーザーの変量効果をサンプリング
  #モデルの応答変数
  u_er <- U - theta_v_vec1 - theta_c_vec1 - uv - uc - vc - uvc   
  
  #ユーザーの変量効果の事後分布のパラメータ
  u_mu <- rep(0, hh)
  for(i in 1:hh){
    u_mu[i] <- mean(u_er[user_list[[i]]])
  }
  weights <- Cov_u1^2 / (Sigma^2/freq_user + Cov_u1^2)   #重み係数
  mu_par <- weights*u_mu + (1-weights)*user_mu1   #事後分布の平均
  
  #正規分布より事後分布をサンプリング
  theta_u1 <- rnorm(hh, mu_par, sqrt(1 / (1/Cov_u1^2 + freq_user/Sigma^2)))
  theta_u_vec1 <- theta_u1[user_id]
  
  ##アイテムの変量効果をサンプリング
  #モデルの応答変数
  v_er <- U - theta_u_vec1 - theta_c_vec1 - uv - uc - vc - uvc
  
  #アイテムの変量効果の事後分布のパラメータ
  v_mu <- rep(0, item)
  for(j in 1:item){
    v_mu[j] <- mean(v_er[item_list[[j]]])
  }
  weights <- Cov_v1^2 / (Sigma^2/freq_item + Cov_v1^2)   #重み係数
  mu_par <- weights*v_mu + (1-weights)*item_mu1   #事後分布の平均
  
  #正規分布より事後分布をサンプリング
  theta_v1 <- rnorm(item, mu_par, sqrt(1 / (1/Cov_v1^2 + freq_item/Sigma^2)))
  theta_v_vec1 <- theta_v1[item_id]
  
  ##コンテキストの変量効果をサンプリング
  #モデルの応答変数
  c_er <- U - theta_u_vec1 - theta_v_vec1 - uv - uc - vc - uvc
  
  #アイテムの変量効果の事後分布のパラメータ
  c_mu <- rep(0, context)
  for(j in 1:context){
    c_mu[j] <- mean(c_er[context_list[[j]]])
  }
  weights <- Cov_c1^2 / (Sigma^2/freq_context + Cov_c1^2)   #重み係数
  mu_par <- weights*c_mu + (1-weights)*alpha_c1   #事後分布の平均
  
  #正規分布より事後分布をサンプリング
  theta_c1 <- rnorm(context, mu_par, sqrt(1 / (1/Cov_c1^2 + freq_context/Sigma^2)))
  theta_c_vec1 <- theta_c1[context_id]
  
  
  ##ユーザーの特徴行列をサンプリング
  #モデルの応答変数
  u_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - vc - uvc
  
  for(i in 1:hh){
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_v2[ui_id[[i]], ] + theta_c2[uc_id[[i]], ]
    Xy <- t(X) %*% u_er[user_list[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u2)
    mu <- inv_XXV %*% (Xy + inv_Cov_u2 %*% user_mu2[i, ])   #事後分布の平均
    
    #多変量正規分布からユーザー特徴行列をサンプリング
    theta_u2[i, ] <- mvrnorm(1, mu, Sigma^2*inv_XXV)
  }
  #行列分解のパラメータを更新
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  
  ##アイテムの特徴行列をサンプリング
  #モデルの応答変数
  v_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uc - uvc
  
  for(j in 1:item){
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_u2[iu_id[[j]], ] + theta_c2[ic_id[[j]], ]
    Xy <- t(X) %*% v_er[item_list[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v2)
    mu <- inv_XXV %*% (Xy + inv_Cov_v2 %*% item_mu2[j, ])   #事後分布の平均
    
    #多変量正規分布からユーザー特徴行列をサンプリング
    theta_v2[j, ] <- mvrnorm(1, mu, Sigma^2*inv_XXV)
  }
  #行列分解のパラメータを更新
  uv <- as.numeric((theta_u2[user_id, ] * theta_v2[item_id, ]) %*% vec)
  
  ##コンテキストの特徴行列をサンプリング
  #モデルの応答変数
  c_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uv - uvc
  
  for(j in 1:context){
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_u2[cu_id[[j]], ] + theta_v2[ci_id[[j]], ]
    Xy <- t(X) %*% c_er[context_list[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_c2)
    mu <- inv_XXV %*% (Xy + inv_Cov_c2 %*% alpha_c2)   #事後分布の平均
    
    #多変量正規分布からユーザー特徴行列をサンプリング
    theta_c2[j, ] <- mvrnorm(1, mu, Sigma^2*inv_XXV)
  }
  #行列分解のパラメータを更新
  uc <- as.numeric((theta_u2[user_id, ] * theta_c2[context_id, ]) %*% vec)
  vc <- as.numeric((theta_v2[item_id, ] * theta_c2[context_id, ]) %*% vec)
  
  
  ##ユーザーのテンソルのパラメータをサンプリング
  #モデルの応答変数
  y_er <- U - theta_u_vec1 - theta_v_vec1 - theta_c_vec1 - uv - uc - vc
  
  for(i in 1:hh){
    #特徴ベクトルのパラメータ
    X <- theta_v3[ui_id[[i]], ] * theta_c3[uc_id[[i]], ]
    Xy <- t(X) %*% y_er[user_list[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_u3 %*% user_mu3[i, ])   #多変量正規分布の平均ベクトル
    
    #多変量正規分布からパラメータをサンプリング
    theta_u3[i, ] <- mvrnorm(1, beta_mu, Sigma^2*inv_XXV)
  }
  
  ##アイテムのテンソルのパラメータをサンプリング
  for(j in 1:item){
    #特徴ベクトルのパラメータ
    X <- theta_u3[iu_id[[j]], ] * theta_c3[ic_id[[j]], ]
    Xy <- t(X) %*% y_er[item_list[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_v3 %*% item_mu3[j, ])   #多変量正規分布の平均ベクトル
    
    #多変量正規分布からパラメータをサンプリング
    theta_v3[j, ] <- mvrnorm(1, beta_mu, Sigma^2*inv_XXV)
  }
  
  ##コンテキストのテンソルのパラメータをサンプリング
  for(j in 1:context){
    #特徴ベクトルのパラメータ
    X <- theta_u3[cu_id[[j]], ] * theta_v3[ci_id[[j]], ]
    Xy <- t(X) %*% y_er[context_list[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_c3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov_c3 %*% alpha_c3)   #多変量正規分布の平均ベクトル
    
    #多変量正規分布からパラメータをサンプリング
    theta_c3[j, ] <- mvrnorm(1, beta_mu, Sigma^2*inv_XXV)
  }
  
  
  ##ユーザーの変量効果の階層モデルのパラメータをサンプリング
  #ユーザーの変量効果の階層モデルの回帰ベクトルを更新
  beta_mu <- inv_uu %*% (t(u) %*% theta_u1 + inv_tau_u %*% gamma_u)
  alpha_u1 <- mvrnorm(1, beta_mu, Cov_u1^2*inv_uu)   #多変量正規分布から回帰ベクトルをサンプリング
  user_mu1 <- as.numeric(u %*% alpha_u1)
  
  #ユーザーの階層モデルの標準偏差を更新
  er <- theta_u1 - user_mu1
  s1 <- t(er) %*% er + s0; v1 <- hh + v0   #ガンマ分布のパラメータ
  Cov_u1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #ガンマ分布から標準偏差をサンプリング
  inv_Cov_u1 <- 1 / Cov_u1
  
  ##アイテムの変量効果の階層モデルのパラメータをサンプリング
  #アイテムの変量効果の階層モデルの回帰ベクトルを更新
  beta_mu <- inv_vv %*% (t(v) %*% theta_v1 + inv_tau_v %*% gamma_v)
  alpha_v1 <- mvrnorm(1, beta_mu, Cov_v1^2*inv_vv)   #多変量正規分布から回帰ベクトルをサンプリング
  item_mu1 <- as.numeric(v %*% alpha_v1)
  
  #ユーザーの階層モデルの標準偏差を更新
  er <- theta_v1 - item_mu1
  s1 <- t(er) %*% er + s0; v1 <- item + v0   #ガンマ分布のパラメータ
  Cov_v1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #ガンマ分布から標準偏差をサンプリング
  inv_Cov_v1 <- 1 / Cov_v1
  
  ##コンテキストの変量効果の階層モデルのパラメータをサンプリング
  #コンテキストの変量効果の階層モデルの平均ベクトルを更新
  mu <- mean(theta_c1)
  weights <- tau_c^2 / (Cov_c1^2/context + tau_c^2)   #重み係数
  mu_par <- weights*mu   #事後分布の平均
  alpha_c1 <- rnorm(1, mu_par, sqrt(1 / (1/tau_c^2 + context/Cov_c1^2)))   #正規分布から変量効果をサンプリング
  
  #コンテキストの階層モデルの標準偏差を更新
  er <- theta_c1 - alpha_c1
  s1 <- t(er) %*% er + s0; v1 <- context + v0   #ガンマ分布のパラメータ
  Cov_c1 <- sqrt(1/rgamma(1, v1/2, s1/2))   #ガンマ分布から標準偏差をサンプリング
  inv_Cov_c1 <- 1 / Cov_c1
  
  
  ##ユーザー特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからユーザーの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_u2, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u2 <- out$B
  user_mu2 <- u %*% alpha_u2   #ユーザー特徴行列の平均構造
  Cov_u2 <- diag(diag(out$Sigma))
  inv_cov_u2 <- solve(Cov_u2)
  
  ##アイテム特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからアイテムの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v2, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v2 <- out$B
  item_mu2 <- v %*% alpha_v2   #アイテム特徴行列の平均構造
  Cov_v2 <- diag(diag(out$Sigma))
  inv_cov_v2 <- solve(Cov_v2)
  
  ##コンテキスト特徴行列の階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  IW <- t(theta_c2) %*% theta_c2 + solve(V3)
  Sn <- nu3 + context
  Cov_c2 <- diag(diag(rwishart(Sn, solve(IW))$IW))   #逆ウィシャート分布からパラメータをサンプリング
  inv_Cov_c2 <- solve(Cov_c2)
  
  
  ##ユーザーテンソルの階層モデルのパラメータをサンプリング
  #多変量回帰モデルからユーザーの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_u3, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u3 <- out$B
  user_mu3 <- u %*% alpha_u3   #ユーザーテンソルの平均構造
  Cov_u3 <- diag(diag(out$Sigma))
  inv_cov_u3 <- solve(Cov_u3)
  
  ##アイテムテンソルの階層モデルのパラメータをサンプリング
  #多変量回帰モデルからアイテムの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v3, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v3 <- out$B
  item_mu3 <- v %*% alpha_v3   #アイテム特徴行列の平均構造
  Cov_v3 <- diag(diag(out$Sigma))
  inv_cov_v3 <- solve(Cov_v3)
  
  ##コンテキストテンソルの階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  IW <- t(theta_c3) %*% theta_c3 + solve(V3)
  Sn <- nu3 + context
  Cov_c3 <- diag(diag(rwishart(Sn, solve(IW))$IW))   #逆ウィシャート分布からパラメータをサンプリング
  inv_Cov_c3 <- solve(Cov_c3)
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    
    #モデルパラメータの格納
    THETA_U1[mkeep, ] <- theta_u1; THETA_U2[, , mkeep] <- theta_u2; THETA_U3[, , mkeep] <- theta_u3
    THETA_V1[mkeep, ] <- theta_v1; THETA_V2[, , mkeep] <- theta_v2; THETA_V3[, , mkeep] <- theta_v3
    THETA_C1[mkeep, ] <- theta_c1; THETA_C2[, , mkeep] <- theta_c2; THETA_C3[, , mkeep] <- theta_c3
    
    #階層モデルのパラメータを格納
    ALPHA_U1[mkeep, ] <- alpha_u1; ALPHA_U2[, , mkeep] <- alpha_u2; ALPHA_U3[, , mkeep] <- alpha_u3
    ALPHA_V1[mkeep, ] <- alpha_v1; ALPHA_V2[, , mkeep] <- alpha_v2; ALPHA_V3[, , mkeep] <- alpha_v3
    COV_U1[mkeep] <- Cov_u1; COV_U2[, , mkeep] <- Cov_u2; COV_U3[, , mkeep] <- Cov_u3
    COV_V1[mkeep] <- Cov_v1; COV_V2[, , mkeep] <- Cov_v2; COV_V3[, , mkeep] <- Cov_v3
    COV_C1[mkeep] <- Cov_c1; COV_C2[, , mkeep] <- Cov_c2; COV_C3[, , mkeep] <- Cov_c3
  }

  if(rp%%disp==0){
    #対数尤度を計算
    mu <- theta_u_vec1 + theta_v_vec1 + theta_c_vec1 + uv + uc + vc + uvc   #潜在効用の期待値
    prob <- pnorm(mu, 0, Sigma)   #購買確率
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #対数尤度
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLst, LLbest))
    print(round(rbind(diag(Cov_u2), diag(Cov_ut2)), 3))
  }
}
