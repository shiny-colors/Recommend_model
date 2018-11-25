#####Bayesian Hierarchical Field Aware Factorization Machines#####
library(MASS)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(extraDistr)
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
hh <- 5000   #ユーザー数
item <- 2000   #アイテム数
tag <- 150   #タグ数
pt <- rtpois(hh, rgamma(hh, 30.0, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)
n <- rtpois(item, 0.75, a=0, b=Inf)   #タグ数
k <- 10   #基底数
vec_k <- rep(1, k)

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

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


##アイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(2.5, topic))) %*% 1:topic)

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
item_data <- colSums(sparseMatrix(1:hhpt, item_id, x=rep(1, hhpt), dims=c(hhpt, item)))

##タグを生成
#パラメータの設定
omega <- extraDistr::rdirichlet(topic, rep(0.5, tag))
z <- as.numeric(rmnom(item, 1, t(phi) / rowSums(t(phi))) %*% 1:topic)

#多項分布からタグを生成
max_n <- max(n)
tag_list <- list()
tag_id <- matrix(0, nrow=hhpt, ncol=max_n)
for(j in 1:item){
  repeat { 
    x <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(x)==1){
      tag_list[[j]] <- (x * 1:tag)[x > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
  
}
tag_data <- sparseMatrix(matrix(rep(1:hhpt, max_n), nrow=hhpt, ncol=max_n)[tag_id!=0], tag_id[tag_id!=0], 
                         x=rep(1, sum(tag_id!=0)), dims=c(hhpt, tag))



####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  ##階層モデルのパラメータを生成
  ##ユーザーベースの階層モデルのパラメータ
  #分散共分散行列を設定
  Cov_ut <- Cov_u <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.0025, 0.25)$covariance
  
  #回帰係数を設定
  alpha_u <- matrix(0, nrow=ncol(u), ncol=k)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u[j, ] <- runif(k, -1.3, -0.5)
    } else {
      alpha_u[j, ] <- runif(k, -0.6, 0.7)
    }
  }
  alpha_ut <- alpha_u
  
  #多変量回帰モデルからユーザー個別の回帰パラメータを生成
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, k), Cov_u)   #テンソル分解のパラメータ
  
  
  ##アイテムベースの階層モデルのパラメータ
  #分散共分散行列を設定
  Cov_vt <- Cov_v <- covmatrix(k, corrM(k, -0.6, 0.8, 0.05, 0.2), 0.025, 0.25)$covariance
  
  #回帰係数を設定
  alpha_v <- matrix(0, nrow=ncol(v), ncol=k)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(k, -1.2, -0.4)
    } else {
      alpha_v[j, ] <- runif(k, -0.7, 0.8)
    }
  }
  alpha_vt <- alpha_v
  
  #多変量回帰モデルからアイテム個別の回帰パラメータを生成
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  
  ##時間ベースの階層モデルのパラメータ
  alpha_t <- alpha_tt <- rep(-0.25, k)
  Cov_t <- Cov_tt <- diag(0.25, k)
  theta_t <- theta_tt <- mvrnorm(time, alpha_t, Cov_t)
  
  
  ##正規分布から効用と購買ベクトルを生成
  #潜在効用を生成
  WHC0 <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC0 <- WHC0 + theta_u[, j] %o% t(theta_v)[j, ] %o% t(theta_t)[j, ]
  }
  #whc0 <- as.numeric((W0[user_id0, ] * t(H0)[item_id0, ] * t(C0)[time_id0, ]) %*% vec)   #こちらでもok
  whc0 <- as.numeric(WHC0)   #テンソルをベクトルに変換
  u_vec0 <- whc0 + rnorm(hh*item*time, 0, 1)   #誤差を生成
  
  #購買ベクトルに変換
  y0 <- ifelse(u_vec0 > 0, 1, 0)
  if(mean(y0) > 0.3 & mean(y0) < 0.4) break   #break条件
}

##欠損ベクトルを生成
#欠損確率を生成
user_prob <- rbeta(hh, 10, 50)
item_prob <- rbeta(item, 15, 55)
time_prob <- rbeta(time, 60, 140)
prob <- user_prob[user_id0]*item_id0[item_id0]*time_prob[time_id0]

#ベルヌーイ分布から欠損ベクトルを生成
z_vec <- rbinom(N0, 1, prob)
N <- sum(z_vec)
y <- y0[z_vec==1]; u_vec <- u_vec0[z_vec==1]; whc <- whc0[z_vec==1]
hist(u_vec, breaks=25, col="grey", main="潜在的なスコア分布", xlab="スコア")

#欠損ベクトルからidを再構成
user_id <- user_id0[z_vec==1]
item_id <- item_id0[z_vec==1]
time_id <- time_id0[z_vec==1]


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

##インデックスを設定
user_index <- item_index <- time_index <- list()
ui_id <- ut_id <- iu_id <- it_id <- tu_id <- ti_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
  ut_id[[i]] <- time_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
  it_id[[j]] <- time_id[item_index[[j]]]
}
for(j in 1:time){
  time_index[[j]] <- which(time_id==j)
  tu_id[[j]] <- user_id[time_index[[j]]]
  ti_id[[j]] <- item_id[time_index[[j]]]
}
vec <- rep(1, k)
const1 <- hh / 2.0  #正規化定数
const2 <- item / 2.0

##事前分布の設定
#ユーザーの階層モデルの事前分布
Deltabar1 <- matrix(rep(0, ncol(u)*k), nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta1 <- 0.01 * diag(rep(1, ncol(u)))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#アイテムの階層モデルの事前分布
Deltabar2 <- matrix(rep(0, ncol(v)*k), nrow=ncol(v), ncol=k)   #階層モデルの回帰係数の事前分布の分散
ADelta2 <- 0.01 * diag(rep(1, ncol(v)))   #階層モデルの回帰係数の事前分布の分散
nu2 <- k + 1   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#時間の階層モデルの事前分布
Aelta <- rep(0, k)   #階層モデルの回帰係数の事前分布の分散
ADelta3 <- 0.01 * diag(k)   #階層モデルの回帰係数の事前分布の分散
nu3 <- k + 1   #逆ウィシャート分布の自由度
V3 <- nu3 * diag(rep(1, k))   #逆ウィシャート分布のパラメータ

##パラメータの真値
alpha_u <- alpha_ut; user_mu <- u %*% alpha_u
alpha_v <- alpha_vt; item_mu <- v %*% alpha_v
alpha_t <- alpha_tt
Cov_u <- Cov_ut; inv_Cov1 <- solve(Cov_u)
Cov_v <- Cov_vt; inv_Cov2 <- solve(Cov_v)
Cov_t <- Cov_tt; inv_Cov3 <- solve(Cov_t)
theta_u <- theta_ut
theta_v <- theta_vt
theta_t <- theta_tt
sigma <- 1

##パラメータの初期値
#テンソル分解のパラメータ
theta_u <- mvrnorm(hh, rep(0, k), 0.25 * diag(k))
theta_v <- mvrnorm(item, rep(0, k), 0.25 * diag(k))
theta_t <- mvrnorm(time, rep(0, k), 0.25 * diag(k))
sigma <- 1

#階層モデルのパラメータ
alpha_u <- matrix(0, nrow=ncol(u), ncol=k); user_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=k); item_mu <- v %*% alpha_v
alpha_t <- rep(0, k)
Cov_u <- diag(0.2, k); inv_Cov1 <- solve(Cov_u)
Cov_v <- diag(0.2, k); inv_Cov2 <- solve(Cov_v)
Cov_t <- diag(0.2, k); inv_Cov3 <- solve(Cov_t)

##サンプリング結果のパラメータの保存用配列
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
THETA_T <- array(0, dim=c(time, k, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), k, R/keep))
ALPHA_T <- matrix(0, nrow=R/keep, ncol=k)
COV_U <- array(0, dim=c(k, k, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
COV_T <- array(0, dim=c(k, k, R/keep))

##切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##対数尤度の基準値
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #対数尤度


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  mu <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ] * theta_t[time_id, ]) %*% vec)   #潜在効用の期待値
  U <- rtnorm(mu, sigma, a, b)   #潜在効用を生成
  
  ##テンソル分解の特徴行列のパラメータをサンプリング
  ##ユーザー特徴行列をサンプリング
  for(i in 1:hh){
    #特徴ベクトルのパラメータ
    X <- theta_v[ui_id[[i]], ] * theta_t[ut_id[[i]], ]
    Xy <- t(X) %*% U[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov1)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov1 %*% user_mu[i, ])   #多変量正規分布の平均ベクトル
    cov <- inv_XXV   #多変量正規分布の分散共分散行列
    
    #多変量正規分布からパラメータをサンプリング
    theta_u[i, ] <- mvrnorm(1, beta_mu, cov)
  }
  #theta_u <- theta_u / matrix(colSums(theta_u), nrow=hh, ncol=k, byrow=T) * const1
  
  ##アイテム特徴行列をサンプリング
  for(j in 1:item){
    #特徴ベクトルのパラメータ
    X <- theta_u[iu_id[[j]], ] * theta_t[it_id[[j]], ]
    Xy <- t(X) %*% U[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X +inv_Cov2)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov2 %*% item_mu[j, ])   #多変量正規分布の平均ベクトル
    cov <- inv_XXV   #多変量正規分布の分散共分散行列
    
    #多変量正規分布からパラメータをサンプリング
    theta_v[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  #theta_v <- theta_v / matrix(colSums(theta_v), nrow=item, ncol=k, byrow=T) * const2
  
  ##時間の特徴行列をサンプリング
  for(j in 1:time){
    #特徴ベクトルのパラメータ
    X <- theta_u[tu_id[[j]], ] * theta_v[ti_id[[j]], ]
    Xy <- t(X) %*% U[time_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov3)
    beta_mu <- inv_XXV %*% (Xy + inv_Cov3 %*% alpha_t)   #多変量正規分布の平均ベクトル
    cov <- inv_XXV   #多変量正規分布の分散共分散行列
    
    #多変量正規分布からパラメータをサンプリング
    theta_t[j, ] <- mvrnorm(1, beta_mu, cov)
  }
  
  
  ##階層モデルのパラメータをサンプリング
  ##多変量回帰モデルからユーザーの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_u, u, Deltabar1, ADelta1, nu1, V1)
  alpha_u <- out$B
  user_mu <- u %*% alpha_u   #ユーザー特徴行列の平均構造
  Cov_u <- out$Sigma
  inv_cov1 <- solve(Cov_u)
  
  ##多変量回帰モデルからアイテムの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar2, ADelta2, nu2, V2)
  alpha_v <- out$B
  item_mu <- v %*% alpha_v   #アイテム特徴行列の平均構造
  Cov_v <- out$Sigma
  inv_cov2 <- solve(Cov_v)
  
  ##時間の階層モデルのパラメータをサンプリング
  #多変量正規分布から平均ベクトルをサンプリング
  alpha_mu <- solve(ADelta3 + time*inv_Cov3) %*% (time*inv_Cov3 %*% alpha_t)
  alpha_t <- mvrnorm(1, alpha_mu, Cov_t/time)
  
  #逆ウィシャート分布から分散共分散行列をサンプリング
  er <- theta_t - matrix(alpha_t, nrow=time, ncol=k, byrow=T)   #誤差を算出
  IW <- solve(V3) + t(er) %*% er  
  Sn <- nu3 + time 
  Cov_t <- rwishart(Sn, solve(IW))$IW   #逆ウィシャーと分布からcovをサンプリング
  inv_Cov3 <- solve(Cov_t)
  
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    THETA_T[, , mkeep] <- theta_t
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    ALPHA_T[mkeep, ] <- alpha_t
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    COV_T[, , mkeep] <- Cov_t
  }
  
  #対数尤度を計算
  if(rp%%disp==0){
    mu <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ] * theta_t[time_id, ]) %*% vec)   #完全データの期待値
    prob <- pnorm(mu, 0, sigma)   #購買確率
    LL <- sum(y[index_y1]*log(prob[index_y1])) + sum((1-y[index_y0])*log(1-prob[index_y0]))   #対数尤度
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLst))
    print(round(alpha_u, 3))
  }
}
mean(prob[y==0])

matplot(t(THETA_U[1, , ]), type="l")
matplot(t(THETA_V[1, , ]), type="l")
matplot(t(THETA_T[1, , ]), type="l")

round(prob, 3)


