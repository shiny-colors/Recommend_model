#####Bayesian Hierarchical Factorization Machines#####
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

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
s1 <- 10; vec1 <- rep(1, s1)   #行列分解の基底数
s2 <- 5; vec2 <- rep(1, s2)   #交互作用の基底数
hh <- 5000   #ユーザー数
item <- 3000   #アイテム数
pt <- rtpois(hh, rgamma(hh, 35.0, 0.2), a=1, b=Inf)   #購買接触数
f <- sum(pt)
vec_s1 <- rep(1, s1)
vec_s2 <- rep(1, s2)

#IDを設定
user_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:f, user_id, rank)))
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

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

#スパース行列を作成
user_data <- sparseMatrix(1:f, user_id, x=rep(1, f), dims=c(f, hh))
user_data_T <- t(user_data)
item_data <- sparseMatrix(1:f, item_id, x=rep(1, f), dims=c(f, item))
item_data_T <- t(item_data)

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


##説明変数の生成
#モデルの説明変数を生成
k1 <- 3; k2 <- 4; k3 <- 4
k <- k1 + k2 + k3
x1 <- matrix(0, nrow=f, ncol=k1)
x2 <- matrix(0, nrow=f, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  x1[, j] <- rbeta(f, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(f, 1, pr)
}
x3 <- rmnom(f, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
z <- cbind(x1, x2 ,x3)

#ユーザーの説明変数を生成
k1 <- 3; k2 <- 3; k3 <- 4
u1 <- matrix(0, nrow=hh, ncol=k1)
u2 <- matrix(0, nrow=hh, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  u1[, j] <- rbeta(hh, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  u2[, j] <- rbinom(hh, 1, pr)
}
u3 <- rmnom(hh, 1, runif(k3, 0.2, 1.25)); u3 <- u3[, -which.min(colSums(u3))]
u <- cbind(1, u1, u2, u3)   #データを結合


#アイテムの説明変数を生成
k1 <- 3; k2 <- 2; k3 <- 4
v1 <- matrix(0, nrow=item, ncol=k1)
v2 <- matrix(0, nrow=item, ncol=k2)
for(j in 1:k1){
  par <- runif(2, 1.0, 2.5)
  v1[, j] <- rbeta(item, par[1], par[2])
}
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  v2[, j] <- rbinom(item, 1, pr)
}
v3 <- rmnom(item, 1, runif(k3, 0.2, 1.25)); v3 <- v3[, -which.min(colSums(v3))]
v <- cbind(1, v1, v2, v3)   #データを結合


##交差項の設定
#組み合わせを作成
index_combine <- t(combn(c(1:ncol(z)), m=2))
combine_list <- list()
combine_n <- rep(0, max(index_combine[, 1]))
for(j in 1:max(index_combine[, 1])){
  combine_list[[j]] <- index_combine[which(index_combine[, 1]==j), 2]
  combine_n[j] <- length(combine_list[[j]])
}

#交差項の配列を設定
z_list <- list()
for(j in 1:length(combine_list)){
  z_list[[j]] <- z[, j] * z[, combine_list[[j]]]
}
zz <- do.call(cbind, z_list)


####応答変数を生成####
rp <- 0
repeat { 
  print(rp <- rp + 1)
  
  ##パラメータと応答変数を生成
  #モデルの標準偏差
  sigma <- sigmat <- 1.0
  
  #階層モデルの分散パラメータ
  Cov_x <- Cov_xt <- runif(ncol(x), 0.05, 0.25) * diag(ncol(x))
  Cov_u <- Cov_ut <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_v <- Cov_vt <- runif(s1, 0.05, 0.20) * diag(s1)
  Cov_z <- Cov_zt <- runif(s2, 0.05, 0.15) * diag(s2)
  
  #階層モデルの回帰係数を設定
  alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x))
  alpha_u <- matrix(0, nrow=ncol(u), ncol=s1)
  alpha_v <- matrix(0, nrow=ncol(v), ncol=s1)
  alpha_z <- array(0, dim=c(ncol(u), s2, k-1))
  
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.4, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.225), nrow=s2, ncol=k-1)
    } else {
      alpha_x[j, ] <- runif(ncol(x), -0.4, 0.3)
      alpha_u[j, ] <- runif(s1, -0.35, 0.2)
      alpha_z[j, , ] <- matrix(rnorm(s2*(k-1), 0, 0.225), nrow=s2, ncol=k-1)
    }
  }
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(s2, -0.5, 0.4)
    } else {
      alpha_v[j, ] <- runif(s2, -0.4, 0.4)
    }
  }
  alpha_xt <- alpha_x; alpha_ut <- alpha_u; alpha_vt <- alpha_v; alpha_zt <- alpha_z   #真値を格納
  
  #多変量回帰モデルからユーザー個別の回帰パラメータを生成
  theta_x <- theta_xt <- u %*% alpha_x + mvrnorm(hh, rep(0, ncol(x)), Cov_x)   #変量効果のパラメータ
  theta_u <- theta_ut <- u %*% alpha_u + mvrnorm(hh, rep(0, s1), Cov_u)   #ユーザーの行列分解のパラメータ
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, s1), Cov_v)   #アイテムの行列分解のパラメータ
  
  theta_z <- array(0, c(hh, s2, k-1))
  for(j in 1:(k-1)){
    theta_z[, , j] <- u %*% alpha_z[, , j] + mvrnorm(hh, rep(0, s2), Cov_z)   #交互作用のパラメータ
  }
  theta_xt <- theta_x; theta_ut <- theta_u; theta_vt <- theta_v; theta_zt <- theta_z
  
  ##正規分布から効用と購買ベクトルを生成
  #変量効果のパラメータ
  x_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  #行列分解のパラメータ
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

  #交互作用のパラメータ
  theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
  for(i in 1:length(combine_n)){
    vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
    for(j in 1:combine_n[[i]]){
      vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
    }
    z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
  }
  
  #潜在効用と応答変数を生成
  mu <- mut <- x_mu + z_vec + uv   #期待値
  U <- mu + rnorm(f, 0, sigma)   #潜在効用を生成
  
  #購買ベクトルに変換
  y <- ifelse(U > 0, 1, 0)
  if(mean(y) > 0.25 & mean(y) < 0.4) break   #break条件
}

#生成した応答変数を確認
mean(y)   #購買確率
prob <- pnorm(U, 0, sigma)   #応答確率
mean(prob[y==1]); mean(prob[y==0])   #購買有無別の応答確率
hist(U, col="grey", main="潜在効用の分布", xlab="潜在効用", breaks=25)


####マルコフ連鎖モンテカルロ法で階層ベイズFactorization Machinesを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  par <- qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma)
  return(par)
}

##アルゴリズムの設定
LL1 <- -100000000   #対数尤度の初期値
R <- 1000
keep <- 2  
iter <- 0
burnin <- 100
disp <- 10

##インデックスとデータの定数を設定
#行列分解のインデックスを作成
user_index <- item_index <- list()
ui_id <- iu_id <- list()
for(i in 1:hh){
  user_index[[i]] <- which(user_id==i)
  ui_id[[i]] <- item_id[user_index[[i]]]
}
for(j in 1:item){
  item_index[[j]] <- which(item_id==j)
  iu_id[[j]] <- user_id[item_index[[j]]]
}

#交互作用項のインデックスを作成
index_allocation1 <- rep(1:(ncol(z)-1), combine_n)
index_allocation2 <- unlist(combine_list)
index_allocation <- matrix(1:ncol(z), nrow=ncol(z), ncol=ncol(z), byrow=T); diag(index_allocation) <- 0

index_z <- index_list1 <- index_list2 <- list()
for(j in 1:(k-1)){
  index_z[[j]] <- sort(c(which(index_allocation1==j), which(index_allocation2==j)))
  index_list1[[j]] <- cbind(index_allocation1[index_z[[j]]], index_allocation2[index_z[[j]]])
  index_list2[[j]] <- cbind(index_allocation1[-index_z[[j]]], index_allocation2[-index_z[[j]]])
}

#交互作用項のデータの設定
zz_array1 <- array(0, dim=c(f, k-2, k-1))
zz_array2 <- array(0, dim=c(f, ncol(zz)-(k-2), k-1))
for(j in 1:(k-1)){
  zz_array1[, , j] <- zz[, index_z[[j]]]
  zz_array2[, , j] <- zz[, -index_z[[j]]]
}

#入力変数の定数を設定
xx_list <- list()
for(i in 1:hh){
  xx_list[[i]] <- t(x[user_index[[i]], ]) %*% x[user_index[[i]], ]
}


##事前分布の設定
#変量効果の階層モデルの事前分布
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の平均
ADelta1 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ

#ユーザーの行列分解の階層モデルの事前分布
Deltabar2 <- matrix(0, nrow=ncol(u), ncol=s1)   #階層モデルの回帰係数の事前分布の平均
ADelta2 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu2 <- s1 + 1   #逆ウィシャート分布の自由度
V2 <- nu2 * diag(rep(1, s1)) #逆ウィシャート分布のパラメータ

#アイテムの行列分解の階層モデルの事前分布
Deltabar3 <- matrix(0, nrow=ncol(v), ncol=s1)   #階層モデルの回帰係数の事前分布の平均
ADelta3 <- 0.01 * diag(1, ncol(v))   #階層モデルの回帰係数の事前分布の分散
nu3 <- s1 + 1   #逆ウィシャート分布の自由度
V3 <- nu3 * diag(rep(1, s1)) #逆ウィシャート分布のパラメータ

#交互作用項の階層モデルの事前分布
Deltabar4 <- matrix(0, nrow=ncol(u), ncol=s2)   #階層モデルの回帰係数の事前分布の平均
ADelta4 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu4 <- s2 + 1   #逆ウィシャート分布の自由度
V4 <- nu4 * diag(rep(1, s2)) #逆ウィシャート分布のパラメータ

#変量効果の階層モデルの事前分布
Deltabar1 <- matrix(0, nrow=ncol(u), ncol=k)   #階層モデルの回帰係数の事前分布の平均
ADelta1 <- 0.01 * diag(1, ncol(u))   #階層モデルの回帰係数の事前分布の分散
nu1 <- k + 1   #逆ウィシャート分布の自由度
V1 <- nu1 * diag(rep(1, k)) #逆ウィシャート分布のパラメータ


##パラメータの真値
#階層モデルのパラメータの真値
Cov_x <- Cov_xt; Cov_u <- Cov_ut; Cov_v <- Cov_vt; Cov_z <- array(Cov_zt, dim=c(s2, s2, k-1))
inv_Cov_x <- solve(Cov_x); inv_Cov_u <- solve(Cov_u); inv_Cov_v <- solve(Cov_v)
inv_Cov_z <- array(0, dim=c(s2, s2, k-1))
for(j in 1:(k-1)){
  inv_Cov_z[, , j] <- solve(Cov_z[, , j])
}
alpha_x <- alpha_xt; alpha_u <- alpha_ut; alpha_v <- alpha_vt; alpha_z <- alpha_zt
x_mu <- u %*% alpha_x; u_mu <- u %*% alpha_u; v_mu <- v %*% alpha_v
z_mu <- array(0, c(hh, s2, k-1))
for(j in 1:(k-1)){
  z_mu[, , j] <- u %*% alpha_z[, , j]
}

#モデルパラメータの真値
sigma <- sigmat
theta_x <- theta_xt; theta_u <- theta_ut; theta_v <- theta_vt; theta_z <- theta_zt
user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))

#行列分解のパラメータ
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

#交互作用項のパラメータ
theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
for(i in 1:length(combine_n)){
  vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
  for(j in 1:combine_n[[i]]){
    vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
  }
  z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
}


##パラメータの初期値
#階層モデルの初期値
Cov_x <- diag(0.01, ncol(x)); Cov_u <- Cov_v <- diag(0.01, s1); Cov_z <- array(diag(0.01, s2), dim=c(s2, s2, k-1))
alpha_x <- matrix(0, nrow=ncol(u), ncol=ncol(x)); x_mu <- u %*% alpha_x
alpha_u <- matrix(0, nrow=ncol(u), ncol=s1); u_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=ncol(v), ncol=s1); v_mu <- v %*% alpha_v
alpha_z <- array(0, dim=c(ncol(u), s2, k-1)); z_mu <- array(0, dim=c(hh, s2, k-1))

#モデルパラメータの初期値
sigma <- 1
theta_x <- mvrnorm(hh, as.numeric(solve(t(x) %*% x) %*% t(x) %*% y), Cov_x) 
theta_u <- mvrnorm(hh, rep(0, s1), Cov_u)
theta_v <- mvrnorm(item, rep(0, s1), Cov_v)
theta_z <- array(0, dim=c(hh, s2, k-1))
for(j in 1:(k-1)){
  theta_z[, , j] <- mvrnorm(hh, rep(0, s2), Cov_z[, , j])
}

#行列分解のパラメータ
uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)

#交互作用項のパラメータ
theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
for(i in 1:length(combine_n)){
  vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
  for(j in 1:combine_n[[i]]){
    vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
  }
  z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
}


##パラメータの格納用配列
#モデルパラメータの格納用配列
d <- 0
THETA_X <- array(0, dim=c(hh, k, R/keep))
THETA_U <- array(0, dim=c(hh, s1, R/keep))
THETA_V <- array(0, dim=c(item, s1, R/keep))
THETA_Z <- array(0, dim=c(hh, s2, k-1))

#階層モデルの格納用配列
ALPHA_X <- array(0, dim=c(ncol(u), k, R/keep))
ALPHA_U <- array(0, dim=c(ncol(u), s1, R/keep))
ALPHA_V <- array(0, dim=c(ncol(v), s1, R/keep))
ALPHA_Z <- array(0, dim=c(ncol(u), s2, k-1, R/keep))
COV_X <- array(0, dim=c(k, k, R/keep))
COV_U <- COV_V <- array(0, dim=c(s1, s1, R/keep))
COV_Z <- array(0, dim=c(s2, s2, k-1, R/keep))


##切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##対数尤度の基準値
#1パラメータモデルの対数尤度
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #対数尤度

#ベストモデルの対数尤度
prob <- pnorm(mut, 0, sigmat)   #購買確率
prob[prob==1] <- 0.9999999; prob[prob==0] <- 0.0000001
LLbest <- sum(y*log(prob) + (1-y)*log(1-prob))   #対数尤度


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
    
  ##切断正規分布より潜在効用を生成
  mu <- user_mu + uv + z_vec   #潜在効用の期待値
  U <- extraDistr::rtnorm(f, mu, sigma, a, b)   #潜在効用を生成
  
  ##ユーザーの回帰ベクトルをサンプリング
  #モデルの応答変数
  y_er <- U - uv - z_vec
  
  for(i in 1:hh){
    #回帰ベクトルの事後分布のパラメータ
    XX <- xx_list[[i]]
    Xy <- t(x[user_index[[i]], ]) %*% y_er[user_index[[i]]]
    inv_XXV <- solve(XX + inv_Cov_x)
    mu <- inv_XXV %*% (Xy + inv_Cov_x %*% x_mu[i, ])   #事後分布の平均
    
    #多変量正規分布から回帰ベクトルをサンプリング
    theta_x[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  user_mu <- as.numeric((x * theta_x[user_id, ]) %*% rep(1, ncol(x)))
  
  
  ##ユーザーの特徴行列をサンプリング
  #モデルの応答変数
  u_er <- U - user_mu - z_vec
  
  for(i in 1:hh){
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_v[ui_id[[i]], ]
    Xy <- t(X) %*% u_er[user_index[[i]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_u)
    mu <- inv_XXV %*% (Xy + inv_Cov_u %*% u_mu[i, ])
    
    #多変量正規分布から特徴ベクトルをサンプリング
    theta_u[i, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  
  ##アイテムの特徴行列をサンプリング
  for(j in 1:item){
    
    #特徴ベクトルの事後分布のパラメータ
    X <- theta_u[iu_id[[j]], ]
    Xy <- t(X) %*% u_er[item_index[[j]]]
    inv_XXV <- solve(t(X) %*% X + inv_Cov_v)
    mu <- inv_XXV %*% (Xy + inv_Cov_v %*% v_mu[j, ])
    
    #多変量正規分布から特徴ベクトルをサンプリング
    theta_v[j, ] <- mvrnorm(1, mu, sigma^2*inv_XXV)
  }
  #行列分解のパラメータを更新
  uv <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec1)
  
  
  ##交互作用項の特徴ベクトルをサンプリング
  #モデルの応答変数
  z_er <- U - user_mu - uv

  for(i in 1:hh){
    #データの抽出
    zz1 <- zz_array1[user_index[[i]], , ]; zz2 <- zz_array2[user_index[[i]], , ]
    theta <- t(theta_z[i, , ])
    er_vec <- z_er[user_index[[i]]]
  
    for(j in 1:(k-1))
      #応答変数の設定
      er <- er_vec - as.numeric(zz2[, , j] %*% (theta[index_list2[[j]][, 1], ] * theta[index_list2[[j]][, 2], ]) %*% vec_s2)
      
      #交互作用項の事後分布のパラメータ
      X <- zz1[, , j] %*% theta[index_list1[[j]][, 2], ]
      Xy <- t(X) %*% er
      inv_XXV <- solve(t(X) %*% X + inv_Cov_z[, , j])
      mu <- as.numeric(inv_XXV %*% (Xy + inv_Cov_z[, , j] %*% z_mu[i, , j]))   #事後分布の平均
      
      #多変量正規分布から交互作用項をサンプリング
      theta_z[i, , j] <- mvrnorm(1, mu, sigma^2*inv_XXV)
      theta[j, ] <- theta_z[i, , j]
  }
  
  #交互作用項のパラメータを更新
  theta_vec <- theta_z[user_id, , ]; z_vec <- rep(0, f)
  for(i in 1:length(combine_n)){
    vv <- matrix(0, nrow=f, ncol=combine_n[[i]])
    for(j in 1:combine_n[[i]]){
      vv[, j] <- (theta_vec[, , i] * theta_vec[, , combine_list[[i]][j]]) %*% vec_s2
    }
    z_vec <- z_vec + as.numeric((z_list[[i]] * vv) %*% rep(1, combine_n[[i]]))
  }
  
  
  ##ユーザーの回帰ベクトルの階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_x, u, Deltabar1, ADelta1, nu1, V1)
  alpha_x <- out$B; x_mu <- u %*% alpha_x   
  Cov_x <- out$Sigma; inv_Cov_x <- solve(Cov_x)
  
  ##ユーザー特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_u, u, Deltabar2, ADelta2, nu2, V2)
  alpha_u <- out$B; u_mu <- u %*% alpha_u   
  Cov_u <- out$Sigma; inv_Cov_u <- solve(Cov_u)
  
  ##アイテムの特徴行列の階層モデルのパラメータをサンプリング
  #多変量回帰モデルからパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar3, ADelta3, nu3, V3)
  alpha_v <- out$B; v_mu <- v %*% alpha_v   
  Cov_v <- out$Sigma; inv_Cov_v <- solve(Cov_v)
  
  ##交互作用項の階層モデルのパラメータをサンプリング
  #多変量回帰モデルパラメータをサンプリング
  for(j in 1:(k-1)){
    out <- rmultireg(theta_z[, , j], u, Deltabar4, ADelta4, nu4, V4)
    alpha_z[, , j] <- out$B; z_mu[, , j] <- u %*% alpha_z[, , j]
    Cov_z[, , j] <- out$Sigma; inv_Cov_z[, , j] <- solve(Cov_z[, , j])
  }
  
  ##パラメータの格納とサンプリング結果の格納
  if(rp%%keep==0){
    mkeep <- rp/keep

    #モデルパラメータの格納
    THETA_X[, , mkeep] <- theta_x
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
    
    #階層モデルの格納用配列
    ALPHA_X[, , mkeep] <- alpha_x
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    ALPHA_Z[, , , mkeep] <- alpha_z
    COV_X[, , mkeep] <- Cov_x
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    COV_Z[, , , mkeep] <- Cov_z
    
    if(rp >= burnin){
      d <- d + 1
      THETA_Z <- THETA_Z + theta_z
    }
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    mu <- user_mu + uv + z_vec   #潜在効用の期待値
    prob <- pnorm(mu, 0, sigma)   #購買確率
    prob[prob==1] <- 0.9999999; prob[prob==0] <- 0.0000001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))   #対数尤度
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLbest, LLst))
  }
}


