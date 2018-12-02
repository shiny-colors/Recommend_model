#####Bayesian Hierarchical Field Aware Factorization Machines#####
library(MASS)
library(Matrix)
library(matrixStats)
library(data.table)
library(FAdist)
library(bayesm)
library(float)
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
pt <- rtpois(hh, rgamma(hh, 37.5, 0.25), a=1, b=Inf)   #購買接触数
hhpt <- sum(pt)
n <- rtpois(item, 1.75, a=0, b=7)   #タグ数
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

##素性ベクトルを生成
k1 <- 2; k2 <- 3; k3 <- 4
x1 <- matrix(runif(hhpt*k1, 0, 1), nrow=hhpt, ncol=k1)
x2 <- matrix(0, nrow=hhpt, ncol=k2)
for(j in 1:k2){
  pr <- runif(1, 0.25, 0.55)
  x2[, j] <- rbinom(hhpt, 1, pr)
}
x3 <- rmnom(hhpt, 1, runif(k3, 0.2, 1.25)); x3 <- x3[, -which.min(colSums(x3))]
x <- cbind(1, x1, x2, x3)   #データを結合
column <- ncol(x)

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
u_col <- ncol(u)

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
v_col <- ncol(v)


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
item_data <- sparseMatrix(1:hhpt, item_id, x=rep(1, hhpt), dims=c(hhpt, item))
item_n <- colSums(item_data)

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
    zeta <- as.numeric(rmnom(1, n[j], omega[z[j], ]))
    if(max(zeta)==1){
      tag_list[[j]] <- (zeta * 1:tag)[zeta > 0]
      break
    }
  }
  tag_id[item_list[[j]], 1:n[j]] <- matrix(tag_list[[j]], nrow=length(item_list[[j]]), ncol=n[j], byrow=T)
}
tag_data <- sparseMatrix(matrix(rep(1:hhpt, max_n), nrow=hhpt, ncol=max_n)[tag_id!=0], tag_id[tag_id!=0], 
                         x=rep(1, sum(tag_id!=0)), dims=c(hhpt, tag))
tag_id0 <- tag_id; tag_id0[tag_id0==0] <- tag + 1
tag_n <- colSums(tag_data)


#タグを縦持ちのIDを設定
no_vec <- rep(1:hhpt, max_n)[as.numeric(tag_id) > 0]
tag_vec <- as.numeric(tag_id)[as.numeric(tag_id) > 0]
user_vec <- rep(user_id, max_n)[as.numeric(tag_id) > 0]
f <- length(user_vec)
v_dt <- sparseMatrix(rep(1:hhpt, max_n)[as.numeric(tag_id) > 0], 1:f, x=rep(1, f), dims=c(hhpt, f))

#インデックスの設定
ug_id2 <- ug_id3 <- g_id <- g_dt <- list()
for(j in 1:tag){
  unique_no <- unique(no_vec[which(tag_vec==j)])
  index <- which(no_vec %in% unique_no)
  ug_id2[[j]] <- user_id[which(rowSums(tag_id==j)==1)]
  ug_id3[[j]] <- user_vec[index]
  g_id[[j]] <- tag_vec[index]
  g_dt[[j]] <- v_dt[sort(unique_no), index]
}


#組み合わせを作成
index_combine <- t(combn(1:max_n, m=2))
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
  
  ##素性ベクトルの回帰パラメ－タ
  beta <- betat <- c(-0.7, runif(ncol(x)-1, -1.25, 1.25))
  Sigma <- Sigmat <- 1   #モデルの標準偏差
  
  ##ユーザーベースの特徴行列のパラメータ
  #分散共分散行列を設定
  Cov_ut1 <- Cov_u1 <- diag(runif(k, 0.01, 0.25))
  Cov_ut2 <- Cov_u2 <- diag(runif(k, 0.01, 0.25))
  
  #回帰係数を設定
  alpha_u1 <- alpha_u2 <- matrix(0, nrow=ncol(u), ncol=k)
  for(j in 1:ncol(u)){
    if(j==1){
      alpha_u1[j, ] <- runif(k, -0.8, 0.2)
      alpha_u2[j, ] <- runif(k, -0.7, 0.2)
    } else {
      alpha_u1[j, ] <- runif(k, -0.7, 0.7)
      alpha_u2[j, ] <- runif(k, -0.7, 0.7)
    }
  }
  alpha_ut1 <- alpha_u1; alpha_ut2 <- alpha_u2
  
  #多変量回帰モデルから特徴行列を生成
  theta_u1 <- theta_ut1 <- u %*% alpha_u1 + mvrnorm(hh, rep(0, k), Cov_u1)   #アイテムとの交互作用特徴行列
  theta_u2 <- theta_ut2 <- u %*% alpha_u2 + mvrnorm(hh, rep(0, k), Cov_u2)   #タグとの交互作用特徴行列
  
  
  ##アイテムベースの特徴行列のパラメータ
  #分散共分散行列を設定
  Cov_vt <- Cov_v <- diag(runif(k, 0.01, 0.25))
  
  #回帰係数を設定
  alpha_v <- matrix(0, nrow=ncol(v), ncol=k)
  for(j in 1:ncol(v)){
    if(j==1){
      alpha_v[j, ] <- runif(k, -0.8, 0.2)
    } else {
      alpha_v[j, ] <- runif(k, -0.8, 0.8)
    }
  }
  alpha_vt <- alpha_v
  
  #多変量回帰モデルから特徴行列を生成
  theta_v <- theta_vt <- v %*% alpha_v + mvrnorm(item, rep(0, k), Cov_v)
  
  
  ##タグベースの特徴行列のパラメータ
  #分散共分散行列の設定
  Cov_g1 <- Cov_gt1 <- diag(runif(k, 0.01, 0.3))
  Cov_g2 <- Cov_gt2 <- diag(runif(k, 0.01, 0.25))
  
  #多変量正規分布から特徴行列を生成
  theta_g1 <- theta_gt1 <- mvrnorm(tag, rep(0, k), Cov_g1)
  theta_g2 <- theta_gt2 <- mvrnorm(tag, rep(0, k), Cov_g2)
  
  
  ##行列分解のパラメータを設定
  #ユーザー、アイテム、タグの行列分解のパラメータ
  UV <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec_k)
  UG <- as.numeric(v_dt %*% ((theta_u2[user_vec, ] * rbind(theta_g1, 0)[tag_vec, ]) %*% vec_k))
  
  #タグ間の行列分解のパラメータ
  WH <- rep(0, hhpt)
  for(j in 1:length(combine_n)){
    W <- rbind(theta_g2, 0)[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% rbind(theta_g2, 0)[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((W * H) %*% vec_k)
  }
  
  ##正規分布から効用と購買ベクトルを生成
  #潜在効用を生成
  mu <- mut <- as.numeric(x %*% beta) + UV + UG + WH
  U <- mu + rnorm(hhpt, 0, Sigma)
  
  #購買ベクトルに変換
  y <- ifelse(U > 0, 1, 0)
  Prob <- pnorm(mu, 0, 1)

  #break条件
  print(mean(y))
  if(mean(y) > 0.2 & mean(y) < 0.4) break   #break条件
}

#要約統計量
mean(y)   #全体の購買確率
as.numeric(tapply(y, user_id, mean))   #ユーザーごとの購買確率
as.numeric(tapply(y, item_id, mean))   #アイテムごとの購買確率
hist(Prob, breaks=25, col="grey", xlab="購買確率", main="購買確率の理論値")


####マルコフ連鎖モンテカルロ法でField Aware Factorization Machinesを推定####
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


##データとインデックスを設定
#ユーザーとアイテムのインデックスを設定
uv_id1 <- uv_id2 <- ug_id1 <- list() 
for(i in 1:hh){
  index <- user_list[[i]]
  uv_id1[[i]] <- item_id[index]
  ug_id1[[i]] <- tag_id[index, ]
}
for(j in 1:item){
  index <- item_list[[j]]
  uv_id2[[j]] <- user_id[index]
}

#タグのインデックスを設定
tag_allocation <- list()
for(i in 1:hh){
  tag_no <- as.numeric((tag_id[user_list[[i]], ] > 0) * matrix(1:pt[i], nrow=pt[i], nco=max_n))
  tag_allocation[[i]] <- sparseMatrix(tag_no[tag_no > 0], 1:sum(tag_no > 0),
                                      x=rep(1, sum(tag_no > 0)), dims=c(pt[i], sum(tag_no > 0)))
}
tag_list <- dt_list <- list()
tag_pt <- rep(0, tag)
for(i in 1:tag){
  index1 <- index2 <- c()
  for(j in 1:max_n){
    index1 <- c(index1, which(tag_id[, j]==i))
    index2 <- c(index2, index_n[which(tag_id[index_n, j]==i)])
  }
  tag_list[[i]] <- sort(index1)
  tag_pt[i] <- length(tag_list[[i]])
  dt_list[[i]] <- sparseMatrix(rep(1:tag_pt[i], max_n), 1:(tag_pt[i]*max_n), x=rep(1, tag_pt[i]*max_n), 
                               dims=c(tag_pt[i], tag_pt[i]*max_n))
}


##事前分布の設定
#素性ベクトルの事前分布
alpha <- rep(0, column)
inv_tau <- solve(diag(100, column))

#ユーザーの階層モデルの事前分布
Deltabar_u <- matrix(0, nrow=u_col, 2*k) 
ADelta_u <- 0.01 * diag(rep(1, u_col))
nu1 <- k + 1 
V1 <- diag(rep(1, 2*k))

#アイテムの階層モデルの事前分布
Deltabar_v <- matrix(0, nrow=v_col, k) 
ADelta_v <- 0.01 * diag(rep(1, v_col))
nu2 <- k + 1 
V2 <- diag(rep(1, k))

#タグの階層モデルの事前分布
g_mu <- rep(0, 2*k)
nu3 <- k + 1
V3 <- diag(rep(1, 2*k))


##パラメータの真値
#素性ベクトルの真値
beta <- betat
Sigma <- Sigmat

#階層モデルの真値
alpha_u <- cbind(alpha_ut1, alpha_ut2); u_mu <- u %*% alpha_u
alpha_v <- alpha_vt; v_mu <- v %*% alpha_v
Cov_u <- as.matrix(bdiag(list(Cov_ut1, Cov_ut2))); inv_Cov_u <- solve(Cov_u)
Cov_v <- Cov_vt; inv_Cov_v <- solve(Cov_v)
Cov_g <- as.matrix(bdiag(list(Cov_gt1, Cov_gt2))); inv_Cov_g <- solve(Cov_g)

#特徴行列のパラメータの真値
theta_u1 <- theta_ut1; theta_u2 <- theta_ut2
theta_v <- theta_vt
theta_g1 <- theta_gt1; theta_g2 <- theta_gt2


##パラメータの初期値
#素性ベクトルの初期値
beta <- rep(0, column)
Sigma <- 1

#階層モデルの初期値
alpha_u <- matrix(0, nrow=u_col, ncol=2*k); u_mu <- u %*% alpha_u
alpha_v <- matrix(0, nrow=v_col, ncol=k); v_mu <- v %*% alpha_v
Cov_u <- diag(rep(0.01, 2*k)); inv_Cov_u <- solve(Cov_u)
Cov_v <- diag(rep(0.01, k)); inv_Cov_v <- solve(Cov_v)
Cov_g1 <- Cov_g2 <- diag(rep(0.01, k)); inv_Cov_g1 <- solve(Cov_g1); inv_Cov_g2 <- solve(Cov_g2)

#特徴行列のパラメータの初期値
theta_u1 <- mvrnorm(hh, rep(0, k), Cov_u1)
theta_u2 <- mvrnorm(hh, rep(0, k), Cov_u2)
theta_v <- mvrnorm(item, rep(0, k), Cov_v)
theta_g1 <- mvrnorm(tag, rep(0, k), Cov_g1)
theta_g2 <- mvrnorm(tag, rep(0, k), Cov_g2)


##行列分解のパラメータを設定
#ユーザー、アイテム、タグの行列分解のパラメータ
UV <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec_k)
UG <- as.numeric(v_dt %*% ((theta_u2[user_vec, ] * theta_g1[tag_vec, ]) %*% vec_k))

#タグ間の行列分解のパラメータ
WH <- rep(0, hhpt)
for(j in 1:length(combine_n)){
  W <- rbind(theta_g2, 0)[tag_id0[index_n, j], ]
  H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% rbind(theta_g2, 0)[tag_id0[index_n, combine_list[[j]]], ])
  WH[index_n] <- WH[index_n] + as.numeric((W * H) %*% vec_k)
}

##サンプリング結果のパラメータの保存用配列
BETA <- matrix(0, nrow=R/keep, ncol=column)
THETA_U1 <- array(0, dim=c(hh, k, R/keep))
THETA_U2 <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
THETA_G1 <- array(0, dim=c(tag, k, R/keep))
THETA_G2 <- array(0, dim=c(tag, k, R/keep))
ALPHA_U <- array(0, dim=c(u_col, 2*k, R/keep))
ALPHA_V <- array(0, dim=c(v_col, k, R/keep))
COV_U <- array(0, dim=c(2*k, 2*k, R/keep))
COV_V <- array(0, dim=c(k, k, R/keep))
COV_G <- array(0, dim=c(2*k, 2*k, R/keep))

##切断領域を定義
index_y1 <- which(y==1)
index_y0 <- which(y==0)
a <- ifelse(y==0, -100, 0)
b <- ifelse(y==1, 100, 0)

##対数尤度の基準値
#1パラメータモデルの対数尤度
prob <- mean(y)
LLst <- sum(y*log(prob)) + sum((1-y)*log(1-prob))   #対数尤度

#真値での対数尤度
prob <- pnorm(mut)
prob[prob==1] <- 0.999999; prob[prob==0] <- 0.000001
LLbest <- sum(y*log(prob) + (1-y)*log(1-prob))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布から潜在効用を生成
  mu <- as.numeric(x %*% beta) + UV + UG + WH   #潜在効用の期待値
  U <- rtnorm(mu, Sigma, a, b)   #潜在効用を生成

  ##素性ベクトルのパラメータをサンプリング
  #モデル誤差を設定
  er <- U - UV - UG - WH
  
  #回帰ベクトルの事後分布のパラメータ
  Xy <- t(x) %*% er
  inv_XXV <- solve(t(x) %*% x + inv_tau)
  mu_vec <- inv_XXV %*% (Xy + inv_tau %*% alpha)

  #多変量正規分布からbetaをサンプリング
  beta <- mvrnorm(1, mu_vec, inv_XXV)
  beta_mu <- as.numeric(x %*% beta)
  
  
  ##ユーザーの特徴行列をサンプリング
  #モデル誤差を設定
  er <- U - beta_mu - WH
  
  for(i in 1:hh){
    #データの設定
    dt1 <- theta_v[uv_id1[[i]], ]
    dt2 <- as.matrix(tag_allocation[[i]] %*% theta_g1[ug_id1[[i]], ])
    dt <- cbind(dt1, dt2)

    #特徴行列の事後分布のパラメータ
    Xy <- t(dt) %*% er[user_list[[i]]]
    inv_XXV <- solve(t(dt) %*% dt + inv_Cov_u)
    mu_vec <- inv_XXV %*% (Xy + inv_Cov_u %*% u_mu[i, ])

    #多変量正規分布から特徴行列をサンプリング
    theta_u <- as.numeric(mvrnorm(1, mu_vec, inv_XXV))
    theta_u1[i, ] <- theta_u[1:k]; theta_u2[i, ] <- theta_u[(k+1):length(theta_u)]
  }
  
  #ユーザーの行列分解のパラメータを更新
  UV <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec_k)
  UG <- as.numeric(v_dt %*% ((theta_u2[user_vec, ] * theta_g1[tag_vec, ]) %*% vec_k))
  
  
  ##アイテムの特徴行列をサンプリング
  #モデル誤差を設定
  er <- U - beta_mu - UG - WH

  for(j in 1:item){
    #特徴行列の事後分布のパラメータ
    dt <- theta_u1[uv_id2[[j]], ]
    Xy <- t(dt) %*% er[item_list[[j]]]
    inv_XXV <- solve(t(dt) %*% dt + inv_Cov_v) 
    mu_vec <- inv_XXV %*% (Xy + inv_Cov_v %*% v_mu[j, ])
    
    #多変量正規分布から特徴行列をサンプリング
    theta_v[j, ] <- as.numeric(mvrnorm(1, mu_vec, inv_XXV))
  }
  
  #アイテムの行列分解のパラメータを更新
  UV <- as.numeric((theta_u1[user_id, ] * theta_v[item_id, ]) %*% vec_k)
  
  
  ##タグの特徴行列のパラメータをサンプリング
  #モデル誤差を設定
  er <- U - beta_mu - UV
  
  for(i in 1:tag){
    #データの設定
    index <- tag_list[[i]]
    tag_id00 <- tag_id0[index, ]
    theta_g01 <- rbind(theta_g1, 0); theta_g01[i, ] <- 0
    theta_g02 <- rbind(theta_g2, 0); theta_g02[i, ] <- 0
    
    #推定対象外のユーザー依存の交互作用ベクトルのパラメータ
    ug <- as.numeric(g_dt[[i]] %*% (theta_u2[ug_id3[[i]], ] * theta_g01[g_id[[i]], ]) %*% vec_k)
    
    #推定対象外のタグ依存の交互作用ベクトルのパラメータ
    wh <- rep(0, length(index))
    for(j in 1:length(combine_n)){
      w <- theta_g02[tag_id00[, j], ]
      h <- as.matrix(dt_list[[i]][, 1:(length(index)*combine_n[j])] %*% theta_g02[tag_id00[, combine_list[[j]]], ])
      wh <- wh + as.numeric((w * h) %*% vec_k)
    }
    
    #特徴行列の事後分布のパラメータ
    er_y <- er[index] - ug - wh
    dt <- cbind(theta_u2[ug_id2[[i]], ], dt_list[[i]] %*% theta_g02[tag_id00, ])
    Xy <- t(dt) %*% er_y
    inv_XXV <- solve(t(dt) %*% dt + inv_Cov_g) 
    mu_vec <- as.numeric(inv_XXV %*% (Xy + inv_Cov_g %*% g_mu))
    
    #多変量正規分布から特徴行列をサンプリング
    theta_g <- as.numeric(mvrnorm(1, mu_vec, inv_XXV))
    theta_g1[i, ] <- theta_g[1:k]; theta_g2[i, ] <- theta_g[(k+1):(2*k)]
  }

  ##行列分解のパラメータを更新
  #ユーザーとタグの行列分解のパラメータを更新
  UG <- as.numeric(v_dt %*% ((theta_u2[user_vec, ] * theta_g1[tag_vec, ]) %*% vec_k))
  
  #タグ間の行列分解のパラメータを更新
  WH <- rep(0, hhpt)
  for(j in 1:length(combine_n)){
    W <- rbind(theta_g2, 0)[tag_id0[index_n, j], ]
    H <- as.matrix(tag_dt[, 1:(tag_n*combine_n[j])] %*% rbind(theta_g2, 0)[tag_id0[index_n, combine_list[[j]]], ])
    WH[index_n] <- WH[index_n] + as.numeric((W * H) %*% vec_k)
  }
  
  
  ##階層モデルのパラメータをサンプリング
  ##多変量回帰モデルからユーザーの階層モデルのパラメータをサンプリング
  out <- rmultireg(cbind(theta_u1, theta_u2), u, Deltabar_u, ADelta_u, nu1, V1)
  alpha_u <- out$B
  u_mu <- u %*% alpha_u   #ユーザー特徴行列の平均構造
  Cov_u <- diag(diag(out$Sigma))
  inv_cov_u <- solve(Cov_u)

  ##多変量回帰モデルからアイテムの階層モデルのパラメータをサンプリング
  out <- rmultireg(theta_v, v, Deltabar_v, ADelta_v, nu2, V2)
  alpha_v <- out$B
  v_mu <- v %*% alpha_v   #アイテム特徴行列の平均構造
  Cov_v <- diag(diag(out$Sigma))
  inv_Cov_v <- solve(Cov_v)
  
  ##タグの階層モデルのパラメータをサンプリング
  #逆ウィシャート分布から分散共分散行列をサンプリング
  er <- cbind(theta_g1, theta_g2) - matrix(colMeans(cbind(theta_g1, theta_g2)), nrow=tag, ncol=2*k, byrow=T)   #誤差を算出
  IW <- t(er) %*% er + V3
  Sn <- nu3 + tag
  Cov_g <- diag(diag(rwishart(Sn, solve(IW))$IW))   #逆ウィシャーと分布からcovをサンプリング
  inv_Cov_g <- solve(Cov_g)
  

  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    THETA_U1[, , mkeep] <- theta_u1
    THETA_U2[, , mkeep] <- theta_u2
    THETA_V[, , mkeep] <- theta_v
    THETA_G1[, , mkeep] <- theta_g1
    THETA_G2[, , mkeep] <- theta_g2
    ALPHA_U[, , mkeep] <- alpha_u
    ALPHA_V[, , mkeep] <- alpha_v
    COV_U[, , mkeep] <- Cov_u
    COV_V[, , mkeep] <- Cov_v
    COV_G[, , mkeep] <- Cov_g
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    mu <- as.numeric(x %*% beta) + UV + UG + WH   #潜在効用の期待値   
    prob <- pnorm(mu, 0, Sigma)   #購買確率
    prob[prob==1] <- 0.999999; prob[prob==0] <- 0.000001
    LL <- sum(y*log(prob) + (1-y)*log(1-prob))   
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLbest, LLst))
    print(round(rbind(beta, betat), 3))
    print(round(cbind(Cov_v, Cov_vt), 3))
  }
}

mean(prob[y==0])

matplot(t(THETA_U[1, , ]), type="l")
matplot(t(THETA_V[1, , ]), type="l")
matplot(t(THETA_T[1, , ]), type="l")

round(prob, 3)


