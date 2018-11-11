#####Mixed Probit Matrix Factorization#####
library(MASS)
library(matrixStats)
library(FAdist)
library(mnormt)
library(NMF)
library(bayesm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####多変量正規乱数を発生させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
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
  D <- ifelse(val < 0, val + abs(val) + 0.1, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}


####データの発生####
hh <- 3000
item <- 300
k <- 10   #潜在変数数
n1 <- rep(item, hh)
n2 <- rep(hh, item)

##モデルの仮定に従いデータを生成
#相関行列を設定
Cor <- Cort <- corrM(k, -0.7, 1.0)

#ユーザーアイテム特徴行列のパラメータを生成
W <- WT <- mvrnorm(hh, rep(0, k), Cor)
H <- HT <- matrix(0, nrow=k, ncol=item)
for(j in 1:k){
  par <- rbeta(1, 20.0, 20.0)
  H[j, ] <- HT[j, ] <- rbinom(item, 1, par)
}

#ユーザーとアイテムの変量効果パラメータを生成
theta1 <- thetat1 <- -1.25
theta2 <- thetat2 <- -1.25
sigma1 <- sigmat1 <- 1
sigma2 <- sigmat2 <- 1
alpha <- alphat <- rnorm(hh, theta1, sigma1)
beta <- betat <- rnorm(item, theta2, sigma2)


#プロビットモデルからアイテム購買行列を生成
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
Util <- alpha_matrix + beta_matrix + W %*% H   #効用関数
y0 <- rnorm(hh*item, as.numeric(Util), 1)   #正規分布から購買有無を生成
y <- ifelse(y0 > 0, 1, 0)
Data <- matrix(y, nrow=hh, ncol=item)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")   #スパース行列化


####マルコフ連鎖モンテカルロ法でMixed Probit Matrix Factorizationを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b, hh, item){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(matrix(runif(length(mu)), nrow=hh, ncol=item)*(FB-FA)+FA, mu, sigma))
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
R <- 10000
keep <- 2
disp <- 10
iter <- 0

##事前分布の設定
#変量効果の分散の事前分布
alpha01 <- 100
beta01 <- 100
alpha02 <- 0.01
beta02 <- 0.01

#行列分解のパラメータの事前分布
nu <- k   #逆ウィシャーと分布の自由度
V <- nu * diag(rep(1, k))   #逆ウィシャート分布のパラメータ


##パラメータの真値
W <- WT
H <- HT
r <- rowMeans(H)
alpha <- alphat
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta <- betat
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
util_mu <- alpha_matrix + beta_matrix + W %*% H   #効用関数
Cor <- Cort
cov_inv <- solve(Cor)
theta1 <- -1.25
theta2 <- -1.25
sigma1 <- 1
sigma2 <- 1



##初期値の設定
#変量効果の初期値
user_rnorm <- rnorm(hh, 0, 1)
item_rnorm <- rnorm(item, 0, 1)
user_sort <- order(user_rnorm)
item_sort <- order(item_rnorm)
alpha <- user_rnorm[user_sort][ceiling(rank(rowSums(Data)))]
alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
beta <- item_rnorm[item_sort][ceiling(rank(colSums(Data)))]
beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
Cor <- corrM(k, -0.3, 0.3)
cov_inv <- solve(Cor)
theta1 <- -1
theta2 <- -1
sigma1 <- 1
sigma2 <- 1

#行列分解のパラメータの初期値
res <- nmf(Data+0.1, k)
W <- scale(basis(res))
H <- round(coef(res)*k)
H[H > 1] <- 1


##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
COR <- array(0, dim=c(k, k, R/keep))
ALPHA <- matrix(0, nrow=R/keep, ncol=hh)
BETA <- matrix(0, nrow=R/keep, ncol=item)
THETA <- matrix(0, nrow=R/keep, ncol=2)
SIGMA <- matrix(0, nrow=R/keep, ncol=2)

##切断領域を定義
a <- ifelse(Data==0, -100, 0)
b <- ifelse(Data==1, 100, 0)
a_vec <- as.numeric(a)
b_vec <- as.numeric(b)
const <- -1/2*log(2*pi)   #標準正規分布の対数尤度の定数



####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##切断正規分布よりアイテム購買の効用値をサンプリング
  util_mu <- alpha_matrix + beta_matrix + W %*% H   #効用関数
  util <- rtnorm(util_mu, 1, a, b, hh, item)
  util[is.infinite(util)] <- 0
  
  ##ユーザー特徴行列をサンプリング
  #アイテム特徴行列の誤差を設定
  W0 <- matrix(0, nrow=hh, ncol=k)
  util_fearture <- util - alpha_matrix - beta_matrix
  
  #アイテム特徴行列の事後分布のパラメータ
  XX <- H %*% t(H)
  XXV <- solve(XX + cov_inv)
  
  for(i in 1:hh){
    XXb <- H %*% util_fearture[i, ]
    beta_mean <- beta_mean <- XXV %*% XXb
    
    #多変量正規分布からbetaをサンプリング
    W0[i, ] <- mnormt::rmnorm(1, beta_mean, XXV)
  }
  W <- scale(W0)
  
  ##逆ウィシャート分布から分散共分散行列をサンプリング
  #逆ウィシャート分布のパラメータ
  V_par <- V + t(W) %*% W
  Sn <- nu + hh
  
  #逆ウィシャート分布から分散共分散行列を発生
  Cor <- cov2cor(rwishart(Sn, solve(V_par))$IW)
  inv_cov <- solve(Cor)
  
  
  ##アイテム特徴行列をサンプリング
  for(j in 1:k){
    
    #変数パターンを設定
    H1 <- H0 <- H
    H1[j, ] <- 1
    H0[j, ] <- 0
    
    #パターンごとの対数尤度を計算
    WH0 <- W %*% H0
    WH1 <- W %*% H1
    LL_item0 <- colSums(const -1/2*(util_fearture - WH0)^2)   
    LL_item1 <- colSums(const -1/2*(util_fearture - WH1)^2)
    LLi0 <- cbind(LL_item0, LL_item1)
    LLi <- exp(LLi0 - rowMaxs(LLi0)) * matrix(c(1-r[j], r[j]), nrow=item, ncol=2, byrow=T)   #尤度に変換
    
    #潜在変数の割当確率からHをサンプリング
    z_rate <- LLi[, 2] / rowSums(LLi)  #潜在変数の割当確率 
    H[j, ] <- rbinom(item, 1, z_rate)
    
    #混合率を更新
    r[j] <- mean(H[j, ])
  }
  
  
  ##ユーザー間変量効果をサンプリング
  util_user <-  util - beta_matrix - W %*% H
  
  #ユーザーごとの平均を推定
  mu <- rowMeans(util_user)
  
  #正規分布から変量効果をサンプリング
  mu_par <- (1/sigma1)/(1/sigma1 + n1)*theta1 + n1/(1/sigma1 + n1)*mu
  sigma_par <-  1 / (1/sigma1 + n1[1])
  alpha <- rnorm(hh, mu_par, sigma_par)
  alpha_matrix <- matrix(alpha, nrow=hh, ncol=item)
  
  ##階層モデルの平均と分散をサンプリング
  #平均をサンプリング
  mu <- mean(alpha)
  mu_par <- (1/alpha01)/(1/alpha01 + hh)*beta01 + hh/(1/alpha01 + hh)*mu
  sigma_par <- sigma1 / (1/alpha01+hh)
  theta1 <- rnorm(1, mu_par, sigma_par)
  
  #分散をサンプリング
  r0 <- alpha02 + hh
  s0 <- beta02 + (hh-1)*var(alpha)
  sigma1 <- 1/rgamma(1, r0/2, s0/2)
  
  
  ##アイテム間変量効果をサンプリング
  util_item <-  util - alpha_matrix - W %*% H
  
  #ユーザーごとの平均を推定
  mu <- colMeans(util_item)
  
  #正規分布から変量効果をサンプリング
  mu_par <- (1/sigma2)/(1/sigma2 + n2)*theta2 + n2/(1/sigma2 + n2)*mu
  sigma_par <-  1 / (1/sigma2 + n2[1])
  beta <- rnorm(item, mu_par, sigma_par)
  beta_matrix <- matrix(beta, nrow=hh, ncol=item, byrow=T)
  
  ##階層モデルの平均と分散をサンプリング
  #平均をサンプリング
  mu <- mean(beta)
  mu_par <- (1/alpha01)/(1/alpha01 + item)*beta01 + item/(1/alpha01 + item)*mu
  sigma_par <- sigma1 / (1/alpha01+item)
  theta2 <- rnorm(1, mu_par, sigma_par)
  
  #分散をサンプリング
  r0 <- alpha02 + item
  s0 <- beta02 + (item-1)*var(beta)
  sigma2 <- 1/rgamma(1, r0/2, s0/2)
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    COR[, , mkeep] <- Cor
    ALPHA[mkeep, ] <- alpha
    BETA[mkeep, ] <- beta
    THETA[mkeep, ] <- c(theta1, theta2)
    SIGMA[mkeep, ] <- c(sigma1, sigma2)
  
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(c(sum(const -1/2*(util - util_mu)^2), sum(const -1/2*(util - mean(util))^2)))
      #print(round(cbind(W[1:5, ], WT[1:5, ]), 2))
      print(cbind(H[, 1:10], HT[, 1:10]))
      print(round(cbind(Cor, Cort), 2))
      print(round(c(theta1, theta2, sigma1, sigma2, thetat1, thetat2, sigmat1, sigmat2), 3))
    }
  }
}

