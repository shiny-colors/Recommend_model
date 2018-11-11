#####ベイジアン非負値テンソル因子分解#####
library(MASS)
library(matrixStats)
library(FAdist)
library(NMF)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

##データの設定
hh <- 3000   #ユーザー数
item <- 500   #アイテム数
time <- 12   #観測期間数
k <- 10   #基底数

##非負値テンソル因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.3; beta01 <- 0.9
alpha02 <- 0.3; beta02 <- 0.85
alpha03 <- 0.15; beta03 <- 1.5

#ガンマ分布から特徴行列を生成
W0 <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)
H0 <- matrix(rgamma(item*k, alpha02, beta02), nrow=k, ncol=item)
C0 <- matrix(rgamma(time*k, alpha03, beta03), nrow=k, ncol=time)

#ポアソン分布の期待値を計算
WHC0 <- array(0, dim=c(hh, item, time))
for(j in 1:k){
  WHC0 <- WHC0 + W0[, j] %o% H0[j, ] %o% C0[j, ]   #特徴行列の外積
}

#ポアソン分布から購買量の生成
Data <- array(0, dim=c(hh, item, time))
LLbest <- 0
for(rp in 1:time){
  for(j in 1:item){
    Data[, j, rp] <- rpois(hh, WHC0[, j, rp])
  }
}
sum(Data)

#対数尤度の計算
LLbest <- 0
for(j in 1:k){
  LLbest <- LLbest + sum(dpois(Data[, , j], WHC0[, , j], log=TRUE))
}
LLbest


####マルコフ連鎖モンテカルロ法で非負値テンソル分解を推定####
##アルゴリズムの設定
R <- 2000
keep <- 2
disp <- 8
iter <- 0

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01
alpha3 <- 0.01; beta3 <- 0.01

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 1.0), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 1.0), nrow=k, ncol=item)
C <- matrix(rgamma(time*k, 0.1, 1.0), nrow=k, ncol=time)

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
C_array <- array(0, dim=c(k, time, R/keep))
lambda <- array(0, dim=c(hh, item, k))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ガンマ分布よりWをサンプリング
  #ポアソン分布の期待値
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #補助変数を用いてWをサンプリング
  W_old <- W   #1つ前のWの特徴行列を格納
  for(j in 1:k){
    lambda <- W_old[, j] %o% H[j, ] %o% C[j, ] / WHC   #補助変数を更新
    w1 <- alpha1 + rowSums(lambda * Data)
    w2 <- beta1 + sum(H[j, ] %o% C[j, ])
    W[, j] <- rgamma(hh, w1, w2)   #ガンマ分布からWをサンプリング
  }
  W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  ##ガンマ分布よりHをサンプリング
  #ポアソン分布の期待値
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #補助変数を用いてHをサンプリング
  H_old <- H   #1つ前のHの特徴行列を格納
  for(j in 1:k){
    lambda <- W[, j] %o% H_old[j, ] %o% C[j, ] / WHC   #補助変数を更新
    h1 <- alpha2 + rowSums(colSums(lambda * Data))
    h2 <- beta2 + sum(W[, j] %o% C[j, ])
    H[j, ] <- rgamma(item, h1, h2)   #ガンマ分布からHをサンプリング
  }
  
  ##ガンマ分布よりCをサンプリング
  #ポアソン分布の期待値
  WHC <- array(0, dim=c(hh, item, time))
  for(j in 1:k){
    WHC <- WHC + W[, j] %o% H[j, ] %o% C[j, ]
  }
  #補助変数を用いてCをサンプリング
  C_old <- C    #1つ前のCの特徴行列を格納
  for(j in 1:k){
    lambda <- W[, j] %o% H[j, ] %o% C_old[j, ] / WHC   #補助変数を更新
    c1 <- alpha3 + colSums(colSums(lambda * Data))
    c2 <- beta3 + sum(W[, j] %o% H[j, ])
    C[j, ] <- rgamma(time, c1, c2)   #ガンマ分布からHをサンプリング
  }
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep ] <- H
    C_array[, , mkeep] <- C
  }
    
  if(rp%%disp==0){
    #対数尤度を計算
    LL <- 0
    for(j in 1:k){
      LL <- LL + sum(dpois(Data[, , j], WHC[, , j], log=TRUE))
    }
    print(rp)
    print(c(LL, LLbest))
  }
}
