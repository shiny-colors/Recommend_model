#####非負値多重行列因子分解#####
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

####データの発生####
#データの設定
hh <- 3000   #ユーザー数
item <- 300   #アイテム数
group <- 20   #グループ数
category <- 15   #カテゴリー数
k <- 10   #潜在変数数

#対応行列の設定
pr1 <- runif(group, 0.25, 1.0)
pr2 <- runif(category, 0.25, 1.0)
V <- rmnom(hh, 1, pr1)   #ユーザーグループ対応行列
W <- rmnom(item, 1, pr2)   #商品カテゴリ対応行列


##非負値行列因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
A0 <- matrix(0, nrow=hh, ncol=k)
B0 <- matrix(0, nrow=k, ncol=item)
alpha01 <- matrix(0, nrow=k, ncol=group)
beta01 <- matrix(0, nrow=k, ncol=group)
alpha02 <- matrix(0, nrow=k, ncol=category)
beta02 <- matrix(0, nrow=k, ncol=category)

for(i in 1:k){
  #ユーザー特徴行列を生成
  for(j1 in 1:group){
    alpha01[i, j1] <- runif(1, 0.01, 0.4)
    beta01[i, j1] <- runif(1, 0.6, 1.8)
    A0[V[, j1]==1, i] <- rgamma(sum(V[, j1]), alpha01[i, j1], beta01[i, j1])
  }
  #アイテム特徴行列を生成
  for(j2 in 1:category){
    alpha02[i, j2] <- runif(1, 0.01, 0.3)
    beta02[i, j2] <- runif(1, 0.6, 1.5)
    B0[i, W[, j2]==1] <- rgamma(sum(W[, j2]), alpha02[i, j2], beta02[i, j2])
  }
}

#ポアソン分布の期待値を設定
AB <- A0 %*% B0

#グループ特徴行列およびカテゴリー特徴行列を生成
C0 <- t(t(W) %*% t(B0))
D0 <- t(V) %*% A0 



##カルバックライブラーダイバージェンス距離に基づき応答変数を生成
#ユーザー商品購買行列を生成
X <- matrix(0, nrow=hh, ncol=item)
for(j in 1:item){
  X[, j] <- rpois(hh, AB[, j])
}
rowSums(X)
sum(X)

#ユーザーカテゴリ購買行列を設定
Y <- matrix(0, nrow=hh, ncol=category)
for(j in 1:category){
  Y[, j] <- rowSums(X[, W[, j]==1])
}
colSums(Y)

#グループ商品購買行列を設定
Z <- matrix(0, nrow=group, ncol=item)
for(j in 1:group){
  Z[j, ] <- colSums(X[V[, j]==1, ])
}

#真のパラメータでの対数尤度
LL <- sum(dpois(X, A0 %*% B0, log=TRUE))


####マルコフ連鎖モンテカルロ法でMMNFを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
iter <- 0
disp <- 10

##事前分布の設定
alpha1 <- 0.01; beta1 <- 0.01
alpha2 <- 0.01; beta2 <- 0.01

##初期値の設定
A <- matrix(rgamma(hh*k, 0.25, 0.5), nrow=hh, ncol=k)
B <- matrix(rgamma(item*k, 0.25, 0.5), nrow=k, ncol=item)
C <- matrix(rgamma(category*k, 0.25, 0.5), nrow=k, ncol=category)
D <- matrix(rgamma(group*k, 0.25, 0.5), nrow=group, ncol=k)


##サンプリング結果の保存用配列
A_array <- array(0, dim=c(hh, k, R/keep))
B_array <- array(0, dim=c(k, item, R/keep))
C_array <- array(0, dim=c(k, category, R/keep))
D_array <- array(0, dim=c(group, k, R/keep))


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##補助変数lambdaを更新
  lambda1 <- array(0, dim=c(hh, item, k))
  lambda2 <- array(0, dim=c(hh, category, k))
  AB <- A %*% B
  AC <- A %*% C
  for(j in 1:k){
    lambda1[, , j] <- A[, j] %*% t(B[j, ]) / AB
    lambda2[, , j] <- A[, j] %*% t(C[j, ]) / AC
  }
  
  ##ガンマ分布よりユーザー特徴行列をサンプリング
  for(j in 1:k){
    a1 <- alpha1 + (rowSums(lambda1[, , j] * X) + rowSums(lambda2[, , j] * Y))/2
    a2 <- sum(B[j, ])
    A[, j] <- rgamma(hh, a1, a2)   
  }
  A <- A / matrix(colSums(A), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  ##補助変数lambdaを更新
  lambda1 <- array(0, dim=c(hh, item, k))
  lambda2 <- array(0, dim=c(group, item, k))
  AB <- A %*% B
  BD <- D %*% B
  for(j in 1:k){
    lambda1[, , j] <- A[, j] %*% t(B[j, ]) / AB
    lambda2[, , j] <- D[, j] %*% t(B[j, ]) / BD
  }
  
  ##ガンマ分布よりアイテム特徴行列をサンプリング
  for(j in 1:k){
    b1 <- alpha2 + (colSums(lambda1[, , j] * X) + colSums(lambda2[, , j] * Z))/2
    b2 <- beta2 + sum(A[, j])
    B[j, ] <- rgamma(item, b1, b2)  
  }
  
  ##線形制約をつけたグループ特徴行列およびカテゴリー特徴行列を更新
  C <- t(t(W) %*% t(B))   #カテゴリー特徴行列の更新
  D <- t(V) %*% A   #グループ特徴行列の更新
  
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    A_array[, , mkeep] <- A
    B_array[, , mkeep ] <- B
    C_array[, , mkeep] <- C
    D_array[, , mkeep ] <- D
    if(rp%%disp==0){
      print(c(sum(dpois(X, A %*% B, log=T)), LL))
      print(rp)
    }
  }
}

####サンプリング結果の要約と可視化####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果の可視化
#基底ベクトルAのパラメータの可視化
matplot(t(A_array[1:10, 1, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[1:10, 2, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[1:10, 3, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[1:10, 4, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[1:10, 5, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[11:20, 6, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[11:20, 7, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[11:20, 8, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[11:20, 9, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")
matplot(t(A_array[11:20, 10, ]), type="l", xlab="サンプリング回数", ylab="基底ベクトルのパラメータ")

#結合係数Bのパラメータの可視化
matplot(t(B_array[1, 1:10, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[2, 1:10, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[3, 1:10, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[4, 1:10, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[5, 1:10, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[6, 11:20, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[7, 11:20, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[8, 11:20, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[9, 11:20, ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(B_array[10, 11:20, ]), type="l", xlab="サンプリング回数", ylab="結合係数")

#結合係数Bのパラメータの可視化
matplot(t(C_array[1, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[2, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[3, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[4, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[5, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[6, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[7, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[8, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[9, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(C_array[10, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")

#結合係数Bのパラメータの可視化
matplot(t(D_array[1, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[2, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[3, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[4, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[5, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[6, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[7, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[8, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[9, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")
matplot(t(D_array[10, , ]), type="l", xlab="サンプリング回数", ylab="結合係数")