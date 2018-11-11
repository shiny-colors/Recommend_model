#####Incomplete data Non Negative Matrix Factorization#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(NMF)
library(stringr)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(78594)

####データの発生####
#データの設定
k <- 10   #基底数
hh <- 5000   #ユーザー数
item <- 2000  #アイテム数
pt <- rtpois(hh, rgamma(hh, 27.5, 0.25), a=1, b=Inf)   #接触数
hhpt <- sum(pt)   #総接触数
vec_k <- rep(1, k)

#IDを設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))
ID <- data.frame(no=1:hhpt, id=user_id, t=pt_id)   #データの結合
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

#生成したデータを可視化
freq_item <- plyr::count(item_id); freq_item$x <- as.character(freq_item$x)
hist(freq_item$freq, breaks=25, col="grey", xlab="アイテムの購買頻度", main="アイテムの購買頻度分布")
gc(); gc()


####応答変数を生成####
rp <- 0
repeat {
  rp <- rp + 1
  
  ##NMFのパラメータを生成
  #ガンマ分布の尺度パラメータを設定
  lambda_u <- abs(rnorm(1, 1, 0.5))
  lambda_v <- abs(rnorm(1, 1, 0.5))
  
  #ガンマ分布の形状パラメータ
  beta_u <- beta_ut <- abs(rnorm(1, 1, 0.5))
  beta_v <- beta_vt <- abs(rnorm(1, 1, 0.5))
  
  #ガンマ分布から行列分解のパラメータを生成
  theta_u <- theta_ut <- matrix(rgamma(hh*k, lambda_u, beta_u), nrow=hh, ncol=k)
  theta_v <- theta_vt <- matrix(rgamma(item*k, lambda_v, beta_v), nrow=item, ncol=k)
  
  ##ポアソン分布から応答変数を生成
  WH <- (theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k   #期待値
  y　<- rpois(hhpt, WH)   
  
  #break条件
  print(rp)
  print(max(y))
  if(max(y) < 75 & max(y) > 25 & sd(y) > 2.5){
    break
  }
}

#生成したデータのヒストグラム
hist(y, main="購買頻度の分布", xlab="購買頻度", col="grey", breaks=50)


####マルコフ連鎖モンテカルロ法でBayesian Hierarchical NMFを推定####
##アルゴリズムの設定
R <- 3000
keep <- 4
burnin <- 500/keep
iter <- 0
disp <- 10

##事前分布の設定
alpha1 <- 0.1; beta1 <- 1
alpha2 <- 0.1; beta2 <- 1

##パラメータの真値
theta_u <- theta_ut
theta_v <- theta_vt
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #期待値

##初期値の設定
theta_u <- matrix(rgamma(hh*k, 0.5, 1.0), nrow=hh, ncol=k)
theta_v <- matrix(rgamma(item*k, 0.5, 1.0), nrow=item, ncol=k)
WH <- as.numeric((theta_u[user_id, ] * theta_v[item_id, ]) %*% vec_k)   #期待値


##サンプリング結果の保存用配列
THETA_U <- array(0, dim=c(hh, k, R/keep))
THETA_V <- array(0, dim=c(item, k, R/keep))
LAMBDA <- array(0, dim=c(hh, item, k))

##ユーザーおよびアイテムのインデックスを作成
#個別に和を取るためのスパース行列
user_dt <- sparseMatrix(sort(user_id), unlist(user_list), x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(sort(item_id), unlist(item_list), x=rep(1, hhpt), dims=c(item, hhpt))

##対数尤度の基準値
LLst <- sum(dpois(y, mean(y), log=TRUE))
LLbest <- sum(dpois(y, as.numeric((theta_ut[user_id, ] * theta_vt[item_id, ]) %*% vec_k), log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザー特徴行列をサンプリング
  #補助変数lambdaを更新
  theta_vec2 <- theta_v[item_id, ]
  lambda <- (theta_u[user_id, ] * theta_vec2) / WH
  
  #ユーザーごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y   #要素ごとの期待値
  W1 <- as.matrix(user_dt %*% lambda_y + alpha1)
  W2 <- as.matrix(user_dt %*% theta_vec2 + beta1)
  
  #ガンマ分布よりパラメータをサンプリング
  theta_u <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  theta_u <- theta_u / matrix(colSums(theta_u), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  
  ##アイテム特徴行列をサンプリング
  #補助変数lambdaを更新
  theta_vec1 <- theta_u[user_id, ]
  WH <- as.numeric((theta_vec1 * theta_vec2) %*% vec_k)
  lambda <- (theta_vec1 * theta_vec2) / WH
  
  #ユーザーごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y   #要素ごとの期待値
  H1 <- as.matrix(item_dt %*% lambda_y + alpha1)
  H2 <- as.matrix(item_dt %*% theta_vec1 + beta1)
  
  #ガンマ分布よりパラメータをサンプリング
  theta_v <- matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k)
  WH <- as.numeric((theta_vec1 * theta_v[item_id, ]) %*% vec_k)
  
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA_U[, , mkeep] <- theta_u
    THETA_V[, , mkeep] <- theta_v
  }
  
  if(rp%%disp==0){
    #対数尤度の更新
    LL <- sum(dpois(y, WH, log=TRUE))
    
    #サンプリング結果を表示
    print(rp)
    print(c(LL, LLbest, LLst))
  }
}

####サンプリング結果の確認####
##サンプリング結果をプロット
matplot(t(THETA_U[1, , ]), type="l")
matplot(t(THETA_V[1, , ]), type="l")

