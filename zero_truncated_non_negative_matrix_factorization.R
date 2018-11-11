#####Zero Truncated Non Negative Matrix Factorization#####
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
hh <- 5000   #ユーザー数
item <- 2000  #カテゴリー数
hhpt <- hh*item
k <- 10   #潜在変数数
vec_k <- rep(1, k)

##IDの設定
user_id0 <- rep(1:hh, rep(item, hh))
item_id0 <- rep(1:item, hh)

##非負値行列因子分解の仮定に従いデータを生成
#ガンマ分布よりパラメータを設定
alpha01 <- 0.25; beta01 <- 1.0
alpha02 <- 0.15; beta02 <- 0.85
W <- WT <- matrix(rgamma(hh*k, alpha01, beta01), nrow=hh, ncol=k)   #ユーザー特徴行列
H <- HT <- matrix(rgamma(item*k, alpha02, beta02), nrow=item, ncol=k)   #アイテム特徴行列
WH <- as.numeric(t(W %*% t(H)))

#欠損有無のベータ分布のパラメータを設定
beta1 <- rbeta(hh, 9.5, 10.0)   #ユーザ-購買確率
beta2 <- rbeta(item, 7.5, 8.0)   #アイテム購買確率

#ポアソン分布よりデータを生成
y_comp <- rpois(hhpt, WH)


##欠損がある購買データを生成
#欠損ベクトルを生成
z_vec0 <- rbinom(hhpt, 1, beta1[user_id0] * beta2[item_id0])

#欠損インデックス
z_vec <- z_vec0 * y_comp > 0
index_z <- which(z_vec0==1)
index_z1 <- which(z_vec==1)
index_z0 <- which(z_vec==0)
N <- length(index_z1)

#欠損のある購買ベクトル
user_id <- user_id0[index_z1]
item_id <- item_id0[index_z1]
y_vec <- y_comp[index_z1]

#購買ベクトルに変換
y <- z_vec0 * y_comp   #欠損データを0に変換した購買ベクトル(観測された購買ベクトル)

##ベストなパラメータに対する対数尤度
LLc <- sum(dpois(y_comp, as.numeric(t(W %*% t(H))), log=TRUE))   #完全データに対する対数尤度
LLc1 <- sum(dpois(y_comp, as.numeric(t(W %*% t(H))), log=TRUE)[index_z1])   #非ゼロのデータに対する対数尤度
LLc2 <- sum(dpois(y_comp[index_z], as.numeric(t(W %*% t(H)))[index_z], log=TRUE))   #真の観測に対する対数尤度


####マルコフ連鎖モンテカルロ法でNMFを推定####
##アルゴリズムの設定
R <- 5000
keep <- 4
burnin <- 1000/keep
iter <- 0
disp <- 10

##事前分布の設定
alpha1 <- 0.1; beta1 <- 1
alpha2 <- 0.1; beta2 <- 1

##パラメータの真値
W <- WT
H <- HT; H_t <- t(H)
r <- rowMeans(matrix(z_vec, nrow=hh, ncol=item, byrow=T))
z_vec <- rep(0, hhpt); z_vec[index_z1] <- 1

##初期値の設定
W <- matrix(rgamma(hh*k, 0.1, 0.25), nrow=hh, ncol=k)
H <- matrix(rgamma(item*k, 0.1, 0.25), nrow=item, ncol=k); H_t <- t(H)
r <- rowMeans(matrix(z_vec, nrow=hh, ncol=item, byrow=T))
z_vec <- rep(0, hhpt); z_vec[index_z1] <- 1

##サンプリング結果の保存用配列
W_array <- array(0, dim=c(hh, k, R/keep))
H_array <- array(0, dim=c(k, item, R/keep))
Z_data <- matrix(0, nrow=hh, ncol=item)


##ユーザーおよびアイテムのインデックスを作成
#個別に和を取るためのスパース行列
user_dt <- sparseMatrix(user_id, 1:N, x=rep(1, N), dims=c(hh, N))
user_dt_full <- sparseMatrix(user_id0, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
item_dt <- sparseMatrix(item_id, 1:N, x=rep(1, N), dims=c(item, N))
item_dt_full <- sparseMatrix(item_id0, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))

#欠損した値のインデックス
user_vec_full <- rep(1, hh)
item_vec_full <- rep(1, item)
user_z0 <- user_id0[index_z0]



####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##欠損有無の潜在変数Zをサンプリング
  WH_comp <- as.numeric(t(W %*% t(H)))   #完全データの行列分解の期待値
  
  #潜在変数zの割当確率のパラメータ
  r_vec <- r[user_z0]   #混合率のベクトル
  Li_zeros <- exp(-WH_comp[index_z0])   #データがゼロの時の尤度
  Posterior_zeros <- r_vec * Li_zeros   #z=1の事後分布のパラメータ
  z_rate <- Posterior_zeros / (Posterior_zeros + (1-r_vec))   #潜在変数の割当確率
  
  #二項分布から潜在変数zをサンプリング
  z_vec[index_z0] <- rbinom(length(index_z0), 1, z_rate)
  Zi <- matrix(z_vec, nrow=hh, ncol=item, byrow=T)
  r <- rowMeans(Zi)   #混合率を更新
  z_comp <- which(z_vec==1); N_comp <- length(z_comp)
  

  ##ガンマ分布よりユーザー特徴行列Wをサンプリング
  #補助変数lambdaを更新
  H_vec <- H[item_id, ]
  WH <- as.numeric((W[user_id, ] * H_vec) %*% vec_k)    #観測データの行列分解の期待値
  lambda <- (W[user_id, ] * H_vec) / WH

  #ユーザーごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y_vec   #要素ごとの期待値
  W1 <- as.matrix(user_dt %*% lambda_y + alpha1)
  W2 <- as.matrix((user_dt_full %*% (H[item_id0, ] * z_vec)) + beta1)

  #ガンマ分布よりユーザー特徴行列Wをサンプリング
  W <- matrix(rgamma(hh*k, W1, W2), nrow=hh, ncol=k)
  #W <- W / matrix(colSums(W), nrow=hh, ncol=k, byrow=T) * hh/5   #各列ベクトルを正規化
  
  
  ##ガンマ分布よりアイテム特徴行列Hをサンプリング
  #補助変数lambdaを更新
  W_vec <- W[user_id, ]
  WH <- as.numeric((W_vec * H_vec) %*% vec_k)   #観測データの行列分解の期待値
  lambda <- (W_vec * H_vec) / WH
  
  #ユーザーごとのガンマ分布のパラメータを設定
  lambda_y <- lambda * y_vec   #要素ごとの期待値
  H1 <- as.matrix(item_dt %*% lambda_y + alpha2)
  H2 <- as.matrix((item_dt_full %*% (W[user_id0, ] * z_vec)) + beta2)
  
  #ガンマ分布よりユーザー特徴行列Wをサンプリング
  H <- matrix(rgamma(item*k, H1, H2), nrow=item, ncol=k)
  

  ##サンプリング結果の保存と表示
  #サンプリング結果の格納
  if(rp%%keep==0){
    mkeep <- rp/keep
    W_array[, , mkeep] <- W
    H_array[, , mkeep] <- H
    if(rp > burnin){
      Z_data <- Z_data + Zi
    }
  }
  
  #サンプリング結果の表示
  if(rp%%disp==0){
    #対数尤度の計算
    LL <- sum(dpois(y_comp[index_z], as.numeric(t(W %*% t(H)))[index_z], log=TRUE))
    
    #パラメータの表示
    print(rp)
    print(c(mean(z_vec), mean(z_vec0)))
    print(c(LL, LLc2))
  }
}


####サンプリング結果の要約と適合度####
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W %*% H))[-index_z1], log=TRUE))
sum(dpois(as.numeric(t(Data0))[-index_z1], as.numeric(t(W0 %*% H0))[-index_z1], log=TRUE))


