#####非負値行列因子分解#####
library(MASS)
library(NMF)
library(reshape2)
library(plyr)
library(lattice)
library(ggplot2)

####データの発生####
seg <- 15   #セグメント数
n <- 50   #セグメントごとのサンプル数
N <- seg*n   #全サンプル数
category <- 100   #カテゴリー数
lambda <- 70   #購買数のポアソン平均
y <- rpois(N, lambda)   #購買発生数
Y <- data.frame(seg=rep(1:seg, rep(n, seg)), y)

#購買頻度行列を発生させる
X.freq <- matrix(0, 0, ncol=category)
for(i in 1:seg){
  p <- runif(category, 0, 5)
  freq <- t(apply(Y[Y$seg==i, ], 1, function(x) rmultinom(1, x[2], p)))
  X.freq <- rbind(X.freq, freq)
}
w.freq <- colSums(X.freq)

####非負値行列因子分解を推定####
#パラメータの初期値を設定
H <- matrix(runif(N*seg, 0, 3), nrow=N, ncol=seg)
U <- matrix(runif(N*seg, 0, 3), nrow=seg, ncol=category) 
HU <- H %*% U

#アルゴリズムの設定
tol <- 0.1
diff <- 100
SE1 <- sum((X.freq - HU)^2)

##パラメータを更新(二乗誤差基準)
while(diff >= tol){
  #アルゴリズム
  U1 <- U * t(H) %*% X.freq / t(H) %*% H %*% U
  H1  <- H * X.freq %*% t(U1) / H %*% U1 %*% t(U1)
  
  #パラメータを更新
  HU <- H1 %*% U1
  SE <- sum((X.freq - HU)^2)
  U <- U1
  H <- H1
  diff <- abs(SE1-SE) 
  SE1 <- SE
  print(diff)
}

#推定されたパラメータ
d <- cbind(as.numeric(X.freq), round(as.numeric(HU), 3))
round(H, 2)
round(U, 2)
round(colSums(HU), 0)
colSums(X.freq)

#関数で推定
res <- nmf(X.freq, seg)
H_fun <- basis(res)
U_fun <- coef(res)
HU_fun <- H_fun %*% U_fun
f <- cbind(as.numeric(X.freq), round(as.numeric(HU), 3), round(as.numeric(HU_fun), 3))
