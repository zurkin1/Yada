SVR <- function(sig_matrix, mixture){
  library("e1071")
  library("preprocessCore")
  X <- sig_matrix
  Y <- mixture
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  tmpc <- colnames(Y)
  tmpr <- rownames(Y)
  Y <- normalize.quantiles(Y)
  colnames(Y) <- tmpc
  rownames(Y) <- tmpr
  
  inter_gene <- intersect(rownames(X), rownames(Y))
  X <- X[inter_gene, ]
  Y <- Y[inter_gene, ]
  
  X <- (X-mean(X))/sd(as.vector(X))
  
  rmse_svm <- rep(0,3)
  corr_svm <- rep(0,3)
  rmse_mix <- c()
  corr_mix <- c()
  for(i in 1:ncol(Y)){
    t <- 1
    y <- Y[,i]
    y <- (y-mean(y))/sd(y)
    for(j in seq(0.25, 0.75, 0.25)){
      model <- svm(X, y, scale = F, type = "nu-regression", kernel = "linear", nu = j)
      weights <- t(model$coefs) %*% model$SV
      weights[which(weights<0)] <- 0
      w <- weights/sum(weights)
      u <- sweep(X, MARGIN=2, w, '*')
      k <- apply(u, 1, sum)
      rmse_svm[t] <- sqrt((mean((k-y)^2)))
      t <- t+1
    }
    mn <- which.min(rmse_svm)
    if(mn==1) model <- svm(X, y, scale=F, type="nu-regression", kernel="linear", nu=0.25)
    if(mn==2) model <- svm(X, y, scale=F, type="nu-regression", kernel="linear", nu=0.5)
    if(mn==3) model <- svm(X, y, scale=F, type="nu-regression", kernel="linear", nu=0.75)
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)] <- 0
    if(i==1) f <- (q/sum(q))
    else f <- rbind(f, (q/sum(q)))
  }
  f=f*100
  result <- f
  colnames(result) <- c(colnames(X))
  rownames(result) <- colnames(Y)
  result
}