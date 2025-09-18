lambda.bound <- function(dt, eps = 1e-6, K = 100) {
  q <- dim(dt[[1]])[2]
  Y_bind <- do.call(rbind, dt)
  
  sigma <- set_init(dt, 1)$weight
  weight <- outer(sqrt(sigma), sqrt(1 / sigma))
  
  Y2 <- matrix(NA, q, q)   # V matrix 기준
  for (j in 1:q) {
    for (k in j:q) {
      tmp <- sum(Y_bind[, j] * Y_bind[, k])
      Y2[j, k] <- Y2[k, j] <- tmp
    }
  }
  
  lambda.upp <- max(abs((weight + t(weight)) * Y2) / q)
  lambda.low <- eps * lambda.upp
  lambda <- exp(seq(log(lambda.low), log(lambda.upp), length.out = K))
  
  return(list(lambda.upp = lambda.upp, lambda.low = lambda.low, lambda = lambda))
}


lambda.function <- function(dt, flag, K) {
  if (flag == "custom") {
    source("./lambda.R")   # Load lambda making function
    
    lambda <- lambda.bound(dt, K = K)$lambda
  } else if (flag == "auto") {
    lambda <- exp(seq(log(1e-6), log(10), length.out = K))
  }
  
  return(list(lambda = lambda))
}



source("./lambda.R")


n_lambda <- 30
target_folder <- c("hub", "random", "cluster")

Y.lambda <- Y_t.lambda <- list()
for (target in target_folder) {
  Y.lambda[[target]] <- vector("list", 50)
  Y_t.lambda[[target]] <- vector("list", 50)
  
  for (i in 1:50) {
    load(sprintf("./%s_data/data%02d.RData", target, i))
    
    Y.lambda[[target]][[i]] <- lambda.bound(result$Y, K = n_lambda)$lambda
    Y_t.lambda[[target]][[i]] <- lambda.bound(result$Y_t, K = n_lambda)$lambda
  }
}

lambda_list.Y <- list(    # Custom lambda / Y (V matrix ONLY)
  hub     = matrix(unlist(Y.lambda$hub), n_lambda),
  random  = matrix(unlist(Y.lambda$hub), n_lambda),
  cluster = matrix(unlist(Y.lambda$cluster), n_lambda)
)

lambda_list.Y_t <- list(    # Custom lambda / Y_t (U matrix ONLY)
  hub     = matrix(unlist(Y_t.lambda$hub), n_lambda),
  random  = matrix(unlist(Y_t.lambda$random), n_lambda),
  cluster = matrix(unlist(Y_t.lambda$cluster), n_lambda)
)


summary(lambda_list.Y$hub)
