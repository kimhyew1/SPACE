# lambda.flag <- "custom"
# lambda.flag <- "auto"

lambda.function <- function(dt, flag, K) {
  if (flag == "custom") {
    source("./lambda.R")   # Load lambda making function
    
    lambda <- lambda.bound(dt, K = K)$lambda
  } else if (flag == "auto") {
    lambda <- exp(seq(log(1e-6), log(10), length.out = K))
  }
  
  return(list(lambda = lambda))
}

