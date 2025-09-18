# L1 Norm
l1_calculate <- function(prec, prec_hat) {
  delta <- abs(prec - prec_hat)
  return(l1 = sum(delta))
}

# Frobenius Norm
frob_calculate <- function(prec, prec_hat) {
  delta <- sum((prec - prec_hat)**2)
  return(frob = sqrt(delta))
}

# Spectral Norm
spec_calculate <- function(prec, prec_hat) {
 delta <- prec - prec_hat
 ans <- max(eigen(t(delta) %*% delta)$values)
 return(spec = sqrt(ans))
}


# V_l1 <- V_frob <- V_spec <- rep(NA, length(x.lambda))
# for (i in 1:length(x.lambda)) {
#   prec_hat <- result_beta[[i]]
#   V_l1[i] <- l1_calculate(dt$V_inv, prec_hat)
#   V_frob[i] <- frob_calculate(dt$V_inv, prec_hat)
#   V_spec[i] <- spec_calculate(dt$V_inv, prec_hat)
# }

# plot(x.lambda, V_l1, type = 'l',
#      xlab = 'lambda', ylab = 'L1 norm', main = 'L1 norm by l1-penalty')
# point(x.lambda[which.min(BIC)], V_l1[which.min(BIC)], col = 'blue', pch = 19)
# 
# plot(x.lambda, V_frob, type = 'l',
#      xlab = 'lambda', ylab = 'Frobenius norm', main = 'Frobenius norm by l1-penalty')
# point(x.lambda[which.min(BIC)], V_frob[which.min(BIC)], col = 'blue', pch = 19)
# 
# plot(x.lambda, V_spec, type = 'l',
#      xlab = 'lambda', ylab = 'Spectral norm', main = 'Spectral norm by l1-penalty')
# point(x.lambda[which.min(BIC)], V_spec[which.min(BIC)], col = 'blue', pch = 19)