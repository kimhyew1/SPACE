soft_thresholding <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

set_init <- function(dt, lambda) {
  n <- length(dt); p <- nrow(dt[[1]]); q <- ncol(dt[[1]])
  
  dt.rbind <- do.call(rbind, dt)
  weight_init <- 1 / apply(dt.rbind, 2, var)
  
  
  numer <- matrix(0, q, q)
  for (h in 1:n) {
    numer <- numer + t(dt[[h]]) %*% dt[[h]]
  }
  numer <- sqrt(outer(weight_init, 1 / weight_init)) * numer
  numer <- numer + t(numer)
  
  denom <- matrix(0, 1, q)
  for (h in 1:n) {
    denom <- denom + colSums(dt[[h]]**2)
  }  
  denom <- matrix(rep(denom, each = q), q, q, byrow = T)
  denom <- outer(weight_init, 1 / weight_init) * denom
  denom <- denom + t(denom)
  
  rho_init <- numer / denom
  rho_init <- soft_thresholding(rho_init, lambda)
  
  
  # --- residual 저장: p × q 행렬 ---
  rss_init <- vector("list", n)
  for (h in 1:n) {
    y <- dt[[h]]
    rss_mat <- matrix(0, p, q)   # 각 j에 대해 residual 열 저장
    for (j in 1:q) {
      idx <- setdiff(1:q, j)
      weights <- sqrt(weight_init[idx] / weight_init[j]) * rho_init[j, idx]
      pred <- y[, idx, drop = F] %*% weights
      rss_mat[, j] <- y[, j] - pred
    }
    rss_init[[h]] <- rss_mat
  }
  
  return(list(rho = rho_init, weight = weight_init, resid = rss_init))
}

rho_update <- function(j, k, dt, resid, weight, rho, lambda) {
  n <- length(dt)
  numer <- 0
  denom <- 0
  
  for (h in 1:n) {
    Y <- dt[[h]]
    
    idx_j <- setdiff(1:ncol(Y), c(j, k))
    idx_k <- idx_j
    
    res_jk <- Y[, j] - (Y[, idx_j, drop = F] %*% 
                          (sqrt(weight[idx_j] / weight[j]) * as.vector(rho[j, idx_j])))
    
    res_kj <- Y[, k] - (Y[, idx_k, drop = F] %*% 
                          (sqrt(weight[idx_k] / weight[k]) * as.vector(rho[k, idx_j])))
    
    numer <- numer +
      sum(res_jk * sqrt(weight[k]/weight[j]) * Y[, k]) +
      sum(res_kj * sqrt(weight[j]/weight[k]) * Y[, j])
    
    denom <- denom +
      sum((weight[j]/weight[k]) * (res_jk**2)) +
      sum((weight[k]/weight[j]) * (res_kj**2))
  }
  
  rho_new <- soft_thresholding(numer / denom, lambda)
  delta <- rho[j, k] - rho_new
  rho[k, j] <- rho[j, k] <- rho_new
  
  return(list(rho = rho, delta = delta))
}

update_sigma_resid <- function(dt, rho_new, weight_old) {
  n <- length(dt); p <- nrow(dt[[1]]); q <- ncol(dt[[1]])
  
  weight_new <- weight_old
  rss_new <- vector("list", n)
  for (h in 1:n) {
    rss_new[[h]] <- matrix(0, p, q) 
  }
  
  for (j in 1:q) {
    rss_sum <- 0
    
    for (h in 1:n) {
      Y <- dt[[h]]
      idx <- setdiff(1:q, j)
      
      pred <- Y[, idx, drop = F] %*% 
        (sqrt(weight_old[idx] / weight_old[j]) * rho_new[j, idx])
      
      resid <- Y[, j] - pred
      rss_sum <- rss_sum + sum(resid**2)
      
      rss_new[[h]][, j] <- resid
    }
    
    weight_inv <- rss_sum / (n * p)
    weight_new[j] <- 1 / weight_inv
  }
  
  return(list(weight = weight_new, resid = rss_new))
}


global_update <- function(dt, rho, resid, weight, lambda, maxdiff) {
  n <- length(dt); p <- nrow(dt[[1]]); q <- ncol(dt[[1]])
  
  tmp <- update_sigma_resid(dt, rho, weight)
  weight_new <- tmp$weight
  resid_new <- tmp$resid
  
  for (j in 1:(q-1)) {
    for (k in (j+1):q) {
      upd <- rho_update(j, k, dt, resid_new, weight_new, rho, lambda)
      rho <- upd$rho
      maxdiff <- max(maxdiff, abs(upd$delta))
    }
  }
  return(list(rho = rho, resid = resid_new, maxdiff = maxdiff, weight = weight_new))
}

active_set_update <- function(dt, weight, rho, lambda, tol) {
  n <- length(dt); p <- nrow(dt[[1]]); q <- ncol(dt[[1]])
  
  maxdiff <- 0
  tmp <- update_sigma_resid(dt, rho, weight)
  weight_new <- tmp$weight
  resid_new <- tmp$resid
  
  active_set <- which(abs(rho) > tol, arr.ind = T)
  if (length(active_set) == 0) {
    return(list(rho = rho, resid = resid_new, weight = weight_new, maxdiff = maxdiff))
  }
  
  active_set <- active_set[active_set[, 1] > active_set[, 2], , drop = FALSE]
  
  if (length(active_set) == 0) {
    return(list(rho = rho, resid = resid_new, weight = weight_new, maxdiff = maxdiff))
  }
  
  for (idx in 1:nrow(active_set)) {
    idx.1 <- active_set[idx, 1]
    idx.2 <- active_set[idx, 2]
    
    upd <- rho_update(idx.1, idx.2, dt, resid_new, weight_new, rho, lambda)
    rho <- upd$rho
    maxdiff <- max(maxdiff, abs(upd$delta))
  }
  
  return(list(rho = rho, resid = resid_new, weight = weight_new, maxdiff = maxdiff))
}

space <- function(dt, lambda, max_iter = 10000, tol = 1e-6) {
  n <- length(dt)
  init <- set_init(dt, lambda)
  rho    <- init$rho
  weight <- init$weight
  resid  <- init$resid
  
  for (iter in 1:max_iter) {
    step <- active_set_update(dt, weight, rho, lambda, tol)
    rho    <- step$rho
    weight <- step$weight
    resid  <- step$resid
    
    cat("Process:", iter, " Maxdiff(active):", step$maxdiff, "\n")
    cat("Nonzero beta count:", sum(abs(rho) > tol), "\n")
    
    if (step$maxdiff < tol) {
      gstep <- global_update(dt, rho, resid, weight, lambda, maxdiff = 0)
      rho    <- gstep$rho
      weight <- gstep$weight
      resid  <- gstep$resid
      
      cat("  Global maxdiff:", gstep$maxdiff, "\n")
      if (gstep$maxdiff < tol) {
        cat("Converged at iteration:", iter, "\n")
        break
      } 
    }
  }
  final <- update_sigma_resid(dt, rho, weight)
  
  return(list(rho = rho, resid = final$resid, weight = final$weight))
}





#######################################################
compute_BIC <- function(rho, E, tol = 1e-8) {
  n <- length(E)         # 시점 개수
  p <- dim(E[[1]])[1]      # gene
  q <- dim(E[[1]])[2]     # time (혹은 환자)
  N <- n*p
  
  RSS <- rep(0, q)
  for (h in 1:n) {
    RSS <- RSS + colSums(E[[h]]**2)  # 각 열(j)에 대한 오차 제곱합 누적
  }
  tmp <- apply(abs(rho) > tol, 2, sum)
  # E <- do.call(rbind, E)
  # RSS_centered <- scale(E, center = T, scale = F)
  # RSS <- (sum(RSS_centered ** 2))
  # BIC_vec <- N*log(RSS)
  BIC_vec <- N*log(RSS) + log(N) * tmp
  return(sum(BIC_vec))
}

