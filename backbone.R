# normalize_data <- function(data) {
#   n <- length(data)
#   p <- nrow(data[[1]])
#   q <- ncol(data[[1]])
#   
#   z <- matrix(NA, n*p, q)
#   for (j in 1:q) {
#     idx <- 1
#     for (h in 1:n) {
#       z[idx:(idx + p -1), j] <- data[[h]][, j]
#       idx <- idx + p
#     }
#   }
#   
#   # z <- z - apply(z, 2, mean)
#   # norm.z <- apply(z**2, 2, sum)
#   # z <- z / sqrt(norm.z)
#   
#   data <- list()
#   for (h in 1:n) {
#     idx <- h*p
#     data[[h]] <- z[((h-1) * p + 1):(h * p), , drop = F]
#   }
#   
#   return(list(data = data))
# }

set_init <- function(data, lambda) {
  n <- length(data)
  q <- ncol(data[[1]])
  
  # z <- do.call(rbind, data)
  sigma <- sigma_update(data)
  
  xi <- rep(1, q) 
  for (j in 1:q) {
    tmp <- 0
    for (h in 1:n) {
      tmp <- tmp + sum(data[[h]][, j]**2)
    }
    xi[j] <- tmp
  }
  
  B <- outer(sqrt(sigma), 1 / sqrt(sigma))
  XTX <- B_s <- matrix(NA, q, q)
  for (j in 1:q) {
    for (k in 1:q) {
      XTX[j, k] <- xi[k] * B[j, k] + xi[j] * B[k, j]
      B_s[j, k] <- B[j, k]**2 * xi[j] + B[k, j]**2 * xi[k]
    }
  }
  A <- matrix(NA, q, q)
  for (j in 1:q) {
    for (k in 1:q) {
      tmp_sum <- 0
      for (h in 1:n) {
        tmp_sum <- tmp_sum + sum(data[[h]][, j] * data[[h]][, k])
      }
      A[j, k] <- tmp_sum * B[k, j]
    }
  }
  
  rho <- matrix(NA, q, q)  # 초기 상관계수 세팅
  # for (j in 1:q) {
  #   for (k in 1:q) {
  #     tmp1 <- A[j, k] + A[k, j]
  #     tmp2 <- ifelse(tmp1 > 0, tmp1 - lambda, -tmp1 - lambda)
  #     if (tmp2 < 0) {
  #       rho[j, k] <- 0
  #     } else {
  #       tmp3 <- tmp2 / B_s[j, k]
  #       rho[j, k] <- ifelse(tmp1 < 0, -tmp3, tmp3)
  #     }
  #   }
  # }
  for (j in 1:(q-1)) {
    for (k in (j+1):q) {
      tmp1 <- A[j, k] + A[k, j]
      tmp2 <- (tmp1 - lambda) * (tmp1 > lambda) +
        (tmp1 + lambda) * (tmp1 < -lambda) +
        0 * (abs(tmp1) <= lambda)
      rho[j, k] <- rho[k, j] <- tmp2
    }
  }
  diag(rho) <- 0  # 대각선 삭제 
  
  E <- list()
  for (h in 1:n) {
    E[[h]] <- data[[h]] - data[[h]] %*% (rho * B)  # 어짜피 j=k 인 대각선 항의 rho = 0
  }
  
  return(list(sigma = sigma, xi = xi, B = B, XTX = XTX, B_s = B_s, rho = rho, E = E))
}

rho_update <- function(j, k, data, E, B, B_s, rho, lambda) {  # 인덱스 받아와서 원하는 rho&잔차만 업데이트 
  n <- length(data)
  
  # Ajk <- Akj <- 0
  # for (h in 1:n) {
  #   Ajk <- Ajk + sum(E[[h]][, k] * data[[h]][, j])
  #   Akj <- Akj + sum(E[[h]][, j] * data[[h]][, k])
  # }
  # Ajk <- Ajk * B[j, k]
  # Akj <- Akj * B[k, j]
  Ajk <- Akj <- 0
  for (h in 1:n) {
    Y_tilde_jk <- E[[h]][, j] + rho[j, k] * data[[h]][, k] * B[k, j]
    Y_tilde_kj <- E[[h]][, k] + rho[k, j] * data[[h]][, j] * B[j, k]
    
    Ajk <- Ajk + sum(Y_tilde_jk * data[[h]][, k] * B[j, k])
    Akj <- Akj + sum(Y_tilde_kj * data[[h]][, j] * B[k, j])
  } 
  
  rho_next <- (Ajk + Akj) / B_s[j, k]
  # rho_next <- (Ajk + Akj) / B_s[j, k] + rho[j,k]
  tmp <- sign(rho_next) * max(abs(rho_next) - lambda, 0)
  
  delta <- rho[j, k] - tmp
  
  if (delta != 0) {
    for (h in 1:n) {
      E[[h]][, j] <- E[[h]][, j] + delta * data[[h]][, k] * B[k, j]
      E[[h]][, k] <- E[[h]][, k] + delta * data[[h]][, j] * B[j, k]
    }
    rho[j, k] <- rho[k, j] <- tmp
  }
  return(list(rho = rho, E = E, delta = delta))
}

sigma_update <- function(E) {  # 시그마 업데이트
  n <- length(E)
  p <- nrow(E[[1]])
  q <- ncol(E[[1]])
  
  W <- rep(NA, q)
  
  for (j in 1:q) {
    rss_j <- 0
    for (h in 1:n) {
      rss_j <- rss_j + sum(E[[h]][, j]**2)
    }
    W[j] <- (n*p) / rss_j
  }
  return(sigma = W)
}

global_update <- function(dt, rho, E,
                          lambda, maxdiff) {
  n <- length(dt); p <- dim(dt[[1]])[1]; q <- dim(dt[[1]])[2]
  
  sigma <- sigma_update(E)
  B <- outer(sqrt(sigma), 1 / sqrt(sigma))
  xi <- rep(NA, q)
  for (j in 1:q) {
    tmp <- 0
    for (h in 1:n) {
      tmp <- tmp + sum(dt[[h]][, j]**2)
    }
    xi[j] <- tmp
  }
  B_s <- matrix(NA, q, q)
  for (j in 1:q) {
    for (k in 1:q) {
      B_s[j, k] <- B[j, k]**2 * xi[j] + B[k, j]**2 * xi[k]
    }
  }
  
  for (j in 1:(q-1)) {
    for (k in (j+1):q) {
      rho_update_result <- rho_update(j, k, dt, E, B, B_s, rho, lambda)
      rho <- rho_update_result$rho
      E <- rho_update_result$E
      maxdiff <- max(maxdiff, abs(rho_update_result$delta))
    }
  }
  diag(rho) <- 0
  
  return(list(rho = rho, E = E, maxdiff = maxdiff))
}

active_set_update <- function(dt, rho, E,
                              lambda, tol) {
  n <- length(dt); p <- dim(dt[[1]])[1]; q <- dim(dt[[1]])[2]
  
  maxdiff <- 0
  sigma <- sigma_update(E)
  B <- outer(sqrt(sigma), 1 / sqrt(sigma))
  xi <- rep(NA, q)
  for (j in 1:q) {
    tmp <- 0
    for (h in 1:n) {
      tmp <- tmp + sum(dt[[h]][, j]**2)
    }
    xi[j] <- tmp
  }
  B_s <- matrix(NA, q, q)
  for (j in 1:q) {
    for (k in 1:q) {
      B_s[j, k] <- B[j, k]**2 * xi[j] + B[k, j]**2 * xi[k]
    }
  }
  
  active_set <- which(abs(rho) > tol, arr.ind = T)
  if (length(active_set) == 0) {
    return(list(rho = rho, E = E, sigma = sigma, maxdiff = maxdiff))
  }
  
  active_set <- active_set[active_set[, 1] > active_set[, 2], , drop = FALSE]
  
  if (nrow(active_set) == 0) {
    return(list(rho = rho, E = E, sigma = sigma, maxdiff = maxdiff))
  }
  
  for (idx in 1:nrow(active_set)) {
    idx.1 <- active_set[idx, 1]
    idx.2 <- active_set[idx, 2]
    
    rho_update_result <- rho_update(idx.1, idx.2, dt, E, B, B_s, rho, lambda)
    rho <- rho_update_result$rho
    E <- rho_update_result$E
    maxdiff <- max(maxdiff, abs(rho_update_result$delta))
  }
  diag(rho) <- 0
  
  return(list(rho = rho, E = E, sigma = sigma, maxdiff = maxdiff))
}

space <- function(dt, lambda, max_iter = 10000, tol = 1e-6) {
  n <- length(dt)
  data <- normalize_data(dt)$data
  initial_values <- set_init(data, lambda)
  rho <- initial_values$rho; E <- initial_values$E
  
  active_rho_update <- active_set_update(data, initial_values$rho, initial_values$E,
                                         lambda, tol)
  
  for (iter in 1:max_iter) {
    active_rho_update <- active_set_update(data, 
                                           active_rho_update$rho, active_rho_update$E,
                                           lambda, tol)
    
    if (active_rho_update$maxdiff < tol) {
      global_rho_update <- global_update(data, 
                                         active_rho_update$rho, active_rho_update$E,
                                         lambda, active_rho_update$maxdiff)
      
      if (global_rho_update$maxdiff < tol) {
        cat('Converge at iteration:', iter, '\n')
        rho <- global_rho_update$rho
        E <- global_rho_update$E
        break
      } else {
        active_rho_update <- global_rho_update
      } 
    }
    cat("Process:", iter, ' Maxdiff:', active_rho_update$maxdiff, '\n')
    cat("Nonzero beta count:", sum(abs(rho) > tol), "\n")
  }
  
  return(list(rho = rho, E = E, sigma = sigma_update(E)))
}