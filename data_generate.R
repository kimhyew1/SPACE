create_graph <- function(i, p, type) {
  set.seed(i)
  res_huge <- huge::huge.generator(n = 2, d = p, graph = type,
                                   verbose = F)
  return(
    list(sigma = res_huge$sigma,
         omega = res_huge$omega)
  )
}

GenerateData <- function(i, n, p, q, type_U, type_V){
  res_U <- create_graph(i, p, type = type_U)
  res_V <- create_graph(i, q, type = type_V)
  
  # res_U <- cov_band(p)
  U_chol <- chol(res_U$sigma)
  # res_V <- cov_band(q)
  V_chol <- chol(res_V$sigma)
  
  Y <- list()
  set.seed(i)
  for(h in 1:n){
    Y[[h]] <- t(U_chol) %*% matrix(rnorm(p * q), p, q) %*% V_chol
  }
  Y_t <- lapply(Y, t)
  
  return(list(Y = Y,
              Y_t = Y_t,
              U = res_U$sigma,
              U_inv = res_U$omega,
              V = res_V$sigma,
              V_inv = res_V$omega))
}
set.seed(1)
dt <- GenerateData(1, 50, 90, 100, 'cluster', 'hub')

length(dt$Y)  # [50] 90 100
length(dt$Y_t)  # [50] 100 90
dim(dt$U)  # 90 90
dim(dt$V)   # 100 100
sum(dt$V_inv != 0)   # 2000
sum(dt$U_inv != 0)   # 1601

Y_data <- dt$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.Y.RData")
Y_data <- dt$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.Y_t.RData")
Y_data <- dt$U
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.U.RData")
Y_data <- dt$V
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.V.RData")

# dt.1 <- GenerateData(50, 90, 100,
#                      type_U = list(type = 'hub', v = 5),
#                      type_V = list(type = 'random', v = 0.0001))
# 
# length(dt.1$Y)  # 50
# length(dt.1$Y_t)  # 50 보통 column-wise 추정에 사용
# dim(dt.1$U)  # 90 90
# dim(dt.1$V)   # 100 100
# sum(dt.1$V_inv != 0)

GenerateData_Cluster <- function(i, n, p, q, sigma.U, omega.U, type_V){
  # res_U <- create_graph(p, type = type_U)
  res_U <- sigma.U
  res_V <- create_graph(i, q, type = type_V)
  
  # res_U <- cov_band(p)
  U_chol <- chol(res_U)
  # res_V <- cov_band(q)
  V_chol <- chol(res_V$sigma)
  
  Y <- list()
  for(h in 1:n){
    Y[[h]] <- t(U_chol) %*% matrix(rnorm(p * q), p, q) %*% V_chol
  }
  Y_t <- lapply(Y, t)
  
  return(list(Y = Y,
              Y_t = Y_t,
              U = sigma.U,
              U_inv = omega.U,
              V = res_V$sigma,
              V_inv = res_V$omega))
}

# dt.2 <- GenerateData(50, 90, 100, 'cluster', 'scale-free')
# length(dt.2$Y)  # 50
# length(dt.2$Y_t)  # 50 보통 column-wise 추정에 사용
# dim(dt.2$U)  # 90 90
# dim(dt.2$V)   # 100 100
# sum(dt.2$V_inv != 0)

# dt.3 <- GenerateData(50, 90, 100, 'cluster', 'cluster')
# length(dt.3$Y)  # 50
# length(dt.3$Y_t)  # 50 보통 column-wise 추정에 사용
# dim(dt.3$U)  # 90 90
# dim(dt.3$V)   # 100 100
# sum(dt.3$V_inv != 0)

# Clsuter Cluster
set.seed(1)
dt.3 <- GenerateData_Cluster(50, 90, 100, dt$U, dt$U_inv, 'cluster')
length(dt.3$Y)  # 50
length(dt.3$Y_t)  # 50 보통 column-wise 추정에 사용
dim(dt.3$U)  # 90 90
dim(dt.3$V)   # 100 100
sum(dt.3$V_inv != 0)

Y_data <- dt.3$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.3.Y.RData")
Y_data <- dt.3$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.3.Y_t.RData")
Y_data <- dt.3$U
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.3.U.RData")
Y_data <- dt.3$V
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.3.V.RData")




########### n 증가
random_data <- function(n, p, q, sigma.U, sigma.V) {
  res_U <- sigma.U
  res_V <- sigma.V

  U_chol <- chol(res_U)
  V_chol <- chol(res_V)
  
  Y <- list()
  for(h in 1:n){
    Y[[h]] <- t(U_chol) %*% matrix(rnorm(p * q), p, q) %*% V_chol
  }
  Y_t <- lapply(Y, t)
  
  return(list(Y = Y,
              Y_t = Y_t,
              V = sigma.V, V_inv = V_chol,
              U = sigma.U, U_inv = U_chol))
}

set.seed(1)
dt.200 <- random_data(200, 90, 100, dt$U, dt$V)
length(dt.200$Y)  # 200
length(dt.200$Y_t)  # 200 
dim(dt.200$U)  # 90 90
dim(dt.200$V)   # 100 100
# sum(dt$U != dt.200$U)  # 0
Y_data <- dt.200$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.200.Y.RData")
Y_data <- dt.200$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.200.Y_t.RData")

set.seed(1)
dt.500 <- random_data(500, 90, 100, dt$U, dt$V)
length(dt.500$Y)  # 200
length(dt.500$Y_t)  # 200 
dim(dt.500$U)  # 90 90
dim(dt.500$V)   # 100 100
Y_data <- dt.500$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.500.Y.RData")
Y_data <- dt.500$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.500.Y_t.RData")

set.seed(1)
dt.700 <- random_data(700, 90, 100, dt$U, dt$V)
length(dt.700$Y)  # 200
length(dt.700$Y_t)  # 200 
dim(dt.700$U)  # 90 90
dim(dt.700$V)   # 100 100
Y_data <- dt.700$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.700.Y.RData")
Y_data <- dt.700$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.700.Y_t.RData")

set.seed(1)
dt.1000 <- random_data(1000, 90, 100, dt$U, dt$V)
length(dt.1000$Y)  # 200
length(dt.1000$Y_t)  # 200 
dim(dt.1000$U)  # 90 90
dim(dt.1000$V)   # 100 100
Y_data <- dt.1000$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1000.Y.RData")
Y_data <- dt.1000$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1000.Y_t.RData")





##### Random Graph 
generate_sparse_omega <- function(i, p, edge_prob = 0.05, min_val = 0.2, max_val = 0.5) {
  adj <- matrix(0, p, p)
  
  set.seed(i)
  adj[upper.tri(adj)] <- rbinom(p * (p - 1) / 2, 1, edge_prob)
  adj <- adj + t(adj)  # 대칭
  
  edge_count <- sum(adj == 1, na.rm = TRUE)  # na.rm 추가
  
  if (edge_count == 0) {
    warning(sprintf("No edges created for i = %d", i))
    edge_count <- 1  # 최소한 1개라도 만들어야 아래에서 에러 안 남
    adj[sample(which(upper.tri(adj)), 1)] <- 1
    adj <- adj + t(adj)
  }
  
  weights <- runif(edge_count, min = min_val, max = max_val)
  signs <- sample(c(-1, 1), edge_count, replace = TRUE)
  adj[adj == 1] <- weights * signs
  
  diag(adj) <- rowSums(abs(adj)) + 0.1
  
  return(adj)
}





################## 
# 이미 생성되어 있는 U matrix 활용
GenerateData_custom <- function(i, n, p, q,
                                sigma.U, omega.U,
                                edge_prob_V = 0.01,
                                epsilon = 1e-4) {
  
  omega_V <- generate_sparse_omega(i, q, edge_prob_V)
  
  sigma_V <- solve(omega_V + diag(epsilon, q))
  
  U_chol <- chol(sigma.U)
  V_chol <- chol(sigma_V)
  
  Y <- vector("list", n)
  set.seed(i)
  for (h in 1:n) {
    Y[[h]] <- t(U_chol) %*% matrix(rnorm(p * q), p, q) %*% V_chol
  }
  Y_t <- lapply(Y, t)
  
  return(list(Y = Y,
              Y_t = Y_t,
              U = sigma.U,
              U_inv = omega.U,
              V = sigma_V,
              V_inv = omega_V))
}

set.seed(1)
dt.1 <- GenerateData_custom(     # Cluster, Random
  n = 50, p = 90, q = 100,
  sigma.V = dt$V, omega.V = dt$V_inv,
  edge_prob_U = 0.2
)


length(dt.1$Y)  # 50
str(dt.1$Y)
length(dt.1$Y_t)  # 50 보통 column-wise 추정에 사용
dim(dt.1$U)  # 90 90
dim(dt.1$V)   # 100 100
sum(dt.1$V_inv != 0)   # 2000
sum(dt.1$U_inv != 0)   # 2172

dt.1 <- list()
load('test_data.1.Y.RData')
dt.1$Y <- Y_data
load('test_data.1.Y_t.RData')
dt.1$Y_t <- Y_data
load('test_data.1.U.RData')
dt.1$U <- Y_data
load('test_data.1.V.RData')
dt.1$V <- Y_data

sum(dt$U != dt.1$U)

class(dt$Y[[1]])
class(dt.1$Y[[1]])
sapply(dt.1$Y, class)
which(sapply(dt.1$Y, function(x) !is.matrix(x)))


Y_data <- dt.1$Y
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1.Y.RData")
Y_data <- dt.1$Y_t
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1.Y_t.RData")
Y_data <- dt.1$U
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1.U.RData")
Y_data <- dt.1$V
save(Y_data, file = "C:/Users/User/Desktop/matrix SPACE/250708_simulation/test_data.1.V.RData")




