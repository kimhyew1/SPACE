library(parallel)
library(foreach)
library(doParallel)


GenerateData <- function(n, p, q, type_U, type_V){
  res_U <- create_graph(p, type = type_U)
  res_V <- create_graph(q, type = type_V)
  
  # res_U <- cov_band(p)
  U_chol <- chol(res_U$sigma)
  # res_V <- cov_band(q)
  V_chol <- chol(res_V$sigma)
  
  Y <- list()
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
GenerateData_custom <- function(n, p, q,
                                sigma.U, omega.U,
                                edge_prob_V = 0.01,
                                epsilon = 1e-4) {
  
  omega_V <- generate_sparse_omega(q, edge_prob_V)
  
  sigma_V <- solve(omega_V + diag(epsilon, q))
  
  U_chol <- chol(sigma.U)
  V_chol <- chol(sigma_V)
  
  Y <- vector("list", n)
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

# numCores <- detectCores() - 1
numCores <- 15
myCluster <- makeCluster(numCores)
registerDoParallel(myCluster)

stopCluster(myCluster)




################ 데이터 생성 및 시각화
# Hub 생성된 U matrix 전부 고정
n <- 15; p <- 20; q <- 30

generate_cluster.hub_data <- function(i, n, p, q) {
  dt <- GenerateData(i, n, p, q, 'cluster', 'hub')
  
  return (list(U = dt$U,
               U_inv = dt$U_inv,
               V = dt$V,
               V_inv = dt$V_inv,
               Y = dt$Y,
               Y_t = dt$Y_t))
}

generate_cluster.random_data <- function(i, n, p, q, sigma.U, omega.U, prob = 0.2, eps = 1e-4) {
  if (!is.matrix(sigma.U)) stop("sigma.U is not a matrix")
  if (!is.matrix(omega.U)) stop("omega.U is not a matrix")
  if (nrow(sigma.U) != ncol(sigma.U)) stop("sigma.U is not square")
  
  dt <- GenerateData_custom(i, n, p, q,
                            sigma.U, omega.U,
                            edge_prob_V = prob, epsilon = eps)
  
  return (list(U = dt$U,
               U_inv = dt$U_inv,
               V = dt$V,
               V_inv = dt$V_inv,
               Y = dt$Y,
               Y_t = dt$Y_t))
}

generate_cluster.cluster_data <- function(i, n, p, q, sigma.U, omega.U) {
  dt <- GenerateData_Cluster(i, n, p, q, sigma.U, omega.U, 'cluster')
  
  return (list(U = dt$U,
               U_inv = dt$U_inv,
               V = dt$V,
               V_inv = dt$V_inv,
               Y = dt$Y,
               Y_t = dt$Y_t))
}



foreach(i = 1:50, .packages = c()) %dopar% {
  result <- generate_cluster.hub_data(i, n, p, q)
  
  save(result, file = sprintf("./hub_data/data%02d.RData", i))
}
# for (i in 1:50) {
#   load(sprintf("hub_data/data%02d.RData", i))
#   cat(sprintf("[generated omega] i = %d, U_edges = %d, V_edges = %d \n", i, sum(result$U_inv != 0), sum(result$V_inv != 0)))
# 
#   pdf(sprintf("hub_data/graph_visualise/plot_data%02d.pdf", i), width = 6, height = 6)
#   plot_graph(result$V_inv)
#   plot_graph(result$U_inv)
#   dev.off()
# }


foreach(i = 1:50, .packages = c()) %dopar% {
  load(sprintf("./hub_data/data%02d.RData", i))
  U.matrix <- result$U
  U_inv.matrix <- result$U_inv
  
  result <- generate_cluster.random_data(i, n, p, q, U.matrix, U_inv.matrix)
  
  save(result, file = sprintf("./random_data/data%02d.RData", i))
}
# for (i in 1:50) {
#   load(sprintf("random_data/data%02d.RData", i))
#   cat(sprintf("[generated omega] i = %d, U_edges = %d, V_edges = %d \n", i, sum(result$U_inv != 0), sum(result$V_inv != 0)))
#   
#   pdf(sprintf("random_data/graph_visualise/plot_data%02d.pdf", i), width = 6, height = 6)
#   plot_graph(result$V_inv)
#   plot_graph(result$U_inv)
#   dev.off()
# }


foreach(i = 1:50, .packages = c()) %dopar% {
  load(sprintf("./hub_data/data%02d.RData", i))
  U.matrix <- result$U
  U_inv.matrix <- result$U_inv
  
  result <- generate_cluster.cluster_data(i, n, p, q, U.matrix, U_inv.matrix)
  
  save(result, file = sprintf("./cluster_data/data%02d.RData", i))
}
# Graph Visualise
for (i in 1:50) {
  load(sprintf("./cluster_data/data%02d.RData", i))
  cat(sprintf("[generated omega] i = %d, U_edges = %d, V_edges = %d \n", i, sum(result$U_inv != 0), sum(result$V_inv != 0)))
  
  pdf(sprintf("cluster_data/graph_visualise/plot_data%02d.pdf", i), width = 6, height = 6)
  plot_graph(result$V_inv)
  plot_graph(result$U_inv)
  dev.off()
}



