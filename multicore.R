library(parallel)
library(foreach)
library(doParallel)

# numCores <- detectCores() - 1
numCores <- 15
myCluster <- makeCluster(numCores)
registerDoParallel(myCluster)

stopCluster(myCluster)




################ 데이터 생성
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






########################
# SPACE
n_cores <- detectCores() - 1
myCluster <- makeCluster(n_cores)
registerDoParallel(n_cores)

data_types <- c("hub", "random", "cluster")
lambda_list.Y <- list(    # Custom lambda / Y (V matrix ONLY)
  hub = lambda.hub,
  random = lambda.random,
  cluster = lambda.cluster
)
lambda_list.Y_t <- list(    # Custom lambda / Y_t (U matrix ONLY)
  hub = lambda.hub,
  random = lambda.random,
  cluster = lambda.cluster
)
lambda_list <- list(    # User-specified lambda
  hub = exp(seq(log(1e-6), log(10), length.out = 20)),
  random = exp(seq(log(1e-6), log(10), length.out = 20)),
  cluster = exp(seq(log(1e-6), log(10), length.out = 20))
)
vals <- exp(seq(log(1e-6), log(10), length.out = 20))
mat <- matrix(rep(vals, each = 50), nrow = 20, ncol = 50, byrow = T)

lambda_list.Y <- list(    # Custom lambda / Y (V matrix ONLY)
  hub = mat,
  random = mat,
  cluster = mat
)
lambda_list.Y_t <- list(    # Custom lambda / Y (V matrix ONLY)
  hub = mat,
  random = mat,
  cluster = mat
)
task_list <- expand.grid(
  type = data_types,
  index = 1:50,
  stringsAsFactors = FALSE
)

#########################################
if (!dir.exists("log")) dir.create("log")

results <- foreach(i = 1:nrow(task_list), .packages = c()) %dopar% {
  type <- task_list$type[i]
  idx  <- task_list$index[i]
  
  log_file <- sprintf("log/%s_%02d.log", type, idx)
  sink(log_file)
  on.exit({
    sink()  # 로그 닫기
  }, add=TRUE)
  
  cat(sprintf("Start: %s %02d\n", type, idx))
  
  # ---- 데이터 로드 (공통) ----
  file_path <- sprintf("./%s_data/data%02d.RData", type, idx)
  load(file_path)  # result 객체 로드됨
  
  # ---- U와 V를 모두 처리 ----
  sim_result_U <- list()
  sim_result_V <- list()
  
  for (label in c("U", "V")) {
    # 데이터 선택
    test_data <- if (label == "U") result$Y_t else result$Y
    lambda_seq <- if (label == "U") lambda_list.Y_t[[type]][, idx] else lambda_list.Y[[type]][, idx]
    
    sim_result <- list()
    for (l in lambda_seq) {
      sim <- space(test_data, l)
      sim_result[[as.character(l)]] <- list(
        rho    = sim$rho,
        prec   = - sim$rho * sqrt(outer(sim$weight, sim$weight)),
        E      = sim$resid,
        weight = sim$weight
      )
    }
    
    if (label == "U") {
      sim_result_U <- sim_result
    } else {
      sim_result_V <- sim_result
    }
  }
  
  # ---- 결과 저장 ----
  save_path_U <- sprintf("./2_Auto_lambda/%s_result/data%02d_result.U.RData", type, idx)
  save(sim_result_U, file = save_path_U)
  
  save_path_V <- sprintf("./2_Auto_lambda/%s_result/data%02d_result.V.RData", type, idx)
  save(sim_result_V, file = save_path_V)
  
  cat(sprintf("End: %s %02d\n", type, idx))
  
  list(type=type, idx=idx, ok = TRUE)  
}


stopCluster(myCluster)                                                                                                                              # stopCluster(n_cores)


