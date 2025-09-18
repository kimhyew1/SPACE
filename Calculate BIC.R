input_dir <- "1_Specified_lambda"
input_dir <- "2_Auto_lambda"

# 람다 별 BIC 그래프 그리기 
lambda_seq <- lambda_list$hub
exclude_idx <- c(33, 38, 44)

for (type in unique(data_types)) {
  for (idx in setdiff(1:50, exclude_idx)) {
    
    test_bic.U <- test_bic.V <- rep(NA, length(lambda_seq))
    
    for (label in c("U", "V")) {
      # lambda 불러오기
      lambda_seq <- if (label == "U") lambda_list.Y_t[[type]][, idx] else lambda_list.Y[[type]][, idx]
      
      # 결과 벡터 초기화
      bic_name <- sprintf("test_bic.%s", label)
      assign(bic_name, numeric(length(lambda_seq)))
      
      # 데이터 로드
      load(sprintf("./%s/%s_result/data%02d_result.%s.RData", input_dir, type, idx, label))
      
      # sim_result 객체 이름 만들기 (sim_result_U / sim_result_V)
      sim_name <- sprintf("sim_result_%s", label)
      sim_result <- get(sim_name)
      
      # BIC 계산
      for (i in seq_along(lambda_seq)) {
        l <- lambda_seq[i]
        tmp <- compute_BIC(
          sim_result[[as.character(l)]]$rho,
          sim_result[[as.character(l)]]$E
        )
        # test_bic.U[i] or test_bic.V[i] 에 저장
        bic_vec <- get(bic_name)
        bic_vec[i] <- tmp
        assign(bic_name, bic_vec)
      }
      
      # 저장
      save_path <- sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label)
      save(list = bic_name, file = save_path)
    }
    cat(sprintf("Finish calculating BIC: %s %s\n", type, idx))
  }
}





##############################################################
# lambda 별 BIC 계산 및 최소점 찾기

v <- list()
bic_min_idx.lambda2 <- list()

for (type in unique(data_types)) {
  bic_min_idx[[type]] <- data.frame(data = integer(),
                                    U_min_idx = integer(),
                                    V_min_idx = integer(),
                                    U_min_lambda = numeric(),
                                    V_min_lambda = numeric())
  
  for (idx in setdiff(1:50, exclude_idx)) {
    min_idx_list <- list()
    min_lambda_list <- list()
    
    for (label in c("U", "V")) {
      # lambda 불러오기
      lambda_seq <- if (label == "U") lambda_list.Y_t[[type]][, idx] else lambda_list.Y[[type]][, idx]
      
      load(sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label))  # loads test_bic.U
      bic_vec <- get(sprintf("test_bic.%s", label))
      
      # 최소값 찾기
      min_idx <- which.min(bic_vec)
      min_lambda <- lambda_seq[min_idx]
      min_idx_list[[label]] <- min_idx
      min_lambda_list[[label]] <- min_lambda
      
      png(sprintf("./%s/%s_bic/data%02d.%s.png", input_dir, type, idx, label),
          width = 1200, height = 800, res = 150)
      plot(log(lambda_seq), bic_vec, type = "l", lwd = 2,
           xlab = "lambda", ylab = sprintf("BIC (%s)", label),
           main = sprintf("%s | data%02d (%s)", type, idx, label))
      points(log(lambda_seq[min_idx]), min(bic_vec), col = "blue", pch = 19)
      dev.off()
    }
    bic_min_idx.lambda2[[type]] <- rbind(
      bic_min_idx.lambda2[[type]],
      data.frame(data = idx,
                 U_min_idx = min_idx_list$U,
                 V_min_idx = min_idx_list$V,
                 U_min_lambda = min_lambda_list$U,
                 V_min_lambda = min_lambda_list$V)
    )
      
    # 저장
    save_path <- sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label)
    save(list = bic_name, file = save_path)
  }
  cat(sprintf("Finish calculating BIC: %s %s\n", type, idx))
}





##############################
# Draw ROC curve

library(pROC)

plot_roc <- function(sim_result, rho_true, lambda_seq, main_title, out_file, test_bic) {
  upper_idx <- upper.tri(rho_true)
  edge_true <- as.numeric(abs(rho_true[upper_idx]) > 1e-12)  # 0/1 벡터
  
  auc_val <- NA
  
  # --- BIC 최소 지점 ---
  min_idx <- which.min(test_bic)
  best_lambda <- lambda_seq[min_idx]
  rho_best <- sim_result[[as.character(best_lambda)]]$rho
  
  # predictor: 절대값을 score로 사용 (0일수록 연결 없음, 클수록 edge 강함)
  predictor <- abs(rho_best[upper_idx])
  
  roc_obj <- roc(response = edge_true, predictor = predictor, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  # --- 그리기 ---
  png(out_file, width=1200, height=800, res=150)
  plot(roc_obj, main = sprintf("%s (AUC=%.3f)", main_title, auc_val),
       col = "black", lwd = 2)
  dev.off()
  
  return(auc_val)
}



exclude_idx <- c(33, 38, 44)
auc_list.lambda1 <- list(
  hub     = list(U = numeric(), V = numeric()),
  cluster = list(U = numeric(), V = numeric()),
  random  = list(U = numeric(), V = numeric())
)
auc_list.lambda2 <- list(
  hub     = list(U = numeric(), V = numeric()),
  cluster = list(U = numeric(), V = numeric()),
  random  = list(U = numeric(), V = numeric())
)

load(sprintf("./%s/%s_result/data%02d_result.%s.RData", input_dir, "hub", 1, "U"))
load(sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label)) 

for (type in unique(task_list$type)) {
  for (idx in setdiff(1:50, exclude_idx)) {
    
    # true 값값
    load(sprintf("./%s_data/data%02d.RData", type, idx))
    for (label in c("U", "V")) {
      # lambda_seq <- if (label == "U") lambda_list.Y_t[[type]][, idx] else lambda_list.Y[[type]][, idx]
      
      load(sprintf("./%s/%s_result/data%02d_result.%s.RData", input_dir, type, idx, label))
      load(sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label)) 
      # sim_result 객체 선택 (sim_result_U / sim_result_V)
      sim_result <- get(sprintf("sim_result_%s", label))
      lambda_seq <- names(sim_result)
      
      # test_bic 객체 선택 (test_bic.U / test_bic.V)
      test_bic <- get(sprintf("test_bic.%s", label))
      
      # true inverse matrix 선택
      rho_true <- if (label == "U") result$U_inv else result$V_inv
      
      out_file <- sprintf("./%s/plot.ROC_curve/%s_data%02d.%s.png", input_dir, type, idx, label)
      auc <- plot_roc(sim_result,
                      rho_true,
                      lambda_seq,
                      main_title = sprintf("%s | data%02d (%s)", type, idx, label),
                      out_file   = out_file,
                      test_bic   = test_bic)
      
      # === AUC 저장 ===
      auc_list.lambda2[[type]][[label]] <- c(auc_list.lambda2[[type]][[label]], auc)
      
      cat(sprintf("Saved ROC plots: %s data%02d | AUC(%s)=%.3f\n",
                  type, idx, label, auc))
    }
  }
}
lambda_list.Y$random[, 1]
load(sprintf("./%s/%s_result/data%02d_result.%s.RData", input_dir, "random", 1, "V"))
names(sim_result_V)

for (type in unique(task_list$type)) {
  for (idx in setdiff(1:50, exclude_idx)) {
    
    # true 값값
    load(sprintf("./%s_data/data%02d.RData", type, idx))
    for (label in c("U", "V")) {
      lambda_seq <- if (label == "U") lambda_list.Y_t[[type]][, idx] else lambda_list.Y[[type]][, idx]
      
      load(sprintf("./%s/%s_result/data%02d_result.%s.RData", input_dir, type, idx, label))
      load(sprintf("./%s/%s_bic/data%02d.%s.RData", input_dir, type, idx, label)) 
      # sim_result 객체 선택 (sim_result_U / sim_result_V)
      sim_result <- get(sprintf("sim_result_%s", label))
      
      # test_bic 객체 선택 (test_bic.U / test_bic.V)
      test_bic <- get(sprintf("test_bic.%s", label))
      
      # true inverse matrix 선택
      rho_true <- if (label == "U") result$U_inv else result$V_inv
      
      out_file <- sprintf("./%s/plot.ROC_curve/%s_data%02d.%s.png", input_dir, type, idx, label)
      auc <- plot_roc(sim_result,
                      rho_true,
                      lambda_seq,
                      main_title = sprintf("%s | data%02d (%s)", type, idx, label),
                      out_file   = out_file,
                      test_bic   = test_bic)
      
      cat(sprintf("Saved ROC plots: %s data%02d | AUC(%s)=%.3f\n",
                  type, idx, label, auc))
    }
  }
}






which.min(auc_list.lambda1$cluster$V)
which.max(auc_list.lambda1$cluster$V)
