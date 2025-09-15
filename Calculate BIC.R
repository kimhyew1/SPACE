lambda_seq <- lambda_list$hub
exclude_idx <- c(33, 38, 44)

for (type in unique(data_types)) {
  for (idx in setdiff(1:50, exclude_idx)) {
    
    test_bic.U <- test_bic.V <- rep(NA, length(lambda_seq))
    
    # ---- U 처리 ----
    load(sprintf("./%s_result/data%02d_result.U.RData", type, idx))
    for (i in seq_along(lambda_seq)) {
      l <- lambda_seq[i]
      test_bic.U[i] <- compute_BIC(
        sim_result_U[[as.character(l)]]$rho,
        sim_result_U[[as.character(l)]]$E
      )
    }
    
    # ---- V 처리 ----
    load(sprintf("./%s_result/data%02d_result.V.RData", type, idx))
    for (i in seq_along(lambda_seq)) {
      l <- lambda_seq[i]
      test_bic.V[i] <- compute_BIC(
        sim_result_V[[as.character(l)]]$rho,
        sim_result_V[[as.character(l)]]$E
      )
    }
    
    # ---- 결과 저장 ----
    save_path_U <- sprintf("./%s_bic/data%02d.U.RData", type, idx)
    save(test_bic.U, file = save_path_U)
    
    save_path_V <- sprintf("./%s_bic/data%02d.V.RData", type, idx)
    save(test_bic.V, file = save_path_V)
    
    # ---- 로그 ----
    cat(sprintf("Finish calculating: idx=%02d, type=%s\n", idx, type))
  }
}




bic_min_idx <- list()

for (type in unique(task_list$type)) {
  out_dir <- sprintf("./%s_bic", type)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  bic_min_idx[[type]] <- data.frame(data = integer(),
                                    U_min_idx = integer(),
                                    V_min_idx = integer(),
                                    U_min_lambda = numeric(),
                                    V_min_lambda = numeric())
  
  for (idx in setdiff(1:50, exclude_idx)) {
    
    ## --- U ---
    load(sprintf("./%s_bic/data%02d.U.RData", type, idx))  # loads test_bic.U
    min_idx_U <- which.min(test_bic.U)
    min_lambda_U <- lambda_seq[min_idx_U]
    
    png(sprintf("./%s_bic/data%02d.U.png", type, idx),
        width = 1200, height = 800, res = 150)
    plot(log(lambda_seq), test_bic.U, type = "l", lwd = 2,
         xlab = "lambda", ylab = "BIC (U)",
         main = sprintf("%s | data%02d (U)", type, idx))
    points(log(lambda_seq[which.min(test_bic.U)]), min(test_bic.U), col = "blue", pch = 19)
    dev.off()
    
    ## --- V ---
    load(sprintf("./%s_bic/data%02d.V.RData", type, idx))  # loads test_bic.V
    min_idx_V <- which.min(test_bic.V)
    min_lambda_V <- lambda_seq[min_idx_V]
    
    png(sprintf("./%s_bic/data%02d.V.png", type, idx),
        width = 1200, height = 800, res = 150)
    plot(log(lambda_seq), test_bic.V, type = "l", lwd = 2,
         xlab = "lambda", ylab = "BIC (V)",
         main = sprintf("%s | data%02d (V)", type, idx))
    points(log(lambda_seq[which.min(test_bic.V)]), min(test_bic.V), col = "blue", pch = 19)
    dev.off()
    
    # --- 결과 저장 ---
    bic_min_idx[[type]] <- rbind(
      bic_min_idx[[type]],
      data.frame(data = idx,
                 U_min_idx = min_idx_U,
                 V_min_idx = min_idx_V,
                 U_min_lambda = min_lambda_U,
                 V_min_lambda = min_lambda_V)
    )
    
    cat(sprintf("Saved plots: %s data%02d\n", type, idx))
  }
}








##############################
# Draw ROC curve
compute_tpr_fpr <- function(rho_hat, rho_true) {
  upper_idx <- upper.tri(rho_true)
  
  edge_true <- abs(rho_true[upper_idx]) > 1e-12
  edge_hat  <- abs(rho_hat[upper_idx])  > 0   # 완전 0만 제거
  
  TP <- sum(edge_true & edge_hat)
  FP <- sum(!edge_true & edge_hat)
  TN <- sum(!edge_true & !edge_hat)
  FN <- sum(edge_true & !edge_hat)
  
  TPR <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  FPR <- ifelse(FP + TN > 0, FP / (FP + TN), 0)
  
  return(c(TPR = TPR, FPR = FPR))
}



load(sprintf("./%s_result/data%02d_result.V.RData", 'hub', 23))
load(sprintf("./%s_bic/data%02d.V.RData", 'hub', 23))  # test_bic.U
for (i in 1:20) {
  print(compute_tpr_fpr(sim_result_V[[i]]$rho, result$V))
}



library(pROC)   # AUC 계산 편리

plot_roc <- function(sim_result, rho_true, lambda_seq, main_title, out_file, test_bic) {
  tpr <- fpr <- rep(NA, length(lambda_seq))
  
  for (i in seq_along(lambda_seq)) {
    l <- lambda_seq[i]
    rho_hat <- sim_result[[as.character(l)]]$rho
    vals <- compute_tpr_fpr(rho_hat, rho_true)
    tpr[i] <- vals["TPR"]
    fpr[i] <- vals["FPR"]
  }
  
  # --- 정렬용 (플롯/auc 계산) ---
  ord <- order(fpr, tpr)
  fpr_plot <- c(0, fpr[ord], 1)
  tpr_plot <- c(0, tpr[ord], 1)
  
  # --- trapezoidal AUC ---
  auc_val <- sum(diff(fpr_plot) * (head(tpr_plot, -1) + tail(tpr_plot, -1)) / 2)
  
  # --- BIC 최소 지점 ---
  min_idx <- which.min(test_bic)
  best_lambda <- lambda_seq[min_idx]
  rho_best <- sim_result[[as.character(best_lambda)]]$rho
  best_vals <- compute_tpr_fpr(rho_best, rho_true)
  
  # --- 그리기 ---
  png(out_file, width=1200, height=800, res=150)
  plot(fpr_plot, tpr_plot, type="b", lwd=2, pch=19,
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=sprintf("%s (AUC=%.3f)", main_title, auc_val))
  abline(0,1,col="grey",lty=2)
  points(best_vals["FPR"], best_vals["TPR"], col="blue", pch=19, cex=1.5)
  dev.off()
  
  return(auc_val)
}


exclude_idx <- c(33, 38, 44)

for (type in unique(task_list$type)) {
  for (idx in setdiff(1:50, exclude_idx)) {
    
    ## true 값
    load(sprintf("./%s_data/data%02d.RData", type, idx))
    
    ## U
    load(sprintf("./%s_result/data%02d_result.U.RData", type, idx))
    load(sprintf("./%s_bic/data%02d.U.RData", type, idx))  # test_bic.U
    out_file_U <- sprintf("./plot1.ROC_curve/%s_data%02d.U.png", type, idx)
    auc_U <- plot_roc(sim_result_U, result$U_inv, lambda_seq,
                      main_title = sprintf("%s | data%02d (U)", type, idx),
                      out_file = out_file_U,
                      test_bic = test_bic.U)
    
    ## V
    load(sprintf("./%s_result/data%02d_result.V.RData", type, idx))
    load(sprintf("./%s_bic/data%02d.V.RData", type, idx))  # test_bic.V
    out_file_V <- sprintf("./plot1.ROC_curve/%s_data%02d.V.png", type, idx)
    auc_V <- plot_roc(sim_result_V, result$V_inv, lambda_seq,
                      main_title = sprintf("%s | data%02d (V)", type, idx),
                      out_file = out_file_V,
                      test_bic = test_bic.V)
    
    cat(sprintf("Saved ROC plots: %s data%02d | AUC(U)=%.3f, AUC(V)=%.3f\n",
                type, idx, auc_U, auc_V))
  }
}
