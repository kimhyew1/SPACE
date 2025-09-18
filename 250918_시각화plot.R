###############################################################################
# AUC_boxplot

library(ggplot2)


df_auc <- data.frame(
  value = c(auc_list.lambda1$hub$U,
            auc_list.lambda1$hub$V,
            auc_list.lambda1$random$U,
            auc_list.lambda1$random$V,
            auc_list.lambda1$cluster$U,
            auc_list.lambda1$cluster$V),
  type = rep(c("hub U", "hub V",
               "random U", "random V",
               "cluster U", "cluster V"),
             times = c(length(auc_list.lambda1$hub$U),
                       length(auc_list.lambda1$hub$V),
                       length(auc_list.lambda1$random$U),
                       length(auc_list.lambda1$random$V),
                       length(auc_list.lambda1$cluster$U),
                       length(auc_list.lambda1$cluster$V)))
)

png("./1_Specified_lambda/AUC_boxplot.png", width=1200, height=800, res=150)
ggplot(df_auc, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "AUC by Type (Lambda 1)", y = "AUC", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


df_auc <- data.frame(
  value = c(auc_list.lambda2$hub$U,
            auc_list.lambda2$hub$V,
            auc_list.lambda2$random$U,
            auc_list.lambda2$random$V,
            auc_list.lambda2$cluster$U,
            auc_list.lambda2$cluster$V),
  type = rep(c("hub U", "hub V",
               "random U", "random V",
               "cluster U", "cluster V"),
             times = c(length(auc_list.lambda2$hub$U),
                       length(auc_list.lambda2$hub$V),
                       length(auc_list.lambda2$random$U),
                       length(auc_list.lambda2$random$V),
                       length(auc_list.lambda2$cluster$U),
                       length(auc_list.lambda2$cluster$V)))
)

png("./2_Auto_lambda/AUC_boxplot.png", width=1200, height=800, res=150)
ggplot(df_auc, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "AUC by Type (Lambda 2)", y = "AUC", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


which.max(auc_list.lambda2$random$U)
which.min(auc_list.lambda2$random$U)
which.max(auc_list.lambda2$random$V)
which.min(auc_list.lambda2$random$V)



###############################################################################
# precision matrix Heatmap


library(pheatmap)

# heatmap 그리는 함수 (U/V 공용)
pheatmap_custom <- function(mat, filename, title) {
  p <- nrow(mat)
  rownames(mat) <- 1:p
  colnames(mat) <- 1:p
  
  # === 모든 값이 0이면 저장하지 않고 넘어감 ===
  if (all(mat == 0)) {
    message(sprintf("Skip heatmap: %s (all zeros)", title))
    return(invisible(NULL))
  }
  
  val_range <- range(mat, na.rm = TRUE)
  max_abs <- max(abs(val_range))   # 절댓값 기준 최대
  
  # 항상 -max_abs ~ +max_abs로 스케일 고정 (0=흰색)
  breaks <- seq(-max_abs, max_abs, length.out = 51)
  colors <- colorRampPalette(c("blue", "white", "red"))(50)
  
  pheatmap(mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           color = colors,
           breaks = breaks,      # 색상 구간 고정
           main = title,
           filename = filename,
           width = 8, height = 6,
           legend = TRUE,
           fontsize_row = 6,
           fontsize_col = 6)
}


for (type in names(bic_min_idx.lambda1)) {
  for (row in 1:nrow(bic_min_idx.lambda1[[type]])) {
    idx <- bic_min_idx.lambda1[[type]]$data[row]
    
    # --- U, V 결과 불러오기 ---
    load(sprintf("./%s/%s_result/data%02d_result.U.RData", input_dir, type, idx))
    load(sprintf("./%s/%s_result/data%02d_result.V.RData", input_dir, type, idx))
    
    sim_result_U <- get("sim_result_U")
    sim_result_V <- get("sim_result_V")
    
    # --- 최소 index 가져오기 ---
    U_min_idx <- bic_min_idx.lambda1[[type]]$U_min_idx[row]
    V_min_idx <- bic_min_idx.lambda1[[type]]$V_min_idx[row]
    
    # --- precision matrix 선택 ---
    rho_hat_U <- sim_result_U[[U_min_idx]]$rho
    # diag(rho_hat_U) <- sim_result_U[[U_min_idx]]$weight
    # diag(rho_hat_U) <- 0
    
    rho_hat_V <- sim_result_V[[V_min_idx]]$rho
    # diag(rho_hat_V) <- sim_result_V[[V_min_idx]]$weight
    # diag(rho_hat_V) <- 0
    
    # --- Heatmap 저장 (U, V 공용 코드) ---
    for (label in c("U", "V")) {
      rho_mat <- if (label == "U") rho_hat_U else rho_hat_V
      filename <- sprintf("./%s/plot.heatmap/%s_data%02d.%s.png", input_dir, type, idx, label)
      title <- sprintf("%s | data%02d (%s)", type, idx, label)
      
      pheatmap_custom(rho_mat, filename, title)
    }
    
    # diag(sim_result_V[[V_min_idx]]$rho) <- 0
    prec_V <- - sim_result_V[[V_min_idx]]$rho * sqrt(outer(sim_result_V[[V_min_idx]]$weight, sim_result_V[[V_min_idx]]$weight))
    # diag(prec_V) <- 0
    # diag(sim_result_U[[U_min_idx]]$rho) <- 0
    prec_U <- - sim_result_U[[U_min_idx]]$rho * sqrt(outer(sim_result_U[[U_min_idx]]$weight, sim_result_U[[U_min_idx]]$weight))
    # diag(prec_U) <- 0
    # --- U matrix rescale for identification ---
    u11 <- prec_U[1,1]
    if (u11 != 0) {
      prec_U <- prec_U / u11   # (1,1) = 1 로 고정
    }
    
    
    # --- V ⊗ U (Kronecker product) ---
    rho_hat_kron <- kronecker(prec_V, prec_U)
    filename_kron <- sprintf("./%s/plot.heatmap/%s_data%02d.kron.png", input_dir, type, idx)
    title_kron <- sprintf("%s | data%02d (V ⊗ U)", type, idx)
    pheatmap_custom(rho_hat_kron, filename_kron, title_kron)
    
    cat(sprintf("Saved heatmaps: %s data%02d (U, V)\n", type, idx))
  }
}




for (type in c("random")) {
  for (row in 1:nrow(bic_min_idx.lambda1[[type]])) {
    idx <- bic_min_idx.lambda1[[type]]$data[row]
    
    # --- U, V 결과 불러오기 ---
    load(sprintf("./%s_data/data%02d.RData", type, idx))
    result_U <- get("result")
    
    load(sprintf("./%s_data/data%02d.RData", type, idx))
    result_V <- get("result")
    
    # --- precision matrix 선택 ---
    rho_U <- result_U$U
    rho_V <- result_V$V
    
    # --- Heatmap 저장 (U, V 공용 코드) ---
    for (label in c("U", "V")) {
      rho_mat <- if (label == "U") rho_U else rho_V
      filename <- sprintf("./plot.heatmap/%s_data%02d.%s.png", type, idx, label)
      title <- sprintf("%s | data%02d (%s)", type, idx, label)
      
      pheatmap_custom(rho_mat, filename, title)
    }
    
    # --- U matrix rescale for identification ---
    u11 <- prec_U[1,1]
    if (u11 != 0) {
      prec_U <- prec_U / u11   # (1,1) = 1 로 고정
    }
    
    # --- V ⊗ U (Kronecker product) ---
    rho_hat_kron <- kronecker(rho_V, rho_U)
    filename_kron <- sprintf("./plot.heatmap/%s_data%02d.kron.png", type, idx)
    title_kron <- sprintf("%s | data%02d (V ⊗ U)", type, idx)
    pheatmap_custom(rho_hat_kron, filename_kron, title_kron)
    
    cat(sprintf("Saved heatmaps: %s data%02d (U, V)\n", type, idx))
  }
}






####################################################
# Calculate norm (L1, Frobenius)
for (type in names(bic_min_idx.lambda1)) {
  for (row in 1:nrow(bic_min_idx.lambda1[[type]])) {
    idx <- bic_min_idx.lambda1[[type]]$data[row]
    
    # --- True 값 불러오기 ---
    load(sprintf("./%s_data/data%02d.RData", type, idx))  # result 객체 불러옴
    
    # --- U, V 결과 불러오기 ---
    load(sprintf("./%s/%s_result/data%02d_result.U.RData", input_dir, type, idx))
    load(sprintf("./%s/%s_result/data%02d_result.V.RData", input_dir, type, idx))
    
    sim_result_U <- get("sim_result_U")
    sim_result_V <- get("sim_result_V")
    
    # --- 최소 index 가져오기 ---
    U_min_idx <- bic_min_idx.lambda1[[type]]$U_min_idx[row]
    V_min_idx <- bic_min_idx.lambda1[[type]]$V_min_idx[row]
    
    # --- 람다 시퀀스 불러오기 ---
    lambda_seq_U <- as.numeric(names(sim_result_U))
    lambda_seq_V <- as.numeric(names(sim_result_V))
    
    # --- True precision matrices ---
    true_U <- result$U
    true_V <- result$V
    
    for (label in c("U", "V")) {
      sim_result <- if (label == "U") sim_result_U else sim_result_V
      true_mat   <- if (label == "U") true_U else true_V
      lambda_seq <- if (label == "U") lambda_seq_U else lambda_seq_V
      min_idx    <- if (label == "U") U_min_idx else V_min_idx
      
      # --- Norm 계산 ---
      l1_vals <- sapply(seq_along(sim_result), function(i) {
        l1_calculate(true_mat, sim_result[[i]]$rho)
      })
      frob_vals <- sapply(seq_along(sim_result), function(i) {
        frob_calculate(true_mat, sim_result[[i]]$rho)
      })
      
      # --- Plot L1 ---
      png(sprintf("./%s/plot.norm/%s_data%02d.%s_L1.png", input_dir, type, idx, label),
          width = 1200, height = 800, res = 150)
      plot(log(lambda_seq), l1_vals, type = "l", lwd = 2,
           xlab = "log(lambda)", ylab = "L1 Norm",
           main = sprintf("%s | data%02d (%s)", type, idx, label))
      points(log(lambda_seq[min_idx]), l1_vals[min_idx], col = "blue", pch = 19)
      dev.off()
      
      # --- Plot Frobenius ---
      png(sprintf("./%s/plot.norm/%s_data%02d.%s_Frob.png", input_dir, type, idx, label),
          width = 1200, height = 800, res = 150)
      plot(log(lambda_seq), frob_vals, type = "l", lwd = 2,
           xlab = "log(lambda)", ylab = "Frobenius Norm",
           main = sprintf("%s | data%02d (%s)", type, idx, label))
      points(log(lambda_seq[min_idx]), frob_vals[min_idx], col = "blue", pch = 19)
      dev.off()
    }
    
    cat(sprintf("Saved Norm plots: %s data%02d (U, V)\n", type, idx))
  }
}
