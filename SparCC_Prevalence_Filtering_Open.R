#==================================================#
# 0. 사전 세팅: 경로 및 데이터 불러오기
#==================================================#

# 작업 디렉토리 설정
setwd("Set your working directory")

# 데이터 불러오기
asv_table <- read.csv("your table output", row.names = 1, check.names = FALSE)
meta <- read.csv("your data", stringsAsFactors = FALSE)

# 공통 샘플 기준 정렬
common_samples <- intersect(colnames(asv_table), meta$SRA)
asv_table <- asv_table[, common_samples]
meta <- meta[meta$SRA %in% common_samples, ]

# Cluster별 데이터 추출
cluster_list <- split(meta$SRA, meta$Cluster)
asv_by_cluster <- lapply(cluster_list, function(sra_ids) {
  asv_table[, sra_ids, drop = FALSE]
})

# ✅ Open 데이터 선택
asv_open <- asv_by_cluster$open


#==================================================#
# 1. 다양한 prevalence 값에 대한 SparCC 반복 분석
#==================================================#

# 패키지 로드
library(SpiecEasi)

# prevalence 값 목록
prevalence_values <- c(your specified value)

# 결과 저장용 객체
sparcc_results_list <- list()
asv_counts_summary <- data.frame()

# 반복 분석 시작
for (prev in prevalence_values) {
  
  cat("\n==============================\n")
  cat("🔹 Processing min_prevalence =", prev, "\n")
  cat("==============================\n")
  
  # [1] 희귀 ASV 필터링
  initial_asv_count <- nrow(asv_open)
  asv_filtered <- asv_open[rowSums(asv_open > 0) >= prev * ncol(asv_open), ]
  filtered_asv_count <- nrow(asv_filtered)
  removed_asv_count <- initial_asv_count - filtered_asv_count
  
  # [2] 샘플 정리 및 정규화
  asv_filtered <- asv_filtered[, colSums(asv_filtered) > 0]
  asv_rel <- sweep(asv_filtered, 2, colSums(asv_filtered), FUN = "/")
  asv_rel[asv_rel == 0] <- 1e-6  # pseudocount
  
  # [3] SparCC 수행
  asv_mat_t <- t(as.matrix(asv_rel))
  set.seed(121)
  sparcc_result <- sparcc(asv_mat_t)
  
  # [4] 결과 저장
  key <- paste0("prev_", prev*100, "pct")
  sparcc_results_list[[key]] <- sparcc_result
  
  # 필터링 요약 저장
  asv_counts_summary <- rbind(asv_counts_summary, data.frame(
    Prevalence = prev,
    Initial_ASV = initial_asv_count,
    Filtered_ASV = filtered_asv_count,
    Removed_ASV = removed_asv_count
  ))
  
  # 정규화된 ASV 테이블 저장
  write.csv(asv_rel, file = paste0("asv_open_rel_121_", key, ".csv"))
  
  # SparCC correlation matrix 저장
  write.csv(sparcc_result$Cor, file = paste0("sparcc_cor_open_121_", key, ".csv"))
}

# 전체 결과 저장
save(sparcc_results_list, file = "sparcc_results_open_all_prevalence_121.RData")
write.csv(asv_counts_summary, file = "ASV_filtering_summary_open_121.csv", row.names = FALSE)
