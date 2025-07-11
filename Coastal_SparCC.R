#================================================================================#
# Working Directory
#================================================================================#

setwd <- ("your working directory")

#================================================================================#
# 2.3.1 Coastal: For Seed(121)
#================================================================================#

# [1] 기본 파일 불러오기 + 이름치환

# 예시: 7% prevalence 기준 파일 불러오기
rel_abund <- read.csv("your relative abundance table", row.names = 1, check.names = FALSE)
cor_mat_raw <- read.csv("your correlation table", row.names = 1, check.names = FALSE)

# 행/열 이름 치환: relative abundance 파일의 행 이름 (ASV ID)을 사용
asv_ids <- rownames(rel_abund)
colnames(cor_mat_raw) <- asv_ids
rownames(cor_mat_raw) <- asv_ids

# 결과 확인
round(cor_mat_raw[1:5, 1:5], 3)


# [2] SparCC 상관계수 → 엣지 리스트 변환

# 상삼각행렬 인덱스 추출
edge_idx <- which(upper.tri(cor_mat_raw), arr.ind = TRUE)

# 엣지리스트 생성
edges <- data.frame(
  Source = rownames(cor_mat_raw)[edge_idx[, 1]],
  Target = colnames(cor_mat_raw)[edge_idx[, 2]],
  Weight = cor_mat_raw[edge_idx]
)

head(edges)


# [3] Threshold 설정 및 엣지 필터링

threshold <- XXXXX
abs_weights <- abs(edges$Weight)
quantile_value <- ecdf(abs_weights)(threshold)
percentile <- round((1 - quantile_value) * 100, 2)

cat("절대 상관계수", threshold, "은 상위", percentile, "%\n")

edges_filtered <- subset(edges, abs(Weight) >= threshold)
cat("유의미한 엣지 수:", nrow(edges_filtered), "\n")


# [4] Positive vs. Negative 엣지 분류

pos_edges <- subset(edges_filtered, Weight > 0)
neg_edges <- subset(edges_filtered, Weight < 0)

option <- 1  # 1 = Positive, 2 = Negative, 3 = All

if (option == 1) {
  selected_edges <- pos_edges
} else if (option == 2) {
  selected_edges <- neg_edges
} else {
  selected_edges <- edges_filtered
}


# [5] 노드 Degree 계산 및 그룹화

all_nodes <- c(selected_edges$Source, selected_edges$Target)
node_degree <- as.data.frame(table(all_nodes))
colnames(node_degree) <- c("Node", "Degree")

# ND vs D 분류 (예시: 이름 내 "-ND_" / "-D" 포함 여부 기준)
node_degree$Group <- ifelse(grepl("-ND_", node_degree$Node), "ND", 
                            ifelse(grepl("-D", node_degree$Node), "D", "Other"))

node_degree_filtered <- subset(node_degree, Group %in% c("ND", "D"))


# [6] 통계 검정 및 시각화

library(dplyr)
library(ggplot2)

# 그룹별 평균
summary_stats <- node_degree_filtered %>%
  group_by(Group) %>%
  summarise(
    Mean_Degree = mean(Degree),
    SD_Degree = sd(Degree),
    Count = n()
  )
print(summary_stats)

# 정규성 테스트
shapiro_nd <- shapiro.test(node_degree_filtered$Degree[node_degree_filtered$Group == "ND"])
shapiro_d  <- shapiro.test(node_degree_filtered$Degree[node_degree_filtered$Group == "D"])

# 적절한 검정 수행
if (shapiro_nd$p.value >= 0.05 & shapiro_d$p.value >= 0.05) {
  result <- t.test(Degree ~ Group, data = node_degree_filtered)
  test_method <- "Welch's t-test"
} else {
  result <- wilcox.test(Degree ~ Group, data = node_degree_filtered)
  test_method <- "Wilcoxon test"
}
print(result)

# 박스플롯 시각화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 📊 그래프 객체 생성: Coastal (Seed121, 7%)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
p_coastal_121_7 <- ggplot(node_degree_filtered, aes(x = Group, y = Degree, fill = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.7,
    color = "black",
    size = 0.7  # 수염 및 박스 테두리 강조
  ) +
  labs(
    x = "Group",
    y = "Degree"
  ) +
  scale_fill_manual(values = c("ND" = "#F8766D", "D" = "#00BFC4")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_line(size = 0.6, color = "black"),
    axis.ticks = element_line(size = 0.4, color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_line(color = "grey90", size = 0.1),
    panel.border = element_rect(color = "black", fill = NA, size = 1.0),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 💾 TIFF 저장 (정사각형 비율 통일: 6 × 6 inch)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ggsave(
  filename = "Your Chemosynthetic Plot",
  plot = p_coastal_121_7,
  width = 6,
  height = 6,      # 정사각형 비율로 강제
  dpi = 600,
  units = "in",
  device = "tiff",
  compression = "lzw"
)
