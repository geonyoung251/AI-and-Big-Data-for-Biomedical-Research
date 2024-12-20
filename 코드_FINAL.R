#0 패키지 다운로드
install.packages('readxl')
library(readxl)

install.packages('dplyr')
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

install.packages("brunnermunzel")
library(brunnermunzel)

install.packages("pROC")
library(pROC)

install.packages("caret")
library(caret)

install.packages("glmnet")
library(glmnet)

install.packages("ROSE")
library(ROSE)

install.packages("PRROC")
library(PRROC)

#1 데이터 불러오기 - drug_response
data_drug_response <- read.table("C:\\Users\\geon0\\OneDrive\\문서\\대학\\M4-2_이주상 교수님 수업\\뉴프로젝트\\drug_response.txt", header = TRUE, sep = "\t")

#2 데이터 불러오기 - LunitScope
data_LunitScope <- read_excel("C:\\Users\\geon0\\OneDrive\\문서\\대학\\M4-2_이주상 교수님 수업\\뉴프로젝트\\LunitScope.xlsx", sheet = 2)

#3 데이터 불러오기 - PANCANCER demographic
data_demographics <- read.table("C:\\Users\\geon0\\OneDrive\\문서\\대학\\M4-2_이주상 교수님 수업\\뉴프로젝트\\PANCANCER.txt", header = TRUE, sep = "\t")

#4 QC; False 제외하기 / blank, [Not Evaluated], [Unknown] 제거하기
data_LunitScope_QC <- data_LunitScope[data_LunitScope$`QC` != FALSE, ]

data_demographics <- data_demographics %>% select(-sample)
data_demographics <- data_demographics %>% distinct()
data_demographics_QC <- subset(data_demographics, !(race %in% c("", "[Not Evaluated]", "[Unknown]")))


#5 drug response responder / nonresponder로 구분하기
data_drug_response$RNR <- ifelse(data_drug_response$response %in% c("Partial Response", "Complete Response"), 
                                 "Responder", 
                                 ifelse(data_drug_response$response %in% c("Clinical Progressive Disease", "Stable Disease"), 
                                        "Nonresponder", 
                                        NA))

#6 R / NR 수 확인하기
table(data_drug_response$RNR)

#7 sample_id 형식 맞추기
data_LunitScope_QC$sample_id <- sub("^([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1", data_LunitScope_QC$sample_id)

#8 sample_id, patient.arr 기준으로 병합하기

# 병합 수행
merged_data_tot <- data_drug_response %>%
  inner_join(data_LunitScope_QC, by = c("patient.arr" = "sample_id"))


#추가 가공
merged_data_tot$RNR_binary <- as.numeric(factor(merged_data_tot$RNR, levels = c("Nonresponder", "Responder"))) - 1

#9 병합 후 cancers, drug.name으로 분류한 후 count 수 세기

count_df_tot <- merged_data_tot %>%
  group_by(cancers, drug.name) %>%
  summarise(count = n()) %>%
  ungroup() 


#10 분류 후 세부 통계량 확인 함수
count_summary <- function(cancer_name, drug_name) {
  # 데이터 필터링
  filtered_data <- merged_data_tot %>%
    filter(cancers == cancer_name & drug.name == drug_name)
  
  filtered_data$response <- factor(
    filtered_data$response,
    levels = c("Clinical Progressive Disease", "Stable Disease", "Partial Response", "Complete Response")
  )
  
  # 각 변수의 빈도수 계산
  rnr_count <- table(filtered_data$RNR)
  response_count <- table(filtered_data$response)
  immune_phenotype_count <- table(filtered_data$immune_phenotype)
  
  # 결과 출력
  cat("\nSummary for:", cancer_name, "-", drug_name, "\n")
  cat("\nRNR Count:\n")
  print(rnr_count)
  
  cat("\nResponse Count:\n")
  print(response_count)
  
  cat("\nImmune Phenotype Count:\n")
  print(immune_phenotype_count)
  
  return(list(
    RNR = rnr_count,
    Response = response_count,
    Immune_Phenotype = immune_phenotype_count
  ))
}

#11 Plot function
graph_plot <- function(cancer_name, drug_name, number) {
  #데이터 필터링 후 데이터프레임에 저장
  new_data <- merged_data_tot %>%
    filter(cancers == cancer_name & drug.name == drug_name)
  dataset_name <- paste0("tot_data_check_", number)
  
  assign(dataset_name, new_data, envir = .GlobalEnv)
  
  dynamic_data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # 비율 플롯 데이터 처리
  data_prop_1 <- dynamic_data_check %>%
    group_by(immune_phenotype, RNR) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(immune_phenotype) %>%
    mutate(prop = count / sum(count))  # 그룹별 비율 계산
  
  data_prop_2 <- dynamic_data_check %>%
    group_by(RNR, immune_phenotype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(RNR) %>%
    mutate(prop = count / sum(count))  # 그룹별 비율 계산
  
  # 플롯 생성 및 리스트에 저장
  plots <- list()
  
  plots$p1 <- ggplot(dynamic_data_check, aes(x = immune_phenotype, fill = RNR)) +
    geom_bar(position = "dodge") +
    labs(title = paste("Distribution of RNR by Immune Phenotype_", cancer_name, " - ", drug_name),
         x = "Immune Phenotype",
         y = "Count") +
    theme_minimal()
  
  plots$p2 <- ggplot(dynamic_data_check, aes(x = RNR, y = IS, fill = RNR)) +
    geom_boxplot() +
    labs(title = paste("IS Distribution by RNR_", cancer_name, " - ", drug_name),
         x = "RNR",
         y = "IS Score") +
    theme_minimal()
  
  plots$p3 <- ggplot(data_prop_1, aes(x = immune_phenotype, y = prop, fill = RNR)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = paste("Proportion of RNR by Immune Phenotype_", cancer_name, " - ", drug_name),
      x = "Immune Phenotype",
      y = "Proportion") +
    scale_y_continuous(labels = scales::percent) +  # y축을 백분율로 표시
    theme_minimal() +
    scale_fill_manual(values = c("Nonresponder" = "red", "Responder" = "blue"))
  
  plots$p4 <- ggplot(data_prop_2, aes(x = RNR, y = prop, fill = immune_phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = paste("Proportion of Immune Phenotype by RNR_", cancer_name, " - ", drug_name),
      x = "RNR",
      y = "Proportion") +
    scale_y_continuous(labels = scales::percent) +  # y축을 백분율로 표시
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")  # 보기 좋은 색상 팔레트
  
  output_dir <- paste0(number, "_plots_", cancer_name, " - ", drug_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  ggsave(filename = file.path(output_dir, paste0("plot_", number, "_p1.png")), plot = plots$p1, width = 8, height = 6)
  ggsave(filename = file.path(output_dir, paste0("plot_", number, "_p2.png")), plot = plots$p2, width = 8, height = 6)
  ggsave(filename = file.path(output_dir, paste0("plot_", number, "_p3.png")), plot = plots$p3, width = 8, height = 6)
  ggsave(filename = file.path(output_dir, paste0("plot_", number, "_p4.png")), plot = plots$p4, width = 8, height = 6)
  # 플롯 리스트 반환
  return(plots)
}


#12 통계 검정 function
statistical_test <- function(cancer_name, drug_name, number) {
  #데이터 필터링 후 데이터프레임에 저장
  new_data_1 <- merged_data_tot %>%
    filter(cancers == cancer_name & drug.name == drug_name)
  dataset_name_1 <- paste0("tot_data_check_", number)
  
  assign(dataset_name_1, new_data_1, envir = .GlobalEnv)
  
  dynamic_data_check <- get(dataset_name_1, envir = .GlobalEnv)
  
  #immune_phenotype - RNR 검정
  if (length(unique(dynamic_data_check$immune_phenotype)) >= 2) {
    # Kruskal-Wallis Test
    new_data_2 <- kruskal.test(RNR_binary ~ immune_phenotype, data = dynamic_data_check)
    dataset_name_2 <- paste0("kruskal_result_ip_rnr_", number)
    assign(dataset_name_2, new_data_2, envir = .GlobalEnv)
    dynamic_kruskal_result_ip_rnr <- get(dataset_name_2, envir = .GlobalEnv)
    print(dynamic_kruskal_result_ip_rnr)
    
    # Pairwise Wilcoxon test
    new_data_3 <- pairwise.wilcox.test(
      x = dynamic_data_check$RNR_binary,
      g = dynamic_data_check$immune_phenotype
    )
    dataset_name_3 <- paste0("pairwise_result_ip_rnr_", number)
    assign(dataset_name_3, new_data_3, envir = .GlobalEnv)
    dynamic_pairwise_result_ip_rnr <- get(dataset_name_3, envir = .GlobalEnv)
    print(dynamic_pairwise_result_ip_rnr)
  }
  
  #IS - RNR 검정
  if (length(unique(dynamic_data_check$RNR)) == 2) {
    new_data_5 <- wilcox.test(IS ~ RNR, data = dynamic_data_check)
    dataset_name_5 <- paste0("wilcoxon_result_rnr_is_", number)
    assign(dataset_name_5, new_data_5, envir = .GlobalEnv)
    dynamic_wilcoxon_result_rnr_is <- get(dataset_name_5, envir = .GlobalEnv)
    print(dynamic_wilcoxon_result_rnr_is)
  }
}

#13 시행
{
plot_list_1 <- graph_plot("BRCA", "Cyclophosphamide", 1)
print(plot_list_1$p1)
print(plot_list_1$p2)
print(plot_list_1$p3)
print(plot_list_1$p4)
print(plot_list_1$p5)
print(plot_list_1$p6)
print(plot_list_1$p7)
print(plot_list_1$p8)

statistical_test("BRCA", "Cyclophosphamide", 1)


plot_list_2 <- graph_plot("STAD", "Fluorouracil", 2)
print(plot_list_2$p1)
print(plot_list_2$p2)
print(plot_list_2$p3)
print(plot_list_2$p4)
print(plot_list_2$p5)
print(plot_list_2$p6)
print(plot_list_2$p7)
print(plot_list_2$p8)

statistical_test("STAD", "Fluorouracil", 2)


plot_list_3 <- graph_plot("PAAD", "Gemcitabine", 3)
print(plot_list_3$p1)
print(plot_list_3$p2)
print(plot_list_3$p3)
print(plot_list_3$p4)
print(plot_list_3$p5)
print(plot_list_3$p6)
print(plot_list_3$p7)
print(plot_list_3$p8)

statistical_test("PAAD", "Gemcitabine", 3)


plot_list_4 <- graph_plot("CESC", "Cisplatin", 4)
print(plot_list_4$p1)
print(plot_list_4$p2)
print(plot_list_4$p3)
print(plot_list_4$p4)
print(plot_list_4$p5)
print(plot_list_4$p6)
print(plot_list_4$p7)
print(plot_list_4$p8)

statistical_test("CESC", "Cisplatin", 4)


plot_list_5 <- graph_plot("BRCA", "Doxorubicin", 5)
print(plot_list_5$p1)
print(plot_list_5$p2)
print(plot_list_5$p3)
print(plot_list_5$p4)
print(plot_list_5$p5)
print(plot_list_5$p6)
print(plot_list_5$p7)
print(plot_list_5$p8)

statistical_test("BRCA", "Doxorubicin", 5)


plot_list_6 <- graph_plot("TGCT", "Bleomycin", 6)
print(plot_list_6$p1)
print(plot_list_6$p2)
print(plot_list_6$p3)
print(plot_list_6$p4)
print(plot_list_6$p5)
print(plot_list_6$p6)
print(plot_list_6$p7)
print(plot_list_6$p8)

statistical_test("TGCT", "Bleomycin", 6)


plot_list_7 <- graph_plot("TGCT", "Etoposide", 7)
print(plot_list_7$p1)
print(plot_list_7$p2)
print(plot_list_7$p3)
print(plot_list_7$p4)
print(plot_list_7$p5)
print(plot_list_7$p6)
print(plot_list_7$p7)
print(plot_list_7$p8)

statistical_test("TGCT", "Etoposide", 7)


plot_list_8 <- graph_plot("BLCA", "Gemcitabine", 8)
print(plot_list_8$p1)
print(plot_list_8$p2)
print(plot_list_8$p3)
print(plot_list_8$p4)
print(plot_list_8$p5)
print(plot_list_8$p6)
print(plot_list_8$p7)
print(plot_list_8$p8)

statistical_test("BLCA", "Gemcitabine", 8)


plot_list_9 <- graph_plot("TGCT", "Cisplatin", 9)
print(plot_list_9$p1)
print(plot_list_9$p2)
print(plot_list_9$p3)
print(plot_list_9$p4)
print(plot_list_9$p5)
print(plot_list_9$p6)
print(plot_list_9$p7)
print(plot_list_9$p8)

statistical_test("TGCT", "Cisplatin", 9)


plot_list_10 <- graph_plot("BRCA", "Docetaxel", 10)
print(plot_list_10$p1)
print(plot_list_10$p2)
print(plot_list_10$p3)
print(plot_list_10$p4)
print(plot_list_10$p5)
print(plot_list_10$p6)
print(plot_list_10$p7)
print(plot_list_10$p8)

statistical_test("BRCA", "Docetaxel", 10)


plot_list_11 <- graph_plot("COAD", "Fluorouracil", 11)
print(plot_list_11$p1)
print(plot_list_11$p2)
print(plot_list_11$p3)
print(plot_list_11$p4)
print(plot_list_11$p5)
print(plot_list_11$p6)
print(plot_list_11$p7)
print(plot_list_11$p8)


statistical_test("COAD", "Fluorouracil", 11)


plot_list_12 <- graph_plot("BRCA", "Paclitaxel", 12)
print(plot_list_12$p1)
print(plot_list_12$p2)
print(plot_list_12$p3)
print(plot_list_12$p4)
print(plot_list_12$p5)
print(plot_list_12$p6)
print(plot_list_12$p7)
print(plot_list_12$p8)

statistical_test("BRCA", "Paclitaxel", 12)


plot_list_13 <- graph_plot("BLCA", "Cisplatin", 13)
print(plot_list_13$p1)
print(plot_list_13$p2)
print(plot_list_13$p3)
print(plot_list_13$p4)
print(plot_list_13$p5)
print(plot_list_13$p6)
print(plot_list_13$p7)
print(plot_list_13$p8)

statistical_test("BLCA", "Cisplatin", 13)


plot_list_14 <- graph_plot("UCEC", "Carboplatin", 14)
print(plot_list_14$p1)
print(plot_list_14$p2)
print(plot_list_14$p3)
print(plot_list_14$p4)
print(plot_list_14$p5)
print(plot_list_14$p6)
print(plot_list_14$p7)
print(plot_list_14$p8)

statistical_test("UCEC", "Carboplatin", 14)


plot_list_15 <- graph_plot("COAD", "Leucovorin", 15)
print(plot_list_15$p1)
print(plot_list_15$p2)
print(plot_list_15$p3)
print(plot_list_15$p4)
print(plot_list_15$p5)
print(plot_list_15$p6)
print(plot_list_15$p7)
print(plot_list_15$p8)

statistical_test("COAD", "Leucovorin", 15)


plot_list_16 <- graph_plot("HNSC", "Cisplatin", 16)
print(plot_list_16$p1)
print(plot_list_16$p2)
print(plot_list_16$p3)
print(plot_list_16$p4)
print(plot_list_16$p5)
print(plot_list_16$p6)
print(plot_list_16$p7)
print(plot_list_16$p8)

statistical_test("HNSC", "Cisplatin", 16)
}

#14 결과 확인
# 반복문을 통해 p-value 출력
#14-1 Kruskal Wallis RNR - immune_phenotype
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # 그룹 순서 고정 (Desert -> Excluded -> Inflamed)
  data_check$immune_phenotype <- factor(
    data_check$immune_phenotype,
    levels = c("Desert", "Excluded", "Inflamed")
  )
  
  kruskal_name <- paste0("kruskal_result_ip_rnr_", i)
  kruskal_result <- get(kruskal_name, envir = .GlobalEnv)
  
  pairwise_name <- paste0("pairwise_result_ip_rnr_", i)
  pairwise_result <- get(pairwise_name, envir = .GlobalEnv)
  
  
  # 현재 그룹 확인
  remaining_groups <- unique(data_check$immune_phenotype)
  num_groups <- length(remaining_groups)  # 그룹 개수 계산
  
  # 그룹 개수에 따른 처리
  if (num_groups == 3) {
    # Kruskal-Wallis 결과 출력
    cat("Kruskal-Wallis Test Result for Dataset", i, ":\n")
    cat("  p-value:", kruskal_result$p.value, "\n")
    
    # p-value가 0.05 미만인 경우 Pairwise Wilcoxon 검정 수행
    if (kruskal_result$p.value < 0.05) {
      cat("  Significant! Performing Pairwise Wilcoxon tests...\n")
      cat("  Pairwise p-values for", pairwise_name, ":\n")
      print(pairwise_result$p.value)
      cat("\n")
      
      # 그룹 조합 생성
      groups <- combn(levels(data_check$immune_phenotype), 2, simplify = FALSE)
      
      # Pairwise Wilcoxon 테스트 수행 및 출력
      for (pair in groups) {
        # 데이터 필터링
        subset_data <- subset(data_check, immune_phenotype %in% pair)
        
        # Greater test
        wilcox_greater <- wilcox.test(
          RNR_binary ~ immune_phenotype, 
          data = subset_data, 
          alternative = "greater"
        )
        
        # Lesser test
        wilcox_lesser <- wilcox.test(
          RNR_binary ~ immune_phenotype, 
          data = subset_data, 
          alternative = "less"
        )
        
        # 결과 출력
        cat("    Pairwise Wilcoxon test for groups:", pair[1], "vs", pair[2], "\n")
        cat("      Greater p-value:", wilcox_greater$p.value, "\n")
        cat("      Lesser p-value:", wilcox_lesser$p.value, "\n")
      }
    } 
    else {
      cat("  Not significant. Skipping pairwise tests.\n")
    }
    
  } 
  else if (num_groups == 2) {
    # Kruskal-Wallis 결과 출력
    cat("Kruskal-Wallis Test Result for Dataset", i, ":\n")
    cat("  p-value:", kruskal_result$p.value, "\n")
    if (kruskal_result$p.value < 0.05) {
      cat("  Significant! Performing Pairwise Wilcoxon tests...\n")
      cat("  Pairwise p-values for", pairwise_name, ":\n")
      print(pairwise_result$p.value)
      cat("\n")
      
      # 2개 그룹이 있는 경우: 단순 Wilcoxon 테스트 수행
      subset_data <- subset(data_check, immune_phenotype %in% remaining_groups)
    
      # Greater test
      wilcox_greater <- wilcox.test(
        RNR_binary ~ immune_phenotype, 
        data = subset_data, 
        alternative = "greater"
      )
    
      # Lesser test
      wilcox_lesser <- wilcox.test(
        RNR_binary ~ immune_phenotype, 
        data = subset_data, 
        alternative = "less"
      )
    
      # 결과 출력
      cat("Simple Wilcoxon test for groups:", levels(data_check$immune_phenotype)[1], 
          "vs", levels(data_check$immune_phenotype)[2], "\n")
      cat("  Greater p-value:", wilcox_greater$p.value, "\n")
      cat("  Lesser p-value:", wilcox_lesser$p.value, "\n")
    }
  } else if (num_groups == 1) {
    # 1개 그룹만 있는 경우
    cat("Not enough groups for Wilcoxon or Kruskal-Wallis test. Only one group present:", remaining_groups, "\n")
  }
  
  cat("\n") # 데이터셋 간 구분을 위한 빈 줄
}

#14-2 Wilcoxon Rank Sum IS ~ RNR
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  data_check$RNR <- factor(
    data_check$RNR,
    levels = c("Nonresponder", "Responder")
  )
  
  # 그룹 순서 고정 (Desert -> Excluded -> Inflamed)
  wilcoxon_name <- paste0("wilcoxon_result_rnr_is_", i)
  wilcoxon_result <- get(wilcoxon_name, envir = .GlobalEnv)
  cat("IS ~ RNR Wilcoxon Test Result for Dataset", i, ":\n")
  cat("  p-value:", wilcoxon_result$p.value, "\n")
  
  # Greater test
  wilcox_greater <- wilcox.test(
    IS ~ RNR, 
    data = data_check, 
    alternative = "greater"
  )
  
  # Lesser test
  wilcox_lesser <- wilcox.test(
    IS ~ RNR, 
    data = data_check, 
    alternative = "less"
  )
  
  # 결과 출력
  cat("Simple Wilcoxon test for groups:", levels(data_check$RNR)[1], 
      "vs", levels(data_check$RNR)[2], "\n")
  cat("  Greater p-value:", wilcox_greater$p.value, "\n")
  cat("  Lesser p-value:", wilcox_lesser$p.value, "\n")
}

#15 반복문으로 RNR에 대해 여러 숫자형 변수 유의성 확인하기 + Logistic Regression
#15-1 significant
auc_list_significant <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  significant_vars <- data.frame(
    Variable = character(0),
    Statistic = numeric(0),
    p_value = numeric(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  # Wilcoxon 테스트 수행
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    # Wilcoxon 검정 수행
    wilcox_result <- wilcox.test(
      data_check[[variable_name]] ~ data_check$RNR,
      exact = FALSE  # large datasets may need exact = FALSE
    )
    
    
    # p-value가 0.05 미만인 경우에만 출력
    if (wilcox_result$p.value < 0.05) {
      significant_vars <- rbind(
        significant_vars,
        data.frame(
          Variable = variable_name,
          Statistic = wilcox_result$statistic,
          p_value = wilcox_result$p.value
        )
      )
      cat("Wilcoxon test for variable:", variable_name, "\n")
      cat("  p-value:", wilcox_result$p.value, "\n\n")
    }
  }
  
  significant_variable_names <- significant_vars$Variable
  
  if(length(significant_variable_names) > 0) {
    formula <- as.formula(
      paste("RNR_binary ~", paste(significant_variable_names, collapse = " + "))
    )
    print(formula)
    
    logistic_model <- glm(formula, data = data_check, family = binomial)
    summary_obj <- summary(logistic_model)
    summary_df <- as.data.frame(summary_obj$coefficients)
    summary_df$Variable <- rownames(summary_df)
    rownames(summary_df) <- NULL
    
    summary_name <- paste0("logistic_summary_significant_", i)
    assign(summary_name, summary_df, envir = .GlobalEnv)
    
    predicted_probs <- predict(logistic_model, newdata = data_check, type = "response")
    predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
    
    roc_curve <- roc(data_check$RNR_binary, predicted_probs)
    auc_value <- auc(roc_curve)
    
    pr <- pr.curve(
      scores.class0 = predicted_probs[data_check$RNR_binary == 1],
      scores.class1 = predicted_probs[data_check$RNR_binary == 0],
      curve = TRUE
    )
    auprc_value <- pr$auc.integral
    
    confusion_matrix <- table(Predicted = predicted_classes, Actual = data_check$RNR_binary)
    
    # 2x2 형태로 보장
    all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
    confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
    
    # 원본 confusion_matrix 값 복사
    for (pred in rownames(confusion_matrix)) {
      for (act in colnames(confusion_matrix)) {
        confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
      }
    }
    
    print(confusion_matrix_full)
    
    precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
    recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
    f1_score <- 2*(precision*recall)/ (precision + recall)
    
    output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
    if(!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    plot_path <- file.path(output_folder, paste0("roc_curve_significant_", i, ".png"))
    png(filename = plot_path, width = 800, height = 600)
    plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i), col = "blue")
    text(
      x = 0.6, y = 0.2,
      labels = paste("AUROC =", round(auc_value, 3)),
      col = "red", cex = 1.5
    )
    dev.off()
    
    pr_plot_path <- file.path(output_folder, paste0("pr_curve_significant_", i, ".png"))
    png(filename = pr_plot_path, width = 800, height = 600)
    plot(pr, main = paste("Precision-Recall Curve - Dataset", i), col = "blue")
    text(
      x = 0.6, y = 0.2,
      labels = paste("AUPRC =", round(auprc_value, 3)),
      col = "red", cex = 1.5
    )
    dev.off()
    
    auc_list_significant[[paste0("Dataset_", i)]] <- list(
      AUROC = auc_value,
      AUPRC = auprc_value,
      F1_Score = f1_score
    )
    
    cat("AUROC:", auc_value, "AUPRC:", auprc_value, "F1 Score:", f1_score, "\n\n")
  }
  else {
    cat("No significant variables found for logistic regression.\n\n")
  }
}
print(auc_list_significant)

#15-2 all
auc_list_all <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  all_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    all_vars <- rbind(
      all_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  all_variable_names <- all_vars$Variable
  
  formula <- as.formula(
     paste("RNR_binary ~", paste(all_variable_names, collapse = " + "))
  )
  print(formula)
    
  logistic_model <- glm(formula, data = data_check, family = binomial)
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
    
  summary_name <- paste0("logistic_summary_all_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
    
  predicted_probs <- predict(logistic_model, newdata = data_check, type = "response")
  predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
    
  roc_curve <- roc(data_check$RNR_binary, predicted_probs)
  auc_value <- auc(roc_curve)
  
  pr <- pr.curve(
    scores.class0 = predicted_probs[data_check$RNR_binary == 1],
    scores.class1 = predicted_probs[data_check$RNR_binary == 0],
    curve = TRUE
  )
  auprc_value <- pr$auc.integral
  
  confusion_matrix <- table(Predicted = predicted_classes, Actual = data_check$RNR_binary)
  
  # 2x2 형태로 보장
  all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
  confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
  
  # 원본 confusion_matrix 값 복사
  for (pred in rownames(confusion_matrix)) {
    for (act in colnames(confusion_matrix)) {
      confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
    }
  }
  
  print(confusion_matrix_full)
  
  precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
  recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
  f1_score <- 2*(precision*recall)/ (precision + recall)
    
  output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
  if(!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
    
  plot_path <- file.path(output_folder, paste0("roc_curve_all_", i, ".png"))
  png(filename = plot_path, width = 800, height = 600)
  plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i), col = "blue")
  text(
    x = 0.6, y = 0.2,
    labels = paste("AUROC =", round(auc_value, 3)),
    col = "red", cex = 1.5
  )
  dev.off()
  
  pr_plot_path <- file.path(output_folder, paste0("pr_curve_all_", i, ".png"))
  png(filename = pr_plot_path, width = 800, height = 600)
  plot(pr, main = paste("Precision-Recall Curve - Dataset", i), col = "blue")
  text(
    x = 0.6, y = 0.2,
    labels = paste("AUPRC =", round(auprc_value, 3)),
    col = "red", cex = 1.5
  )
  dev.off()
  
  auc_list_all[[paste0("Dataset_", i)]] <- list(
    AUROC = auc_value,
    AUPRC = auprc_value,
    F1_Score = f1_score
  )
  
  cat("AUROC:", auc_value, "AUPRC:", auprc_value, "F1 Score:", f1_score, "\n\n")
}
print(auc_list_all)

#15-3 IS, IES, IDS
auc_list_three <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 20
  
  three_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")

  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    three_vars <- rbind(
      three_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  three_variable_names <- three_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(three_variable_names, collapse = " + "))
  )
  print(formula)
  
  logistic_model <- glm(formula, data = data_check, family = binomial)
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  summary_name <- paste0("logistic_summary_three_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
  
  predicted_probs <- predict(logistic_model, newdata = data_check, type = "response")
  predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
  
  roc_curve <- roc(data_check$RNR_binary, predicted_probs)
  auc_value <- auc(roc_curve)
  
  pr <- pr.curve(
    scores.class0 = predicted_probs[data_check$RNR_binary == 1],
    scores.class1 = predicted_probs[data_check$RNR_binary == 0],
    curve = TRUE
  )
  auprc_value <- pr$auc.integral
  
  confusion_matrix <- table(Predicted = predicted_classes, Actual = data_check$RNR_binary)
  
  # 2x2 형태로 보장
  all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
  confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
  
  # 원본 confusion_matrix 값 복사
  for (pred in rownames(confusion_matrix)) {
    for (act in colnames(confusion_matrix)) {
      confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
    }
  }
  
  print(confusion_matrix_full)
  
  precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
  recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
  f1_score <- 2*(precision*recall)/ (precision + recall)
  
  output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
  if(!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  plot_path <- file.path(output_folder, paste0("roc_curve_three_", i, ".png"))
  png(filename = plot_path, width = 800, height = 600)
  plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i), col = "blue")
  text(
    x = 0.6, y = 0.2,
    labels = paste("AUROC =", round(auc_value, 3)),
    col = "red", cex = 1.5
  )
  dev.off()
  
  pr_plot_path <- file.path(output_folder, paste0("pr_curve_three_", i, ".png"))
  png(filename = pr_plot_path, width = 800, height = 600)
  plot(pr, main = paste("Precision-Recall Curve - Dataset", i), col = "blue")
  text(
    x = 0.6, y = 0.2,
    labels = paste("AUPRC =", round(auprc_value, 3)),
    col = "red", cex = 1.5
  )
  dev.off()
  
  auc_list_three[[paste0("Dataset_", i)]] <- list(
    AUROC = auc_value,
    AUPRC = auprc_value,
    F1_Score = f1_score
  )
  
  cat("AUROC:", auc_value, "AUPRC:", auprc_value, "F1 Score:", f1_score, "\n\n")
}
print(auc_list_three)

#16 반복문으로 RNR logistic regression 0.8 / 0.2
#16-1 significant
auc_list_significant_mean <- list()
auc_list_significant_stdev <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  significant_vars <- data.frame(
    Variable = character(0),
    Statistic = numeric(0),
    p_value = numeric(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  # Wilcoxon 테스트 수행
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    # Wilcoxon 검정 수행
    wilcox_result <- wilcox.test(
      data_check[[variable_name]] ~ data_check$RNR,
      exact = FALSE  # large datasets may need exact = FALSE
    )
    
    
    # p-value가 0.05 미만인 경우에만 출력
    if (wilcox_result$p.value < 0.05) {
      significant_vars <- rbind(
        significant_vars,
        data.frame(
          Variable = variable_name,
          Statistic = wilcox_result$statistic,
          p_value = wilcox_result$p.value
        )
      )
      cat("Wilcoxon test for variable:", variable_name, "\n")
      cat("  p-value:", wilcox_result$p.value, "\n\n")
    }
  }
  
  significant_variable_names <- significant_vars$Variable
  
  if(length(significant_variable_names) > 0) {
    formula <- as.formula(
      paste("RNR_binary ~", paste(significant_variable_names, collapse = " + "))
    )
    print(formula)
    
    set.seed(123)
    table_rnr <- table(data_check$RNR_binary)
    if (any(table_rnr <= 1)) {
      cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
      next
    }
    
    num_iterations <- 10
    auc_values <- numeric(num_iterations)
    auprc_values <- numeric(num_iterations)
    F1_values <- numeric(num_iterations)
    
    
    for(j in 1:num_iterations) {
      repeat{
       tryCatch({
         repeat {
           # 데이터 0.8 / 0.2로 분할
           train_index <- createDataPartition(data_check$RNR_binary, p = 0.8, list = FALSE)
           train_data <- data_check[train_index, ]
           test_data <- data_check[-train_index, ]
           
           # RNR_binary가 0과 1을 모두 포함하는지 확인
           if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
             break
           }
         }
         
         logistic_model <- glm(formula, data = train_data, family = binomial)
         summary_obj <- summary(logistic_model)
         summary_df <- as.data.frame(summary_obj$coefficients)
         summary_df$Variable <- rownames(summary_df)
         rownames(summary_df) <- NULL
         
         summary_name <- paste0("logistic_summary_significant_", i, "_iteration_", j)
         assign(summary_name, summary_df, envir = .GlobalEnv)
         
         predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
         predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
         
         roc_curve <- roc(test_data$RNR_binary, predicted_probs)
         auc_values[j] <- auc(roc_curve)
         
         pr <- pr.curve(
           scores.class0 = predicted_probs[test_data$RNR_binary == 1],
           scores.class1 = predicted_probs[test_data$RNR_binary == 0],
           curve = TRUE
         )
         
         auprc_values[j] <- pr$auc.integral
         
         
         confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
         
         # 2x2 형태로 보장
         all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
         confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
         
         # 원본 confusion_matrix 값 복사
         for (pred in rownames(confusion_matrix)) {
           for (act in colnames(confusion_matrix)) {
             confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
           }
         }
         
         print(confusion_matrix_full)
         
         precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
         recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
         F1_values[j] <- 2*(precision*recall)/ (precision + recall)
         
         output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
         if(!dir.exists(output_folder)) {
           dir.create(output_folder, recursive = TRUE)
         }
         
         plot_path <- file.path(output_folder, paste0("roc_curve_significant_", i, "_iteration_", j, ".png"))
         png(filename = plot_path, width = 800, height = 600)
         plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
         text(
           x = 0.6, y = 0.2,
           labels = paste("AUROC =", round(auc_values[j], 3)),
           col = "red", cex = 1.5
         )
         dev.off()
         
         pr_plot_path <- file.path(output_folder, paste0("pr_curve_significant_", i, "_iteration_", j, ".png"))
         png(filename = pr_plot_path, width = 800, height = 600)
         plot(pr, main = paste("Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
         text(
           x = 0.6, y = 0.2,
           labels = paste("AUPRC =", round(auprc_values[j], 3)),
           col = "red", cex = 1.5
         )
         dev.off()
         
         cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
         
         break
       }, error = function(e) {
         cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
         cat("Retrying iteration", j, "...\n")
       }) 
      }
      
    }

    auc_list_significant_mean[[paste0("Dataset_", i)]] <- list(
      AUROC = mean(auc_values),
      AUPRC = mean(auprc_values),
      F1_score = mean(F1_values)
    )
    
    auc_list_significant_stdev[[paste0("Dataset_", i)]] <- list(
      AUROC = sd(auc_values),
      AUPRC = sd(auprc_values),
      F1_score = sd(F1_values)
    )
  }
  else {
    cat("No significant variables found for logistic regression.\n\n")
  }
}
print(auc_list_significant_mean)
print(auc_list_significant_stdev)

#16-2 all
auc_list_all_mean <- list()
auc_list_all_stdev <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  all_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    all_vars <- rbind(
      all_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  all_variable_names <- all_vars$Variable
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(all_variable_names, collapse = " + "))
  )
  print(formula)
  
  set.seed(123)
  table_rnr <- table(data_check$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check[train_index, ]
          test_data <- data_check[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(formula, data = train_data, family = binomial)
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("logistic_summary_all_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("roc_curve_all_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("pr_curve_all_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      }) 
    }
    
  }
  
  auc_list_all_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  auc_list_all_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(auc_list_all_mean)
print(auc_list_all_stdev)

#16-3 IS, IES, IDS
auc_list_three_mean <- list()
auc_list_three_stdev <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 20
  
  three_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    three_vars <- rbind(
      three_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  three_variable_names <- three_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(three_variable_names, collapse = " + "))
  )
  print(formula)
  
  set.seed(123)
  # RNR_binary 데이터 확인
  table_rnr <- table(data_check$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check[train_index, ]
          test_data <- data_check[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(formula, data = train_data, family = binomial)
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("logistic_summary_three_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("roc_curve_three_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("pr_curve_three_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      }) 
    }
    
  }
  
  auc_list_three_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  auc_list_three_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(auc_list_three_mean)
print(auc_list_three_stdev)

#R/(R+NR) 확인
for(i in 1:16) {
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  class_distribution <- table(data_check$RNR_binary)
  
  proportion_of_1 <- class_distribution[as.character(1)] / sum(class_distribution)
  
  cat("Proportion of 1 in RNR_binary:", proportion_of_1, "\n")
}

#17 race, age, gender 병합 수행
merged_data_demoadded <- merged_data_tot %>%
  inner_join(data_demographics_QC, by = c("patient.arr" = "X_PATIENT"))

count_df_demoadded <- merged_data_demoadded %>%
  group_by(cancers, drug.name) %>%
  summarise(count = n()) %>%
  ungroup() 



grouping <- function(cancer_name, drug_name, number) {
  new_data <- merged_data_demoadded %>%
    filter(cancers == cancer_name & drug.name == drug_name)
  dataset_name <- paste0("demo_data_check_", number)
  
  assign(dataset_name, new_data, envir = .GlobalEnv)
}

{
grouping("BRCA", "Cyclophosphamide", 1)
grouping("PAAD", "Gemcitabine", 2)
grouping("STAD", "Fluorouracil", 3)
grouping("TGCT", "Etoposide", 4)
grouping("BLCA", "Gemcitabine", 5)
grouping("TGCT", "Bleomycin", 6)
grouping("CESC", "Cisplatin", 7)
grouping("TGCT", "Cisplatin", 8)
grouping("BRCA", "Doxorubicin", 9)
grouping("BLCA", "Cisplatin", 10)
grouping("BRCA", "Paclitaxel", 11)
grouping("BRCA", "Docetaxel", 12)
grouping("HNSC", "Cisplatin", 13)
grouping("LUAD", "Cisplatin", 14)
grouping("UCEC", "Carboplatin", 15)
grouping("UCEC", "Paclitaxel", 16)
}

#R/(R+NR) 확인
for(i in 1:16) {
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  class_distribution <- table(data_check$RNR_binary)
  
  proportion_of_1 <- class_distribution[as.character(1)] / sum(class_distribution)
  
  cat("Proportion of 1 in RNR_binary:", proportion_of_1, "\n")
}

#18 반복문으로 RNR logistic regression 0.8 / 0.2
#18-0 race, age, gender
demo_auc_list_confounder_mean <- list()
demo_auc_list_confounder_stdev <- list()
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  #진짜 formula는 아님
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis", 
      "RNR_binary"
    )
  }
  
  set.seed(111)
  table_rnr <- table(data_check_encoded$RNR_binary)
  
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check_encoded$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check_encoded[train_index, ]
          test_data <- data_check_encoded[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(
          RNR_binary ~ .,
          data = train_data[, selected_columns],
          family = binomial
        )
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("demo_logistic_summary_confounder_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0("demo_", i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("demo_roc_curve_confounder_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("demo_ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("demo_pr_curve_confounder_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("demo_Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      })
    }
  }
  
  demo_auc_list_confounder_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  demo_auc_list_confounder_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(demo_auc_list_confounder_mean)
print(demo_auc_list_confounder_stdev)

#18-1 significant + confounder
demo_auc_list_significant_mean <- list()
demo_auc_list_significant_stdev <- list()
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  significant_vars <- data.frame(
    Variable = character(0),
    Statistic = numeric(0),
    p_value = numeric(0)
  )
  
  # Wilcoxon 테스트 수행
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    # Wilcoxon 검정 수행
    wilcox_result <- wilcox.test(
      data_check[[variable_name]] ~ data_check$RNR,
      exact = FALSE  # large datasets may need exact = FALSE
    )
    
    
    # p-value가 0.05 미만인 경우에만 출력
    if (wilcox_result$p.value < 0.05) {
      significant_vars <- rbind(
        significant_vars,
        data.frame(
          Variable = variable_name,
          Statistic = wilcox_result$statistic,
          p_value = wilcox_result$p.value
        )
      )
      cat("Wilcoxon test for variable:", variable_name, "\n")
      cat("  p-value:", wilcox_result$p.value, "\n\n")
    }
  }
  
  significant_variable_names <- significant_vars$Variable
  variable_names = c(variable_names, significant_variable_names)
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      significant_variable_names,
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis",
      significant_variable_names,
      "RNR_binary"
    )
  }
  
  if(length(significant_variable_names) > 0) {
    set.seed(111)
    table_rnr <- table(data_check_encoded$RNR_binary)
    if (any(table_rnr <= 1)) {
      cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
      next
    }

    num_iterations <- 10
    auc_values <- numeric(num_iterations)
    auprc_values <- numeric(num_iterations)
    F1_values <- numeric(num_iterations)
    
    for(j in 1:num_iterations) {
      repeat{
        tryCatch({
          repeat {
            # 데이터 0.8 / 0.2로 분할
            train_index <- createDataPartition(data_check_encoded$RNR_binary, p = 0.8, list = FALSE)
            train_data <- data_check_encoded[train_index, ]
            test_data <- data_check_encoded[-train_index, ]
            
            # RNR_binary가 0과 1을 모두 포함하는지 확인
            if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
              break
            }
          }
          
          logistic_model <- glm(
            RNR_binary ~ .,
            data = train_data[, selected_columns],
            family = binomial
          )
          summary_obj <- summary(logistic_model)
          summary_df <- as.data.frame(summary_obj$coefficients)
          summary_df$Variable <- rownames(summary_df)
          rownames(summary_df) <- NULL
          
          summary_name <- paste0("demo_logistic_summary_significant_", i, "_iteration_", j)
          assign(summary_name, summary_df, envir = .GlobalEnv)
          
          predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
          predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
          
          roc_curve <- roc(test_data$RNR_binary, predicted_probs)
          auc_values[j] <- auc(roc_curve)
          
          pr <- pr.curve(
            scores.class0 = predicted_probs[test_data$RNR_binary == 1],
            scores.class1 = predicted_probs[test_data$RNR_binary == 0],
            curve = TRUE
          )
          
          auprc_values[j] <- pr$auc.integral
          
          
          confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
          
          # 2x2 형태로 보장
          all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
          confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
          
          # 원본 confusion_matrix 값 복사
          for (pred in rownames(confusion_matrix)) {
            for (act in colnames(confusion_matrix)) {
              confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
            }
          }
          
          print(confusion_matrix_full)
          
          precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
          recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
          F1_values[j] <- 2*(precision*recall)/ (precision + recall)
          
          output_folder <- paste0("demo_", i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
          if(!dir.exists(output_folder)) {
            dir.create(output_folder, recursive = TRUE)
          }
          
          plot_path <- file.path(output_folder, paste0("demo_roc_curve_significant_", i, "_iteration_", j, ".png"))
          png(filename = plot_path, width = 800, height = 600)
          plot.roc(roc_curve, main = paste("demo_ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
          text(
            x = 0.6, y = 0.2,
            labels = paste("AUROC =", round(auc_values[j], 3)),
            col = "red", cex = 1.5
          )
          dev.off()
          
          pr_plot_path <- file.path(output_folder, paste0("demo_pr_curve_significant_", i, "_iteration_", j, ".png"))
          png(filename = pr_plot_path, width = 800, height = 600)
          plot(pr, main = paste("demo_Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
          text(
            x = 0.6, y = 0.2,
            labels = paste("AUPRC =", round(auprc_values[j], 3)),
            col = "red", cex = 1.5
          )
          dev.off()
          
          cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
          
          break
        }, error = function(e) {
          cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
          cat("Retrying iteration", j, "...\n")
        })
      }
    }
    
    demo_auc_list_significant_mean[[paste0("Dataset_", i)]] <- list(
      AUROC = mean(auc_values),
      AUPRC = mean(auprc_values),
      F1_score = mean(F1_values, na.rm = TRUE)
    )
    
    demo_auc_list_significant_stdev[[paste0("Dataset_", i)]] <- list(
      AUROC = sd(auc_values),
      AUPRC = sd(auprc_values),
      F1_score = sd(F1_values, na.rm = TRUE)
    )
  }
  else {
    cat("No significant variables found for logistic regression.\n\n")
  }
}
print(demo_auc_list_significant_mean)
print(demo_auc_list_significant_stdev)

#18-2 all  + confounder
demo_auc_list_all_mean <- list()
demo_auc_list_all_stdev <- list()
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 45
  
  all_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    all_vars <- rbind(
      all_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  all_variable_names <- all_vars$Variable
  
  variable_names = c(variable_names, all_variable_names)
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      all_variable_names,
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis",
      all_variable_names,
      "RNR_binary"
    )
  }
  
  set.seed(111)
  table_rnr <- table(data_check_encoded$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check_encoded$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check_encoded[train_index, ]
          test_data <- data_check_encoded[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(
          RNR_binary ~ .,
          data = train_data[, selected_columns],
          family = binomial
        )
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("demo_logistic_summary_all_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0("demo_", i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("demo_roc_curve_all_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("demo_ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("demo_pr_curve_all_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("demo_Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      })
    }
  }
  
  demo_auc_list_all_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  demo_auc_list_all_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(demo_auc_list_all_mean)
print(demo_auc_list_all_stdev)

#18-3 IS, IES, IDS + confounder
demo_auc_list_three_mean <- list()
demo_auc_list_three_stdev <- list()
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 20
  
  three_vars <- data.frame(
    Variable = character(0)
  )
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    three_vars <- rbind(
      three_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  three_variable_names <- three_vars$Variable
  
  variable_names = c(variable_names, three_variable_names)
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      three_variable_names,
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis",
      three_variable_names,
      "RNR_binary"
    )
  }
  
  set.seed(111)
  table_rnr <- table(data_check_encoded$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check_encoded$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check_encoded[train_index, ]
          test_data <- data_check_encoded[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(
          RNR_binary ~ .,
          data = train_data[, selected_columns],
          family = binomial
        )
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("demo_logistic_summary_three_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0("demo_", i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("demo_roc_curve_three_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("demo_ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("demo_pr_curve_three_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("demo_Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      })
    }
  }
  
  demo_auc_list_three_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  demo_auc_list_three_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(demo_auc_list_three_mean)
print(demo_auc_list_three_stdev)

#19 교란변수 보정 방법
#19-1 IS, IES만 - merged_data_tot에 대해서
auc_list_two_mean <- list()
auc_list_two_stdev <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 19
  
  two_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    two_vars <- rbind(
      two_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  two_variable_names <- two_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(two_variable_names, collapse = " + "))
  )
  print(formula)
  
  set.seed(123)
  # RNR_binary 데이터 확인
  table_rnr <- table(data_check$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check[train_index, ]
          test_data <- data_check[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(formula, data = train_data, family = binomial)
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("logistic_summary_two_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0(i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("roc_curve_two_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("pr_curve_two_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      }) 
    }
    
  }
  
  auc_list_two_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  auc_list_two_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(auc_list_two_mean)
print(auc_list_two_stdev)

#전체 데이터셋을 이용한 logistic model 만들기
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("tot_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 19
  
  two_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    two_vars <- rbind(
      two_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  two_variable_names <- two_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(two_variable_names, collapse = " + "))
  )
  print(formula)
  
  logistic_model <- glm(formula, data = data_check, family = binomial)
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  summary_name <- paste0("logistic_summary_two_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
  
  cat("ASSIGN COMPLETE\n")
}

tot_regression_result_two <- data.frame(
  Dataset = paste0("Dataset", 1:16),
  IS = NA,                           
  IES = NA                           
)

# 반복문을 통해 각 데이터셋에서 p-value 추출
for (i in 1:16) {
  # 데이터셋 이름 생성
  dataset_name <- paste0("logistic_summary_two_", i)
  
  # 데이터셋 가져오기
  logistic_summary <- get(dataset_name, envir = .GlobalEnv)
  
  # IS와 IES의 p-value 추출
  tot_regression_result_two$IS[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IS"]
  tot_regression_result_two$IES[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IES"]
}

# 결과 출력
print(tot_regression_result_two)


#19-2 IS, IES만 - merged_data_final에 대해서
demo_auc_list_onlytwo_mean <- list()
demo_auc_list_onlytwo_stdev <- list()
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 19
  
  two_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    two_vars <- rbind(
      two_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  two_variable_names <- two_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(two_variable_names, collapse = " + "))
  )
  print(formula)
  
  set.seed(123)
  # RNR_binary 데이터 확인
  table_rnr <- table(data_check$RNR_binary)
  if (any(table_rnr <= 1)) {
    cat("Error: RNR_binary contains fewer than 2 observations for one or more categories. Analysis cannot proceed.\n")
    next
  }
  
  num_iterations <- 10
  auc_values <- numeric(num_iterations)
  auprc_values <- numeric(num_iterations)
  F1_values <- numeric(num_iterations)
  
  
  for(j in 1:num_iterations) {
    repeat{
      tryCatch({
        repeat {
          # 데이터 0.8 / 0.2로 분할
          train_index <- createDataPartition(data_check$RNR_binary, p = 0.8, list = FALSE)
          train_data <- data_check[train_index, ]
          test_data <- data_check[-train_index, ]
          
          # RNR_binary가 0과 1을 모두 포함하는지 확인
          if (length(unique(train_data$RNR_binary)) == 2 && length(unique(test_data$RNR_binary)) == 2) {
            break
          }
        }
        
        logistic_model <- glm(formula, data = train_data, family = binomial)
        summary_obj <- summary(logistic_model)
        summary_df <- as.data.frame(summary_obj$coefficients)
        summary_df$Variable <- rownames(summary_df)
        rownames(summary_df) <- NULL
        
        summary_name <- paste0("demo_logistic_summary_onlytwo_", i, "_iteration_", j)
        assign(summary_name, summary_df, envir = .GlobalEnv)
        
        predicted_probs <- predict(logistic_model, newdata = test_data, type = "response")
        predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
        
        roc_curve <- roc(test_data$RNR_binary, predicted_probs)
        auc_values[j] <- auc(roc_curve)
        
        pr <- pr.curve(
          scores.class0 = predicted_probs[test_data$RNR_binary == 1],
          scores.class1 = predicted_probs[test_data$RNR_binary == 0],
          curve = TRUE
        )
        
        auprc_values[j] <- pr$auc.integral
        
        
        confusion_matrix <- table(Predicted = predicted_classes, Actual = test_data$RNR_binary)
        
        # 2x2 형태로 보장
        all_classes <- c(0, 1)  # 포함해야 하는 모든 클래스
        confusion_matrix_full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(Predicted = all_classes, Actual = all_classes))
        
        # 원본 confusion_matrix 값 복사
        for (pred in rownames(confusion_matrix)) {
          for (act in colnames(confusion_matrix)) {
            confusion_matrix_full[pred, act] <- confusion_matrix[pred, act]
          }
        }
        
        print(confusion_matrix_full)
        
        precision <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[2, ])
        recall <- confusion_matrix_full[2, 2] / sum(confusion_matrix_full[, 2])
        F1_values[j] <- 2*(precision*recall)/ (precision + recall)
        
        output_folder <- paste0("demo_", i, "_plots_", data_check$cancers[1], " - ", data_check$drug.name[1])
        if(!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        plot_path <- file.path(output_folder, paste0("demo_roc_curve_onlytwo_", i, "_iteration_", j, ".png"))
        png(filename = plot_path, width = 800, height = 600)
        plot.roc(roc_curve, main = paste("ROC Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUROC =", round(auc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        pr_plot_path <- file.path(output_folder, paste0("demo_pr_curve_onlytwo_", i, "_iteration_", j, ".png"))
        png(filename = pr_plot_path, width = 800, height = 600)
        plot(pr, main = paste("Precision-Recall Curve - Dataset", i, "_iteration_", j), col = "blue")
        text(
          x = 0.6, y = 0.2,
          labels = paste("AUPRC =", round(auprc_values[j], 3)),
          col = "red", cex = 1.5
        )
        dev.off()
        
        cat("Iteration", j, "AUROC:", auc_values[j], "AUPRC:", auprc_values[j], "F1 Score:", F1_values[j], "\n\n")
        
        break
      }, error = function(e) {
        cat("Error during PR Curve calculation for iteration", j, ":", conditionMessage(e), "\n")
        cat("Retrying iteration", j, "...\n")
      }) 
    }
    
  }
  
  demo_auc_list_onlytwo_mean[[paste0("Dataset_", i)]] <- list(
    AUROC = mean(auc_values),
    AUPRC = mean(auprc_values),
    F1_score = mean(F1_values, na.rm = TRUE)
  )
  
  demo_auc_list_onlytwo_stdev[[paste0("Dataset_", i)]] <- list(
    AUROC = sd(auc_values),
    AUPRC = sd(auprc_values),
    F1_score = sd(F1_values, na.rm = TRUE)
  )
}
print(demo_auc_list_onlytwo_mean)
print(demo_auc_list_onlytwo_stdev)

#전체 데이터셋을 이용한 logistic model 만들기
for (i in 1:16) {
  # 데이터셋 동적으로 가져오기
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 19
  
  two_vars <- data.frame(
    Variable = character(0)
  )
  
  cat("Test Result for Dataset", i, ":\n")
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    two_vars <- rbind(
      two_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  two_variable_names <- two_vars$Variable
  
  
  formula <- as.formula(
    paste("RNR_binary ~", paste(two_variable_names, collapse = " + "))
  )
  print(formula)
  
  logistic_model <- glm(formula, data = data_check, family = binomial)
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  summary_name <- paste0("demo_logistic_summary_onlytwo_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
  
  cat("ASSIGN COMPLETE\n")
}

demo_regression_result_onlytwo <- data.frame(
  Dataset = paste0("Dataset", 1:16),
  IS = NA,                           
  IES = NA                           
)

# 반복문을 통해 각 데이터셋에서 p-value 추출
for (i in 1:16) {
  # 데이터셋 이름 생성
  dataset_name <- paste0("demo_logistic_summary_onlytwo_", i)
  
  # 데이터셋 가져오기
  logistic_summary <- get(dataset_name, envir = .GlobalEnv)
  
  # IS와 IES의 p-value 추출
  demo_regression_result_onlytwo$IS[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IS"]
  demo_regression_result_onlytwo$IES[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IES"]
}

print(demo_regression_result_onlytwo)

#19-3 IS + race, age, gender - merged_data_final에 대해서
#전체 데이터셋을 이용한 logistic model 만들기
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  # Column 범위 지정
  start_col <- 18
  end_col <- 18
  
  one_vars <- data.frame(
    Variable = character(0)
  )
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    one_vars <- rbind(
      one_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  one_variable_names <- one_vars$Variable
  
  variable_names = c(variable_names, one_variable_names)
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      one_variable_names,
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis",
      one_variable_names,
      "RNR_binary"
    )
  }
  
  logistic_model <- glm(
    RNR_binary ~ .,
    data = data_check_encoded[, selected_columns],
    family = binomial
  )
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  summary_name <- paste0("demo_logistic_summary_IS_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
}

demo_regression_result_IS <- data.frame(
  Dataset = paste0("Dataset", 1:16),
  IS = NA                         
)

# 반복문을 통해 각 데이터셋에서 p-value 추출
for (i in 1:16) {
  # 데이터셋 이름 생성
  dataset_name <- paste0("demo_logistic_summary_IS_", i)
  
  # 데이터셋 가져오기
  logistic_summary <- get(dataset_name, envir = .GlobalEnv)
  
  # IS와 IES의 p-value 추출
  demo_regression_result_IS$IS[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IS"]
}

# 결과 출력
print(demo_regression_result_IS)

#19-4 IES + race, age, gender - merged_data_final에 대해서
for (i in 1:16) {
  cat("Test Result for Dataset", i, ":\n")
  
  dataset_name <- paste0("demo_data_check_", i)
  data_check <- get(dataset_name, envir = .GlobalEnv)
  
  valid_factors <- c()  # 유효한 범주형 변수 이름을 저장할 리스트
  
  if ("gender" %in% names(data_check)) {
    data_check$gender <- factor(data_check$gender)
    if (nlevels(data_check$gender) >= 2) {
      data_check$gender <- droplevels(data_check$gender)
      valid_factors <- c(valid_factors, "gender")
    }
  }
  
  if ("race" %in% names(data_check)) {
    data_check$race <- factor(data_check$race)
    if (nlevels(data_check$race) >= 2) {
      data_check$race <- droplevels(data_check$race)
      valid_factors <- c(valid_factors, "race")
    }
  }
  
  variable_names <- c(valid_factors, "age_at_initial_pathologic_diagnosis")
  
  # Column 범위 지정
  start_col <- 19
  end_col <- 19
  
  one_vars <- data.frame(
    Variable = character(0)
  )
  
  for (col in start_col:end_col) {
    # 현재 열 이름 가져오기
    variable_name <- colnames(data_check)[col]
    
    # 변수 변환 시도 (character -> numeric)
    if (!is.numeric(data_check[[variable_name]])) {
      #cat("Attempting to convert variable", variable_name, "to numeric...\n")
      # 숫자로 변환
      data_check[[variable_name]] <- suppressWarnings(as.numeric(data_check[[variable_name]]))
      
      # 변환 결과 확인
      if (all(is.na(data_check[[variable_name]]))) {
        cat("Warning: Variable", variable_name, "could not be converted to numeric. Skipping...\n\n")
        next  # 변환 실패 시 다음 변수로 넘어감
      }
    }
    
    one_vars <- rbind(
      one_vars,
      data.frame(
        Variable = variable_name
      )      
    )
  }
  
  one_variable_names <- one_vars$Variable
  
  variable_names = c(variable_names, one_variable_names)
  formula <- as.formula(
    paste("RNR_binary ~", paste(variable_names, collapse = " + "))
  )
  print(formula)
  
  if(length(valid_factors) > 0) {
    encoding_formula <- paste(" ~ ", paste(valid_factors, collapse = " + "))
    encoded_vars <- dummyVars(encoding_formula, data = data_check, fullRank = TRUE)
    one_hot_encoded <- as.data.frame(predict(encoded_vars, newdata = data_check))
    
    data_check_encoded <- cbind(data_check, one_hot_encoded)
    selected_columns <- c(
      names(one_hot_encoded), # One-Hot Encoding된 열
      "age_at_initial_pathologic_diagnosis", 
      one_variable_names,
      "RNR_binary"
    )
    
  } else {
    data_check_encoded <- data_check
    selected_columns <- c(
      "age_at_initial_pathologic_diagnosis",
      one_variable_names,
      "RNR_binary"
    )
  }
  
  logistic_model <- glm(
    RNR_binary ~ .,
    data = data_check_encoded[, selected_columns],
    family = binomial
  )
  
  summary_obj <- summary(logistic_model)
  summary_df <- as.data.frame(summary_obj$coefficients)
  summary_df$Variable <- rownames(summary_df)
  rownames(summary_df) <- NULL
  
  summary_name <- paste0("demo_logistic_summary_IES_", i)
  assign(summary_name, summary_df, envir = .GlobalEnv)
}

demo_regression_result_IES <- data.frame(
  Dataset = paste0("Dataset", 1:16),
  IES = NA                         
)

# 반복문을 통해 각 데이터셋에서 p-value 추출
for (i in 1:16) {
  # 데이터셋 이름 생성
  dataset_name <- paste0("demo_logistic_summary_IES_", i)
  
  # 데이터셋 가져오기
  logistic_summary <- get(dataset_name, envir = .GlobalEnv)
  
  # IS와 IES의 p-value 추출
  demo_regression_result_IES$IES[i] <- logistic_summary$'Pr(>|z|)'[logistic_summary$Variable == "IES"]
}

# 결과 출력
print(demo_regression_result_IES)