# =============================================================================
#                    직장암 완전관해(pCR) 예측 모델
#                    Group LASSO 기반 변수 선택 분석
# =============================================================================
#
# [프로젝트 개요]
# 본 코드는 직장암 환자의 치료 후 완전관해(pCR: pathologic Complete Response)를
# 예측하기 위한 통계 모델을 구축합니다. 다중 판독자(Multi-reader) MRI 영상 데이터를
# 활용하여 판독자 간 일치도(Consensus)와 변동성(Dispersion)을 동시에 고려합니다.
#
# [분석 방법론]
# - Consensus-Dispersion (CD) 접근법: 각 환자에 대해 여러 판독자의 평가를
#   요약하여 평균(Consensus)과 표준편차(Dispersion)를 계산
# - Group LASSO: 관련 변수들을 그룹으로 묶어 변수 선택을 수행하는 정규화 기법
#   → 순서형 변수의 경우 Consensus와 Dispersion을 하나의 그룹으로 처리
#
# [변수 처리 전략]
# - 이진형 변수 (4개): Consensus만 사용 (0/1 변수의 SD는 의미 제한적)
#   예) Pre_EMVI, normalized_wall, ulcer, split_scar
# - 순서형 변수 (10개): Consensus + Dispersion 쌍으로 그룹화
#   예) Pre_cT, Post_cT, MERCURY, ESGAR 등
#
# [코드 구성]
# ┌─────────────────────────────────────────────────────────────────────────┐
# │ Part 1. Group LASSO 분석 (Lines ~166)                                   │
# │   1. 데이터 로딩 및 전처리                                              │
# │   2. 변수 분류 (이진형/순서형)                                          │
# │   3. CD 데이터 생성 (환자별 Consensus, Dispersion 계산)                 │
# │   4. Group LASSO 모델 적합 (10-fold CV)                                 │
# │   5. 변수 선택 결과 및 예측 성능(AUC) 평가                              │
# ├─────────────────────────────────────────────────────────────────────────┤
# │ Part 2. 이분산 구조 검증 (Lines ~357)                                   │
# │   1. Consensus-Dispersion 관계 시각화                                   │
# │   2. 패턴 분류 (선형/역U자형/무관)                                      │
# │   3. 이분산 구조 요약 테이블 생성                                       │
# └─────────────────────────────────────────────────────────────────────────┘
#
# [데이터 구조]
# - 입력: data.xlsx (다중 판독자 MRI 평가 데이터)
# - Training: Center 1 (내부 검증)
# - External Validation: Center 2, 3 (외부 검증)
#
# [출력물]
# - 선택된 변수 목록 및 계수
# - Training/External AUC
# - 이분산 구조 시각화 (heteroscedasticity_no_pvalue.png)
#
# =============================================================================

# =============================================================================
# Part 1. Group LASSO 분석
# =============================================================================

# --- 패키지 로딩 ---
library(tidyverse)
library(readxl)
library(grpreg)
library(pROC)

# --- 1.1 데이터 로딩 및 전처리 ---
data <- read_excel("data.xlsx")

data_clean <- data %>%
  slice(-(1:2)) %>%                              # 헤더 행 제거
  rename(
    Pre_DWI_ADC  = `DWI-ADC...12`,
    Post_DWI_ADC = `DWI-ADC...23`,
    pCR          = `...28`
  ) %>%
  mutate(
    across(c(Pre_EMVI, normalized_wall, ulcer, split_scar,
             Pre_cT, Pre_DWI, Pre_ADC, Post_cT, Fibrosis_grading,
             Post_DWI, Post_ADC, Post_restriction, MERCURY, ESGAR,
             pCR, Center), as.numeric),
    reviewer = as.numeric(reviewer)
  )

# --- 1.2 변수 분류 ---
# 이진형: 있음/없음 (0/1)
binary_vars <- c("Pre_EMVI", "normalized_wall", "ulcer", "split_scar")

# 순서형: 등급 척도 (예: 1-4점)
ordinal_vars <- c("Pre_cT", "Pre_DWI", "Pre_ADC", "Post_cT", "Fibrosis_grading",
                  "Post_DWI", "Post_ADC", "Post_restriction", "MERCURY", "ESGAR")

# --- 1.3 환자별 CD(Consensus-Dispersion) 데이터 생성 ---

cd_data <- data_clean %>%
  group_by(SubjectNo) %>%
  summarise(
    pCR    = first(pCR),
    Center = first(Center),

    # 이진형: Consensus(평균)만 계산
    across(all_of(binary_vars), list(cons = ~ mean(.x, na.rm = TRUE))),

    # 순서형: Consensus(평균) + Dispersion(표준편차) 계산
    across(all_of(ordinal_vars), list(
      cons = ~ mean(.x, na.rm = TRUE),
      disp = ~ sd(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_disp"), ~ replace_na(.x, 0)))  # 단일 판독자 → SD=0

# --- 1.4 Group LASSO 변수 및 그룹 구성 ---
binary_cons    <- paste0(binary_vars, "_cons")                        # 4개 (단독 그룹)
ordinal_cons   <- paste0(ordinal_vars, "_cons")
ordinal_disp   <- paste0(ordinal_vars, "_disp")
ordinal_paired <- as.vector(rbind(ordinal_cons, ordinal_disp))        # 20개 (10그룹 × 2)

all_vars <- c(binary_cons, ordinal_paired)                            # 총 24개 변수
groups   <- c(1:4, rep(5:14, each = 2))                               # 총 14개 그룹

# 구성 확인 출력
cat("변수 구성:\n")
cat(sprintf("  이진형: %d 변수 (각각 단독, 그룹 1-4)\n", length(binary_cons)))
cat(sprintf("  순서형: %d 변수 (%d 그룹 × 2, 그룹 5-14)\n",
            length(ordinal_paired), length(ordinal_vars)))
cat(sprintf("  총 %d 변수, %d 그룹\n\n", length(all_vars), max(groups)))

tibble(position = 1:length(all_vars), variable = all_vars, group = groups) %>%
  print(n = 30)

# --- 1.5 Training/Test 분할 및 모델 적합 ---

# 데이터 분할
train_data <- cd_data %>% filter(Center == 1)           # 내부 검증용
test_data  <- cd_data %>% filter(Center %in% c(2, 3))   # 외부 검증용

cat(sprintf("\nTraining: n=%d, pCR=%.1f%%\n", nrow(train_data), 100*mean(train_data$pCR)))
cat(sprintf("External: n=%d, pCR=%.1f%%\n",  nrow(test_data),  100*mean(test_data$pCR)))

# 설계 행렬 구성
X_train <- as.matrix(train_data[, all_vars])
y_train <- train_data$pCR

# 10-fold CV를 통한 Group LASSO 적합
set.seed(42)
cv_fit <- cv.grpreg(X_train, y_train, group = groups, family = "binomial",
                    penalty = "grLasso", nfolds = 10)

# --- 1.6 결과 추출 ---
coefs <- coef(cv_fit, s = "lambda.min")[-1]    # 절편 제외
names(coefs) <- all_vars

# 이진형 결과 (Consensus 계수만)
binary_result <- tibble(
  feature    = binary_vars,
  type       = "Binary",
  beta_cons  = coefs[binary_cons],
  gamma_disp = NA_real_,
  L2         = abs(beta_cons)
)

# 순서형 결과 (Consensus + Dispersion 계수)
ordinal_result <- tibble(
  feature    = ordinal_vars,
  type       = "Ordinal",
  beta_cons  = coefs[ordinal_cons],
  gamma_disp = coefs[ordinal_disp],
  L2         = sqrt(beta_cons^2 + gamma_disp^2)    # 그룹 L2 norm
)

# 결과 통합 및 선택 여부 판단
result <- bind_rows(binary_result, ordinal_result) %>%
  mutate(selected = L2 > 1e-6) %>%
  arrange(desc(L2))

cat("\n", strrep("=", 60), "\n선택된 변수\n", strrep("=", 60), "\n\n")
result %>% filter(selected) %>% print(n = Inf)

# --- 1.7 예측 성능 평가 (AUC) ---
pred_train <- predict(cv_fit, X_train, type = "response", lambda = cv_fit$lambda.min)
auc_train  <- auc(roc(y_train, as.vector(pred_train), quiet = TRUE))

X_test    <- as.matrix(test_data[, all_vars])
y_test    <- test_data$pCR
pred_test <- predict(cv_fit, X_test, type = "response", lambda = cv_fit$lambda.min)
auc_test  <- auc(roc(y_test, as.vector(pred_test), quiet = TRUE))

cat("\n", strrep("=", 60), "\n")
cat(sprintf("Training AUC:            %.3f\n", auc_train))
cat(sprintf("External Validation AUC: %.3f\n", auc_test))
cat(strrep("=", 60), "\n")

# --- 1.8 Leave-One-Center-Out (LOCO) 검증 ---
cat("\n", strrep("=", 60), "\n")
cat("Leave-One-Center-Out (LOCO) Validation\n")
cat(strrep("=", 60), "\n\n")

centers <- unique(cd_data$Center)
loco_results <- list()

for (test_center in centers) {
  # Train on all other centers, test on left-out center
  train_loco <- cd_data %>% filter(Center != test_center)
  test_loco  <- cd_data %>% filter(Center == test_center)

  X_train_loco <- as.matrix(train_loco[, all_vars])
  y_train_loco <- train_loco$pCR
  X_test_loco  <- as.matrix(test_loco[, all_vars])
  y_test_loco  <- test_loco$pCR

  # Fit Group LASSO
  set.seed(42)
  cv_loco <- cv.grpreg(X_train_loco, y_train_loco, group = groups,
                       family = "binomial", penalty = "grLasso", nfolds = 10)

  # Predictions
  pred_train_loco <- predict(cv_loco, X_train_loco, type = "response", lambda = cv_loco$lambda.min)
  pred_test_loco  <- predict(cv_loco, X_test_loco, type = "response", lambda = cv_loco$lambda.min)

  # AUC calculation
  auc_train_loco <- auc(roc(y_train_loco, as.vector(pred_train_loco), quiet = TRUE))
  auc_test_loco  <- auc(roc(y_test_loco, as.vector(pred_test_loco), quiet = TRUE))

  loco_results[[as.character(test_center)]] <- list(
    test_center = test_center,
    n_train = nrow(train_loco),
    n_test = nrow(test_loco),
    auc_train = auc_train_loco,
    auc_test = auc_test_loco
  )

  cat(sprintf("Test Center %d: n_train=%d, n_test=%d, Train AUC=%.3f, Test AUC=%.3f\n",
              test_center, nrow(train_loco), nrow(test_loco), auc_train_loco, auc_test_loco))
}

# LOCO Summary
loco_summary <- map_dfr(loco_results, ~ tibble(
  Center = .x$test_center,
  n_train = .x$n_train,
  n_test = .x$n_test,
  AUC_train = as.numeric(.x$auc_train),
  AUC_test = as.numeric(.x$auc_test)
))

cat("\n[LOCO Summary]\n")
print(loco_summary)

cat(sprintf("\nLOCO Test AUC Range: %.3f - %.3f\n",
            min(loco_summary$AUC_test), max(loco_summary$AUC_test)))
cat(sprintf("LOCO Test AUC Mean:  %.3f\n", mean(loco_summary$AUC_test)))



# =============================================================================
# Part 2. 이분산 구조 검증 (Heteroscedasticity Analysis)
# =============================================================================
# 목적: Consensus(평균 점수)와 Dispersion(판독자 간 변동성) 사이의 관계를 분석하여
#       CD 접근법의 통계적 타당성을 검증합니다.
#
# 패턴 유형:
#   - 양의 선형: 점수가 높을수록 판독자 간 불일치 증가
#   - 음의 선형: 점수가 높을수록 판독자 간 불일치 감소
#   - 역U자형: 중간 점수에서 불일치가 최대 (경계 효과)
#   - 무관: 점수와 불일치 간 유의한 관계 없음
# =============================================================================

# --- 2.1 추가 패키지 로딩 및 한글 폰트 설정 ---
if (!require("showtext")) install.packages("showtext")
library(patchwork)
library(showtext)

font_add("apple", "/System/Library/Fonts/Supplemental/AppleGothic.ttf")
showtext_auto()

# --- 2.2 데이터 재로딩 (독립 실행 가능하도록) ---
data <- read_excel("data.xlsx")

data_clean <- data %>%
  slice(-(1:2)) %>%
  rename(Pre_DWI_ADC = `DWI-ADC...12`, Post_DWI_ADC = `DWI-ADC...23`, pCR = `...28`) %>%
  mutate(
    across(c(Pre_EMVI, normalized_wall, ulcer, split_scar,
             Pre_cT, Pre_DWI, Pre_ADC, Post_cT, Fibrosis_grading,
             Post_DWI, Post_ADC, Post_restriction, MERCURY, ESGAR, pCR, Center),
           as.numeric),
    reviewer = as.numeric(reviewer)
  )

# --- 2.3 Training 데이터 필터링 ---
train_clean <- data_clean %>% filter(Center == 1)

# --- 2.4 환자별 Consensus/Dispersion 계산 ---
ordinal_vars <- c("Pre_cT", "Pre_DWI", "Pre_ADC", "Post_cT", "Fibrosis_grading",
                  "Post_DWI", "Post_ADC", "Post_restriction", "MERCURY", "ESGAR")

variance_check <- train_clean %>%
  group_by(SubjectNo) %>%
  summarise(
    across(all_of(ordinal_vars),
           list(mean = ~ mean(.x, na.rm = TRUE),
                sd   = ~ sd(.x, na.rm = TRUE))),
    .groups = "drop"
  )

# --- 2.5 이분산 시각화 함수 ---

plot_heteroscedasticity <- function(var) {
  df <- tibble(
    mean = variance_check[[paste0(var, "_mean")]],
    sd = variance_check[[paste0(var, "_sd")]]
  ) %>% filter(!is.na(sd))
  
  # 통계 계산 (내부 판단용)
  r <- cor(df$mean, df$sd, use = "complete.obs")
  linear_pval <- cor.test(df$mean, df$sd)$p.value
  fit_quad <- lm(sd ~ mean + I(mean^2), data = df)
  quad_pval <- summary(fit_quad)$coefficients[3, 4]
  peak <- -coef(fit_quad)[2] / (2 * coef(fit_quad)[3])
  
  linear_sig <- linear_pval < 0.05
  quad_sig <- quad_pval < 0.05
  
  # 패턴 결정
  pattern <- case_when(
    linear_sig & r > 0 ~ "양의 선형",
    linear_sig & r < 0 ~ "음의 선형",
    !linear_sig & quad_sig ~ "역U자형",
    TRUE ~ "무관"
  )
  
  # 자막에서 p-value 제거, 상관계수만 표시
  subtitle_text <- sprintf("상관계수(r) = %+.2f", r)
  
  p <- ggplot(df, aes(x = mean, y = sd)) +
    geom_point(alpha = 0.6, size = 2.5) +
    theme_minimal(base_size = 10, base_family = "apple") +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      axis.title = element_text(size = 9)
    ) +
    labs(
      title = sprintf("%s [%s]", var, pattern),
      subtitle = subtitle_text,
      x = "Consensus (Mean)", 
      y = "Within-subject SD"
    )
  
  # 패턴별 회귀선 추가
  if (linear_sig) {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1, alpha = 0.2)
  } else if (quad_sig) {
    p <- p + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "blue", linewidth = 1, alpha = 0.2) +
      geom_vline(xintercept = peak, linetype = "dashed", color = "darkblue", linewidth = 0.7) +
      annotate("text", x = peak, y = max(df$sd, na.rm = TRUE) * 0.9, 
               label = sprintf("Peak: %.1f", peak), family = "apple", size = 3, color = "darkblue")
  } else {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "gray50", linewidth = 1, alpha = 0.2)
  }
  
  return(p)
}

# --- 2.6 플롯 생성 및 저장 ---
all_plots <- map(ordinal_vars, plot_heteroscedasticity)

all_combined <- wrap_plots(all_plots, ncol = 4) +
  plot_annotation(
    title    = "이분산 구조 검증: 점수대별 일치도 분석 (Training Data, n = 50)",
    subtitle = "빨간선: 선형 경향 | 파란선: 역U자형 경향 | 회색선: 뚜렷한 경향성 없음",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", family = "apple"),
      plot.subtitle = element_text(size = 11, color = "gray40", family = "apple")
    )
  )

print(all_combined)

ggsave("heteroscedasticity_no_pvalue.png", all_combined,
       width = 14, height = 10, dpi = 300)

cat("\n[완료] 이분산 구조 시각화 저장: heteroscedasticity_no_pvalue.png\n")

# --- 2.7 이분산 구조 요약 테이블 ---
hetero_table <- map_dfr(ordinal_vars, function(var) {
  df <- tibble(
    mean = variance_check[[paste0(var, "_mean")]],
    sd   = variance_check[[paste0(var, "_sd")]]
  ) %>% filter(!is.na(sd))
  
  cor_test <- cor.test(df$mean, df$sd)
  fit_quad <- lm(sd ~ mean + I(mean^2), data = df)
  
  r <- cor_test$estimate
  linear_p <- cor_test$p.value
  quad_p <- summary(fit_quad)$coefficients[3, 4]
  
  pattern <- case_when(
    linear_p < 0.05 & r > 0 ~ "양의 선형",
    linear_p < 0.05 & r < 0 ~ "음의 선형",
    linear_p >= 0.05 & quad_p < 0.05 ~ "역U자형",
    TRUE ~ "무관"
  )
  
  tibble(
    Variable      = var,
    Correlation_r = round(r, 2),
    Linear_p      = signif(linear_p, 3),
    Quadratic_p   = signif(quad_p, 3),
    Pattern       = pattern
  )
})

cat("\n[이분산 구조 요약 테이블]\n")
print(hetero_table)

# =============================================================================
# End of Script
# =============================================================================
