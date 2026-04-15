# =============================================================================
# PRIMED 실제 데이터 적용 (Data A: 직장암 pCR 예측)
# CD Group LASSO, Consensus LASSO와 비교
# =============================================================================

library(tidyverse)
library(readxl)
library(grpreg)
library(pROC)
library(glmnet)

# =============================================================================
# 1. 데이터 로딩 및 전처리 (main.R과 동일)
# =============================================================================

data <- read_excel("/Users/kwon-yonghan/Documents/graduate school/project/han/RectalCR_GroupPenalty/data/data.xlsx")

data_clean <- data %>%
  slice(-(1:2)) %>%
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

binary_vars <- c("Pre_EMVI", "normalized_wall", "ulcer", "split_scar")
ordinal_vars <- c("Pre_cT", "Pre_DWI", "Pre_ADC", "Post_cT", "Fibrosis_grading",
                  "Post_DWI", "Post_ADC", "Post_restriction", "MERCURY", "ESGAR")

cd_data <- data_clean %>%
  group_by(SubjectNo) %>%
  summarise(
    pCR    = first(pCR),
    Center = first(Center),
    across(all_of(binary_vars), list(cons = ~ mean(.x, na.rm = TRUE))),
    across(all_of(ordinal_vars), list(
      cons = ~ mean(.x, na.rm = TRUE),
      disp = ~ sd(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_disp"), ~ replace_na(.x, 0)))

# 데이터 분할
train_data <- cd_data %>% filter(Center == 1)
test_data  <- cd_data %>% filter(Center %in% c(2, 3))

cat(sprintf("Training: n=%d, pCR=%.1f%%\n", nrow(train_data), 100*mean(train_data$pCR)))
cat(sprintf("External: n=%d, pCR=%.1f%%\n",  nrow(test_data),  100*mean(test_data$pCR)))

# =============================================================================
# 2. PRIMED 구현 (R version)
# =============================================================================

primed_fit <- function(X_bar, S, Y, X_bin = NULL, lambda,
                       eta = 0.005, max_iter = 10000, tol = 1e-6) {
  n <- length(Y)
  K <- ncol(X_bar)   # ordinal features
  J <- if (!is.null(X_bin)) ncol(X_bin) else 0

  # Null losses (Eq. in Sec 2.3)
  L_med_null <- sum(apply(S, 2, var)) / 2.0
  y_bar <- mean(Y)
  y_bar <- max(min(y_bar, 1 - 1e-10), 1e-10)
  L_out_null <- -(y_bar * log(y_bar) + (1 - y_bar) * log(1 - y_bar))

  # Initialize
  alpha  <- rep(0, K)
  alpha0 <- colMeans(S)
  beta   <- rep(0, K)
  gamma  <- rep(0, K)
  delta  <- if (J > 0) rep(0, J) else NULL
  beta0  <- log(y_bar / (1 - y_bar))

  sqrt2 <- sqrt(2)  # group size = 2: (beta_k, gamma_k)

  for (iter in 1:max_iter) {
    old_params <- c(alpha, alpha0, beta, gamma, delta)

    # --- Residual disagreement: S_tilde = S - alpha0 - alpha * X_bar ---
    S_tilde <- S - rep(1, n) %o% alpha0 - X_bar * rep(1, n) %o% alpha

    # --- Dispersion loss gradients (w.r.t. alpha, alpha0) ---
    grad_alpha  <- -colMeans(S_tilde * X_bar) / L_med_null
    grad_alpha0 <- -colMeans(S_tilde) / L_med_null

    # --- Outcome model using residual S_tilde (paper Sec 2.3) ---
    lp <- beta0 + X_bar %*% beta + S_tilde %*% gamma
    if (J > 0) lp <- lp + X_bin %*% delta
    p <- 1 / (1 + exp(-lp))
    r <- as.vector(p - Y)

    grad_beta0 <- mean(r) / L_out_null
    grad_beta  <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_gamma <- as.vector(t(S_tilde) %*% r) / (n * L_out_null)
    if (J > 0) grad_delta <- as.vector(t(X_bin) %*% r) / (n * L_out_null)

    # --- Joint gradient: outcome loss -> alpha (through S_tilde) ---
    #   dL_out/d(alpha_k)  = -gamma_k * mean(r * X_bar_k) / L_out_null
    #   dL_out/d(alpha0_k) = -gamma_k * mean(r) / L_out_null
    grad_alpha  <- grad_alpha  - gamma * as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_alpha0 <- grad_alpha0 - gamma * mean(r) / L_out_null

    # --- Update intercepts (unpenalized, plain gradient descent) ---
    beta0  <- beta0  - eta * grad_beta0
    alpha0 <- alpha0 - eta * grad_alpha0

    # --- Update alpha by plain gradient descent (unpenalized, paper Eq.2) ---
    alpha <- alpha - eta * grad_alpha

    # --- Proximal step: group LASSO on (beta_k, gamma_k) only (paper Eq.2) ---
    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k]  - eta * grad_beta[k]
      gt <- gamma[k] - eta * grad_gamma[k]
      nt <- sqrt(bt^2 + gt^2)

      if (nt <= threshold) {
        beta[k] <- 0; gamma[k] <- 0
      } else {
        sh <- 1 - threshold / nt
        beta[k] <- sh * bt; gamma[k] <- sh * gt
      }
    }

    # --- Proximal step for binary variables (LASSO) ---
    if (J > 0) {
      threshold_l1 <- eta * lambda
      for (j in 1:J) {
        dt <- delta[j] - eta * grad_delta[j]
        delta[j] <- sign(dt) * max(abs(dt) - threshold_l1, 0)
      }
    }

    # Convergence check
    new_params <- c(alpha, alpha0, beta, gamma, delta)
    rel_change <- sqrt(sum((new_params - old_params)^2)) /
                  (sqrt(sum(old_params^2)) + 1e-10)
    if (rel_change < tol) break
  }

  list(alpha = alpha, beta = beta, gamma = gamma, delta = delta,
       alpha0 = alpha0, beta0 = beta0, iter = iter, converged = iter < max_iter)
}

primed_predict <- function(fit, X_bar, S, X_bin = NULL) {
  n <- nrow(X_bar); K <- ncol(X_bar)
  S_tilde <- S - rep(1, n) %o% fit$alpha0 - X_bar * rep(1, n) %o% fit$alpha
  lp <- fit$beta0 + X_bar %*% fit$beta + S_tilde %*% fit$gamma
  if (!is.null(X_bin) && !is.null(fit$delta)) lp <- lp + X_bin %*% fit$delta
  1 / (1 + exp(-lp))
}

# =============================================================================
# 3. CV를 통한 lambda 선택
# =============================================================================

primed_cv <- function(X_bar, S, Y, X_bin = NULL,
                      lambda_grid = seq(0.01, 3.0, length.out = 100),
                      n_folds = 5, seed = 42) {
  set.seed(seed)
  n <- length(Y)

  # Stratified folds
  idx1 <- which(Y == 1)
  idx0 <- which(Y == 0)
  fold_ids <- rep(0, n)
  fold_ids[idx1] <- rep(1:n_folds, length.out = length(idx1))[sample(length(idx1))]
  fold_ids[idx0] <- rep(1:n_folds, length.out = length(idx0))[sample(length(idx0))]

  cv_dev <- rep(0, length(lambda_grid))

  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]
    total_dev <- 0

    for (fold in 1:n_folds) {
      test_mask  <- fold_ids == fold
      train_mask <- !test_mask

      X_bar_tr <- X_bar[train_mask, , drop = FALSE]
      S_tr     <- S[train_mask, , drop = FALSE]
      Y_tr     <- Y[train_mask]
      X_bin_tr <- if (!is.null(X_bin)) X_bin[train_mask, , drop = FALSE] else NULL

      X_bar_te <- X_bar[test_mask, , drop = FALSE]
      S_te     <- S[test_mask, , drop = FALSE]
      Y_te     <- Y[test_mask]
      X_bin_te <- if (!is.null(X_bin)) X_bin[test_mask, , drop = FALSE] else NULL

      fit <- primed_fit(X_bar_tr, S_tr, Y_tr, X_bin_tr, lam)
      pred <- as.vector(primed_predict(fit, X_bar_te, S_te, X_bin_te))
      pred <- pmax(pmin(pred, 1 - 1e-10), 1e-10)

      dev <- -2 * mean(Y_te * log(pred) + (1 - Y_te) * log(1 - pred))
      total_dev <- total_dev + dev
    }
    cv_dev[j] <- total_dev / n_folds
  }

  best_idx <- which.min(cv_dev)
  best_fit <- primed_fit(X_bar, S, Y, X_bin, lambda_grid[best_idx])

  list(fit = best_fit, lambda = lambda_grid[best_idx], cv_deviance = cv_dev)
}

# =============================================================================
# 4. 데이터 준비 및 PRIMED 적용
# =============================================================================

ordinal_cons_vars <- paste0(ordinal_vars, "_cons")
ordinal_disp_vars <- paste0(ordinal_vars, "_disp")
binary_cons_vars  <- paste0(binary_vars, "_cons")

# Training data
X_bar_train <- as.matrix(train_data[, ordinal_cons_vars])
S_train     <- as.matrix(train_data[, ordinal_disp_vars])
X_bin_train <- as.matrix(train_data[, binary_cons_vars])
y_train     <- train_data$pCR

# Test data
X_bar_test <- as.matrix(test_data[, ordinal_cons_vars])
S_test     <- as.matrix(test_data[, ordinal_disp_vars])
X_bin_test <- as.matrix(test_data[, binary_cons_vars])
y_test     <- test_data$pCR

# PRIMED with CV
cat("\n=== PRIMED 적용 중... ===\n")
primed_res <- primed_cv(X_bar_train, S_train, y_train, X_bin_train,
                        lambda_grid = seq(0.01, 3.0, length.out = 100),
                        n_folds = 5, seed = 42)

cat(sprintf("Best lambda: %.3f\n", primed_res$lambda))

# =============================================================================
# 5. 결과 정리
# =============================================================================

fit <- primed_res$fit

cat("\n", strrep("=", 70), "\n")
cat("PRIMED 결과\n")
cat(strrep("=", 70), "\n\n")

# Ordinal features
primed_result <- tibble(
  feature = ordinal_vars,
  alpha   = round(fit$alpha, 4),
  beta    = round(fit$beta, 4),
  gamma   = round(fit$gamma, 4),
  L2_norm = round(sqrt(fit$alpha^2 + fit$beta^2 + fit$gamma^2), 4),
  selected = L2_norm > 1e-6,
  indirect = round(fit$alpha * fit$gamma, 4)
)

cat("[순서형 변수]\n")
print(primed_result, n = Inf)

# Binary features
if (!is.null(fit$delta)) {
  binary_result <- tibble(
    feature  = binary_vars,
    delta    = round(fit$delta, 4),
    selected = abs(fit$delta) > 1e-6
  )
  cat("\n[이진형 변수]\n")
  print(binary_result, n = Inf)
}

# 선택된 변수
cat("\n[선택된 순서형 변수]\n")
primed_result %>% filter(selected) %>% print(n = Inf)

if (!is.null(fit$delta)) {
  cat("\n[선택된 이진형 변수]\n")
  binary_result %>% filter(selected) %>% print(n = Inf)
}

# =============================================================================
# 6. 예측 성능 비교
# =============================================================================

# PRIMED
pred_primed_train <- as.vector(primed_predict(fit, X_bar_train, S_train, X_bin_train))
pred_primed_test  <- as.vector(primed_predict(fit, X_bar_test, S_test, X_bin_test))
auc_primed_train  <- auc(roc(y_train, pred_primed_train, quiet = TRUE))
auc_primed_test   <- auc(roc(y_test, pred_primed_test, quiet = TRUE))

# CD Group LASSO (재현)
all_vars <- c(binary_cons_vars, as.vector(rbind(ordinal_cons_vars, ordinal_disp_vars)))
groups   <- c(1:4, rep(5:14, each = 2))

X_train_gl <- as.matrix(train_data[, all_vars])
X_test_gl  <- as.matrix(test_data[, all_vars])

set.seed(42)
cv_gl <- cv.grpreg(X_train_gl, y_train, group = groups, family = "binomial",
                   penalty = "grLasso", nfolds = 10)
pred_gl_train <- predict(cv_gl, X_train_gl, type = "response", lambda = cv_gl$lambda.min)
pred_gl_test  <- predict(cv_gl, X_test_gl, type = "response", lambda = cv_gl$lambda.min)
auc_gl_train  <- auc(roc(y_train, as.vector(pred_gl_train), quiet = TRUE))
auc_gl_test   <- auc(roc(y_test, as.vector(pred_gl_test), quiet = TRUE))

# Consensus LASSO
cons_vars <- c(binary_cons_vars, ordinal_cons_vars)
X_train_cl <- as.matrix(train_data[, cons_vars])
X_test_cl  <- as.matrix(test_data[, cons_vars])

set.seed(42)
cv_cl <- cv.glmnet(X_train_cl, y_train, family = "binomial", alpha = 1, nfolds = 10)
pred_cl_train <- predict(cv_cl, X_train_cl, type = "response", s = "lambda.min")
pred_cl_test  <- predict(cv_cl, X_test_cl, type = "response", s = "lambda.min")
auc_cl_train  <- auc(roc(y_train, as.vector(pred_cl_train), quiet = TRUE))
auc_cl_test   <- auc(roc(y_test, as.vector(pred_cl_test), quiet = TRUE))

# 비교 테이블
cat("\n", strrep("=", 70), "\n")
cat("예측 성능 비교\n")
cat(strrep("=", 70), "\n\n")

comparison <- tibble(
  Method = c("PRIMED", "CD + Group LASSO", "Consensus LASSO"),
  Training_AUC = c(auc_primed_train, auc_gl_train, auc_cl_train),
  External_AUC = c(auc_primed_test, auc_gl_test, auc_cl_test)
) %>%
  mutate(across(where(is.numeric), ~ round(as.numeric(.x), 3)))

print(comparison)

# =============================================================================
# 7. PRIMED 계수 상세 (CD GL과 비교)
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PRIMED vs CD Group LASSO: 선택된 변수 비교\n")
cat(strrep("=", 70), "\n\n")

coefs_gl <- coef(cv_gl, s = "lambda.min")[-1]
names(coefs_gl) <- all_vars

gl_result <- tibble(
  feature   = ordinal_vars,
  gl_beta   = coefs_gl[ordinal_cons_vars],
  gl_gamma  = coefs_gl[ordinal_disp_vars],
  gl_L2     = sqrt(gl_beta^2 + gl_gamma^2),
  gl_selected = gl_L2 > 1e-6
)

combined <- primed_result %>%
  left_join(gl_result, by = "feature") %>%
  mutate(
    primed_sel = selected,
    gl_sel     = gl_selected,
    change     = case_when(
      primed_sel & !gl_sel ~ "NEW (PRIMED only)",
      !primed_sel & gl_sel ~ "DROPPED",
      primed_sel & gl_sel  ~ "BOTH",
      TRUE ~ "-"
    )
  ) %>%
  dplyr::select(feature, alpha, beta, gamma, indirect, primed_sel,
         gl_beta, gl_gamma, gl_sel, change)

print(combined, n = Inf, width = 120)

# 결과 저장
saveRDS(list(
  primed_fit = fit,
  primed_res = primed_res,
  comparison = comparison,
  combined   = combined,
  primed_result = primed_result,
  binary_result = if (!is.null(fit$delta)) binary_result else NULL
), file = "/Users/kwon-yonghan/Documents/graduate school/project/han/RectalCR_GroupPenalty/data/primed_real_data_results.rds")

cat("\n결과 저장 완료\n")
