# =============================================================================
# Quadratic PRIMED vs Linear PRIMED: Application 비교
# Dispersion model: S_ik = α0k + α1k * X_bar_ik + α2k * X_bar_ik^2 + ε_ik
# =============================================================================

library(tidyverse)
library(readxl)
library(grpreg)
library(pROC)
library(glmnet)

# =============================================================================
# 1. 데이터 로딩 및 전처리 (기존과 동일)
# =============================================================================

data <- read_excel("../../data/data.xlsx")

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

var_labels <- c("Pre_cT" = "E", "Pre_DWI" = "F", "Pre_ADC" = "G",
                "Post_cT" = "H", "Fibrosis_grading" = "B",
                "Post_DWI" = "I", "Post_ADC" = "J",
                "Post_restriction" = "K", "MERCURY" = "L", "ESGAR" = "A")
bin_labels <- c("Pre_EMVI" = "D", "normalized_wall" = "M", "ulcer" = "N", "split_scar" = "C")

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

train_data <- cd_data %>% filter(Center == 1)
test_data  <- cd_data %>% filter(Center %in% c(2, 3))

cat(sprintf("Training: n=%d, pCR=%.1f%%\n", nrow(train_data), 100*mean(train_data$pCR)))
cat(sprintf("External: n=%d, pCR=%.1f%%\n",  nrow(test_data),  100*mean(test_data$pCR)))

# =============================================================================
# 2. Linear PRIMED (기존 - 비교 기준)
# =============================================================================

linear_primed_fit <- function(X_bar, S, Y, X_bin = NULL, lambda,
                              eta = 0.005, max_iter = 10000, tol = 1e-6) {
  n <- length(Y); K <- ncol(X_bar)
  J <- if (!is.null(X_bin)) ncol(X_bin) else 0

  L_disp_null <- sum(apply(S, 2, var)) / 2.0
  y_bar <- max(min(mean(Y), 1 - 1e-10), 1e-10)
  L_out_null <- -(y_bar * log(y_bar) + (1 - y_bar) * log(1 - y_bar))

  alpha <- rep(0, K); alpha0 <- colMeans(S)
  beta <- rep(0, K); gamma <- rep(0, K)
  delta <- if (J > 0) rep(0, J) else NULL
  beta0 <- log(y_bar / (1 - y_bar))
  sqrt2 <- sqrt(2)

  for (iter in 1:max_iter) {
    old_params <- c(alpha, beta, gamma, delta)

    # Dispersion model: S = alpha0 + alpha * X_bar
    resid_S <- matrix(0, n, K)
    for (k in 1:K) {
      resid <- S[, k] - alpha0[k] - alpha[k] * X_bar[, k]
      resid_S[, k] <- resid
      alpha0[k] <- alpha0[k] + eta * mean(resid) / L_disp_null
      alpha[k]  <- alpha[k]  + eta * mean(resid * X_bar[, k]) / L_disp_null
    }

    # Outcome model: logit(p) = beta0 + X_bar*beta + resid_S*gamma
    lp <- beta0 + X_bar %*% beta + resid_S %*% gamma
    if (J > 0) lp <- lp + X_bin %*% delta
    p <- 1 / (1 + exp(-lp)); r <- as.vector(p - Y)

    # Coupling gradient: outcome loss -> alpha
    for (k in 1:K) {
      coupling <- gamma[k] * mean(r * X_bar[, k]) / L_out_null
      alpha[k] <- alpha[k] + eta * coupling  # note: +coupling because d(resid)/d(alpha) = -X_bar
    }

    beta0 <- beta0 - eta * mean(r) / L_out_null
    grad_beta  <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_gamma <- as.vector(t(resid_S) %*% r) / (n * L_out_null)
    if (J > 0) grad_delta <- as.vector(t(X_bin) %*% r) / (n * L_out_null)

    # Group proximal step for (beta, gamma)
    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k] - eta * grad_beta[k]
      gt <- gamma[k] - eta * grad_gamma[k]
      nt <- sqrt(bt^2 + gt^2)
      if (nt <= threshold) { beta[k] <- 0; gamma[k] <- 0
      } else { sh <- 1 - threshold / nt; beta[k] <- sh * bt; gamma[k] <- sh * gt }
    }

    # Binary LASSO
    if (J > 0) {
      threshold_l1 <- eta * lambda
      for (j in 1:J) {
        dt <- delta[j] - eta * grad_delta[j]
        delta[j] <- sign(dt) * max(abs(dt) - threshold_l1, 0)
      }
    }

    new_params <- c(alpha, beta, gamma, delta)
    if (sqrt(sum((new_params - old_params)^2)) / (sqrt(sum(old_params^2)) + 1e-10) < tol) break
  }

  list(alpha0 = alpha0, alpha1 = alpha, beta = beta, gamma = gamma,
       delta = delta, beta0 = beta0, iter = iter, converged = iter < max_iter)
}

linear_primed_predict <- function(fit, X_bar, S, X_bin = NULL) {
  K <- ncol(X_bar); n <- nrow(X_bar)
  resid_S <- matrix(0, n, K)
  for (k in 1:K) resid_S[, k] <- S[, k] - fit$alpha0[k] - fit$alpha1[k] * X_bar[, k]
  lp <- fit$beta0 + X_bar %*% fit$beta + resid_S %*% fit$gamma
  if (!is.null(X_bin) && !is.null(fit$delta)) lp <- lp + X_bin %*% fit$delta
  1 / (1 + exp(-lp))
}

# =============================================================================
# 3. Quadratic PRIMED (새로운 확장)
# Dispersion: S_ik = α0k + α1k * X_bar_ik + α2k * X_bar_ik^2 + ε_ik
# =============================================================================

quadratic_primed_fit <- function(X_bar, S, Y, X_bin = NULL, lambda,
                                  eta = 0.005, max_iter = 10000, tol = 1e-6) {
  n <- length(Y); K <- ncol(X_bar)
  J <- if (!is.null(X_bin)) ncol(X_bin) else 0
  X_bar_sq <- X_bar^2  # precompute

  L_disp_null <- sum(apply(S, 2, var)) / 2.0
  y_bar <- max(min(mean(Y), 1 - 1e-10), 1e-10)
  L_out_null <- -(y_bar * log(y_bar) + (1 - y_bar) * log(1 - y_bar))

  alpha0 <- colMeans(S)
  alpha1 <- rep(0, K)  # linear term
  alpha2 <- rep(0, K)  # quadratic term
  beta <- rep(0, K); gamma <- rep(0, K)
  delta <- if (J > 0) rep(0, J) else NULL
  beta0 <- log(y_bar / (1 - y_bar))
  sqrt2 <- sqrt(2)

  for (iter in 1:max_iter) {
    old_params <- c(alpha1, alpha2, beta, gamma, delta)

    # Dispersion model: S = alpha0 + alpha1*X + alpha2*X^2
    resid_S <- matrix(0, n, K)
    for (k in 1:K) {
      resid <- S[, k] - alpha0[k] - alpha1[k] * X_bar[, k] - alpha2[k] * X_bar_sq[, k]
      resid_S[, k] <- resid
      alpha0[k] <- alpha0[k] + eta * mean(resid) / L_disp_null
      alpha1[k] <- alpha1[k] + eta * mean(resid * X_bar[, k]) / L_disp_null
      alpha2[k] <- alpha2[k] + eta * mean(resid * X_bar_sq[, k]) / L_disp_null
    }

    # Outcome model
    lp <- beta0 + X_bar %*% beta + resid_S %*% gamma
    if (J > 0) lp <- lp + X_bin %*% delta
    p <- 1 / (1 + exp(-lp)); r <- as.vector(p - Y)

    # Coupling gradients: outcome loss -> alpha1, alpha2
    for (k in 1:K) {
      coupling1 <- gamma[k] * mean(r * X_bar[, k]) / L_out_null
      coupling2 <- gamma[k] * mean(r * X_bar_sq[, k]) / L_out_null
      alpha1[k] <- alpha1[k] + eta * coupling1
      alpha2[k] <- alpha2[k] + eta * coupling2
    }

    beta0 <- beta0 - eta * mean(r) / L_out_null
    grad_beta  <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_gamma <- as.vector(t(resid_S) %*% r) / (n * L_out_null)
    if (J > 0) grad_delta <- as.vector(t(X_bin) %*% r) / (n * L_out_null)

    # Group proximal step for (beta, gamma)
    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k] - eta * grad_beta[k]
      gt <- gamma[k] - eta * grad_gamma[k]
      nt <- sqrt(bt^2 + gt^2)
      if (nt <= threshold) { beta[k] <- 0; gamma[k] <- 0
      } else { sh <- 1 - threshold / nt; beta[k] <- sh * bt; gamma[k] <- sh * gt }
    }

    # Binary LASSO
    if (J > 0) {
      threshold_l1 <- eta * lambda
      for (j in 1:J) {
        dt <- delta[j] - eta * grad_delta[j]
        delta[j] <- sign(dt) * max(abs(dt) - threshold_l1, 0)
      }
    }

    new_params <- c(alpha1, alpha2, beta, gamma, delta)
    if (sqrt(sum((new_params - old_params)^2)) / (sqrt(sum(old_params^2)) + 1e-10) < tol) break
  }

  list(alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
       beta = beta, gamma = gamma, delta = delta,
       beta0 = beta0, iter = iter, converged = iter < max_iter)
}

quadratic_primed_predict <- function(fit, X_bar, S, X_bin = NULL) {
  K <- ncol(X_bar); n <- nrow(X_bar)
  resid_S <- matrix(0, n, K)
  for (k in 1:K) resid_S[, k] <- S[, k] - fit$alpha0[k] - fit$alpha1[k] * X_bar[, k] - fit$alpha2[k] * X_bar[, k]^2
  lp <- fit$beta0 + X_bar %*% fit$beta + resid_S %*% fit$gamma
  if (!is.null(X_bin) && !is.null(fit$delta)) lp <- lp + X_bin %*% fit$delta
  1 / (1 + exp(-lp))
}

# =============================================================================
# 4. CV 함수 (linear / quadratic 공통)
# =============================================================================

primed_cv <- function(fit_fn, predict_fn, X_bar, S, Y, X_bin = NULL,
                      lambda_grid = seq(0.01, 3.0, length.out = 100),
                      n_folds = 5, seed = 42) {
  set.seed(seed); n <- length(Y)
  idx1 <- which(Y == 1); idx0 <- which(Y == 0)
  fold_ids <- rep(0, n)
  fold_ids[idx1] <- rep(1:n_folds, length.out = length(idx1))[sample(length(idx1))]
  fold_ids[idx0] <- rep(1:n_folds, length.out = length(idx0))[sample(length(idx0))]

  cv_dev <- rep(0, length(lambda_grid))
  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]; total_dev <- 0
    for (fold in 1:n_folds) {
      te <- fold_ids == fold; tr <- !te
      fit <- fit_fn(X_bar[tr,,drop=F], S[tr,,drop=F], Y[tr],
                    if (!is.null(X_bin)) X_bin[tr,,drop=F] else NULL, lam)
      pred <- as.vector(predict_fn(fit, X_bar[te,,drop=F], S[te,,drop=F],
                                   if (!is.null(X_bin)) X_bin[te,,drop=F] else NULL))
      pred <- pmax(pmin(pred, 1 - 1e-10), 1e-10)
      total_dev <- total_dev - 2 * mean(Y[te] * log(pred) + (1 - Y[te]) * log(1 - pred))
    }
    cv_dev[j] <- total_dev / n_folds
  }
  best_idx <- which.min(cv_dev)
  best_fit <- fit_fn(X_bar, S, Y, X_bin, lambda_grid[best_idx])
  list(fit = best_fit, lambda = lambda_grid[best_idx], cv_deviance = cv_dev)
}

# =============================================================================
# 5. 데이터 준비
# =============================================================================

ordinal_cons_vars <- paste0(ordinal_vars, "_cons")
ordinal_disp_vars <- paste0(ordinal_vars, "_disp")
binary_cons_vars  <- paste0(binary_vars, "_cons")

X_bar_train <- as.matrix(train_data[, ordinal_cons_vars])
S_train     <- as.matrix(train_data[, ordinal_disp_vars])
X_bin_train <- as.matrix(train_data[, binary_cons_vars])
y_train     <- train_data$pCR

X_bar_test <- as.matrix(test_data[, ordinal_cons_vars])
S_test     <- as.matrix(test_data[, ordinal_disp_vars])
X_bin_test <- as.matrix(test_data[, binary_cons_vars])
y_test     <- test_data$pCR

lambda_grid <- seq(0.01, 3.0, length.out = 30)

# =============================================================================
# 6. 두 모델 Fitting
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("  LINEAR PRIMED 적용 중...\n")
cat(strrep("=", 70), "\n")
res_lin <- primed_cv(linear_primed_fit, linear_primed_predict,
                     X_bar_train, S_train, y_train, X_bin_train,
                     lambda_grid, n_folds = 5, seed = 42)
cat(sprintf("Best lambda: %.3f, Converged: %s (iter=%d)\n",
            res_lin$lambda, res_lin$fit$converged, res_lin$fit$iter))

cat("\n", strrep("=", 70), "\n")
cat("  QUADRATIC PRIMED 적용 중...\n")
cat(strrep("=", 70), "\n")
res_quad <- primed_cv(quadratic_primed_fit, quadratic_primed_predict,
                      X_bar_train, S_train, y_train, X_bin_train,
                      lambda_grid, n_folds = 5, seed = 42)
cat(sprintf("Best lambda: %.3f, Converged: %s (iter=%d)\n",
            res_quad$lambda, res_quad$fit$converged, res_quad$fit$iter))

# =============================================================================
# 7. Dispersion model 비교: R² (linear vs quadratic)
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("  DISPERSION MODEL 비교: Linear vs Quadratic R²\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-20s  %8s  %8s  %8s  %10s\n",
            "Feature", "R²(lin)", "R²(quad)", "ΔR²", "α2 (quad)"))
cat(strrep("-", 70), "\n")

# Use full data for dispersion comparison
X_bar_all <- as.matrix(cd_data[, ordinal_cons_vars])
S_all     <- as.matrix(cd_data[, ordinal_disp_vars])

for (k in seq_along(ordinal_vars)) {
  x <- X_bar_all[, k]
  s <- S_all[, k]

  # Linear
  fit_lin <- lm(s ~ x)
  r2_lin <- summary(fit_lin)$r.squared

  # Quadratic
  fit_quad <- lm(s ~ x + I(x^2))
  r2_quad <- summary(fit_quad)$r.squared
  alpha2_coef <- coef(fit_quad)[3]
  p_quad <- summary(fit_quad)$coefficients[3, 4]

  cat(sprintf("%-6s %-13s  %8.3f  %8.3f  %+8.3f  %+10.4f %s\n",
              var_labels[ordinal_vars[k]], ordinal_vars[k],
              r2_lin, r2_quad, r2_quad - r2_lin, alpha2_coef,
              ifelse(p_quad < 0.001, "***",
                     ifelse(p_quad < 0.01, "**",
                            ifelse(p_quad < 0.05, "*", "")))))
}

# =============================================================================
# 8. Variable Selection 비교
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("  VARIABLE SELECTION 비교: Linear vs Quadratic PRIMED\n")
cat(strrep("=", 70), "\n\n")

fit_l <- res_lin$fit
fit_q <- res_quad$fit

result_comparison <- tibble(
  Label = var_labels[ordinal_vars],
  Feature = ordinal_vars,
  # Linear
  lin_alpha1 = round(fit_l$alpha1, 3),
  lin_beta   = round(fit_l$beta, 3),
  lin_gamma  = round(fit_l$gamma, 3),
  lin_sel    = abs(fit_l$beta) > 1e-6 | abs(fit_l$gamma) > 1e-6,
  # Quadratic
  quad_alpha1 = round(fit_q$alpha1, 3),
  quad_alpha2 = round(fit_q$alpha2, 3),
  quad_beta   = round(fit_q$beta, 3),
  quad_gamma  = round(fit_q$gamma, 3),
  quad_sel    = abs(fit_q$beta) > 1e-6 | abs(fit_q$gamma) > 1e-6,
  # Change
  change = case_when(
    lin_sel & quad_sel   ~ "Both",
    !lin_sel & quad_sel  ~ "NEW in Quad",
    lin_sel & !quad_sel  ~ "LOST in Quad",
    TRUE ~ "-"
  )
)

cat("[순서형 변수]\n")
print(result_comparison, n = Inf, width = 130)

# Binary
if (!is.null(fit_l$delta)) {
  bin_comp <- tibble(
    Label = bin_labels[binary_vars],
    Feature = binary_vars,
    lin_delta  = round(fit_l$delta, 3),
    lin_sel    = abs(fit_l$delta) > 1e-6,
    quad_delta = round(fit_q$delta, 3),
    quad_sel   = abs(fit_q$delta) > 1e-6
  )
  cat("\n[이진형 변수]\n")
  print(bin_comp, n = Inf)
}

# =============================================================================
# 9. 예측 성능 비교
# =============================================================================

pred_lin_train  <- as.vector(linear_primed_predict(fit_l, X_bar_train, S_train, X_bin_train))
pred_lin_test   <- as.vector(linear_primed_predict(fit_l, X_bar_test, S_test, X_bin_test))
pred_quad_train <- as.vector(quadratic_primed_predict(fit_q, X_bar_train, S_train, X_bin_train))
pred_quad_test  <- as.vector(quadratic_primed_predict(fit_q, X_bar_test, S_test, X_bin_test))

auc_lin_train  <- auc(roc(y_train, pred_lin_train, quiet = TRUE))
auc_lin_test   <- auc(roc(y_test, pred_lin_test, quiet = TRUE))
auc_quad_train <- auc(roc(y_train, pred_quad_train, quiet = TRUE))
auc_quad_test  <- auc(roc(y_test, pred_quad_test, quiet = TRUE))

# CD Group LASSO (비교 기준)
all_vars <- c(binary_cons_vars, as.vector(rbind(ordinal_cons_vars, ordinal_disp_vars)))
groups   <- c(1:4, rep(5:14, each = 2))
X_train_gl <- as.matrix(train_data[, all_vars])
X_test_gl  <- as.matrix(test_data[, all_vars])
set.seed(42)
cv_gl <- cv.grpreg(X_train_gl, y_train, group = groups, family = "binomial",
                   penalty = "grLasso", nfolds = 10)
auc_gl_train <- auc(roc(y_train, as.vector(predict(cv_gl, X_train_gl, type="response", lambda=cv_gl$lambda.min)), quiet=TRUE))
auc_gl_test  <- auc(roc(y_test, as.vector(predict(cv_gl, X_test_gl, type="response", lambda=cv_gl$lambda.min)), quiet=TRUE))

# Consensus LASSO
cons_vars <- c(binary_cons_vars, ordinal_cons_vars)
X_train_cl <- as.matrix(train_data[, cons_vars])
X_test_cl  <- as.matrix(test_data[, cons_vars])
set.seed(42)
cv_cl <- cv.glmnet(X_train_cl, y_train, family = "binomial", alpha = 1, nfolds = 10)
auc_cl_train <- auc(roc(y_train, as.vector(predict(cv_cl, X_train_cl, type="response", s="lambda.min")), quiet=TRUE))
auc_cl_test  <- auc(roc(y_test, as.vector(predict(cv_cl, X_test_cl, type="response", s="lambda.min")), quiet=TRUE))

cat("\n\n", strrep("=", 70), "\n")
cat("  PREDICTION PERFORMANCE 비교\n")
cat(strrep("=", 70), "\n\n")

n_sel_lin  <- sum(result_comparison$lin_sel)  + sum(abs(fit_l$delta) > 1e-6)
n_sel_quad <- sum(result_comparison$quad_sel) + sum(abs(fit_q$delta) > 1e-6)

final <- tibble(
  Method = c("PRIMED (Linear)", "PRIMED (Quadratic)", "CD Group LASSO", "Consensus LASSO"),
  Lambda = c(res_lin$lambda, res_quad$lambda, NA, NA),
  Vars = c(n_sel_lin, n_sel_quad, NA, NA),
  Train_AUC = as.numeric(c(auc_lin_train, auc_quad_train, auc_gl_train, auc_cl_train)),
  External_AUC = as.numeric(c(auc_lin_test, auc_quad_test, auc_gl_test, auc_cl_test))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(final)

# =============================================================================
# 10. CV deviance 비교 (linear vs quadratic)
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("  CV DEVIANCE 비교\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("Linear  PRIMED: min CV deviance = %.4f (at lambda=%.3f)\n",
            min(res_lin$cv_deviance), res_lin$lambda))
cat(sprintf("Quadratic PRIMED: min CV deviance = %.4f (at lambda=%.3f)\n",
            min(res_quad$cv_deviance), res_quad$lambda))

cat("\n\n완료.\n")
