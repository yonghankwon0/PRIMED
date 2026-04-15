# =============================================================================
# Step 4: PRIMED 적용 (High-Dimensional Radiomics)
#
# 입력: radiomics_cd.csv (03_extract_radiomics.py 출력)
# 비교: PRIMED vs CD Group LASSO vs Consensus LASSO
#
# 사용법: Rscript 04_apply_primed.R [radiomics_cd.csv 경로]
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(grpreg)
library(pROC)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) > 0) args[1] else "output/radiomics_cd.csv"
stopifnot(file.exists(csv_path))

# =============================================================================
# 1. 데이터 로딩
# =============================================================================
cd_data <- read_csv(csv_path, show_col_types = FALSE)
cat(sprintf("Loaded: %d nodules\n", nrow(cd_data)))

# Outcome 이진화
cd_data <- cd_data %>%
  filter(malignancy_mean != 3.0) %>%
  mutate(malignant = as.integer(malignancy_mean > 3))

cat(sprintf("After filtering: %d nodules (mal %.1f%%)\n",
            nrow(cd_data), 100 * mean(cd_data$malignant)))

# Feature 식별
cons_cols <- names(cd_data)[grepl("_cons$", names(cd_data))]
disp_cols <- names(cd_data)[grepl("_disp$", names(cd_data))]
feature_names <- sub("_cons$", "", cons_cols)
K <- length(feature_names)
cat(sprintf("Features: K = %d\n", K))

# NA/Inf 처리
cd_data <- cd_data %>%
  mutate(across(all_of(c(cons_cols, disp_cols)),
                ~ ifelse(is.finite(.x), .x, NA))) %>%
  mutate(across(all_of(disp_cols), ~ ifelse(is.na(.x), 0, .x)))

# 결측 feature 제거 (모든 결절에서 NA인 feature)
valid <- sapply(cons_cols, function(c) mean(is.na(cd_data[[c]])) < 0.1)
cons_cols <- cons_cols[valid]
disp_cols <- disp_cols[valid]
feature_names <- feature_names[valid]
K <- length(feature_names)
cat(sprintf("After NA filter: K = %d features\n", K))

# 남은 NA는 median imputation
for (col in c(cons_cols, disp_cols)) {
  med <- median(cd_data[[col]], na.rm = TRUE)
  cd_data[[col]][is.na(cd_data[[col]])] <- med
}

# =============================================================================
# 2. Train/Test Split
# =============================================================================
set.seed(2024)
train_idx <- cd_data %>%
  mutate(row_id = row_number()) %>%
  group_by(malignant) %>%
  slice_sample(prop = 0.7) %>%
  pull(row_id)

train_data <- cd_data[train_idx, ]
test_data  <- cd_data[-train_idx, ]

X_bar_train <- as.matrix(scale(train_data[, cons_cols]))
S_train     <- as.matrix(scale(train_data[, disp_cols]))
y_train     <- train_data$malignant

# 같은 scaling을 test에 적용
train_center_X <- attr(X_bar_train, "scaled:center")
train_scale_X  <- attr(X_bar_train, "scaled:scale")
train_center_S <- attr(S_train, "scaled:center")
train_scale_S  <- attr(S_train, "scaled:scale")

X_bar_test <- scale(test_data[, cons_cols], center = train_center_X, scale = train_scale_X)
S_test     <- scale(test_data[, disp_cols], center = train_center_S, scale = train_scale_S)
y_test     <- test_data$malignant

cat(sprintf("\nTrain: n=%d, Test: n=%d\n", nrow(train_data), nrow(test_data)))

# =============================================================================
# 3. PRIMED (Quadratic Dispersion, Joint Estimation)
# =============================================================================
primed_fit <- function(X_bar, S, Y, lambda,
                       eta = 0.001, max_iter = 15000, tol = 1e-6) {
  n <- length(Y); K <- ncol(X_bar)

  L_med_null <- sum(apply(S, 2, var)) / 2.0
  y_bar <- max(min(mean(Y), 1 - 1e-10), 1e-10)
  L_out_null <- -(y_bar * log(y_bar) + (1 - y_bar) * log(1 - y_bar))

  X_bar_c <- scale(X_bar, scale = FALSE)
  X_bar_sq <- X_bar_c^2

  alpha0 <- colMeans(S); alpha1 <- rep(0, K); alpha2 <- rep(0, K)
  beta <- rep(0, K); gamma <- rep(0, K)
  beta0 <- log(y_bar / (1 - y_bar))
  sqrt2 <- sqrt(2)

  for (iter in 1:max_iter) {
    old <- c(alpha0, alpha1, alpha2, beta, gamma)

    S_tilde <- S - rep(1, n) %o% alpha0 -
               X_bar_c * rep(1, n) %o% alpha1 -
               X_bar_sq * rep(1, n) %o% alpha2

    ga0 <- -colMeans(S_tilde) / L_med_null
    ga1 <- -colMeans(S_tilde * X_bar_c) / L_med_null
    ga2 <- -colMeans(S_tilde * X_bar_sq) / L_med_null

    lp <- beta0 + X_bar %*% beta + S_tilde %*% gamma
    p <- 1 / (1 + exp(-lp)); r <- as.vector(p - Y)

    gb0 <- mean(r) / L_out_null
    gb <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
    gg <- as.vector(t(S_tilde) %*% r) / (n * L_out_null)

    rXc <- as.vector(t(X_bar_c) %*% r) / (n * L_out_null)
    rXs <- as.vector(t(X_bar_sq) %*% r) / (n * L_out_null)
    rm  <- mean(r) / L_out_null

    ga0 <- ga0 - gamma * rm; ga1 <- ga1 - gamma * rXc; ga2 <- ga2 - gamma * rXs

    beta0  <- beta0  - eta * gb0
    alpha0 <- alpha0 - eta * ga0
    alpha1 <- alpha1 - eta * ga1
    alpha2 <- alpha2 - eta * ga2

    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k] - eta * gb[k]; gt <- gamma[k] - eta * gg[k]
      nt <- sqrt(bt^2 + gt^2)
      if (nt <= threshold) { beta[k] <- 0; gamma[k] <- 0 }
      else { sh <- 1 - threshold/nt; beta[k] <- sh*bt; gamma[k] <- sh*gt }
    }

    new <- c(alpha0, alpha1, alpha2, beta, gamma)
    rc <- sqrt(sum((new - old)^2)) / (sqrt(sum(old^2)) + 1e-10)
    if (is.na(rc) || rc < tol) break
  }

  list(alpha0=alpha0, alpha1=alpha1, alpha2=alpha2,
       beta=beta, gamma=gamma, beta0=beta0,
       X_bar_center=attr(X_bar_c, "scaled:center"),
       iter=iter, converged=iter<max_iter,
       n_selected=sum(abs(beta) > 1e-8 | abs(gamma) > 1e-8))
}

primed_predict <- function(fit, X_bar, S) {
  n <- nrow(X_bar)
  Xc <- sweep(X_bar, 2, fit$X_bar_center)
  St <- S - rep(1,n) %o% fit$alpha0 - Xc * rep(1,n) %o% fit$alpha1 -
        (Xc^2) * rep(1,n) %o% fit$alpha2
  1 / (1 + exp(-(fit$beta0 + X_bar %*% fit$beta + St %*% fit$gamma)))
}

# CV lambda selection
cat("\n=== PRIMED CV ===\n")
lambda_grid <- seq(0.05, 2.0, length.out = 20)
n_folds <- 5
set.seed(42)
idx1 <- which(y_train == 1); idx0 <- which(y_train == 0)
fold_ids <- rep(0, length(y_train))
fold_ids[idx1] <- rep(1:n_folds, length.out=length(idx1))[sample(length(idx1))]
fold_ids[idx0] <- rep(1:n_folds, length.out=length(idx0))[sample(length(idx0))]

cv_dev <- rep(0, length(lambda_grid))
for (j in seq_along(lambda_grid)) {
  cat(sprintf("  lambda %d/%d (%.3f)...\n", j, length(lambda_grid), lambda_grid[j]))
  td <- 0
  for (fold in 1:n_folds) {
    te <- fold_ids == fold; tr <- !te
    fit <- primed_fit(X_bar_train[tr,], S_train[tr,], y_train[tr], lambda_grid[j])
    pred <- as.vector(primed_predict(fit, X_bar_train[te,], S_train[te,]))
    pred <- pmax(pmin(pred, 1-1e-10), 1e-10)
    td <- td - 2*mean(y_train[te]*log(pred) + (1-y_train[te])*log(1-pred))
  }
  cv_dev[j] <- td / n_folds
}

best_lam <- lambda_grid[which.min(cv_dev)]
cat(sprintf("Best lambda: %.3f\n", best_lam))

final_fit <- primed_fit(X_bar_train, S_train, y_train, best_lam)
cat(sprintf("Selected features: %d / %d\n", final_fit$n_selected, K))

# =============================================================================
# 4. 비교 모델
# =============================================================================
pred_primed_test <- as.vector(primed_predict(final_fit, X_bar_test, S_test))
auc_primed <- as.numeric(auc(roc(y_test, pred_primed_test, quiet=TRUE)))

# CD Group LASSO
cd_vars_train <- cbind(X_bar_train, S_train)
cd_vars_test  <- cbind(X_bar_test, S_test)
groups <- rep(1:K, each=2)
set.seed(42)
cv_gl <- cv.grpreg(cd_vars_train, y_train, group=groups, family="binomial",
                   penalty="grLasso", nfolds=10)
pred_gl <- predict(cv_gl, cd_vars_test, type="response", lambda=cv_gl$lambda.min)
auc_gl <- as.numeric(auc(roc(y_test, as.vector(pred_gl), quiet=TRUE)))

# Consensus LASSO
set.seed(42)
cv_cl <- cv.glmnet(X_bar_train, y_train, family="binomial", alpha=1, nfolds=10)
pred_cl <- predict(cv_cl, X_bar_test, type="response", s="lambda.min")
auc_cl <- as.numeric(auc(roc(y_test, as.vector(pred_cl), quiet=TRUE)))

# =============================================================================
# 5. 결과 출력 및 저장
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("Performance Comparison (High-Dimensional Radiomics)\n")
cat(strrep("=", 60), "\n\n")

comparison <- tibble(
  Method = c("PRIMED (quadratic)", "CD + Group LASSO", "Consensus LASSO"),
  Test_AUC = round(c(auc_primed, auc_gl, auc_cl), 3),
  Selected = c(final_fit$n_selected,
               sum(abs(coef(cv_gl, s="lambda.min")[-1]) > 1e-8),
               sum(abs(coef(cv_cl, s="lambda.min")[-1]) > 1e-8))
)
print(comparison)

# 선택된 feature 출력
sel_idx <- which(abs(final_fit$beta) > 1e-8 | abs(final_fit$gamma) > 1e-8)
if (length(sel_idx) > 0) {
  cat("\n=== PRIMED 선택된 Features ===\n")
  sel_df <- tibble(
    feature = feature_names[sel_idx],
    beta = round(final_fit$beta[sel_idx], 4),
    gamma = round(final_fit$gamma[sel_idx], 4),
    alpha1 = round(final_fit$alpha1[sel_idx], 4),
    alpha2 = round(final_fit$alpha2[sel_idx], 4)
  ) %>% arrange(desc(abs(beta) + abs(gamma)))
  print(sel_df, n = 30)
}

# 결과 저장
out_dir <- dirname(csv_path)
saveRDS(list(fit=final_fit, comparison=comparison, selected=sel_df,
             feature_names=feature_names),
        file.path(out_dir, "primed_radiomics_results.rds"))
cat(sprintf("\n결과 저장: %s\n", file.path(out_dir, "primed_radiomics_results.rds")))
