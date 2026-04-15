# =============================================================================
# LIDC-IDRI: Bootstrap Selection Stability
#
# B=200 bootstrap resampling으로 selection frequency + Jaccard stability 계산
#
# 사용법: Rscript 02_lidc_bootstrap_stability.R
# =============================================================================

library(tidyverse)
library(grpreg)
library(pROC)
library(glmnet)

# =============================================================================
# 1. 데이터 로딩
# =============================================================================
raw <- read_csv("data/lidc_annotations.csv", show_col_types = FALSE) %>%
  filter(reader <= 4)

ordinal_vars <- c("subtlety", "internalStructure", "calcification",
                  "sphericity", "margin", "lobulation", "spiculation", "texture")

cd_data <- raw %>%
  group_by(nodule_id, patient_id) %>%
  summarise(
    n_readers = n(),
    mal_consensus = mean(malignancy, na.rm = TRUE),
    across(all_of(ordinal_vars), list(
      cons = ~ mean(.x, na.rm = TRUE),
      disp = ~ sd(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_disp"), ~ replace_na(.x, 0))) %>%
  filter(mal_consensus != 3.0) %>%
  mutate(malignant = as.integer(mal_consensus > 3))

cat(sprintf("Nodules: %d, Malignant: %.1f%%\n",
            nrow(cd_data), 100 * mean(cd_data$malignant)))

# =============================================================================
# 2. PRIMED 구현 (quadratic dispersion, joint estimation)
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

    ga0 <- ga0 - gamma * rm
    ga1 <- ga1 - gamma * rXc
    ga2 <- ga2 - gamma * rXs

    beta0  <- beta0  - eta * gb0
    alpha0 <- alpha0 - eta * ga0
    alpha1 <- alpha1 - eta * ga1
    alpha2 <- alpha2 - eta * ga2

    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k] - eta * gb[k]; gt <- gamma[k] - eta * gg[k]
      nt <- sqrt(bt^2 + gt^2)
      if (nt <= threshold) { beta[k] <- 0; gamma[k] <- 0 }
      else { sh <- 1 - threshold / nt; beta[k] <- sh * bt; gamma[k] <- sh * gt }
    }

    new <- c(alpha0, alpha1, alpha2, beta, gamma)
    rc <- sqrt(sum((new - old)^2)) / (sqrt(sum(old^2)) + 1e-10)
    if (is.na(rc) || rc < tol) break
  }

  list(alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
       beta = beta, gamma = gamma, beta0 = beta0,
       X_bar_center = attr(X_bar_c, "scaled:center"))
}

primed_cv_lambda <- function(X_bar, S, Y,
                             lambda_grid = seq(0.05, 2.0, length.out = 20),
                             n_folds = 5, seed = 42) {
  set.seed(seed); n <- length(Y)
  idx1 <- which(Y == 1); idx0 <- which(Y == 0)
  fold_ids <- rep(0, n)
  fold_ids[idx1] <- rep(1:n_folds, length.out = length(idx1))[sample(length(idx1))]
  fold_ids[idx0] <- rep(1:n_folds, length.out = length(idx0))[sample(length(idx0))]

  cv_dev <- rep(0, length(lambda_grid))
  for (j in seq_along(lambda_grid)) {
    td <- 0
    for (fold in 1:n_folds) {
      te <- fold_ids == fold; tr <- !te
      fit <- primed_fit(X_bar[tr,,drop=F], S[tr,,drop=F], Y[tr], lambda_grid[j])
      Xc_te <- sweep(X_bar[te,,drop=F], 2, fit$X_bar_center)
      St_te <- S[te,,drop=F] - rep(1,sum(te)) %o% fit$alpha0 -
               Xc_te * rep(1,sum(te)) %o% fit$alpha1 -
               (Xc_te^2) * rep(1,sum(te)) %o% fit$alpha2
      lp <- fit$beta0 + X_bar[te,,drop=F] %*% fit$beta + St_te %*% fit$gamma
      pred <- pmax(pmin(1/(1+exp(-lp)), 1-1e-10), 1e-10)
      td <- td - 2*mean(Y[te]*log(pred) + (1-Y[te])*log(1-pred))
    }
    cv_dev[j] <- td / n_folds
  }
  lambda_grid[which.min(cv_dev)]
}

# =============================================================================
# 3. Bootstrap Stability
# =============================================================================
B <- 200
K <- length(ordinal_vars)

ordinal_cons_vars <- paste0(ordinal_vars, "_cons")
ordinal_disp_vars <- paste0(ordinal_vars, "_disp")

X_bar_full <- as.matrix(cd_data[, ordinal_cons_vars])
S_full     <- as.matrix(cd_data[, ordinal_disp_vars])
y_full     <- cd_data$malignant
n_full     <- length(y_full)

primed_sel_mat <- matrix(0, B, K)
cdgl_sel_mat   <- matrix(0, B, K)
cl_sel_mat     <- matrix(0, B, K)
colnames(primed_sel_mat) <- colnames(cdgl_sel_mat) <- colnames(cl_sel_mat) <- ordinal_vars

cat(sprintf("\n=== Bootstrap Stability (B=%d) ===\n", B))

for (b in 1:B) {
  if (b %% 20 == 0) cat(sprintf("  Bootstrap %d/%d\n", b, B))
  set.seed(b)

  idx1 <- which(y_full == 1); idx0 <- which(y_full == 0)
  boot_idx <- c(sample(idx1, length(idx1), replace = TRUE),
                sample(idx0, length(idx0), replace = TRUE))

  X_bar_b <- X_bar_full[boot_idx, ]
  S_b     <- S_full[boot_idx, ]
  y_b     <- y_full[boot_idx]

  # --- PRIMED ---
  best_lam <- tryCatch(
    primed_cv_lambda(X_bar_b, S_b, y_b,
                     lambda_grid = seq(0.05, 2.0, length.out = 15),
                     n_folds = 5, seed = b),
    error = function(e) 0.5)
  fit_p <- tryCatch(
    primed_fit(X_bar_b, S_b, y_b, best_lam),
    error = function(e) NULL)
  if (!is.null(fit_p)) {
    primed_sel_mat[b, ] <- as.integer(abs(fit_p$beta) > 1e-8 | abs(fit_p$gamma) > 1e-8)
  }

  # --- CD Group LASSO ---
  cd_vars <- as.vector(rbind(ordinal_cons_vars, ordinal_disp_vars))
  X_cd <- as.matrix(cd_data[boot_idx, cd_vars])
  groups <- rep(1:K, each = 2)
  cv_gl <- tryCatch(
    cv.grpreg(X_cd, y_b, group = groups, family = "binomial",
              penalty = "grLasso", nfolds = 5),
    error = function(e) NULL)
  if (!is.null(cv_gl)) {
    co <- coef(cv_gl, s = "lambda.min")[-1]
    for (k in 1:K) {
      cdgl_sel_mat[b, k] <- as.integer(abs(co[2*k-1]) > 1e-8 | abs(co[2*k]) > 1e-8)
    }
  }

  # --- Consensus LASSO ---
  X_cl <- X_bar_full[boot_idx, ]
  cv_cl <- tryCatch(
    cv.glmnet(X_cl, y_b, family = "binomial", alpha = 1, nfolds = 5),
    error = function(e) NULL)
  if (!is.null(cv_cl)) {
    co <- as.vector(coef(cv_cl, s = "lambda.min"))[-1]
    cl_sel_mat[b, ] <- as.integer(abs(co) > 1e-8)
  }
}

# =============================================================================
# 4. Selection Frequency
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("Selection Frequency (B=", B, ")\n")
cat(strrep("=", 70), "\n\n")

freq_table <- tibble(
  feature = ordinal_vars,
  PRIMED  = round(colMeans(primed_sel_mat), 3),
  CD_GL   = round(colMeans(cdgl_sel_mat), 3),
  C_LASSO = round(colMeans(cl_sel_mat), 3)
) %>% arrange(desc(PRIMED))
print(freq_table, n = Inf)

# =============================================================================
# 5. Jaccard Stability
# =============================================================================
jaccard <- function(a, b) {
  inter <- sum(a & b); union <- sum(a | b)
  if (union == 0) return(1)
  inter / union
}

compute_stability <- function(sel_mat) {
  B <- nrow(sel_mat); jaccards <- numeric(0)
  for (i in 1:(B-1)) {
    for (j in (i+1):min(i+50, B)) {
      jaccards <- c(jaccards, jaccard(sel_mat[i,], sel_mat[j,]))
    }
  }
  c(mean = mean(jaccards), sd = sd(jaccards))
}

cat("\n", strrep("=", 70), "\n")
cat("Jaccard Stability (1=완벽, 0=불안정)\n")
cat(strrep("=", 70), "\n\n")

stab_p  <- compute_stability(primed_sel_mat)
stab_cd <- compute_stability(cdgl_sel_mat)
stab_cl <- compute_stability(cl_sel_mat)

stab_result <- tibble(
  Method = c("PRIMED", "CD Group LASSO", "Consensus LASSO"),
  Jaccard_mean = round(c(stab_p["mean"], stab_cd["mean"], stab_cl["mean"]), 3),
  Jaccard_sd   = round(c(stab_p["sd"], stab_cd["sd"], stab_cl["sd"]), 3)
)
print(stab_result)

# =============================================================================
# 6. 저장
# =============================================================================
saveRDS(list(
  primed_sel_mat = primed_sel_mat,
  cdgl_sel_mat   = cdgl_sel_mat,
  cl_sel_mat     = cl_sel_mat,
  freq_table     = freq_table,
  stab_result    = stab_result
), "results/lidc_bootstrap_stability.rds")

cat("\n결과 저장: results/lidc_bootstrap_stability.rds\n")
