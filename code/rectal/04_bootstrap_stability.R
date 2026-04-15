# =============================================================================
# Rectal Cancer: Bootstrap Selection Stability
#
# 원본: Center 1 train → Centers 2-3 test (single split, stability 불가)
# 여기서: B=200 bootstrap resampling으로 selection frequency 계산
#
# 사용법: Rscript 01_rectal_bootstrap_stability.R
# =============================================================================

library(tidyverse)
library(readxl)
library(grpreg)
library(pROC)
library(glmnet)

# =============================================================================
# 1. 데이터 로딩
# =============================================================================
data <- read_excel("data/data.xlsx")

data_clean <- data %>%
  slice(-(1:2)) %>%
  rename(Pre_DWI_ADC = `DWI-ADC...12`, Post_DWI_ADC = `DWI-ADC...23`, pCR = `...28`) %>%
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
    pCR = first(pCR), Center = first(Center),
    across(all_of(binary_vars), list(cons = ~ mean(.x, na.rm = TRUE))),
    across(all_of(ordinal_vars), list(
      cons = ~ mean(.x, na.rm = TRUE),
      disp = ~ sd(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_disp"), ~ replace_na(.x, 0)))

cat(sprintf("Total subjects: %d, pCR rate: %.1f%%\n",
            nrow(cd_data), 100 * mean(cd_data$pCR)))

# =============================================================================
# 2. PRIMED 구현 (논문 일치: residual S_tilde, alpha unpenalized, joint grad)
# =============================================================================
primed_fit <- function(X_bar, S, Y, X_bin = NULL, lambda,
                       eta = 0.005, max_iter = 10000, tol = 1e-6) {
  n <- length(Y); K <- ncol(X_bar)
  J <- if (!is.null(X_bin)) ncol(X_bin) else 0

  L_med_null <- sum(apply(S, 2, var)) / 2.0
  y_bar <- max(min(mean(Y), 1 - 1e-10), 1e-10)
  L_out_null <- -(y_bar * log(y_bar) + (1 - y_bar) * log(1 - y_bar))

  alpha <- rep(0, K); alpha0 <- colMeans(S)
  beta <- rep(0, K); gamma <- rep(0, K)
  delta <- if (J > 0) rep(0, J) else NULL
  beta0 <- log(y_bar / (1 - y_bar))
  sqrt2 <- sqrt(2)

  for (iter in 1:max_iter) {
    old_params <- c(alpha, alpha0, beta, gamma, delta)

    S_tilde <- S - rep(1, n) %o% alpha0 - X_bar * rep(1, n) %o% alpha

    grad_alpha  <- -colMeans(S_tilde * X_bar) / L_med_null
    grad_alpha0 <- -colMeans(S_tilde) / L_med_null

    lp <- beta0 + X_bar %*% beta + S_tilde %*% gamma
    if (J > 0) lp <- lp + X_bin %*% delta
    p <- 1 / (1 + exp(-lp)); r <- as.vector(p - Y)

    grad_beta0 <- mean(r) / L_out_null
    grad_beta  <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_gamma <- as.vector(t(S_tilde) %*% r) / (n * L_out_null)
    if (J > 0) grad_delta <- as.vector(t(X_bin) %*% r) / (n * L_out_null)

    grad_alpha  <- grad_alpha  - gamma * as.vector(t(X_bar) %*% r) / (n * L_out_null)
    grad_alpha0 <- grad_alpha0 - gamma * mean(r) / L_out_null

    beta0  <- beta0  - eta * grad_beta0
    alpha0 <- alpha0 - eta * grad_alpha0
    alpha  <- alpha  - eta * grad_alpha

    threshold <- eta * lambda * sqrt2
    for (k in 1:K) {
      bt <- beta[k] - eta * grad_beta[k]
      gt <- gamma[k] - eta * grad_gamma[k]
      nt <- sqrt(bt^2 + gt^2)
      if (nt <= threshold) { beta[k] <- 0; gamma[k] <- 0 }
      else { sh <- 1 - threshold / nt; beta[k] <- sh * bt; gamma[k] <- sh * gt }
    }

    if (J > 0) {
      threshold_l1 <- eta * lambda
      for (j in 1:J) {
        dt <- delta[j] - eta * grad_delta[j]
        delta[j] <- sign(dt) * max(abs(dt) - threshold_l1, 0)
      }
    }

    new_params <- c(alpha, alpha0, beta, gamma, delta)
    rel_change <- sqrt(sum((new_params - old_params)^2)) / (sqrt(sum(old_params^2)) + 1e-10)
    if (is.na(rel_change) || rel_change < tol) break
  }

  list(alpha = alpha, beta = beta, gamma = gamma, delta = delta,
       alpha0 = alpha0, beta0 = beta0)
}

primed_cv_lambda <- function(X_bar, S, Y, X_bin = NULL,
                             lambda_grid = seq(0.01, 3.0, length.out = 50),
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
      fit <- primed_fit(X_bar[tr,,drop=F], S[tr,,drop=F], Y[tr],
                        if(!is.null(X_bin)) X_bin[tr,,drop=F] else NULL, lambda_grid[j])
      S_tilde_te <- S[te,,drop=F] - rep(1,sum(te)) %o% fit$alpha0 -
                    X_bar[te,,drop=F] * rep(1,sum(te)) %o% fit$alpha
      lp <- fit$beta0 + X_bar[te,,drop=F] %*% fit$beta + S_tilde_te %*% fit$gamma
      if (!is.null(X_bin)) lp <- lp + X_bin[te,,drop=F] %*% fit$delta
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
B <- 200  # bootstrap 반복 수
K_ord <- length(ordinal_vars)
K_bin <- length(binary_vars)

ordinal_cons_vars <- paste0(ordinal_vars, "_cons")
ordinal_disp_vars <- paste0(ordinal_vars, "_disp")
binary_cons_vars  <- paste0(binary_vars, "_cons")

# 전체 데이터 사용 (Center 구분 없이)
X_bar_full <- as.matrix(cd_data[, ordinal_cons_vars])
S_full     <- as.matrix(cd_data[, ordinal_disp_vars])
X_bin_full <- as.matrix(cd_data[, binary_cons_vars])
y_full     <- cd_data$pCR
n_full     <- length(y_full)

# 저장 행렬: 각 bootstrap에서 변수별 선택 여부
primed_sel_mat <- matrix(0, B, K_ord)  # ordinal
primed_bin_mat <- matrix(0, B, K_bin)  # binary
cdgl_sel_mat   <- matrix(0, B, K_ord)
cl_sel_mat     <- matrix(0, B, K_ord + K_bin)

colnames(primed_sel_mat) <- ordinal_vars
colnames(primed_bin_mat) <- binary_vars
colnames(cdgl_sel_mat)   <- ordinal_vars

cat(sprintf("\n=== Bootstrap Stability (B=%d) ===\n", B))

for (b in 1:B) {
  if (b %% 20 == 0) cat(sprintf("  Bootstrap %d/%d\n", b, B))
  set.seed(b)

  # Stratified bootstrap
  idx1 <- which(y_full == 1); idx0 <- which(y_full == 0)
  boot_idx <- c(sample(idx1, length(idx1), replace = TRUE),
                sample(idx0, length(idx0), replace = TRUE))

  X_bar_b <- X_bar_full[boot_idx, ]
  S_b     <- S_full[boot_idx, ]
  X_bin_b <- X_bin_full[boot_idx, ]
  y_b     <- y_full[boot_idx]

  # --- PRIMED ---
  best_lam <- primed_cv_lambda(X_bar_b, S_b, y_b, X_bin_b,
                               lambda_grid = seq(0.1, 3.0, length.out = 30),
                               n_folds = 5, seed = b)
  fit_p <- primed_fit(X_bar_b, S_b, y_b, X_bin_b, best_lam)
  primed_sel_mat[b, ] <- as.integer(abs(fit_p$beta) > 1e-8 | abs(fit_p$gamma) > 1e-8)
  primed_bin_mat[b, ] <- as.integer(abs(fit_p$delta) > 1e-8)

  # --- CD Group LASSO ---
  cd_vars <- as.vector(rbind(ordinal_cons_vars, ordinal_disp_vars))
  X_cd <- cbind(X_bin_b, as.matrix(cd_data[boot_idx, cd_vars]))
  groups <- c(1:K_bin, rep((K_bin+1):(K_bin+K_ord), each = 2))
  cv_gl <- tryCatch(
    cv.grpreg(X_cd, y_b, group = groups, family = "binomial",
              penalty = "grLasso", nfolds = 5),
    error = function(e) NULL)
  if (!is.null(cv_gl)) {
    co <- coef(cv_gl, s = "lambda.min")[-1]
    for (k in 1:K_ord) {
      cdgl_sel_mat[b, k] <- as.integer(
        abs(co[K_bin + 2*k - 1]) > 1e-8 | abs(co[K_bin + 2*k]) > 1e-8)
    }
  }

  # --- Consensus LASSO ---
  cons_vars <- c(binary_cons_vars, ordinal_cons_vars)
  X_cl <- as.matrix(cd_data[boot_idx, cons_vars])
  cv_cl <- tryCatch(
    cv.glmnet(X_cl, y_b, family = "binomial", alpha = 1, nfolds = 5),
    error = function(e) NULL)
  if (!is.null(cv_cl)) {
    co <- as.vector(coef(cv_cl, s = "lambda.min"))[-1]
    cl_sel_mat[b, ] <- as.integer(abs(co) > 1e-8)
  }
}

# =============================================================================
# 4. Selection Frequency 계산
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("PRIMED: Ordinal Feature Selection Frequency (B=", B, ")\n")
cat(strrep("=", 70), "\n\n")

primed_freq <- colMeans(primed_sel_mat)
cdgl_freq   <- colMeans(cdgl_sel_mat)

freq_table <- tibble(
  feature = ordinal_vars,
  PRIMED  = round(primed_freq, 3),
  CD_GL   = round(cdgl_freq, 3),
  diff    = round(primed_freq - cdgl_freq, 3)
) %>% arrange(desc(PRIMED))
print(freq_table, n = Inf)

cat("\nPRIMED: Binary Feature Selection Frequency\n")
bin_freq <- tibble(
  feature = binary_vars,
  PRIMED  = round(colMeans(primed_bin_mat), 3)
)
print(bin_freq, n = Inf)

# Consensus LASSO frequency
cat("\nConsensus LASSO: Selection Frequency\n")
cl_freq <- tibble(
  feature = c(binary_vars, ordinal_vars),
  C_LASSO = round(colMeans(cl_sel_mat), 3)
) %>% arrange(desc(C_LASSO))
print(cl_freq, n = Inf)

# =============================================================================
# 5. Pairwise Jaccard Stability
# =============================================================================
jaccard <- function(a, b) {
  inter <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(1)
  inter / union
}

compute_stability <- function(sel_mat) {
  B <- nrow(sel_mat)
  jaccards <- numeric(0)
  for (i in 1:(B-1)) {
    for (j in (i+1):min(i+50, B)) {  # 계산 효율을 위해 근접 쌍만
      jaccards <- c(jaccards, jaccard(sel_mat[i,], sel_mat[j,]))
    }
  }
  c(mean = mean(jaccards), sd = sd(jaccards))
}

cat("\n", strrep("=", 70), "\n")
cat("Pairwise Jaccard Stability (1=완벽 안정, 0=완전 불안정)\n")
cat(strrep("=", 70), "\n\n")

# PRIMED: ordinal + binary 합쳐서
primed_full_sel <- cbind(primed_sel_mat, primed_bin_mat)
stab_primed <- compute_stability(primed_full_sel)
stab_cdgl   <- compute_stability(cdgl_sel_mat)
stab_cl     <- compute_stability(cl_sel_mat)

stab_result <- tibble(
  Method = c("PRIMED", "CD Group LASSO", "Consensus LASSO"),
  Jaccard_mean = round(c(stab_primed["mean"], stab_cdgl["mean"], stab_cl["mean"]), 3),
  Jaccard_sd   = round(c(stab_primed["sd"], stab_cdgl["sd"], stab_cl["sd"]), 3)
)
print(stab_result)

# =============================================================================
# 6. 결과 저장
# =============================================================================
saveRDS(list(
  primed_sel_mat = primed_sel_mat,
  primed_bin_mat = primed_bin_mat,
  cdgl_sel_mat   = cdgl_sel_mat,
  cl_sel_mat     = cl_sel_mat,
  freq_table     = freq_table,
  stab_result    = stab_result
), "results/rectal_bootstrap_stability.rds")

cat("\n결과 저장: results/rectal_bootstrap_stability.rds\n")
