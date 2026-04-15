# =============================================================================
# 10-Fold CV Evaluation v2: 통일된 tuning 설정
#
# 통일 사항:
#   - 외부 CV: 10-fold stratified (동일 fold split)
#   - 내부 CV: 5-fold stratified (lambda 선택)
#   - Lambda: 모두 data-driven, log-spaced 50개
#   - PRIMED: lambda_max 계산 후 log grid
#   - CD GL: grpreg 자동 lambda path (내부 50개)
#   - C LASSO: glmnet 자동 lambda path (내부 50개)
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(grpreg)
library(pROC)
library(glmnet)
library(parallel)

N_CORES <- 20

csv_path <- "/home/yonghankwon0/PRIMED/PRIMED_submission/results/lidc/radiomics_cd.csv"
cd_data <- read_csv(csv_path, show_col_types = FALSE)

cd_data <- cd_data %>%
  filter(malignancy_mean != 3.0) %>%
  mutate(malignant = as.integer(malignancy_mean > 3))

cons_cols <- names(cd_data)[grepl("_cons$", names(cd_data))]
disp_cols <- names(cd_data)[grepl("_disp$", names(cd_data))]
feature_names <- sub("_cons$", "", cons_cols)
K <- length(feature_names)

cd_data <- cd_data %>%
  mutate(across(all_of(c(cons_cols, disp_cols)),
                ~ ifelse(is.finite(.x), .x, NA))) %>%
  mutate(across(all_of(disp_cols), ~ ifelse(is.na(.x), 0, .x)))

valid <- sapply(cons_cols, function(c) mean(is.na(cd_data[[c]])) < 0.1)
cons_cols <- cons_cols[valid]; disp_cols <- disp_cols[valid]
feature_names <- feature_names[valid]; K <- length(feature_names)

for (col in c(cons_cols, disp_cols)) {
  med <- median(cd_data[[col]], na.rm = TRUE)
  cd_data[[col]][is.na(cd_data[[col]])] <- med
}

cat(sprintf("Data: n=%d, K=%d, mal=%.1f%%\n\n", nrow(cd_data), K, 100*mean(cd_data$malignant)))

# =============================================================================
# PRIMED core
# =============================================================================
primed_fit <- function(X_bar, S, Y, lambda, eta=0.001, max_iter=15000, tol=1e-6) {
  n <- length(Y); K <- ncol(X_bar)
  L_med_null <- sum(apply(S,2,var))/2.0
  y_bar <- max(min(mean(Y),1-1e-10),1e-10)
  L_out_null <- -(y_bar*log(y_bar)+(1-y_bar)*log(1-y_bar))
  X_bar_c <- scale(X_bar, scale=FALSE); X_bar_sq <- X_bar_c^2
  alpha0<-colMeans(S); alpha1<-rep(0,K); alpha2<-rep(0,K)
  beta<-rep(0,K); gamma<-rep(0,K)
  beta0<-log(y_bar/(1-y_bar)); sqrt2<-sqrt(2)
  for(iter in 1:max_iter){
    old<-c(alpha0,alpha1,alpha2,beta,gamma)
    S_tilde<-S-rep(1,n)%o%alpha0-X_bar_c*rep(1,n)%o%alpha1-X_bar_sq*rep(1,n)%o%alpha2
    ga0<- -colMeans(S_tilde)/L_med_null; ga1<- -colMeans(S_tilde*X_bar_c)/L_med_null
    ga2<- -colMeans(S_tilde*X_bar_sq)/L_med_null
    lp<-beta0+X_bar%*%beta+S_tilde%*%gamma; p<-1/(1+exp(-lp)); r<-as.vector(p-Y)
    gb0<-mean(r)/L_out_null; gb<-as.vector(t(X_bar)%*%r)/(n*L_out_null)
    gg<-as.vector(t(S_tilde)%*%r)/(n*L_out_null)
    rXc<-as.vector(t(X_bar_c)%*%r)/(n*L_out_null)
    rXs<-as.vector(t(X_bar_sq)%*%r)/(n*L_out_null); rm_val<-mean(r)/L_out_null
    ga0<-ga0-gamma*rm_val; ga1<-ga1-gamma*rXc; ga2<-ga2-gamma*rXs
    beta0<-beta0-eta*gb0; alpha0<-alpha0-eta*ga0; alpha1<-alpha1-eta*ga1; alpha2<-alpha2-eta*ga2
    threshold<-eta*lambda*sqrt2
    for(k in 1:K){bt<-beta[k]-eta*gb[k]; gt<-gamma[k]-eta*gg[k]; nt<-sqrt(bt^2+gt^2)
      if(nt<=threshold){beta[k]<-0;gamma[k]<-0} else{sh<-1-threshold/nt; beta[k]<-sh*bt; gamma[k]<-sh*gt}}
    new<-c(alpha0,alpha1,alpha2,beta,gamma)
    rc<-sqrt(sum((new-old)^2))/(sqrt(sum(old^2))+1e-10)
    if(is.na(rc)||rc<tol) break
  }
  list(alpha0=alpha0,alpha1=alpha1,alpha2=alpha2,beta=beta,gamma=gamma,
       beta0=beta0,X_bar_center=attr(X_bar_c,"scaled:center"),
       selected=which(abs(beta)>1e-8|abs(gamma)>1e-8))
}

primed_predict <- function(fit, X_bar, S) {
  n<-nrow(X_bar); Xc<-sweep(X_bar,2,fit$X_bar_center)
  St<-S-rep(1,n)%o%fit$alpha0-Xc*rep(1,n)%o%fit$alpha1-(Xc^2)*rep(1,n)%o%fit$alpha2
  1/(1+exp(-(fit$beta0+X_bar%*%fit$beta+St%*%fit$gamma)))
}

# Compute data-driven lambda_max for PRIMED
# At null model (beta=gamma=0), gradient of outcome loss w.r.t. (beta_k, gamma_k)
# lambda_max = max_k sqrt(grad_beta_k^2 + grad_gamma_k^2) / sqrt(2)
compute_primed_lambda_max <- function(X_bar, S, Y) {
  n <- length(Y); K <- ncol(X_bar)
  y_bar <- max(min(mean(Y),1-1e-10),1e-10)
  L_out_null <- -(y_bar*log(y_bar)+(1-y_bar)*log(1-y_bar))
  # At null: p = y_bar for all i, r = y_bar - Y
  r <- rep(y_bar, n) - Y
  # S_tilde at null (alpha=0): S_tilde = S - colMeans(S)
  S_tilde <- scale(S, scale=FALSE)
  gb <- as.vector(t(X_bar) %*% r) / (n * L_out_null)
  gg <- as.vector(t(S_tilde) %*% r) / (n * L_out_null)
  group_norms <- sqrt(gb^2 + gg^2)
  max(group_norms) / sqrt(2)
}

# PRIMED CV with data-driven lambda path
primed_cv_fit <- function(X_bar_tr, S_tr, y_tr, n_lambda=50, n_inner_folds=5,
                          lambda_ratio=0.001) {
  # Data-driven lambda path
  lam_max <- compute_primed_lambda_max(X_bar_tr, S_tr, y_tr)
  lam_min <- lam_max * lambda_ratio
  lambda_grid <- exp(seq(log(lam_max), log(lam_min), length.out=n_lambda))

  cat(sprintf("[lambda: %.4f to %.4f, %d values] ", lam_max, lam_min, n_lambda))

  n <- length(y_tr)
  idx1<-which(y_tr==1); idx0<-which(y_tr==0)
  fold_ids<-rep(0,n)
  fold_ids[idx1]<-rep(1:n_inner_folds,length.out=length(idx1))[sample(length(idx1))]
  fold_ids[idx0]<-rep(1:n_inner_folds,length.out=length(idx0))[sample(length(idx0))]

  cv_dev <- mclapply(seq_along(lambda_grid), function(j) {
    lam <- lambda_grid[j]; td <- 0
    for(fold in 1:n_inner_folds) {
      te<-fold_ids==fold; tr<-!te
      fit<-primed_fit(X_bar_tr[tr,],S_tr[tr,],y_tr[tr],lam)
      pred<-as.vector(primed_predict(fit,X_bar_tr[te,],S_tr[te,]))
      pred<-pmax(pmin(pred,1-1e-10),1e-10)
      td<-td-2*mean(y_tr[te]*log(pred)+(1-y_tr[te])*log(1-pred))
    }
    td/n_inner_folds
  }, mc.cores = N_CORES)

  cv_dev <- unlist(cv_dev)
  best_idx <- which.min(cv_dev)
  best_lam <- lambda_grid[best_idx]
  best_fit <- primed_fit(X_bar_tr, S_tr, y_tr, best_lam)

  # Report where in the path the best lambda is
  cat(sprintf("best=%d/%d (%.4f) ", best_idx, n_lambda, best_lam))

  list(fit=best_fit, lambda=best_lam, lambda_idx=best_idx, n_lambda=n_lambda,
       lambda_grid=lambda_grid, cv_dev=cv_dev)
}

# =============================================================================
# 10-Fold CV
# =============================================================================
n_folds <- 10
n_inner_folds <- 5  # 통일: 모든 method 5-fold inner CV
set.seed(2024)

y_all <- cd_data$malignant
idx1 <- which(y_all==1); idx0 <- which(y_all==0)
fold_ids <- rep(0, nrow(cd_data))
fold_ids[idx1] <- rep(1:n_folds, length.out=length(idx1))[sample(length(idx1))]
fold_ids[idx0] <- rep(1:n_folds, length.out=length(idx0))[sample(length(idx0))]

results <- list()
selected_vars <- list(primed=list(), cdgl=list(), classo=list())

cat(sprintf("=== 10-Fold CV v2 (n=%d, K=%d, %d cores) ===\n", nrow(cd_data), K, N_CORES))
cat(sprintf("Inner CV: %d-fold for all methods\n", n_inner_folds))
cat(sprintf("PRIMED lambda: data-driven log-spaced 100 values\n"))
cat(sprintf("CD GL / C LASSO: package default lambda path\n\n"))

for (fold in 1:n_folds) {
  cat(sprintf("--- Fold %d/%d ---\n", fold, n_folds))
  te <- fold_ids == fold; tr <- !te

  train_data <- cd_data[tr, ]; test_data <- cd_data[te, ]
  X_bar_tr <- as.matrix(scale(train_data[, cons_cols]))
  S_tr <- as.matrix(scale(train_data[, disp_cols]))
  y_tr <- train_data$malignant

  center_X <- attr(X_bar_tr, "scaled:center"); scale_X <- attr(X_bar_tr, "scaled:scale")
  center_S <- attr(S_tr, "scaled:center"); scale_S <- attr(S_tr, "scaled:scale")
  X_bar_te <- scale(test_data[, cons_cols], center=center_X, scale=scale_X)
  S_te <- scale(test_data[, disp_cols], center=center_S, scale=scale_S)
  y_te <- test_data$malignant

  # --- PRIMED (data-driven lambda) ---
  cat("  PRIMED: ")
  t0 <- proc.time()
  pr <- primed_cv_fit(X_bar_tr, S_tr, y_tr, n_lambda=100,
                       n_inner_folds=n_inner_folds, lambda_ratio=0.001)
  pred_pr <- as.vector(primed_predict(pr$fit, X_bar_te, S_te))
  auc_pr <- as.numeric(auc(roc(y_te, pred_pr, quiet=TRUE)))
  sel_pr <- pr$fit$selected
  cat(sprintf("AUC=%.3f, sel=%d (%.0fs)\n",
      auc_pr, length(sel_pr), (proc.time()-t0)[3]))

  # --- CD Group LASSO (5-fold inner CV) ---
  cat("  CD GL:  ")
  cd_tr <- cbind(X_bar_tr, S_tr); cd_te <- cbind(X_bar_te, S_te)
  groups <- rep(1:K, each=2)
  cv_gl <- tryCatch(cv.grpreg(cd_tr, y_tr, group=groups, family="binomial",
                               penalty="grLasso", nfolds=n_inner_folds),
                    error=function(e) NULL)
  if (!is.null(cv_gl)) {
    pred_gl <- as.vector(predict(cv_gl, cd_te, type="response", lambda=cv_gl$lambda.min))
    auc_gl <- as.numeric(auc(roc(y_te, pred_gl, quiet=TRUE)))
    co_gl <- coef(cv_gl, s="lambda.min")[-1]
    sel_gl <- unique(groups[which(abs(co_gl) > 1e-8)])
  } else { auc_gl <- NA; sel_gl <- integer(0) }
  cat(sprintf("AUC=%.3f, sel=%d\n", auc_gl, length(sel_gl)))

  # --- Consensus LASSO (5-fold inner CV) ---
  cat("  C LAS:  ")
  cv_cl <- tryCatch(cv.glmnet(X_bar_tr, y_tr, family="binomial", alpha=1,
                               nfolds=n_inner_folds),
                    error=function(e) NULL)
  if (!is.null(cv_cl)) {
    pred_cl <- as.vector(predict(cv_cl, X_bar_te, type="response", s="lambda.min"))
    auc_cl <- as.numeric(auc(roc(y_te, pred_cl, quiet=TRUE)))
    co_cl <- as.vector(coef(cv_cl, s="lambda.min"))[-1]
    sel_cl <- which(abs(co_cl) > 1e-8)
  } else { auc_cl <- NA; sel_cl <- integer(0) }
  cat(sprintf("AUC=%.3f, sel=%d\n\n", auc_cl, length(sel_cl)))

  results[[fold]] <- data.frame(fold=fold,
    primed_auc=auc_pr, cdgl_auc=auc_gl, classo_auc=auc_cl,
    primed_nsel=length(sel_pr), cdgl_nsel=length(sel_gl), classo_nsel=length(sel_cl),
    primed_lambda=pr$lambda, primed_lambda_idx=pr$lambda_idx)

  selected_vars$primed[[fold]] <- sel_pr
  selected_vars$cdgl[[fold]] <- sel_gl
  selected_vars$classo[[fold]] <- sel_cl
}

# =============================================================================
# 결과
# =============================================================================
res_df <- bind_rows(results)

cat("\n", strrep("=", 60), "\n")
cat("  10-Fold CV Results (v2: unified tuning)\n")
cat(strrep("=", 60), "\n\n")

cat("Per-fold:\n")
print(res_df %>% dplyr::select(fold, primed_auc, cdgl_auc, classo_auc,
                         primed_nsel, cdgl_nsel, classo_nsel,
                         primed_lambda, primed_lambda_idx), row.names=FALSE)

cat(sprintf("\n\nMean AUC (SD):\n"))
cat(sprintf("  PRIMED:  %.3f (%.3f)\n", mean(res_df$primed_auc,na.rm=T), sd(res_df$primed_auc,na.rm=T)))
cat(sprintf("  CD GL:   %.3f (%.3f)\n", mean(res_df$cdgl_auc,na.rm=T), sd(res_df$cdgl_auc,na.rm=T)))
cat(sprintf("  C LASSO: %.3f (%.3f)\n", mean(res_df$classo_auc,na.rm=T), sd(res_df$classo_auc,na.rm=T)))

cat(sprintf("\nMean # selected (SD):\n"))
cat(sprintf("  PRIMED:  %.1f (%.1f)\n", mean(res_df$primed_nsel), sd(res_df$primed_nsel)))
cat(sprintf("  CD GL:   %.1f (%.1f)\n", mean(res_df$cdgl_nsel), sd(res_df$cdgl_nsel)))
cat(sprintf("  C LASSO: %.1f (%.1f)\n", mean(res_df$classo_nsel), sd(res_df$classo_nsel)))

cat(sprintf("\nPRIMED lambda: mean=%.4f (%.4f), mean idx=%.1f/%d\n",
    mean(res_df$primed_lambda), sd(res_df$primed_lambda),
    mean(res_df$primed_lambda_idx), 100))

cat(sprintf("\nPaired t-test (AUC):\n"))
tt1 <- t.test(res_df$primed_auc, res_df$cdgl_auc, paired=TRUE)
tt2 <- t.test(res_df$primed_auc, res_df$classo_auc, paired=TRUE)
tt3 <- t.test(res_df$cdgl_auc, res_df$classo_auc, paired=TRUE)
cat(sprintf("  PRIMED vs CD GL:    diff=%.4f, p=%.4f\n", tt1$estimate, tt1$p.value))
cat(sprintf("  PRIMED vs C LASSO:  diff=%.4f, p=%.4f\n", tt2$estimate, tt2$p.value))
cat(sprintf("  CD GL vs C LASSO:   diff=%.4f, p=%.4f\n", tt3$estimate, tt3$p.value))

# Variable selection frequency
pr_freq <- table(unlist(selected_vars$primed))
gl_freq <- table(unlist(selected_vars$cdgl))
cl_freq <- table(unlist(selected_vars$classo))

sel_summary <- data.frame(idx=1:K, feature=feature_names, primed=0, cdgl=0, classo=0)
for (i in names(pr_freq)) sel_summary$primed[as.integer(i)] <- pr_freq[i]
for (i in names(gl_freq)) sel_summary$cdgl[as.integer(i)] <- gl_freq[i]
for (i in names(cl_freq)) sel_summary$classo[as.integer(i)] <- cl_freq[i]

cat("\n\nFeatures selected in >= 5/10 folds:\n")
sel_any5 <- sel_summary %>% filter(primed>=5|cdgl>=5|classo>=5) %>% arrange(desc(primed),desc(cdgl))
print(as.data.frame(sel_any5), row.names=FALSE)

write.csv(res_df, "/home/yonghankwon0/PRIMED/PRIMED_submission/results/lidc/cv10v2_auc_results.csv", row.names=FALSE)
write.csv(sel_summary, "/home/yonghankwon0/PRIMED/PRIMED_submission/results/lidc/cv10v2_selection_frequency.csv", row.names=FALSE)
saveRDS(list(results=res_df, selected_vars=selected_vars, sel_summary=sel_summary),
        "/home/yonghankwon0/PRIMED/PRIMED_submission/results/lidc/cv10v2_full_results.rds")
cat(sprintf("\nSaved: cv10v2_*.csv/rds\n"))

# Jaccard stability
compute_pairwise_jaccard <- function(sl) {
  n <- length(sl); jac <- c()
  for (i in 1:(n-1)) for (j in (i+1):n) {
    inter <- length(intersect(sl[[i]], sl[[j]]))
    uni <- length(union(sl[[i]], sl[[j]]))
    if (uni > 0) jac <- c(jac, inter/uni)
  }
  jac
}

cat("\n\nJaccard Stability:\n")
for (m in c("primed","cdgl","classo")) {
  mn <- c(primed="PRIMED",cdgl="CD GL",classo="C LASSO")[m]
  jac <- compute_pairwise_jaccard(selected_vars[[m]])
  cat(sprintf("  %s: Jaccard = %.3f (%.3f)\n", mn, mean(jac), sd(jac)))
}
