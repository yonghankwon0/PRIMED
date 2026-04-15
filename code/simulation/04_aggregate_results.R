# ============================================================
# 전체 48개 시뮬레이션 조합 — SD 포함 결과 출력
#
# Sim 1: K=20, τ=0.3, DGP{1,2} × S{1,2} × SNR{0.5,1.0}  = 8
# Sim 2: K=20, DGP{1,2} × S{1,2} × τ{0.3,0.6,0.9} × SNR{0.5,1.0} = 24
#        (Sim1은 Sim2의 τ=0.3 부분집합이므로 Sim2만 실행하면 Sim1 포함)
# Sim 3: τ=0.3, DGP{1,2} × S{1,2} × K{50,100} × SNR{0.5,1.0} = 16
# 중복 제거 후 실제: Sim2(24) + Sim3(16) = 40개
# (Sim1의 8개는 Sim2의 τ=0.3에 포함)
# ============================================================
if (!requireNamespace("grpreg", quietly=TRUE)) install.packages("grpreg")
if (!requireNamespace("glmnet", quietly=TRUE)) install.packages("glmnet")
library(grpreg); library(glmnet); library(dplyr, warn.conflicts=FALSE)

WORK_DIR <- getwd()
JULIA_CMD <- "julia"
B <- 50; n <- 200; m <- 5; K_a <- 7

# ============================================================
# 1. Calibration
# ============================================================
calibrate_dgp1 <- function(target_snr, K_a, tau) {
  f <- function(c_val) {
    s2 <- c_val^2 + tau^2
    V_k <- c_val^2 * (1 + exp(s2)*(exp(s2)-1) + 2*c_val*exp(s2/2))
    K_a * V_k / (pi^2/3) - target_snr
  }
  uniroot(f, c(0.01, 3.0), tol=1e-8)$root
}

calibrate_dgp2 <- function(target_snr, K_a) sqrt(target_snr * (pi^2/3) / K_a)

calibrate_beta0 <- function(c_val, tau, K_a, gamma_val, n_mc=200000, seed=42) {
  set.seed(seed)
  mu <- matrix(rnorm(n_mc * K_a), n_mc, K_a)
  sig <- matrix(0, n_mc, K_a)
  for (k in 1:K_a) sig[,k] <- exp(c_val * mu[,k] + rnorm(n_mc, 0, tau))
  lp <- mu %*% rep(c_val, K_a) + sig %*% rep(gamma_val, K_a)
  uniroot(function(b0) mean(1/(1+exp(-(b0+lp)))) - 0.5, c(-20, 20), tol=1e-6)$root
}

# ============================================================
# 2. Helpers
# ============================================================
generate_data <- function(n, K, B, alpha_true, beta_true, gamma_true,
                          beta0, m=5, tau=0.3, seed=2024) {
  set.seed(seed); all_data <- vector("list", B)
  for (b in 1:B) {
    mu <- matrix(rnorm(n*K), n, K)
    X_bar <- S <- sig <- matrix(0, n, K)
    for (k in 1:K) {
      sk <- exp(alpha_true[k]*mu[,k] + rnorm(n,0,tau)); sig[,k] <- sk
      xr <- matrix(rnorm(n*m, rep(mu[,k],m), rep(sk,m)), nrow=n)
      X_bar[,k] <- rowMeans(xr); S[,k] <- apply(xr,1,sd)
    }
    lp <- beta0 + mu %*% beta_true + sig %*% gamma_true
    all_data[[b]] <- list(X_bar=X_bar, S=S, Y=rbinom(n,1,1/(1+exp(-lp))))
  }
  all_data
}

save_csv <- function(all_data, fp) {
  K <- ncol(all_data[[1]]$X_bar); rows <- list()
  for (b in seq_along(all_data)) { d <- all_data[[b]]; nn <- nrow(d$X_bar)
    rows[[b]] <- cbind(rep(b,nn),1:nn,d$X_bar,d$S,d$Y) }
  mat <- do.call(rbind, rows)
  colnames(mat) <- c("rep","i",paste0("X",1:K),paste0("S",1:K),"Y")
  write.csv(mat, fp, row.names=FALSE)
}

fit_cd <- function(dat, nf=5) {
  K <- ncol(dat$X_bar); X <- matrix(0,nrow(dat$X_bar),2*K)
  for(k in 1:K){X[,2*k-1]<-dat$X_bar[,k]; X[,2*k]<-dat$S[,k]}
  colnames(X)<-paste0(rep(c("X","S"),K),rep(1:K,each=2)); g<-rep(1:K,each=2)
  cv<-tryCatch(cv.grpreg(X,dat$Y,g,family="binomial",penalty="grLasso",nfolds=nf),error=function(e)NULL)
  if(is.null(cv))return(NULL)
  co<-coef(cv,s="lambda.min"); bh<-co[seq(2,2*K,2)]; gh<-co[seq(3,2*K+1,2)]
  list(beta0=co[1],beta=bh,gamma=gh,selected=which(bh!=0|gh!=0))
}

fit_cl <- function(dat, nf=5) {
  K<-ncol(dat$X_bar); X<-dat$X_bar; colnames(X)<-paste0("X",1:K)
  cv<-tryCatch(cv.glmnet(X,dat$Y,family="binomial",alpha=1,nfolds=nf),error=function(e)NULL)
  if(is.null(cv))return(NULL)
  co<-as.vector(coef(cv,s="lambda.min")); bh<-co[-1]
  list(beta0=co[1],beta=bh,selected=which(bh!=0))
}

compute_auc <- function(y,p) {
  pos<-p[y==1]; neg<-p[y==0]
  if(length(pos)==0||length(neg)==0)return(NA)
  sum(sapply(pos,function(x)mean(x>neg)+0.5*mean(x==neg)))/length(pos)
}

# ============================================================
# 3. Runner — returns per-replication metrics
# ============================================================
run_sim <- function(tag, K, alpha_true, beta_true, gamma_true,
                    beta0, tau, active, noise_idx=integer(0),
                    julia_core="01_primed_core.jl") {
  inact <- setdiff(1:K, active)
  na_act <- length(active); ni <- length(inact)
  cat(sprintf("\n=== %s (K=%d, tau=%.1f) ===\n", tag, K, tau)); flush.console()

  all_data <- generate_data(n, K, B, alpha_true, beta_true, gamma_true,
                            beta0, m=m, tau=tau, seed=2024)

  # PRIMED via Julia
  dcf  <- file.path(WORK_DIR, sprintf("sim_data_%s.csv", tag))
  mpcf <- file.path(WORK_DIR, sprintf("sim_modprimed_%s.csv", tag))
  save_csv(all_data, dcf)
  cmd <- sprintf('"%s" --threads=auto "%s" "%s" "%s" %d',
                 JULIA_CMD, file.path(WORK_DIR, julia_core), dcf, mpcf, 5)
  cat("  Julia...\n"); flush.console()
  system(cmd, ignore.stdout=TRUE)
  mpdf <- read.csv(mpcf)

  # CD GL + C LASSO
  cat("  CD GL + C LASSO...\n"); flush.console()
  cd_res <- lapply(1:B, function(b) fit_cd(all_data[[b]], 5))
  cl_res <- lapply(1:B, function(b) fit_cl(all_data[[b]], 5))

  # Per-replication metrics
  cat("  Metrics...\n"); flush.console()
  nt <- 1000; set.seed(2024+9999)

  rep_metrics <- data.frame(
    rep=1:B,
    mp_tpr=NA, mp_ppv=NA, mp_f1=NA, mp_auc=NA,
    cd_tpr=NA, cd_ppv=NA, cd_f1=NA, cd_auc=NA,
    cl_tpr=NA, cl_ppv=NA, cl_f1=NA, cl_auc=NA
  )

  f1_fn <- function(tpr, ppv) ifelse(tpr+ppv>0, 2*tpr*ppv/(tpr+ppv), 0)

  for (b in 1:B) {
    # Test data
    mut<-matrix(rnorm(nt*K),nt,K); Xt<-St<-sigt<-matrix(0,nt,K)
    for(k in 1:K){
      sk<-exp(alpha_true[k]*mut[,k]+rnorm(nt,0,tau)); sigt[,k]<-sk
      xr<-matrix(rnorm(nt*m,rep(mut[,k],m),rep(sk,m)),nrow=nt)
      Xt[,k]<-rowMeans(xr); St[,k]<-apply(xr,1,sd)
    }
    lpt<-beta0+mut%*%beta_true+sigt%*%gamma_true; Yt<-rbinom(nt,1,1/(1+exp(-lpt)))

    # PRIMED per-rep
    mp_sel <- as.numeric(mpdf[b, paste0("sel_", 1:K)])
    mp_tpr_b <- mean(mp_sel[active])
    mp_fpr_b <- mean(mp_sel[inact])
    mp_ppv_b <- (mp_tpr_b*na_act) / max(mp_tpr_b*na_act + mp_fpr_b*ni, 1e-10)
    bmp<-as.numeric(mpdf[b,paste0("beta_",1:K)]); gmp<-as.numeric(mpdf[b,paste0("gamma_",1:K)])
    mp_auc_b<-compute_auc(Yt,1/(1+exp(-(mpdf$beta0[b]+Xt%*%bmp+St%*%gmp))))

    rep_metrics$mp_tpr[b] <- mp_tpr_b
    rep_metrics$mp_ppv[b] <- mp_ppv_b
    rep_metrics$mp_f1[b]  <- f1_fn(mp_tpr_b, mp_ppv_b)
    rep_metrics$mp_auc[b] <- mp_auc_b

    # CD GL per-rep
    cd <- cd_res[[b]]
    if (!is.null(cd)) {
      cd_sel <- rep(0, K); cd_sel[cd$selected] <- 1
      cd_tpr_b <- mean(cd_sel[active]); cd_fpr_b <- mean(cd_sel[inact])
      cd_ppv_b <- (cd_tpr_b*na_act) / max(cd_tpr_b*na_act + cd_fpr_b*ni, 1e-10)
      bc<-cd$beta; gc<-cd$gamma
      cd_auc_b <- compute_auc(Yt, 1/(1+exp(-(cd$beta0 + Xt%*%bc + St%*%gc))))
    } else { cd_tpr_b<-NA; cd_ppv_b<-NA; cd_auc_b<-NA }
    rep_metrics$cd_tpr[b] <- cd_tpr_b
    rep_metrics$cd_ppv[b] <- cd_ppv_b
    rep_metrics$cd_f1[b]  <- f1_fn(cd_tpr_b, cd_ppv_b)
    rep_metrics$cd_auc[b] <- cd_auc_b

    # C LASSO per-rep
    cl <- cl_res[[b]]
    if (!is.null(cl)) {
      cl_sel <- rep(0, K); cl_sel[cl$selected] <- 1
      cl_tpr_b <- mean(cl_sel[active]); cl_fpr_b <- mean(cl_sel[inact])
      cl_ppv_b <- (cl_tpr_b*na_act) / max(cl_tpr_b*na_act + cl_fpr_b*ni, 1e-10)
      bcl<-cl$beta
      cl_auc_b <- compute_auc(Yt, 1/(1+exp(-(cl$beta0 + Xt%*%bcl))))
    } else { cl_tpr_b<-NA; cl_ppv_b<-NA; cl_auc_b<-NA }
    rep_metrics$cl_tpr[b] <- cl_tpr_b
    rep_metrics$cl_ppv[b] <- cl_ppv_b
    rep_metrics$cl_f1[b]  <- f1_fn(cl_tpr_b, cl_ppv_b)
    rep_metrics$cl_auc[b] <- cl_auc_b
  }

  unlink(dcf); unlink(mpcf)

  # Summarize: mean (sd)
  summarize_method <- function(prefix) {
    cols <- paste0(prefix, c("_tpr","_ppv","_f1","_auc"))
    vals <- rep_metrics[, cols]
    data.frame(
      TPR_mean = round(mean(vals[,1], na.rm=T), 3),
      TPR_sd   = round(sd(vals[,1], na.rm=T), 3),
      PPV_mean = round(mean(vals[,2], na.rm=T), 3),
      PPV_sd   = round(sd(vals[,2], na.rm=T), 3),
      F1_mean  = round(mean(vals[,3], na.rm=T), 3),
      F1_sd    = round(sd(vals[,3], na.rm=T), 3),
      AUC_mean = round(mean(vals[,4], na.rm=T), 3),
      AUC_sd   = round(sd(vals[,4], na.rm=T), 3)
    )
  }

  res <- rbind(
    cbind(method="PRIMED",  summarize_method("mp")),
    cbind(method="CD GL",   summarize_method("cd")),
    cbind(method="C LASSO", summarize_method("cl"))
  )

  cat(sprintf("  PRIMED: F1=%.3f(%.3f) AUC=%.3f(%.3f)\n",
              res$F1_mean[1], res$F1_sd[1], res$AUC_mean[1], res$AUC_sd[1]))
  cat(sprintf("  CD GL:  F1=%.3f(%.3f) AUC=%.3f(%.3f)\n",
              res$F1_mean[2], res$F1_sd[2], res$AUC_mean[2], res$AUC_sd[2]))
  cat(sprintf("  C LAS:  F1=%.3f(%.3f) AUC=%.3f(%.3f)\n",
              res$F1_mean[3], res$F1_sd[3], res$AUC_mean[3], res$AUC_sd[3]))
  flush.console()
  res
}

# ============================================================
# 4. Run ALL combinations
# ============================================================
all_results <- list(); idx <- 0

# --- Sim 2: τ sensitivity (includes Sim 1 as τ=0.3 subset) ---
for (tau_val in c(0.3, 0.6, 0.9)) {
  for (snr_val in c(0.5, 1.0)) {
    for (dgp in c("DGP1", "DGP2")) {
      if (dgp == "DGP1") {
        c_v <- calibrate_dgp1(snr_val, K_a, tau_val)
        g_v <- c_v
      } else {
        c_v <- calibrate_dgp2(snr_val, K_a)
        g_v <- 0
      }
      b0 <- calibrate_beta0(c_v, tau_val, K_a, gamma_val=g_v)

      for (scen in c("S1", "S2")) {
        alpha_t <- c(rep(c_v, 7),
                     if (scen=="S2") rep(c_v, 3) else rep(0, 3),
                     rep(0, 10))
        beta_t  <- c(rep(c_v, 7), rep(0, 13))
        gamma_t <- c(rep(g_v, 7), rep(0, 13))
        noise   <- if (scen=="S2") 8:10 else integer(0)

        tag <- sprintf("t%.1f_%s_%s_snr%.1f", tau_val, dgp, scen, snr_val)
        idx <- idx + 1
        all_results[[idx]] <- cbind(
          sim="tau", dgp=dgp, scenario=scen, tau=tau_val, snr=snr_val, K=20,
          run_sim(tag, 20, alpha_t, beta_t, gamma_t, b0, tau_val,
                  active=1:7, noise_idx=noise))
      }
    }
  }
}

# --- Sim 3: High-dimensional ---
tau_val <- 0.3
for (K_val in c(50, 100)) {
  jcore <- if (K_val <= 50) "01_primed_core.jl" else "01_primed_core_highdim.jl"
  for (snr_val in c(0.5, 1.0)) {
    for (dgp in c("DGP1", "DGP2")) {
      if (dgp == "DGP1") {
        c_v <- calibrate_dgp1(snr_val, K_a, tau_val)
        g_v <- c_v
      } else {
        c_v <- calibrate_dgp2(snr_val, K_a)
        g_v <- 0
      }
      b0 <- calibrate_beta0(c_v, tau_val, K_a, gamma_val=g_v)

      for (scen in c("S1", "S2")) {
        alpha_t <- c(rep(c_v, 7),
                     if (scen=="S2") rep(c_v, 3) else rep(0, 3),
                     rep(0, K_val - 10))
        beta_t  <- c(rep(c_v, 7), rep(0, K_val - 7))
        gamma_t <- c(rep(g_v, 7), rep(0, K_val - 7))
        noise   <- if (scen=="S2") 8:10 else integer(0)

        tag <- sprintf("K%d_%s_%s_snr%.1f", K_val, dgp, scen, snr_val)
        idx <- idx + 1
        all_results[[idx]] <- cbind(
          sim="hdim", dgp=dgp, scenario=scen, tau=0.3, snr=snr_val, K=K_val,
          run_sim(tag, K_val, alpha_t, beta_t, gamma_t, b0, tau_val,
                  active=1:7, noise_idx=noise, julia_core=jcore))
      }
    }
  }
}

# ============================================================
# 5. 저장
# ============================================================
final_df <- bind_rows(all_results)
write.csv(final_df, file.path(WORK_DIR, "sim_all_with_sd_results.csv"), row.names=FALSE)
cat(sprintf("\n\n=== 전체 %d개 조합 완료. sim_all_with_sd_results.csv 저장 ===\n", idx))
print(final_df, n=Inf, width=150)
