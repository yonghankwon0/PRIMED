# Compute pairwise Jaccard index for Table 5
# Targets: DGP1, S2, SNR=1.0, K=20/50/100
# Runs PRIMED (Julia), CD GL, C LASSO, then computes Jaccard from per-rep selections.

WORK_DIR <- getwd()
JULIA_CMD <- "julia"
JULIA_THREADS <- 8L
library(grpreg); library(glmnet)

B <- 50; n <- 200; m <- 5; K_a <- 7; tau_val <- 0.3; snr_val <- 1.0

# Calibration
calibrate_dgp1 <- function(target_snr, K_a, tau) {
  f <- function(c_val) {
    s2 <- c_val^2 + tau^2
    V_k <- c_val^2 * (1 + exp(s2)*(exp(s2)-1) + 2*c_val*exp(s2/2))
    K_a * V_k / (pi^2/3) - target_snr
  }
  uniroot(f, c(0.01, 3.0), tol=1e-8)$root
}
calibrate_beta0 <- function(c_val, tau, K_a, gamma_val, n_mc=200000, seed=42) {
  set.seed(seed); mu <- matrix(rnorm(n_mc*K_a), n_mc, K_a)
  sig <- matrix(0, n_mc, K_a)
  for (k in 1:K_a) sig[,k] <- exp(c_val*mu[,k]+rnorm(n_mc,0,tau))
  lp <- mu%*%rep(c_val,K_a)+sig%*%rep(gamma_val,K_a)
  uniroot(function(b0) mean(1/(1+exp(-(b0+lp))))-0.5, c(-20,20), tol=1e-6)$root
}

generate_data <- function(n,K,B,alpha_true,beta_true,gamma_true,beta0,m=5,tau=0.3,seed=2024) {
  set.seed(seed); all_data <- vector("list",B)
  for(b in 1:B) { mu<-matrix(rnorm(n*K),n,K); X_bar<-S<-sig<-matrix(0,n,K)
    for(k in 1:K){sk<-exp(alpha_true[k]*mu[,k]+rnorm(n,0,tau)); sig[,k]<-sk
      xr<-matrix(rnorm(n*m,rep(mu[,k],m),rep(sk,m)),nrow=n)
      X_bar[,k]<-rowMeans(xr); S[,k]<-apply(xr,1,sd)}
    lp<-beta0+mu%*%beta_true+sig%*%gamma_true
    all_data[[b]]<-list(X_bar=X_bar,S=S,Y=rbinom(n,1,1/(1+exp(-lp))))}
  all_data
}
save_csv <- function(all_data,fp) {
  K<-ncol(all_data[[1]]$X_bar); rows<-list()
  for(b in seq_along(all_data)){d<-all_data[[b]]; nn<-nrow(d$X_bar)
    rows[[b]]<-cbind(rep(b,nn),1:nn,d$X_bar,d$S,d$Y)}
  mat<-do.call(rbind,rows); colnames(mat)<-c("rep","i",paste0("X",1:K),paste0("S",1:K),"Y")
  write.csv(mat,fp,row.names=FALSE)
}
fit_cd <- function(dat,nf=5) {
  K<-ncol(dat$X_bar); X<-matrix(0,nrow(dat$X_bar),2*K)
  for(k in 1:K){X[,2*k-1]<-dat$X_bar[,k]; X[,2*k]<-dat$S[,k]}
  colnames(X)<-paste0(rep(c("X","S"),K),rep(1:K,each=2)); g<-rep(1:K,each=2)
  cv<-tryCatch(cv.grpreg(X,dat$Y,g,family="binomial",penalty="grLasso",nfolds=nf),error=function(e)NULL)
  if(is.null(cv))return(NULL); co<-coef(cv,s="lambda.min"); bh<-co[seq(2,2*K,2)]; gh<-co[seq(3,2*K+1,2)]
  list(selected=which(bh!=0|gh!=0))
}
fit_cl <- function(dat,nf=5) {
  K<-ncol(dat$X_bar); X<-dat$X_bar; colnames(X)<-paste0("X",1:K)
  cv<-tryCatch(cv.glmnet(X,dat$Y,family="binomial",alpha=1,nfolds=nf),error=function(e)NULL)
  if(is.null(cv))return(NULL); co<-as.vector(coef(cv,s="lambda.min")); bh<-co[-1]
  list(selected=which(bh!=0))
}

# Pairwise Jaccard: mean over all C(B,2) pairs
pairwise_jaccard <- function(sel_mat) {
  # sel_mat: B x K binary matrix
  B <- nrow(sel_mat)
  jaccards <- numeric(B*(B-1)/2)
  idx <- 0
  for (i in 1:(B-1)) {
    for (j in (i+1):B) {
      idx <- idx + 1
      inter <- sum(sel_mat[i,] & sel_mat[j,])
      union <- sum(sel_mat[i,] | sel_mat[j,])
      jaccards[idx] <- if (union > 0) inter / union else 0
    }
  }
  mean(jaccards)
}

# DGP1 S2 SNR=1.0
c_v <- calibrate_dgp1(snr_val, K_a, tau_val)
g_v <- c_v
b0 <- calibrate_beta0(c_v, tau_val, K_a, gamma_val=g_v)

results <- data.frame(K=integer(), PRIMED=numeric(), CD_GL=numeric(), C_LASSO=numeric())

for (K_val in c(20, 50, 100)) {
  cat(sprintf("\n========== K=%d ==========\n", K_val)); flush.console()
  jcore <- if(K_val<=50) "01_primed_core.jl" else "01_primed_core_highdim.jl"

  alpha_t <- c(rep(c_v,7), rep(c_v,3), rep(0,K_val-10))  # S2: noise present
  beta_t  <- c(rep(c_v,7), rep(0,K_val-7))
  gamma_t <- c(rep(g_v,7), rep(0,K_val-7))

  # Generate data
  cat("  Generating data...\n"); flush.console()
  all_data <- generate_data(n, K_val, B, alpha_t, beta_t, gamma_t, b0, m=m, tau=tau_val, seed=2024)

  # PRIMED via Julia
  tag <- sprintf("jaccard_K%d", K_val)
  dcf  <- file.path(WORK_DIR, sprintf("sim_data_%s.csv", tag))
  mpcf <- file.path(WORK_DIR, sprintf("sim_modprimed_%s.csv", tag))
  save_csv(all_data, dcf)
  cmd <- sprintf('"%s" --threads=%d "%s" "%s" "%s" %d',
                 JULIA_CMD, JULIA_THREADS, file.path(WORK_DIR, jcore), dcf, mpcf, 5)
  cat(sprintf("  Running PRIMED (Julia, %d threads)...\n", JULIA_THREADS)); flush.console()
  system(cmd, ignore.stdout=TRUE)
  mpdf <- read.csv(mpcf)

  # Per-rep selection matrices
  mp_sel_mat <- as.matrix(mpdf[, paste0("sel_", 1:K_val)])

  cat("  Running CD GL + C LASSO...\n"); flush.console()
  cd_sel_mat <- matrix(0, B, K_val)
  cl_sel_mat <- matrix(0, B, K_val)
  for (b in 1:B) {
    cd <- fit_cd(all_data[[b]], 5)
    if (!is.null(cd)) cd_sel_mat[b, cd$selected] <- 1
    cl <- fit_cl(all_data[[b]], 5)
    if (!is.null(cl)) cl_sel_mat[b, cl$selected] <- 1
  }

  # Compute Jaccard
  cat("  Computing pairwise Jaccard...\n"); flush.console()
  j_mp <- pairwise_jaccard(mp_sel_mat)
  j_cd <- pairwise_jaccard(cd_sel_mat)
  j_cl <- pairwise_jaccard(cl_sel_mat)

  cat(sprintf("  Jaccard: PRIMED=%.3f, CD GL=%.3f, C LASSO=%.3f\n", j_mp, j_cd, j_cl))
  results <- rbind(results, data.frame(K=K_val, PRIMED=round(j_mp,3), CD_GL=round(j_cd,3), C_LASSO=round(j_cl,3)))

  unlink(dcf); unlink(mpcf)
}

cat("\n\n===== FINAL JACCARD RESULTS =====\n")
print(results)
write.csv(results, file.path(WORK_DIR, "jaccard_results.csv"), row.names=FALSE)
cat("Saved: jaccard_results.csv\n")
