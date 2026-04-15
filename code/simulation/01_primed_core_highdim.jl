# ============================================================
# Modified PRIMED Algorithm — No-omega version
# L_joint = L_med / L_med_null + L_out / L_out_null
#
# Key difference from PRIMED GL:
#   - α is estimated freely (not penalized)
#   - Group penalty on (β_k, γ_k) only
#   - Proximal norm = sqrt(β² + γ²), threshold = η·λ·√2
#   - Selection based on β_k ≠ 0 or γ_k ≠ 0
#
# Usage: julia --threads=auto 01_primed_core_highdim.jl <data_csv> <result_csv> [n_folds]
# ============================================================

using LinearAlgebra, Statistics, Random, Printf, DelimitedFiles

# ============================================================
# 1. Modified PRIMED: Joint Loss + Group Penalty on (β,γ) only
# ============================================================

function modified_primed(X_bar::Matrix{Float64}, S::Matrix{Float64}, Y::Vector{Float64},
                         lambda::Float64;
                         eta::Float64=0.005, max_iter::Int=10000, tol::Float64=1e-6)
    n, K = size(X_bar)

    # Null losses for normalization
    L_med_null = sum(var(S[:, k]; corrected=false) for k in 1:K) / 2.0
    y_bar = clamp(mean(Y), 1e-10, 1.0 - 1e-10)
    L_out_null = -(y_bar * log(y_bar) + (1.0 - y_bar) * log(1.0 - y_bar))

    alpha  = zeros(K)
    alpha0 = vec(mean(S, dims=1))
    beta   = zeros(K)
    gamma  = zeros(K)
    beta0  = log(clamp(mean(Y), 1e-10, 1.0-1e-10) / clamp(1.0 - mean(Y), 1e-10, 1.0))

    sqrt2 = sqrt(2.0)

    for iter in 1:max_iter
        old_alpha = copy(alpha)
        old_beta  = copy(beta)
        old_gamma = copy(gamma)
        old_norm_sq = sum(alpha.^2) + sum(beta.^2) + sum(gamma.^2)

        # --- Mediator gradients: ∂(L_med/L_med_null)/∂α_k ---
        grad_alpha = zeros(K)
        for k in 1:K
            s = 0.0; sx = 0.0
            @inbounds for i in 1:n
                r = S[i, k] - alpha0[k] - alpha[k] * X_bar[i, k]
                s += r; sx += r * X_bar[i, k]
            end
            alpha0[k] += eta * s / (n * L_med_null)
            grad_alpha[k] = -sx / (n * L_med_null)
        end

        # --- Outcome gradients using residual S_tilde ---
        grad_beta  = zeros(K)
        grad_gamma = zeros(K)
        s_resid = 0.0
        @inbounds for i in 1:n
            lp = beta0
            for k in 1:K
                s_tilde_ik = S[i, k] - alpha0[k] - alpha[k] * X_bar[i, k]
                lp += beta[k] * X_bar[i, k] + gamma[k] * s_tilde_ik
            end
            p = 1.0 / (1.0 + exp(-lp))
            r = p - Y[i]
            s_resid += r
            for k in 1:K
                s_tilde_ik = S[i, k] - alpha0[k] - alpha[k] * X_bar[i, k]
                grad_beta[k]  += r * X_bar[i, k]
                grad_gamma[k] += r * s_tilde_ik
            end
        end
        beta0 -= eta * s_resid / (n * L_out_null)
        grad_beta  ./= (n * L_out_null)
        grad_gamma ./= (n * L_out_null)

        # --- Coupling gradient: ∂(L_out/L_out^0)/∂(α_k, α_{0k}) ---
        for k in 1:K
            grad_alpha[k] -= gamma[k] * grad_beta[k]
            alpha0[k] += eta * gamma[k] * s_resid / (n * L_out_null)
        end

        # --- Proximal step ---
        # α: FREE (gradient step only, no shrinkage)
        for k in 1:K
            alpha[k] -= eta * grad_alpha[k]
        end

        # (β,γ): Group LASSO proximal with norm = sqrt(β² + γ²)
        diff_sq = 0.0
        threshold = eta * lambda * sqrt2
        for k in 1:K
            b0 = beta[k]; g0 = gamma[k]
            bt = b0 - eta * grad_beta[k]
            gt = g0 - eta * grad_gamma[k]
            nt = sqrt(bt^2 + gt^2)

            if nt <= threshold
                beta[k] = 0.0; gamma[k] = 0.0
            else
                sh = 1.0 - threshold / nt
                beta[k] = sh * bt; gamma[k] = sh * gt
            end
            diff_sq += (alpha[k]-old_alpha[k])^2 + (beta[k]-b0)^2 + (gamma[k]-g0)^2
        end

        if sqrt(diff_sq) / (sqrt(old_norm_sq) + 1e-10) < tol
            return (alpha=copy(alpha), beta=copy(beta), gamma=copy(gamma),
                    alpha0=copy(alpha0), beta0=beta0, iter=iter, converged=true)
        end
    end
    return (alpha=copy(alpha), beta=copy(beta), gamma=copy(gamma),
            alpha0=copy(alpha0), beta0=beta0, iter=max_iter, converged=false)
end


# ============================================================
# 2. Outcome deviance (for CV evaluation)
# ============================================================

function outcome_deviance(X_bar, S, Y, fit)
    n = length(Y)
    K = size(X_bar, 2)
    dev = 0.0
    @inbounds for i in 1:n
        lp = fit.beta0
        for k in 1:K
            s_tilde_ik = S[i, k] - fit.alpha0[k] - fit.alpha[k] * X_bar[i, k]
            lp += fit.beta[k] * X_bar[i, k] + fit.gamma[k] * s_tilde_ik
        end
        p = clamp(1.0 / (1.0 + exp(-lp)), 1e-10, 1.0 - 1e-10)
        dev -= 2.0 * (Y[i] * log(p) + (1.0 - Y[i]) * log(1.0 - p))
    end
    return dev / n
end


# ============================================================
# 3. Data-driven lambda_max
# ============================================================

function compute_lambda_max(X_bar, S, Y)
    n, K = size(X_bar)
    y_bar = clamp(mean(Y), 1e-10, 1.0 - 1e-10)
    L_out_null = -(y_bar * log(y_bar) + (1.0 - y_bar) * log(1.0 - y_bar))
    S_centered = S .- mean(S, dims=1)
    r = fill(y_bar, n) .- Y
    max_norm = 0.0
    for k in 1:K
        gb = 0.0; gg = 0.0
        @inbounds for i in 1:n
            gb += r[i] * X_bar[i, k]
            gg += r[i] * S_centered[i, k]
        end
        gb /= (n * L_out_null); gg /= (n * L_out_null)
        max_norm = max(max_norm, sqrt(gb^2 + gg^2))
    end
    return max_norm / sqrt(2.0)
end

function make_lambda_grid(X_bar, S, Y; n_lambda=100, lambda_ratio=0.01)
    lam_max = compute_lambda_max(X_bar, S, Y)
    lam_min = lam_max * lambda_ratio
    return exp.(range(log(lam_max), log(lam_min), length=n_lambda))
end


# ============================================================
# 4. CV-based lambda selection
# ============================================================

function cv_modified_primed(X_bar, S, Y, lambda_grid; n_folds=5, seed=42)
    n = length(Y)
    K = size(X_bar, 2)
    n_lam = length(lambda_grid)

    rng = MersenneTwister(seed)
    fold_ids = Vector{Int}(undef, n)
    idx1 = findall(Y .== 1.0)
    idx0 = findall(Y .== 0.0)
    fold_ids[idx1] = mod1.(shuffle(rng, 1:length(idx1)), n_folds)
    fold_ids[idx0] = mod1.(shuffle(rng, 1:length(idx0)), n_folds)

    cv_dev = zeros(n_lam)

    Threads.@threads for j in 1:n_lam
        lam = lambda_grid[j]
        total_dev = 0.0

        for fold in 1:n_folds
            test_mask  = fold_ids .== fold
            train_mask = .!test_mask

            fit = modified_primed(X_bar[train_mask, :], S[train_mask, :], Y[train_mask], lam)
            total_dev += outcome_deviance(X_bar[test_mask, :], S[test_mask, :], Y[test_mask], fit)
        end
        cv_dev[j] = total_dev / n_folds
    end

    best_idx = argmin(cv_dev)
    best_fit = modified_primed(X_bar, S, Y, lambda_grid[best_idx])

    return (fit=best_fit, lambda=lambda_grid[best_idx], cv_deviance=cv_dev)
end


# ============================================================
# 4. Main
# ============================================================

function main()
    data_csv   = ARGS[1]
    result_csv = ARGS[2]
    n_folds    = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 5

    raw = readdlm(data_csv, ','; header=true)
    data = raw[1]
    reps = unique(data[:, 1])
    B = length(reps)

    ncols = size(data, 2)
    K = div(ncols - 3, 2)

    nt = Threads.nthreads()
    @printf("=== Modified PRIMED with %d-fold CV (Julia, %d threads) ===\n", n_folds, nt)
    @printf("B=%d, K=%d, data-driven lambda grid\n", B, K)
    @printf("L_joint = L_med/L_med_null + L_out/L_out_null  [α free, penalty on (β,γ)]\n\n")

    n_out_cols = 3 + K + 4*K
    results = zeros(B, n_out_cols)

    t0 = time()

    for (idx, rep_id) in enumerate(reps)
        rows = data[:, 1] .== rep_id
        block = data[rows, :]

        X_bar = block[:, 3:(2+K)]
        S     = block[:, (3+K):(2+2K)]
        Y     = block[:, end]

        lambda_grid = make_lambda_grid(X_bar, S, Y; n_lambda=100, lambda_ratio=0.01)
        res = cv_modified_primed(X_bar, S, Y, lambda_grid; n_folds=n_folds, seed=42+idx)

        results[idx, 1] = rep_id
        results[idx, 2] = res.lambda
        results[idx, 3] = res.fit.beta0
        for k in 1:K
            # Selection: based on β and γ only (NOT α)
            sel = (res.fit.beta[k] != 0 || res.fit.gamma[k] != 0) ? 1.0 : 0.0
            results[idx, 3+k]      = sel
            results[idx, 3+K+k]    = res.fit.alpha[k]
            results[idx, 3+2K+k]   = res.fit.alpha0[k]
            results[idx, 3+3K+k]   = res.fit.beta[k]
            results[idx, 3+4K+k]   = res.fit.gamma[k]
        end

        if idx % max(1, B÷5) == 0
            @printf("  Rep %d/%d done (lambda=%.3f, selected=%d)\n",
                    idx, B, res.lambda, count(results[idx, 4:(3+K)] .> 0))
        end
    end

    elapsed = time() - t0
    @printf("\nDone in %.1f sec (%.2f sec/rep)\n", elapsed, elapsed/B)

    header = ["rep" "lambda" "beta0" [string("sel_", k) for k in 1:K]... [string("alpha_", k) for k in 1:K]... [string("alpha0_", k) for k in 1:K]... [string("beta_", k) for k in 1:K]... [string("gamma_", k) for k in 1:K]...]

    open(result_csv, "w") do io
        println(io, join(header, ","))
        for i in 1:B
            println(io, join([@sprintf("%.6f", results[i, j]) for j in 1:n_out_cols], ","))
        end
    end

    @printf("Results saved to: %s\n", result_csv)
end

main()
