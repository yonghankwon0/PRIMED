#!/bin/bash
cd "$(dirname "$0")"
TOTAL_JOBS=40
MAX_PARALLEL=8
JULIA_THREADS=8
mkdir -p logs
rm -f result_full_job*.csv logs/full_job*.log

echo "============================================================"
echo "  Full Simulation with SD: $TOTAL_JOBS jobs"
echo "  $(date)"
echo "  Strategy: $MAX_PARALLEL parallel × $JULIA_THREADS Julia threads"
echo "============================================================"

declare -a PIDS
declare -a JOB_IDS
RUNNING=0
COMPLETED=0
FAILED=0
START_TIME=$(date +%s)

launch_job() {
    local job_id=$1
    echo "[$(date +%H:%M:%S)] Launching job $job_id/$TOTAL_JOBS"
    Rscript sim_full_single_job.R "$job_id" "$JULIA_THREADS" \
        > "logs/full_job$(printf '%02d' $job_id).log" 2>&1 &
    PIDS+=($!)
    JOB_IDS+=($job_id)
    RUNNING=$((RUNNING + 1))
}

wait_for_slot() {
    while [ $RUNNING -ge $MAX_PARALLEL ]; do
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                wait "${PIDS[$i]}"
                local exit_code=$?
                local jid=${JOB_IDS[$i]}
                if [ $exit_code -eq 0 ] && [ -f "result_full_job$(printf '%02d' $jid).csv" ]; then
                    echo "[$(date +%H:%M:%S)] Job $jid COMPLETED"
                    COMPLETED=$((COMPLETED + 1))
                else
                    echo "[$(date +%H:%M:%S)] Job $jid FAILED"
                    FAILED=$((FAILED + 1))
                fi
                unset 'PIDS[i]'
                unset 'JOB_IDS[i]'
                PIDS=("${PIDS[@]}")
                JOB_IDS=("${JOB_IDS[@]}")
                RUNNING=$((RUNNING - 1))
                return
            fi
        done
        sleep 2
    done
}

for job_id in $(seq 1 $TOTAL_JOBS); do
    wait_for_slot
    launch_job $job_id
done

echo ""
echo "All jobs launched. Waiting for remaining $RUNNING..."
for i in "${!PIDS[@]}"; do
    wait "${PIDS[$i]}"
    local_exit=$?
    jid=${JOB_IDS[$i]}
    if [ $local_exit -eq 0 ] && [ -f "result_full_job$(printf '%02d' $jid).csv" ]; then
        echo "[$(date +%H:%M:%S)] Job $jid COMPLETED"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "[$(date +%H:%M:%S)] Job $jid FAILED"
        FAILED=$((FAILED + 1))
    fi
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "============================================================"
echo "  Done! Completed: $COMPLETED/$TOTAL_JOBS Failed: $FAILED"
echo "  Time: $((ELAPSED/60))m $((ELAPSED%60))s"
echo "============================================================"

if [ $COMPLETED -gt 0 ]; then
    Rscript -e '
    library(dplyr, warn.conflicts=FALSE)

    # Merge main results
    files <- sort(Sys.glob("result_full_job*.csv"))
    cat(sprintf("Merging %d result files\n", length(files)))
    df <- bind_rows(lapply(files, read.csv))
    num_cols <- c("TPR_mean","TPR_sd","PPV_mean","PPV_sd","F1_mean","F1_sd","AUC_mean","AUC_sd",
                  "mp_nsel_mean","mp_nsel_sd","cd_nsel_mean","cd_nsel_sd","cl_nsel_mean","cl_nsel_sd",
                  "mp_lambda_mean","mp_lambda_sd")
    for(col in num_cols) if(col %in% names(df)) df[[col]] <- as.numeric(df[[col]])
    write.csv(df, "sim_all_with_sd_results.csv", row.names=FALSE)
    cat(sprintf("Saved: sim_all_with_sd_results.csv (%d rows)\n", nrow(df)))

    # Print summary with selection info
    df %>% filter(method=="PRIMED") %>%
      arrange(sim,dgp,scenario,tau,snr,K) %>%
      select(sim,dgp,scenario,tau,snr,K,F1_mean,F1_sd,AUC_mean,AUC_sd,
             mp_nsel_mean,mp_nsel_sd,mp_lambda_mean,mp_lambda_sd) %>%
      as.data.frame() %>% print(row.names=FALSE)

    # Merge selection frequency
    sf_files <- sort(Sys.glob("selfreq_full_job*.csv"))
    if(length(sf_files) > 0) {
      cat(sprintf("\nMerging %d selfreq files\n", length(sf_files)))
      sf_list <- lapply(seq_along(sf_files), function(i) {
        d <- read.csv(sf_files[i])
        d$job_id <- i
        d
      })
      sf_all <- bind_rows(sf_list)
      write.csv(sf_all, "sim_all_selfreq.csv", row.names=FALSE)
      cat("Saved: sim_all_selfreq.csv\n")
    }
    '
fi
