#!/usr/bin/env Rscript
# ==============================================================================
# FOREST SUCCESSION MODEL FITTING - INLAjoint
# Method B: Interval-censored multi-state survival models
# ==============================================================================

library(tidyverse)
library(data.table)
library(INLA)
library(INLAjoint)
library(here)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths
RECON_DATA_PATH <- here("Data-Output", "INLA-realData", "forest_reconstruction_complete.rds")
OUTPUT_DIR <- here("Data-Output", "INLA-realResults")
PROGRESS_FILE <- file.path(OUTPUT_DIR, "fitting_log.txt")

# Model parameters
N_TRANSITIONS <- 3        # How many transitions to fit (top N by event count)
MIN_EVENTS <- 10000        # Minimum events required per transition
INCLUDE_SOIL <- FALSE      # Include soil as covariate? (risky - many levels)
INCLUDE_PERT <- FALSE      # Include perturbation as covariate?

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PROGRESS LOGGING
# ==============================================================================

log_progress <- function(message, progress_file = PROGRESS_FILE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- sprintf("[%s] %s\n", timestamp, message)
  cat(log_message)
  cat(log_message, file = progress_file, append = TRUE)
}

cat("", file = PROGRESS_FILE)  # Clear file
log_progress("=== FOREST SUCCESSION MODEL FITTING - STARTED ===")

# ==============================================================================
# LOAD RECONSTRUCTED DATA
# ==============================================================================

log_progress("Loading reconstructed data...")
start_load <- Sys.time()

model_data <- readRDS(RECON_DATA_PATH)

end_load <- Sys.time()
log_progress(sprintf("Data loaded! Time: %.1f seconds", 
                     as.numeric(difftime(end_load, start_load, units = "secs"))))
log_progress(sprintf("Dimensions: %d rows, %d columns", 
                     nrow(model_data), ncol(model_data)))
log_progress(sprintf("Unique plots: %d", length(unique(model_data$id))))
log_progress(sprintf("Total transitions: %d", sum(model_data$status == 3)))

# ==============================================================================
# SELECT TRANSITIONS TO FIT
# ==============================================================================

log_progress("\n=== SELECTING TRANSITIONS ===")

# Calculate transition counts
transition_summary <- model_data %>%
  filter(status == 3) %>%
  group_by(from, to) %>%
  summarise(n_events = n(), .groups = "drop") %>%
  arrange(desc(n_events)) %>%
  mutate(transition = paste0(from, "→", to))

# Apply filters
transitions_to_fit <- transition_summary %>%
  filter(n_events >= MIN_EVENTS) %>%
  head(N_TRANSITIONS)

log_progress(sprintf("Transitions with ≥%d events: %d", 
                     MIN_EVENTS, 
                     nrow(transition_summary %>% filter(n_events >= MIN_EVENTS))))
log_progress(sprintf("Fitting top %d transitions", nrow(transitions_to_fit)))

# Save transition summary
write.csv(transition_summary, 
          file.path(OUTPUT_DIR, "all_transition_counts.csv"), 
          row.names = FALSE)
write.csv(transitions_to_fit, 
          file.path(OUTPUT_DIR, "transitions_to_fit.csv"), 
          row.names = FALSE)

log_progress("\nTop 10 transitions by event count:")
for(i in 1:min(10, nrow(transitions_to_fit))) {
  log_progress(sprintf("  %s: %d events", 
                       transitions_to_fit$transition[i], 
                       transitions_to_fit$n_events[i]))
}

# ==============================================================================
# PREPARE DATA FOR INLA
# ==============================================================================

log_progress("\n=== PREPARING DATA FOR INLAjoint ===")

# Convert to data.frame if needed (INLAjoint prefers data.frame)
if(is.data.table(model_data)) {
  model_data <- as.data.frame(model_data)
}

# Create list of datasets, one per transition
log_progress("Creating data list for each transition...")
data_list <- list()

for(i in 1:nrow(transitions_to_fit)) {
  from_s <- transitions_to_fit$from[i]
  to_s <- transitions_to_fit$to[i]
  
  data_list[[i]] <- model_data %>%
    filter(from == from_s, to == to_s)
  
  if(i %% 10 == 0) {
    log_progress(sprintf("  Prepared %d/%d datasets", i, nrow(transitions_to_fit)))
  }
}

log_progress(sprintf("Data list created: %d transitions", length(data_list)))

# Check for missing covariates
log_progress("Checking for missing covariate data...")
for(i in 1:length(data_list)) {
  n_missing_cmi <- sum(is.na(data_list[[i]]$cov_CMI_std))
  n_missing_tmean <- sum(is.na(data_list[[i]]$cov_Tmean_std))
  
  if(n_missing_cmi > 0 || n_missing_tmean > 0) {
    log_progress(sprintf("  WARNING: Transition %d has missing covariates (CMI: %d, Tmean: %d)",
                         i, n_missing_cmi, n_missing_tmean))
  }
}

# ==============================================================================
# CREATE SURVIVAL OBJECTS
# ==============================================================================

log_progress("\n=== CREATING SURVIVAL OBJECTS ===")
log_progress(sprintf("Creating %d survival objects...", length(data_list)))

start_surv <- Sys.time()

for(i in 1:length(data_list)) {
  d <- data_list[[i]]
  
  # Time vector: Tstart for interval-censored (status=3), Tstop for right-censored (status=0)
  time_vec <- d$Tstart
  idx_censored <- d$status == 0
  time_vec[idx_censored] <- d$Tstop[idx_censored]
  
  # Create survival object
  assign(paste0("s", i), 
         inla.surv(time = time_vec, 
                   time2 = d$Tstop, 
                   event = d$status, 
                   truncation = d$approx_entry),
         envir = .GlobalEnv)
  
  if(i %% 10 == 0) {
    log_progress(sprintf("  Created %d/%d survival objects", i, length(data_list)))
  }
}

end_surv <- Sys.time()
log_progress(sprintf("Survival objects created! Time: %.1f seconds", 
                     as.numeric(difftime(end_surv, start_surv, units = "secs"))))

# ==============================================================================
# BUILD FORMULAS
# ==============================================================================

log_progress("\n=== BUILDING FORMULAS ===")

# Decide which covariates to include
covariate_formula <- "cov_CMI_std + cov_Tmean_std"

if(INCLUDE_SOIL) {
  covariate_formula <- paste0(covariate_formula, " + cov_soil")
  log_progress("Including SOIL as covariate (WARNING: may cause convergence issues)")
}

if(INCLUDE_PERT) {
  covariate_formula <- paste0(covariate_formula, " + cov_pert_class")
  log_progress("Including PERTURBATION as covariate")
}

log_progress(sprintf("Formula: s ~ %s", covariate_formula))

# Create formula list
formulas <- lapply(1:nrow(transitions_to_fit), function(i) {
  as.formula(paste0("s", i, " ~ ", covariate_formula))
})

log_progress(sprintf("Created %d formulas", length(formulas)))

# ==============================================================================
# FIT MODEL
# ==============================================================================

log_progress("\n=== FITTING MODEL ===")
log_progress(sprintf("Transitions: %d", nrow(transitions_to_fit)))
log_progress(sprintf("Parameters per transition: %d (α + β₀ + covariates)", 
                     1 + str_count(covariate_formula, "\\+") + 1))
log_progress(sprintf("Total parameters: ~%d", 
                     nrow(transitions_to_fit) * (2 + str_count(covariate_formula, "\\+") + 1)))
log_progress("This may take 1-4 hours depending on complexity...")
log_progress(sprintf("Started: %s", Sys.time()))

start_fit <- Sys.time()

# Fit the model
fit_model <- tryCatch({
  joint(
    formSurv = formulas,
    basRisk = rep("weibullsurv", nrow(transitions_to_fit)),
    dataSurv = data_list,
    id = "id",
    control = list(config = TRUE, variant = 1)
  )
}, error = function(e) {
  log_progress(sprintf("ERROR during fitting: %s", e$message))
  return(NULL)
})

end_fit <- Sys.time()
fit_time <- as.numeric(difftime(end_fit, start_fit, units = "mins"))

if(is.null(fit_model)) {
  log_progress("MODEL FITTING FAILED!")
  log_progress("Check error messages above")
  quit(status = 1)
}

log_progress(sprintf("Model fitted successfully! Time: %.1f minutes (%.1f hours)", 
                     fit_time, fit_time / 60))

# ==============================================================================
# CHECK FOR CONVERGENCE ISSUES
# ==============================================================================

log_progress("\n=== CHECKING CONVERGENCE ===")

# Check for warnings in model
if(!is.null(fit_model$mode$convergence)) {
  if(fit_model$mode$convergence != 0) {
    log_progress(sprintf("WARNING: Convergence code = %d", fit_model$mode$convergence))
  } else {
    log_progress("Convergence: OK")
  }
}

# Check parameter uncertainties
n_params_per_trans <- 2 + str_count(covariate_formula, "\\+") + 1
n_fixed <- nrow(transitions_to_fit) * (1 + str_count(covariate_formula, "\\+") + 1)

large_sd_fixed <- sum(fit_model$summary.fixed$sd > 1, na.rm = TRUE)
large_sd_hyper <- sum(fit_model$summary.hyperpar$sd > 0.5, na.rm = TRUE)

log_progress(sprintf("Fixed effects with SD > 1: %d / %d", large_sd_fixed, n_fixed))
log_progress(sprintf("Shape parameters with SD > 0.5: %d / %d", 
                     large_sd_hyper, nrow(transitions_to_fit)))

if(large_sd_fixed > n_fixed * 0.1) {
  log_progress("WARNING: >10% of fixed effects have large uncertainty")
  log_progress("Consider simplifying model (fewer transitions, fewer covariates)")
}

# ==============================================================================
# EXTRACT RESULTS
# ==============================================================================

log_progress("\n=== EXTRACTING RESULTS ===")

results_table <- data.frame()

for(i in 1:nrow(transitions_to_fit)) {
  from_s <- transitions_to_fit$from[i]
  to_s <- transitions_to_fit$to[i]
  n_events <- transitions_to_fit$n_events[i]
  
  # Get parameter indices
  param_start <- (i - 1) * 3 + 1  # 3 params: intercept, CMI, Tmean
  
  # Extract alpha (shape parameter) - FIX: use [[ ]] for scalar extraction
  alpha_mean <- fit_model$summary.hyperpar[i, "mean"]
  alpha_sd <- fit_model$summary.hyperpar[i, "sd"]
  alpha_lower <- fit_model$summary.hyperpar[i, "0.025quant"]
  alpha_upper <- fit_model$summary.hyperpar[i, "0.975quant"]
  
  # Extract intercept
  beta0_mean <- fit_model$summary.fixed[param_start, "mean"]
  beta0_sd <- fit_model$summary.fixed[param_start, "sd"]
  beta0_lower <- fit_model$summary.fixed[param_start, "0.025quant"]
  beta0_upper <- fit_model$summary.fixed[param_start, "0.975quant"]
  
  # Extract CMI effect
  cmi_mean <- fit_model$summary.fixed[param_start + 1, "mean"]
  cmi_sd <- fit_model$summary.fixed[param_start + 1, "sd"]
  cmi_lower <- fit_model$summary.fixed[param_start + 1, "0.025quant"]
  cmi_upper <- fit_model$summary.fixed[param_start + 1, "0.975quant"]
  cmi_sig <- cmi_lower * cmi_upper > 0
  
  # Extract Tmean effect
  tmean_mean <- fit_model$summary.fixed[param_start + 2, "mean"]
  tmean_sd <- fit_model$summary.fixed[param_start + 2, "sd"]
  tmean_lower <- fit_model$summary.fixed[param_start + 2, "0.025quant"]
  tmean_upper <- fit_model$summary.fixed[param_start + 2, "0.975quant"]
  tmean_sig <- tmean_lower * tmean_upper > 0
  
  results_table <- rbind(results_table, data.frame(
    Transition = paste0(from_s, "→", to_s),
    From = from_s,
    To = to_s,
    N_events = n_events,
    
    # Shape parameter
    Alpha = alpha_mean,
    Alpha_SD = alpha_sd,
    Alpha_lower = alpha_lower,
    Alpha_upper = alpha_upper,
    
    # Intercept
    Beta0 = beta0_mean,
    Beta0_SD = beta0_sd,
    Beta0_lower = beta0_lower,
    Beta0_upper = beta0_upper,
    
    # CMI effect
    CMI_effect = cmi_mean,
    CMI_SD = cmi_sd,
    CMI_lower = cmi_lower,
    CMI_upper = cmi_upper,
    CMI_sig = cmi_sig,
    
    # Temperature effect
    Tmean_effect = tmean_mean,
    Tmean_SD = tmean_sd,
    Tmean_lower = tmean_lower,
    Tmean_upper = tmean_upper,
    Tmean_sig = tmean_sig
  ))
}

log_progress(sprintf("Extracted results for %d transitions", nrow(results_table)))

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

log_progress("\n=== SAVING RESULTS ===")

# Save results table
write.csv(results_table, 
          file.path(OUTPUT_DIR, "model_results.csv"), 
          row.names = FALSE)
log_progress("Saved: model_results.csv")

# Save full model object
saveRDS(fit_model, 
        file.path(OUTPUT_DIR, "fitted_model.rds"),
        compress = TRUE)
log_progress("Saved: fitted_model.rds")

# Save summary statistics
summary_stats <- list(
  n_transitions = nrow(transitions_to_fit),
  n_plots = length(unique(model_data$id)),
  n_total_events = sum(transitions_to_fit$n_events),
  fitting_time_mins = fit_time,
  n_significant_CMI = sum(results_table$CMI_sig),
  n_significant_Tmean = sum(results_table$Tmean_sig),
  convergence_issues = large_sd_fixed > n_fixed * 0.1,
  timestamp = Sys.time()
)

saveRDS(summary_stats, file.path(OUTPUT_DIR, "summary_stats.rds"))
log_progress("Saved: summary_stats.rds")

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

log_progress("\n=== SUMMARY STATISTICS ===")
log_progress(sprintf("Transitions fitted: %d", nrow(results_table)))
log_progress(sprintf("Total events: %d", sum(results_table$N_events)))
log_progress(sprintf("Significant CMI effects: %d / %d (%.1f%%)", 
                     sum(results_table$CMI_sig), 
                     nrow(results_table),
                     100 * sum(results_table$CMI_sig) / nrow(results_table)))
log_progress(sprintf("Significant Tmean effects: %d / %d (%.1f%%)", 
                     sum(results_table$Tmean_sig), 
                     nrow(results_table),
                     100 * sum(results_table$Tmean_sig) / nrow(results_table)))

log_progress("\nMean absolute effects:")
log_progress(sprintf("  CMI: %.3f", mean(abs(results_table$CMI_effect))))
log_progress(sprintf("  Tmean: %.3f", mean(abs(results_table$Tmean_effect))))

log_progress("\nTop 5 strongest CMI effects:")
top_cmi <- results_table %>% 
  arrange(desc(abs(CMI_effect))) %>% 
  head(5) %>%
  select(Transition, CMI_effect, CMI_sig)
for(i in 1:nrow(top_cmi)) {
  sig_mark <- ifelse(top_cmi$CMI_sig[i], "*", "")
  log_progress(sprintf("  %s: %.3f%s", 
                       top_cmi$Transition[i], 
                       top_cmi$CMI_effect[i],
                       sig_mark))
}

log_progress("\nTop 5 strongest Temperature effects:")
top_tmean <- results_table %>% 
  arrange(desc(abs(Tmean_effect))) %>% 
  head(5) %>%
  select(Transition, Tmean_effect, Tmean_sig)
for(i in 1:nrow(top_tmean)) {
  sig_mark <- ifelse(top_tmean$Tmean_sig[i], "*", "")
  log_progress(sprintf("  %s: %.3f%s", 
                       top_tmean$Transition[i], 
                       top_tmean$Tmean_effect[i],
                       sig_mark))
}

# ==============================================================================
# FINISH
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_load, units = "mins"))

log_progress("\n=================================================================")
log_progress("MODEL FITTING COMPLETE")
log_progress(sprintf("Total runtime: %.1f minutes (%.1f hours)", total_time, total_time / 60))
log_progress(sprintf("Output directory: %s", OUTPUT_DIR))
log_progress("=================================================================")

cat("\n\n")
cat("##############################################\n")
cat("# MODEL FITTING COMPLETE\n")
cat(sprintf("# Transitions: %d\n", nrow(results_table)))
cat(sprintf("# Significant effects: CMI=%d, Tmean=%d\n", 
            sum(results_table$CMI_sig), sum(results_table$Tmean_sig)))
cat(sprintf("# Time: %.1f hours\n", total_time / 60))
cat(sprintf("# Results: %s/model_results.csv\n", OUTPUT_DIR))
cat("##############################################\n")
