#!/usr/bin/env Rscript
# ==============================================================================
# FOREST SUCCESSION MODEL FITTING - INDIVIDUAL TRANSITIONS
# Method B: Interval-censored multi-state survival models
# Fits each transition separately for efficiency and parallelization
# ==============================================================================

library(tidyverse)
library(data.table)
library(INLA)
library(here)
library(parallel)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths
RECON_DATA_PATH <- here("Data-Output", "INLA-realData", "forest_reconstruction_complete.rds")
OUTPUT_DIR <- here("Data-Output", "INLA-realResults-IND")
PROGRESS_FILE <- file.path(OUTPUT_DIR, "fitting_log_ind.txt")

# Model parameters
N_TRANSITIONS <- 8      # How many transitions to fit (top N by event count)
MIN_EVENTS <- 1000       # Minimum events required per transition
INCLUDE_SOIL <- T     # Include soil as covariate? (risky - many levels)
INCLUDE_PERT <- T     # Include perturbation as covariate?

# Parallel processing
USE_PARALLEL <- TRUE
N_CORES <- 8            # Adjust for your cluster

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
log_progress("Fitting transitions INDIVIDUALLY (not jointly)")

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
# PREPARE COVARIATES
# ==============================================================================

log_progress("\n=== PREPARING COVARIATES ===")

# Filter: Keep only no perturbation (0) or pest outbreak (3)
n_before <- nrow(model_data)
model_data <- model_data %>%
  filter(cov_pert_class %in% c(0, 3))
n_after <- nrow(model_data)

log_progress(sprintf("Filtered perturbation data: %d → %d rows (%.1f%% retained)",
                     n_before, n_after, 100 * n_after / n_before))

# Create binary pest variable
model_data$pert_binary <- ifelse(model_data$cov_pert_class == 3, 1, 0)
log_progress(sprintf("Pest outbreak events: %d (%.2f%%)", 
                     sum(model_data$pert_binary), 
                     100 * mean(model_data$pert_binary)))

# Create 3-level soil classification
model_data$soil_class <- case_when(
  model_data$cov_soil %in% c(0, 1, 2) ~ "Dry",
  model_data$cov_soil %in% c(3, 4, 5, 6) ~ "Optimal",
  model_data$cov_soil %in% c(7, 8, 9) ~ "Wet",
  TRUE ~ NA_character_
)

# Set Optimal as reference level
model_data$soil_class <- factor(model_data$soil_class, 
                                levels = c("Optimal", "Dry", "Wet"))

log_progress("Soil distribution:")
soil_table <- table(model_data$soil_class, useNA = "ifany")
for(level in names(soil_table)) {
  log_progress(sprintf("  %s: %d (%.1f%%)", 
                       level, soil_table[level], 
                       100 * soil_table[level] / sum(soil_table)))
}

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

log_progress("\nTransitions to fit:")
for(i in 1:nrow(transitions_to_fit)) {
  log_progress(sprintf("  %s: %d events", 
                       transitions_to_fit$transition[i], 
                       transitions_to_fit$n_events[i]))
}

# ==============================================================================
# PREPARE DATA FOR INLA
# ==============================================================================

log_progress("\n=== PREPARING DATA FOR INLA ===")

# Convert to data.frame if needed
if(is.data.table(model_data)) {
  model_data <- as.data.frame(model_data)
}

# Build covariate formula
covariate_formula <- "cov_CMI_std + cov_Tmean_std + soil_class + pert_binary"

log_progress(sprintf("Formula: surv ~ %s", covariate_formula))
log_progress("Covariates included:")
log_progress("  - Climate: CMI (standardized), Temperature (standardized)")
log_progress("  - Soil: 3-level (Dry/Optimal/Wet)")
log_progress("  - Perturbation: Binary (Pest vs None)")

log_progress(sprintf("Formula: surv ~ %s", covariate_formula))

# Check for missing covariates
log_progress("Checking for missing covariate data...")
n_missing_cmi <- sum(is.na(model_data$cov_CMI_std))
n_missing_tmean <- sum(is.na(model_data$cov_Tmean_std))

if(n_missing_cmi > 0 || n_missing_tmean > 0) {
  log_progress(sprintf("  WARNING: Missing covariates (CMI: %d, Tmean: %d)",
                       n_missing_cmi, n_missing_tmean))
  log_progress("  Rows with missing data will be removed per transition")
}

# ==============================================================================
# FUNCTION TO FIT SINGLE TRANSITION
# ==============================================================================

fit_single_transition <- function(from_state, to_state, data, 
                                  covariate_formula, 
                                  transition_label,
                                  output_dir,
                                  variant = 1) {
  
  # Extract data for this transition
  trans_data <- data %>%
    filter(from == from_state, to == to_state)
  
  # Remove rows with missing covariates
  # Remove rows with missing covariates
  trans_data <- trans_data %>%
    filter(!is.na(cov_CMI_std), !is.na(cov_Tmean_std), !is.na(soil_class))
  
  n_events <- sum(trans_data$status == 3)
  n_censored <- sum(trans_data$status == 0)
  
  # Check for sufficient events
  if(n_events < 10) {
    warning(sprintf("Transition %s has <10 events after removing NAs - skipping", 
                    transition_label))
    return(NULL)
  }
  
  # Create survival object
  time_vec <- trans_data$Tstart
  time_vec[trans_data$status == 0] <- trans_data$Tstop[trans_data$status == 0]
  
  surv_obj <- inla.surv(
    time = time_vec,
    time2 = trans_data$Tstop,
    event = trans_data$status,
    truncation = trans_data$approx_entry
  )
  
  # Build formula
  formula <- as.formula(paste0("surv_obj ~ ", covariate_formula))
  
  # Fit model
  fit <- tryCatch({
    inla(
      formula = formula,
      data = trans_data,
      family = "weibullsurv",
      control.family = list(variant = variant),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
      control.predictor = list(compute = TRUE)
    )
  }, error = function(e) {
    warning(sprintf("Transition %s failed: %s", transition_label, e$message))
    return(NULL)
  })
  
  if(is.null(fit)) return(NULL)
  
  # Check convergence
  convergence_ok <- TRUE
  if(!is.null(fit$mode$convergence)) {
    if(fit$mode$convergence != 0) {
      warning(sprintf("Transition %s: convergence code = %d", 
                      transition_label, fit$mode$convergence))
      convergence_ok <- FALSE
    }
  }
  
  # Extract results
  results <- list(
    transition = transition_label,
    from = from_state,
    to = to_state,
    n_rows = nrow(trans_data),
    n_events = n_events,
    n_censored = n_censored,
    convergence_ok = convergence_ok,
    
    # Shape parameter (alpha)
    alpha = fit$summary.hyperpar["alpha for weibullsurv", "mean"],
    alpha_sd = fit$summary.hyperpar["alpha for weibullsurv", "sd"],
    alpha_lower = fit$summary.hyperpar["alpha for weibullsurv", "0.025quant"],
    alpha_upper = fit$summary.hyperpar["alpha for weibullsurv", "0.975quant"],
    
    # Intercept
    beta0 = fit$summary.fixed["(Intercept)", "mean"],
    beta0_sd = fit$summary.fixed["(Intercept)", "sd"],
    beta0_lower = fit$summary.fixed["(Intercept)", "0.025quant"],
    beta0_upper = fit$summary.fixed["(Intercept)", "0.975quant"],
    
    # CMI effect
    cmi_effect = fit$summary.fixed["cov_CMI_std", "mean"],
    cmi_sd = fit$summary.fixed["cov_CMI_std", "sd"],
    cmi_lower = fit$summary.fixed["cov_CMI_std", "0.025quant"],
    cmi_upper = fit$summary.fixed["cov_CMI_std", "0.975quant"],
    
    # Temperature effect
    tmean_effect = fit$summary.fixed["cov_Tmean_std", "mean"],
    tmean_sd = fit$summary.fixed["cov_Tmean_std", "sd"],
    tmean_lower = fit$summary.fixed["cov_Tmean_std", "0.025quant"],
    tmean_upper = fit$summary.fixed["cov_Tmean_std", "0.975quant"],
    
    # Soil effects
    soil_dry_effect = fit$summary.fixed["soil_classDry", "mean"],
    soil_dry_sd = fit$summary.fixed["soil_classDry", "sd"],
    soil_dry_lower = fit$summary.fixed["soil_classDry", "0.025quant"],
    soil_dry_upper = fit$summary.fixed["soil_classDry", "0.975quant"],
    
    soil_wet_effect = fit$summary.fixed["soil_classWet", "mean"],
    soil_wet_sd = fit$summary.fixed["soil_classWet", "sd"],
    soil_wet_lower = fit$summary.fixed["soil_classWet", "0.025quant"],
    soil_wet_upper = fit$summary.fixed["soil_classWet", "0.975quant"],
    
    # Pest effect
    pest_effect = fit$summary.fixed["pert_binary", "mean"],
    pest_sd = fit$summary.fixed["pert_binary", "sd"],
    pest_lower = fit$summary.fixed["pert_binary", "0.025quant"],
    pest_upper = fit$summary.fixed["pert_binary", "0.975quant"],
    
    # Model fit
    dic = fit$dic$dic,
    waic = fit$waic$waic,
    
    # Full model object
    model = fit
  )
  
  # Save individual model
  saveRDS(results, 
          file.path(output_dir, sprintf("fit_%s.rds", 
                                        gsub("→", "_to_", transition_label))))
  
  return(results)
}

# ==============================================================================
# FIT ALL TRANSITIONS
# ==============================================================================

log_progress("\n=== FITTING MODELS ===")
log_progress(sprintf("Transitions: %d", nrow(transitions_to_fit)))
log_progress(sprintf("Parameters per transition: %d (α + β₀ + covariates)", 
                     1 + str_count(covariate_formula, "\\+") + 1))
log_progress(sprintf("Mode: %s", ifelse(USE_PARALLEL, 
                                        sprintf("Parallel (%d cores)", N_CORES),
                                        "Sequential")))
log_progress(sprintf("Started: %s", Sys.time()))

start_fit <- Sys.time()

if(USE_PARALLEL) {
  log_progress(sprintf("Initializing parallel cluster with %d cores...", N_CORES))
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, c("model_data", "covariate_formula", "OUTPUT_DIR",
                      "fit_single_transition", "transitions_to_fit"))
  clusterEvalQ(cl, {
    library(INLA)
    library(tidyverse)
  })
  
  results_list <- parLapply(cl, 1:nrow(transitions_to_fit), function(i) {
    fit_single_transition(
      from_state = transitions_to_fit$from[i],
      to_state = transitions_to_fit$to[i],
      data = model_data,
      covariate_formula = covariate_formula,
      transition_label = transitions_to_fit$transition[i],
      output_dir = OUTPUT_DIR,
      variant = 1
    )
  })
  
  stopCluster(cl)
  
} else {
  log_progress("Fitting transitions sequentially...")
  
  results_list <- list()
  for(i in 1:nrow(transitions_to_fit)) {
    log_progress(sprintf("\nFitting transition %d/%d: %s", 
                         i, nrow(transitions_to_fit), 
                         transitions_to_fit$transition[i]))
    
    results_list[[i]] <- fit_single_transition(
      from_state = transitions_to_fit$from[i],
      to_state = transitions_to_fit$to[i],
      data = model_data,
      covariate_formula = covariate_formula,
      transition_label = transitions_to_fit$transition[i],
      output_dir = OUTPUT_DIR,
      variant = 1
    )
    
    if(!is.null(results_list[[i]])) {
      log_progress(sprintf("  Alpha = %.3f (95%% CI: %.3f - %.3f)", 
                           results_list[[i]]$alpha,
                           results_list[[i]]$alpha_lower,
                           results_list[[i]]$alpha_upper))
      log_progress(sprintf("  CMI effect = %.3f (95%% CI: %.3f - %.3f)",
                           results_list[[i]]$cmi_effect,
                           results_list[[i]]$cmi_lower,
                           results_list[[i]]$cmi_upper))
    }
  }
}

end_fit <- Sys.time()
fit_time <- as.numeric(difftime(end_fit, start_fit, units = "mins"))

# Remove NULL results (failed fits)
results_list <- results_list[!sapply(results_list, is.null)]

log_progress(sprintf("\nModel fitting complete! Time: %.1f minutes (%.1f hours)", 
                     fit_time, fit_time / 60))
log_progress(sprintf("Successfully fitted: %d/%d transitions", 
                     length(results_list), nrow(transitions_to_fit)))

if(length(results_list) == 0) {
  log_progress("ERROR: No transitions fitted successfully!")
  quit(status = 1)
}

# ==============================================================================
# CHECK FOR CONVERGENCE ISSUES
# ==============================================================================

log_progress("\n=== CHECKING CONVERGENCE ===")

n_convergence_issues <- sum(!sapply(results_list, function(x) x$convergence_ok))
log_progress(sprintf("Transitions with convergence issues: %d / %d",
                     n_convergence_issues, length(results_list)))

# Check parameter uncertainties
large_sd_alpha <- sum(sapply(results_list, function(x) x$alpha_sd > 0.5))
large_sd_cmi <- sum(sapply(results_list, function(x) x$cmi_sd > 1))
large_sd_tmean <- sum(sapply(results_list, function(x) x$tmean_sd > 1))

log_progress(sprintf("Shape parameters with SD > 0.5: %d / %d", 
                     large_sd_alpha, length(results_list)))
log_progress(sprintf("CMI effects with SD > 1: %d / %d",
                     large_sd_cmi, length(results_list)))
log_progress(sprintf("Tmean effects with SD > 1: %d / %d",
                     large_sd_tmean, length(results_list)))

# ==============================================================================
# EXTRACT RESULTS
# ==============================================================================

log_progress("\n=== EXTRACTING RESULTS ===")

results_table <- do.call(rbind, lapply(results_list, function(res) {
  data.frame(
    Transition = res$transition,
    From = res$from,
    To = res$to,
    N_events = res$n_events,
    
    # Shape parameter
    Alpha = res$alpha,
    Alpha_SD = res$alpha_sd,
    Alpha_lower = res$alpha_lower,
    Alpha_upper = res$alpha_upper,
    
    # Intercept
    Beta0 = res$beta0,
    Beta0_SD = res$beta0_sd,
    Beta0_lower = res$beta0_lower,
    Beta0_upper = res$beta0_upper,
    
    # CMI effect
    CMI_effect = res$cmi_effect,
    CMI_SD = res$cmi_sd,
    CMI_lower = res$cmi_lower,
    CMI_upper = res$cmi_upper,
    CMI_sig = res$cmi_lower * res$cmi_upper > 0,
    
    # Temperature effect
    Tmean_effect = res$tmean_effect,
    Tmean_SD = res$tmean_sd,
    Tmean_lower = res$tmean_lower,
    Tmean_upper = res$tmean_upper,
    Tmean_sig = res$tmean_lower * res$tmean_upper > 0,
    
    # Soil effects
    Soil_Dry_effect = res$soil_dry_effect,
    Soil_Dry_SD = res$soil_dry_sd,
    Soil_Dry_lower = res$soil_dry_lower,
    Soil_Dry_upper = res$soil_dry_upper,
    Soil_Dry_sig = res$soil_dry_lower * res$soil_dry_upper > 0,
    
    Soil_Wet_effect = res$soil_wet_effect,
    Soil_Wet_SD = res$soil_wet_sd,
    Soil_Wet_lower = res$soil_wet_lower,
    Soil_Wet_upper = res$soil_wet_upper,
    Soil_Wet_sig = res$soil_wet_lower * res$soil_wet_upper > 0,
    
    # Pest effect
    Pest_effect = res$pest_effect,
    Pest_SD = res$pest_sd,
    Pest_lower = res$pest_lower,
    Pest_upper = res$pest_upper,
    Pest_sig = res$pest_lower * res$pest_upper > 0,
    
    # Model fit
    DIC = res$dic,
    WAIC = res$waic,
    
    # Convergence
    Convergence_OK = res$convergence_ok
  )
}))

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

# Save summary statistics
summary_stats <- list(
  n_transitions = nrow(results_table),
  n_plots = length(unique(model_data$id)),
  n_total_events = sum(results_table$N_events),
  fitting_time_mins = fit_time,
  n_significant_CMI = sum(results_table$CMI_sig),
  n_significant_Tmean = sum(results_table$Tmean_sig),
  n_convergence_issues = n_convergence_issues,
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

log_progress(sprintf("Significant Dry soil effects: %d / %d (%.1f%%)", 
                     sum(results_table$Soil_Dry_sig), 
                     nrow(results_table),
                     100 * sum(results_table$Soil_Dry_sig) / nrow(results_table)))
log_progress(sprintf("Significant Wet soil effects: %d / %d (%.1f%%)", 
                     sum(results_table$Soil_Wet_sig), 
                     nrow(results_table),
                     100 * sum(results_table$Soil_Wet_sig) / nrow(results_table)))
log_progress(sprintf("Significant Pest effects: %d / %d (%.1f%%)", 
                     sum(results_table$Pest_sig), 
                     nrow(results_table),
                     100 * sum(results_table$Pest_sig) / nrow(results_table)))

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