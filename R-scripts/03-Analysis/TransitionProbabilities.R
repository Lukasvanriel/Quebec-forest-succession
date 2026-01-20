#!/usr/bin/env Rscript
# ==============================================================================
# CALCULATE TRANSITION PROBABILITIES FROM FITTED MODELS
# Generates probability matrices for different covariate scenarios
# ==============================================================================

library(tidyverse)
library(data.table)
library(here)
library(gridExtra)
library(viridis)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths
RESULTS_DIR <- here("Data-Output", "INLA-realResults-IND")
OUTPUT_DIR <- here("Data-Output", "INLA-Probabilities")
PROGRESS_FILE <- file.path(OUTPUT_DIR, "probability_calculation_log.txt")

# Time points for probability calculation
TIME_POINTS <- c(5, 10, 15, 20)  # Years
DEFAULT_TIME <- 10  # Default time for tables/matrices

# Parallel processing
N_CORES <- 8  # Number of cores for parallel probability calculation

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Progress file for detailed monitoring
PROB_PROGRESS <- file.path(OUTPUT_DIR, "calculation_progress.txt")
cat("", file = PROB_PROGRESS)

progress_update <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg), 
      file = PROB_PROGRESS, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ==============================================================================
# PROGRESS LOGGING
# ==============================================================================

log_progress <- function(message, progress_file = PROGRESS_FILE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- sprintf("[%s] %s\n", timestamp, message)
  cat(log_message)
  cat(log_message, file = progress_file, append = TRUE)
}

cat("", file = PROGRESS_FILE)
log_progress("=== PROBABILITY CALCULATION STARTED ===")

# ==============================================================================
# LOAD FITTED MODELS
# ==============================================================================
log_progress("Loading fitted models...")

# Load individual model files
model_files <- list.files(RESULTS_DIR, pattern = "^fit_.*\\.rds$", full.names = TRUE)
log_progress(sprintf("Found %d model files", length(model_files)))

# Create a mapping of transition labels to file paths
# This avoids loading everything into RAM at once
model_map <- list()
for(file in model_files) {
  # This converts "fit_1_to_6.rds" into "1_to_6"
  model_name <- gsub("fit_|\\.rds", "", basename(file))
  model_map[[model_name]] <- file
}

log_progress("Model map created successfully")

# ==============================================================================
# DEFINE COVARIATE SCENARIOS
# ==============================================================================

log_progress("\n=== DEFINING COVARIATE SCENARIOS ===")

scenarios <- list(
  # Baseline scenarios
  baseline_optimal = list(
    name = "Baseline, Optimal soil",
    CMI = 0, Tmean = 0, soil = "Optimal", pest = 0
  ),
  baseline_dry = list(
    name = "Baseline, Dry soil",
    CMI = 0, Tmean = 0, soil = "Dry", pest = 0
  ),
  baseline_wet = list(
    name = "Baseline, Wet soil",
    CMI = 0, Tmean = 0, soil = "Wet", pest = 0
  ),
  
  # Post-pest scenarios
  post_pest_optimal = list(
    name = "Post-Pest, Optimal soil",
    CMI = 0, Tmean = 0, soil = "Optimal", pest = 1
  ),
  post_pest_dry = list(
    name = "Post-Pest, Dry soil",
    CMI = 0, Tmean = 0, soil = "Dry", pest = 1
  ),
  post_pest_wet = list(
    name = "Post-Pest, Wet soil",
    CMI = 0, Tmean = 0, soil = "Wet", pest = 1
  ),
  
  # Climate variation (optimal soil, no pest)
  warm_wet = list(
    name = "Warm & Wet (favorable for SBW)",
    CMI = 1, Tmean = 1, soil = "Optimal", pest = 0
  ),
  warm_dry = list(
    name = "Warm & Dry",
    CMI = -1, Tmean = 1, soil = "Optimal", pest = 0
  ),
  cool_wet = list(
    name = "Cool & Wet",
    CMI = 1, Tmean = -1, soil = "Optimal", pest = 0
  ),
  cool_dry = list(
    name = "Cool & Dry",
    CMI = -1, Tmean = -1, soil = "Optimal", pest = 0
  ),
  
  # Post-pest + favorable climate
  post_pest_favorable = list(
    name = "Post-Pest, Favorable climate",
    CMI = 1, Tmean = 1, soil = "Optimal", pest = 1
  )
)

log_progress(sprintf("Defined %d scenarios:", length(scenarios)))
for(i in seq_along(scenarios)) {
  log_progress(sprintf("  %d. %s", i, scenarios[[i]]$name))
}

# ==============================================================================
# FUNCTION TO CALCULATE TRANSITION PROBABILITY
# ==============================================================================

# Weibull hazard with variant 1: h(t) = α * t^(α-1) * λ^α
# where λ = exp(linear predictor)
# Cumulative hazard: H(t) = (λ * t)^α

calculate_transition_prob <- function(model_result, time, scenario) {
  
  # Extract parameters
  alpha <- model_result$alpha
  beta0 <- model_result$beta0
  beta_cmi <- model_result$cmi_effect
  beta_tmean <- model_result$tmean_effect
  beta_soil_dry <- model_result$soil_dry_effect
  beta_soil_wet <- model_result$soil_wet_effect
  beta_pest <- model_result$pest_effect
  
  # Build linear predictor
  linear_pred <- beta0 + 
    beta_cmi * scenario$CMI + 
    beta_tmean * scenario$Tmean +
    beta_pest * scenario$pest
  
  # Add soil effect
  if(scenario$soil == "Dry") {
    linear_pred <- linear_pred + beta_soil_dry
  } else if(scenario$soil == "Wet") {
    linear_pred <- linear_pred + beta_soil_wet
  }
  # Optimal is reference (no addition)
  
  # Calculate lambda
  lambda <- exp(linear_pred)
  
  # Cumulative hazard at time t
  cum_hazard <- (lambda * time)^alpha
  
  return(list(
    lambda = lambda,
    alpha = alpha,
    cum_hazard = cum_hazard,
    linear_pred = linear_pred
  ))
}

# ==============================================================================
# CALCULATE PROBABILITIES FOR ALL TRANSITIONS AND SCENARIOS
# ==============================================================================

log_progress("\n=== CALCULATING PROBABILITIES (OPTIMIZED) ===")

# 1. PRE-LOAD PARAMETERS (The Secret to Speed)
log_progress("Pre-loading model parameters to avoid disk bottleneck...")
coef_list <- lapply(names(model_map), function(m_name) {
  res <- readRDS(model_map[[m_name]])
  # Extract only the numbers we need
  data.table(
    transition = m_name,
    alpha      = res$alpha,
    beta0      = res$beta0,
    cmi_b      = res$cmi_effect,
    tmean_b    = res$tmean_effect,
    soil_dry_b = res$soil_dry_effect,
    soil_wet_b = res$soil_wet_effect,
    pest_b     = res$pest_effect
  )
})
coef_dt <- rbindlist(coef_list)
setkey(coef_dt, transition) # Enables instant searching

# 2. RUN SCENARIOS
# We use the pre-loaded coef_dt instead of reading files inside the loop
scenario_results <- mclapply(names(scenarios), function(scen_name) {
  
  progress_update(sprintf("Processing scenario: %s", scen_name))
  scenario_vals <- scenarios[[scen_name]]
  scenario_data <- prob_data[prob_data$scenario_name == scen_name, ]
  
  for(i in 1:nrow(scenario_data)) {
    m_name <- paste0(scenario_data$from[i], "_to_", scenario_data$to[i])
    
    # Grab params from RAM (instant)
    m_coefs <- coef_dt[transition == m_name]
    
    if(nrow(m_coefs) > 0) {
      # Calculate linear predictor
      lp <- m_coefs$beta0 + 
        (scenario_vals$CMI * m_coefs$cmi_b) + 
        (scenario_vals$Tmean * m_coefs$tmean_b) +
        (scenario_vals$pest * m_coefs$pest_b)
      
      if(scenario_vals$soil == "Dry") lp <- lp + m_coefs$soil_dry_b
      if(scenario_vals$soil == "Wet") lp <- lp + m_coefs$soil_wet_b
      
      lambda <- exp(lp)
      set(scenario_data, i, "lambda", lambda)
      set(scenario_data, i, "alpha", m_coefs$alpha)
      set(scenario_data, i, "cum_hazard", (lambda * scenario_data$time[i])^m_coefs$alpha)
    }
  }
  
  progress_update(sprintf("Completed scenario: %s", scen_name))
  return(scenario_data)
}, mc.cores = N_CORES)

prob_data <- rbindlist(scenario_results)
log_progress("All probabilities calculated using in-memory parameters.")

# ==============================================================================
# CONVERT TO PROBABILITIES (ACCOUNTING FOR COMPETING RISKS)
# ==============================================================================
log_progress("Converting to probabilities with competing risks...")

# 1. Calculate total hazard per state/time
prob_data[, total_cum_hazard := sum(cum_hazard, na.rm = TRUE), 
          by = .(scenario_name, time, from)]

# 2. Calculate survival probability (staying in same state)
prob_data[, survival_in_state := exp(-total_cum_hazard)]

# 3. Apportion total exit probability (1 - survival) to specific transitions
prob_data[, probability := (cum_hazard / total_cum_hazard) * (1 - survival_in_state)]

# Replace NaN with 0
prob_data[is.nan(probability), probability := 0]

# Add diagonal (probability of staying in same state)
diagonal_probs <- prob_data %>%
  group_by(scenario_name, time, from) %>%
  summarise(
    survival_in_state = first(survival_in_state),
    .groups = "drop"
  ) %>%
  mutate(
    to = from,
    probability = survival_in_state,
    lambda = NA,
    alpha = NA,
    total_cum_hazard = NA
  )

prob_data <- bind_rows(
  prob_data %>% select(scenario_name, time, from, to, probability, lambda, alpha),
  diagonal_probs %>% select(scenario_name, time, from, to, probability, lambda, alpha)
)

log_progress("Probabilities calculated successfully")

# SAVE CHECKPOINT BEFORE VISUALIZATION
log_progress("Saving probability data checkpoint...")
saveRDS(prob_data, file.path(OUTPUT_DIR, "prob_data_checkpoint.rds"))
saveRDS(list(
  prob_data = prob_data,
  scenarios = scenarios,
  fitted_models = names(model_map), # Changed from fitted_models to names(model_map)
  timestamp = Sys.time()
), file.path(OUTPUT_DIR, "full_checkpoint.rds"))
log_progress("Checkpoint saved - visualization can be re-run from here if needed")

# ==============================================================================
# CREATE PROBABILITY MATRICES
# ==============================================================================

log_progress("\n=== CREATING PROBABILITY MATRICES ===")

# Function to create transition matrix for a scenario/time
create_prob_matrix <- function(scenario_name, time_point) {
  
  probs <- prob_data %>%
    filter(scenario_name == !!scenario_name, time == !!time_point) %>%
    select(from, to, probability) %>%
    pivot_wider(names_from = to, values_from = probability, values_fill = 0)
  
  # Convert to matrix
  from_states <- probs$from
  probs <- probs %>% select(-from) %>% as.matrix()
  rownames(probs) <- from_states
  
  return(probs)
}

# Create matrices for default time point
matrices <- list()
for(scenario_name in names(scenarios)) {
  matrices[[scenario_name]] <- create_prob_matrix(scenario_name, DEFAULT_TIME)
}

# Save matrices
saveRDS(matrices, file.path(OUTPUT_DIR, "probability_matrices.rds"))
log_progress(sprintf("Saved %d probability matrices", length(matrices)))

# Save as CSV files
for(scenario_name in names(scenarios)) {
  filename <- paste0("matrix_", gsub("[^A-Za-z0-9]", "_", scenario_name), ".csv")
  write.csv(matrices[[scenario_name]], 
            file.path(OUTPUT_DIR, filename))
}

log_progress("Matrices saved as CSV files")

# ==============================================================================
# CREATE COMPARISON TABLES
# ==============================================================================

log_progress("\n=== CREATING COMPARISON TABLES ===")

# Key transitions for SBW analysis
key_transitions <- list(
  "Balsam Fir mortality" = c("6_to_7", "6_to_1", "6_to_2", "6_to_5"),
  "Spruce mortality" = c("7_to_6", "7_to_1", "7_to_2", "7_to_5"),
  "Pioneer to mature" = c("1_to_6", "2_to_6", "1_to_7", "2_to_7"),
  "Regeneration to BF/Spruce" = c("5_to_6", "5_to_7")
)

# Create comparison table for key transitions
comparison_table <- prob_data %>%
  filter(time == DEFAULT_TIME) %>%
  mutate(transition = paste0(from, "_to_", to)) %>%
  filter(transition %in% unlist(key_transitions)) %>%
  select(scenario_name, transition, probability) %>%
  pivot_wider(names_from = scenario_name, values_from = probability)

write.csv(comparison_table, 
          file.path(OUTPUT_DIR, "key_transitions_comparison.csv"),
          row.names = FALSE)

log_progress("Key transitions comparison table saved")

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

log_progress("\n=== CREATING VISUALIZATIONS ===")

# 1. Probabilities over time for key transitions
log_progress("Creating time-series plots...")

for(trans_group in names(key_transitions)) {
  
  transitions <- key_transitions[[trans_group]]
  
  plot_data <- prob_data %>%
    mutate(transition = paste0(from, "→", to)) %>%
    filter(transition %in% transitions) %>%
    left_join(
      tibble(scenario_name = names(scenarios),
             scenario_label = sapply(scenarios, function(x) x$name)),
      by = "scenario_name"
    )
  
  p <- ggplot(plot_data, aes(x = time, y = probability, 
                             color = transition, linetype = scenario_label)) +
    geom_line(size = 1) +
    facet_wrap(~ transition, scales = "free_y") +
    scale_color_viridis_d(option = "plasma", end = 0.9) +
    labs(
      title = paste("Transition Probabilities Over Time:", trans_group),
      subtitle = sprintf("Comparing scenarios at %d time points", length(TIME_POINTS)),
      x = "Time (years)",
      y = "Probability",
      color = "Transition",
      linetype = "Scenario"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = element_text(face = "bold")
    )
  
  filename <- paste0("timeseries_", gsub("[^A-Za-z0-9]", "_", trans_group), ".png")
  ggsave(file.path(OUTPUT_DIR, filename), p, width = 14, height = 10, dpi = 300)
  
  log_progress(sprintf("  Saved: %s", filename))
}

# 2. Heatmaps of probability matrices
log_progress("Creating heatmap visualizations...")

create_heatmap <- function(scenario_name, time_point = DEFAULT_TIME) {
  
  plot_data <- prob_data %>%
    filter(scenario_name == !!scenario_name, time == !!time_point) %>%
    mutate(
      from = factor(from, levels = sort(unique(from))),
      to = factor(to, levels = sort(unique(to)))
    )
  
  scenario_label <- scenarios[[scenario_name]]$name
  
  p <- ggplot(plot_data, aes(x = to, y = from, fill = probability)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", probability)), 
              size = 3, color = "black") +
    scale_fill_viridis_c(option = "magma", begin = 0.1, limits = c(0, 1)) +
    labs(
      title = scenario_label,
      subtitle = sprintf("Transition probabilities over %d years", time_point),
      x = "To State",
      y = "From State",
      fill = "Probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}

# Create heatmaps for key scenarios
key_scenarios <- c("baseline_optimal", "post_pest_optimal", 
                   "post_pest_favorable", "warm_wet")

for(scenario_name in key_scenarios) {
  p <- create_heatmap(scenario_name)
  filename <- paste0("heatmap_", scenario_name, ".png")
  ggsave(file.path(OUTPUT_DIR, filename), p, width = 10, height = 8, dpi = 300)
  log_progress(sprintf("  Saved: %s", filename))
}

# 3. Comparison plots: Pest effect
log_progress("Creating pest effect comparison plots...")

pest_comparison <- prob_data %>%
  filter(time == DEFAULT_TIME) %>%
  filter(scenario_name %in% c("baseline_optimal", "post_pest_optimal")) %>%
  mutate(transition = paste0(from, "→", to)) %>%
  select(scenario_name, transition, probability) %>%
  pivot_wider(names_from = scenario_name, values_from = probability) %>%
  mutate(
    pest_effect = post_pest_optimal - baseline_optimal,
    relative_change = (post_pest_optimal - baseline_optimal) / baseline_optimal
  ) %>%
  filter(!is.na(pest_effect)) %>%
  arrange(desc(abs(pest_effect))) %>%
  head(20)

p_pest <- ggplot(pest_comparison, 
                 aes(x = reorder(transition, pest_effect), y = pest_effect)) +
  geom_col(aes(fill = pest_effect > 0)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                    labels = c("Decrease", "Increase")) +
  labs(
    title = "Top 20 Transitions Most Affected by Pest Outbreak",
    subtitle = sprintf("Change in %d-year transition probability (Post-Pest - Baseline)", 
                       DEFAULT_TIME),
    x = "Transition",
    y = "Change in Probability",
    fill = "Effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "pest_effect_comparison.png"), p_pest,
       width = 10, height = 8, dpi = 300)

log_progress("  Saved: pest_effect_comparison.png")

# 4. Soil effect comparison
log_progress("Creating soil effect comparison plots...")

soil_comparison <- prob_data %>%
  filter(time == DEFAULT_TIME) %>%
  filter(scenario_name %in% c("baseline_optimal", "baseline_dry", "baseline_wet")) %>%
  mutate(transition = paste0(from, "→", to)) %>%
  select(scenario_name, transition, probability) %>%
  pivot_wider(names_from = scenario_name, values_from = probability) %>%
  mutate(
    dry_effect = baseline_dry - baseline_optimal,
    wet_effect = baseline_wet - baseline_optimal
  ) %>%
  filter(!is.na(dry_effect) | !is.na(wet_effect)) %>%
  pivot_longer(cols = c(dry_effect, wet_effect), 
               names_to = "soil_type", values_to = "effect") %>%
  group_by(soil_type) %>%
  arrange(desc(abs(effect))) %>%
  slice_head(n = 15) %>%
  ungroup()

p_soil <- ggplot(soil_comparison, 
                 aes(x = reorder(transition, effect), y = effect, fill = soil_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(
    values = c("dry_effect" = "#E69F00", "wet_effect" = "#56B4E9"),
    labels = c("Dry soil", "Wet soil")
  ) +
  labs(
    title = "Top Transitions Affected by Suboptimal Soil",
    subtitle = sprintf("Change in %d-year probability vs. Optimal soil", DEFAULT_TIME),
    x = "Transition",
    y = "Change in Probability",
    fill = "Soil Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "soil_effect_comparison.png"), p_soil,
       width = 10, height = 10, dpi = 300)

log_progress("  Saved: soil_effect_comparison.png")

# ==============================================================================
# SAVE FULL PROBABILITY DATA
# ==============================================================================

log_progress("\n=== SAVING COMPLETE RESULTS ===")

# Save full probability dataset
fwrite(prob_data, file.path(OUTPUT_DIR, "all_probabilities.csv"))
log_progress("Saved complete probability data")

# Create summary document
summary_doc <- sprintf("
TRANSITION PROBABILITY CALCULATION SUMMARY
==========================================

Date: %s

INPUT:
- Fitted models directory: %s
- Number of transitions: %d
- Number of scenarios: %d

TIME POINTS: %s years

SCENARIOS:
%s

OUTPUTS:
- Probability matrices: probability_matrices.rds
- Individual CSV matrices: matrix_*.csv
- Complete probability data: all_probabilities.csv
- Key transitions comparison: key_transitions_comparison.csv
- Visualizations: *.png files

NOTES:
- Probabilities calculated using Weibull hazard model (variant 1)
- Competing risks accounted for in probability calculation
- Default time for matrices: %d years
",
                       Sys.time(),
                       RESULTS_DIR,
                       nrow(results_table),
                       length(scenarios),
                       paste(TIME_POINTS, collapse = ", "),
                       paste(sprintf("  - %s: %s", names(scenarios), 
                                     sapply(scenarios, function(x) x$name)), collapse = "\n"),
                       DEFAULT_TIME
)

writeLines(summary_doc, file.path(OUTPUT_DIR, "README.txt"))
log_progress("Created summary document")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), 
                                  as.POSIXct(readLines(PROGRESS_FILE)[1] %>%
                                               str_extract("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
                                             format = "%Y-%m-%d %H:%M:%S"),
                                  units = "mins"))

log_progress("\n=== CALCULATION COMPLETE ===")
log_progress(sprintf("Total scenarios: %d", length(scenarios)))
log_progress(sprintf("Total probabilities calculated: %d", nrow(prob_data)))
log_progress(sprintf("Time points: %s", paste(TIME_POINTS, collapse = ", ")))
log_progress(sprintf("Total time: %.1f minutes", total_time))
log_progress(sprintf("Output directory: %s", OUTPUT_DIR))
log_progress("============================")

cat("\n\n")
cat("##############################################\n")
cat("# PROBABILITY CALCULATION COMPLETE\n")
cat(sprintf("# Scenarios: %d\n", length(scenarios)))
cat(sprintf("# Probabilities: %d\n", nrow(prob_data)))
cat(sprintf("# Output: %s\n", OUTPUT_DIR))
cat("##############################################\n")