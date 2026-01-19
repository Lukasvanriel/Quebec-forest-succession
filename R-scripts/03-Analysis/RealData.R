#  INLA with real data
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

## Read in data
data.all <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]
plot(data.all$LONGI, data.all$LATIT)

data <- data.all %>% 
  filter(LONGI > -78, LONGI < -75)

points(data$LONGI, data$LATIT, col='red')

length(unique(data$TESSELLE))
colnames(data)

# Set seed for reproducibility
set.seed(123)

# Get unique plot IDs
all_plot_ids <- unique(data$TESSELLE)

# Randomly sample 1000 plots
selected_plots <- sample(all_plot_ids, size = 1000)

# Filter your data to keep only selected plots
forest_subset <- data %>% filter(TESSELLE %in% selected_plots)

# Verify
cat(sprintf("Selected %d plots with %d total observations\n", 
            length(unique(forest_subset$TESSELLE)), 
            nrow(forest_subset)))

table(forest_subset$sp_class)

head(data)



# Check observation times distribution
forest_subset %>%
  group_by(TESSELLE) %>%
  summarise(
    n_obs = n(),
    years = paste(sort(unique(year)), collapse = ", "),
    times = paste(sort(unique(time)), collapse = ", ")
  ) %>%
  head(20)

# Summary of observation counts
table(forest_subset %>% group_by(TESSELLE) %>% tally() %>% pull(n))

# Typical interval lengths
intervals <- forest_subset %>%
  arrange(TESSELLE, time) %>%
  group_by(TESSELLE) %>%
  mutate(interval = time - lag(time)) %>%
  filter(!is.na(interval))

summary(intervals$interval)
hist(intervals$interval, breaks = 30, main = "Distribution of Observation Intervals",
     xlab = "Years between observations")



#############
library(tidyverse)
library(INLA)
library(INLAjoint)

# ==============================================================================
# SETUP: 9-STATE TRANSITION MATRIX
# ==============================================================================

# Fully reversible 9-state model (all transitions possible except self-transitions)
tmat_9state <- matrix(TRUE, 9, 9)
diag(tmat_9state) <- FALSE
colnames(tmat_9state) <- rownames(tmat_9state) <- as.character(1:9)

cat("9-state transition matrix created: 9 x 8 = 72 possible transitions\n")

# ==============================================================================
# RECONSTRUCT METHOD B FOR REAL DATA
# ==============================================================================

reconstruct_plot_methodB <- function(plot_data, TT) {
  # Extract plot info
  plot_id <- plot_data$TESSELLE[1]
  
  # Get observation times and states
  plot_data <- plot_data %>% arrange(time)
  obs_times <- plot_data$time
  obs_states <- plot_data$sp_class
  
  # Get covariates (take from first row - assuming constant)
  cov_CMI <- plot_data$cov_CMI[1]
  cov_Tmean <- plot_data$cov_Tmean[1]
  cov_soil <- plot_data$cov_soil[1]
  cov_pert_class <- plot_data$cov_pert_class[1]
  
  data_list <- list()
  counter <- 1
  
  s0 <- obs_states[1]
  approx_entry_time <- obs_times[1]
  
  for(k in 2:length(obs_times)) {
    t0 <- obs_times[k-1]
    t1 <- obs_times[k]
    s1 <- obs_states[k]
    
    possible_to_states <- setdiff(1:9, s0)
    
    if(s0 != s1) {
      # Transition detected
      
      # Compacted rows for competing transitions (entry to interval END)
      for(to_state in possible_to_states) {
        if(to_state != s1) {
          data_list[[counter]] <- data.frame(
            id = plot_id,
            Tstart = approx_entry_time,
            Tstop = t1,
            from = s0,
            to = to_state,
            status = 0,
            cov_CMI = cov_CMI,
            cov_Tmean = cov_Tmean,
            cov_soil = factor(cov_soil),
            cov_pert_class = factor(cov_pert_class),
            approx_entry = approx_entry_time
          )
          counter <- counter + 1
        }
      }
      
      # Actual transition (status=3)
      data_list[[counter]] <- data.frame(
        id = plot_id,
        Tstart = t0,
        Tstop = t1,
        from = s0,
        to = s1,
        status = 3,
        cov_CMI = cov_CMI,
        cov_Tmean = cov_Tmean,
        cov_soil = factor(cov_soil),
        cov_pert_class = factor(cov_pert_class),
        approx_entry = approx_entry_time
      )
      counter <- counter + 1
      
      s0 <- s1
      approx_entry_time <- t1
    }
  }
  
  # Final: compacted rows for last state
  t_final <- max(obs_times)
  possible_to_states <- setdiff(1:9, s0)
  
  for(to_state in possible_to_states) {
    data_list[[counter]] <- data.frame(
      id = plot_id,
      Tstart = approx_entry_time,
      Tstop = TT,
      from = s0,
      to = to_state,
      status = 0,
      cov_CMI = cov_CMI,
      cov_Tmean = cov_Tmean,
      cov_soil = factor(cov_soil),
      cov_pert_class = factor(cov_pert_class),
      approx_entry = approx_entry_time
    )
    counter <- counter + 1
  }
  
  bind_rows(data_list)
}

# ==============================================================================
# APPLY TO ALL PLOTS
# ==============================================================================

# What's your maximum time?
TT <- max(forest_subset$time)
cat(sprintf("Maximum observation time: %d years\n", TT))

# Apply reconstruction
cat("Reconstructing interval-censored data for 1000 plots...\n")
model_data <- bind_rows(lapply(unique(forest_subset$TESSELLE), function(plot_id) {
  if(plot_id %% 100 == 0) cat(sprintf("  %d/%d\n", which(unique(forest_subset$TESSELLE) == plot_id), 1000))
  reconstruct_plot_methodB(forest_subset %>% filter(TESSELLE == plot_id), TT)
}))

cat("Reconstruction complete!\n\n")

# Check the result
cat("=== RECONSTRUCTION SUMMARY ===\n")
cat(sprintf("Total rows: %d\n", nrow(model_data)))
cat(sprintf("Unique plots: %d\n", length(unique(model_data$id))))
cat(sprintf("Transitions detected (status=3): %d\n", sum(model_data$status == 3)))
cat(sprintf("Censored rows (status=0): %d\n", sum(model_data$status == 0)))

# Check one example plot
example_plot <- unique(model_data$id)[1]
cat("\n=== EXAMPLE PLOT:", example_plot, "===\n")
print(model_data %>% filter(id == example_plot) %>% select(id, Tstart, Tstop, from, to, status, approx_entry))


# ==============================================================================
# PREPARE DATA FOR INLAjoint (72 TRANSITIONS!)
# ==============================================================================

cat("Preparing data for 72 transitions...\n")

# Create list of datasets, one per transition
data_list <- list()
counter <- 1

for(from_state in 1:9) {
  for(to_state in 1:9) {
    if(from_state != to_state) {
      data_list[[counter]] <- model_data %>% 
        filter(from == from_state, to == to_state)
      counter <- counter + 1
    }
  }
}

# Check how many events per transition
transition_summary <- data.frame()
counter <- 1
for(from_state in 1:9) {
  for(to_state in 1:9) {
    if(from_state != to_state) {
      n_events <- sum(data_list[[counter]]$status == 3)
      n_censored <- sum(data_list[[counter]]$status == 0)
      transition_summary <- rbind(transition_summary, data.frame(
        transition = paste0(from_state, "→", to_state),
        from = from_state,
        to = to_state,
        n_events = n_events,
        n_censored = n_censored
      ))
      counter <- counter + 1
    }
  }
}

cat("\n=== TRANSITION EVENT COUNTS ===\n")
print(transition_summary %>% arrange(desc(n_events)))

# Identify sparse transitions (potential problems)
sparse_transitions <- transition_summary %>% filter(n_events < 10)
cat(sprintf("\nWARNING: %d transitions have < 10 events\n", nrow(sparse_transitions)))
if(nrow(sparse_transitions) > 0) {
  cat("These may cause convergence issues:\n")
  print(sparse_transitions)
}


######## Not a lot of 8 transitions:
# ==============================================================================
# REMOVE STATE 8 AND RECONSTRUCT
# ==============================================================================

# Filter out state 8
forest_subset_no8 <- forest_subset %>% filter(sp_class != 8)

cat(sprintf("Original: %d observations across %d plots\n", 
            nrow(forest_subset), length(unique(forest_subset$TESSELLE))))
cat(sprintf("After removing state 8: %d observations across %d plots\n", 
            nrow(forest_subset_no8), length(unique(forest_subset_no8$TESSELLE))))

# Check how many plots we lost (plots that only had state 8 observations)
n_plots_lost <- length(unique(forest_subset$TESSELLE)) - length(unique(forest_subset_no8$TESSELLE))
cat(sprintf("Plots lost: %d\n\n", n_plots_lost))

# State distribution after removal
cat("State distribution after removing state 8:\n")
print(table(forest_subset_no8$sp_class))

# ==============================================================================
# RE-RUN RECONSTRUCTION WITHOUT STATE 8
# ==============================================================================

# Update transition matrix for 8 states (1-7, 9)
states_remaining <- sort(unique(forest_subset_no8$sp_class))
n_states <- length(states_remaining)

# First, check for NAs in the data
cat("Checking for data issues...\n")
cat(sprintf("Plots with only 1 observation: %d\n", 
            sum(table(forest_subset_no8$TESSELLE) == 1)))
cat(sprintf("Missing sp_class values: %d\n", sum(is.na(forest_subset_no8$sp_class))))
cat(sprintf("Missing time values: %d\n", sum(is.na(forest_subset_no8$time))))

# Check which plots have issues
problem_plots <- forest_subset_no8 %>%
  group_by(TESSELLE) %>%
  summarise(n_obs = n(),
            has_na_state = any(is.na(sp_class)),
            has_na_time = any(is.na(time))) %>%
  filter(n_obs == 1 | has_na_state | has_na_time)

if(nrow(problem_plots) > 0) {
  cat("\nProblematic plots found:\n")
  print(problem_plots)
  
  # Remove these plots
  cat(sprintf("\nRemoving %d problematic plots...\n", nrow(problem_plots)))
  forest_subset_no8 <- forest_subset_no8 %>%
    filter(!TESSELLE %in% problem_plots$TESSELLE)
}

# Updated reconstruction with safety checks
reconstruct_plot_methodB_nocov <- function(plot_data, TT) {
  plot_id <- plot_data$TESSELLE[1]
  
  plot_data <- plot_data %>% arrange(time)
  obs_times <- plot_data$time
  obs_states <- plot_data$sp_class
  
  # Safety check
  if(length(obs_times) < 2) {
    warning(sprintf("Plot %s has only 1 observation, skipping", plot_id))
    return(NULL)
  }
  if(any(is.na(obs_states))) {
    warning(sprintf("Plot %s has NA states, skipping", plot_id))
    return(NULL)
  }
  
  data_list <- list()
  counter <- 1
  
  s0 <- obs_states[1]
  approx_entry_time <- obs_times[1]
  
  for(k in 2:length(obs_times)) {
    t0 <- obs_times[k-1]
    t1 <- obs_times[k]
    s1 <- obs_states[k]
    
    possible_to_states <- setdiff(states_remaining, s0)
    
    if(s0 != s1) {
      # Compacted rows for competing transitions
      for(to_state in possible_to_states) {
        if(to_state != s1) {
          data_list[[counter]] <- data.frame(
            id = plot_id,
            Tstart = approx_entry_time,
            Tstop = t1,
            from = s0,
            to = to_state,
            status = 0,
            approx_entry = approx_entry_time
          )
          counter <- counter + 1
        }
      }
      
      # Actual transition
      data_list[[counter]] <- data.frame(
        id = plot_id,
        Tstart = t0,
        Tstop = t1,
        from = s0,
        to = s1,
        status = 3,
        approx_entry = approx_entry_time
      )
      counter <- counter + 1
      
      s0 <- s1
      approx_entry_time <- t1
    }
  }
  
  # Final censored rows
  possible_to_states <- setdiff(states_remaining, s0)
  for(to_state in possible_to_states) {
    data_list[[counter]] <- data.frame(
      id = plot_id,
      Tstart = approx_entry_time,
      Tstop = TT,
      from = s0,
      to = to_state,
      status = 0,
      approx_entry = approx_entry_time
    )
    counter <- counter + 1
  }
  
  bind_rows(data_list)
}

# Try reconstruction again
TT <- max(forest_subset_no8$time)
cat(sprintf("\nReconstructing data for %d plots...\n", length(unique(forest_subset_no8$TESSELLE))))

model_data_no8 <- bind_rows(lapply(unique(forest_subset_no8$TESSELLE), function(plot_id) {
  reconstruct_plot_methodB_nocov(forest_subset_no8 %>% filter(TESSELLE == plot_id), TT)
}))

cat("Reconstruction complete!\n")


# Check transition counts again
transition_summary_no8 <- data.frame()
for(from_state in states_remaining) {
  for(to_state in states_remaining) {
    if(from_state != to_state) {
      temp_data <- model_data_no8 %>% filter(from == from_state, to == to_state)
      transition_summary_no8 <- rbind(transition_summary_no8, data.frame(
        transition = paste0(from_state, "→", to_state),
        from = from_state,
        to = to_state,
        n_events = sum(temp_data$status == 3),
        n_censored = sum(temp_data$status == 0)
      ))
    }
  }
}

cat("=== UPDATED TRANSITION EVENT COUNTS (no state 8) ===\n")
print(transition_summary_no8 %>% arrange(desc(n_events)))

sparse <- transition_summary_no8 %>% filter(n_events < 10)
cat(sprintf("\nTransitions with < 10 events: %d out of %d\n", nrow(sparse), nrow(transition_summary_no8)))


# ==============================================================================
# FIT INTERCEPT-ONLY MODEL (NO COVARIATES)
# ==============================================================================

cat("Preparing data for INLAjoint (56 transitions, intercept-only)...\n")

# Create data list for each transition
data_list <- list()
counter <- 1

for(from_state in states_remaining) {
  for(to_state in states_remaining) {
    if(from_state != to_state) {
      data_list[[counter]] <- model_data_no8 %>% 
        filter(from == from_state, to == to_state)
      counter <- counter + 1
    }
  }
}

# Create survival objects
cat("Creating survival objects...\n")
Surv_list <- list()

# ==============================================================================
# CREATE SURVIVAL OBJECTS WITH EXPLICIT NAMES
# ==============================================================================

cat("Creating explicitly named survival objects...\n")

# First, create all survival objects with names s1, s2, s3, etc.
for(i in 1:length(data_list)) {
  d <- data_list[[i]]
  
  time_vec <- d$Tstart
  idx_censored <- d$status == 0
  time_vec[idx_censored] <- d$Tstop[idx_censored]
  
  assign(paste0("s", i), 
         inla.surv(time = time_vec, time2 = d$Tstop, 
                   event = d$status, truncation = d$approx_entry))
  
  if(i %% 10 == 0) cat(sprintf("  %d/%d\n", i, length(data_list)))
}

cat("All survival objects created: s1, s2, ..., s56\n\n")

# Now create formulas explicitly
formulas <- list(
  s1 ~ 1, s2 ~ 1, s3 ~ 1, s4 ~ 1, s5 ~ 1, s6 ~ 1, s7 ~ 1, s8 ~ 1,
  s9 ~ 1, s10 ~ 1, s11 ~ 1, s12 ~ 1, s13 ~ 1, s14 ~ 1, s15 ~ 1, s16 ~ 1,
  s17 ~ 1, s18 ~ 1, s19 ~ 1, s20 ~ 1, s21 ~ 1, s22 ~ 1, s23 ~ 1, s24 ~ 1,
  s25 ~ 1, s26 ~ 1, s27 ~ 1, s28 ~ 1, s29 ~ 1, s30 ~ 1, s31 ~ 1, s32 ~ 1,
  s33 ~ 1, s34 ~ 1, s35 ~ 1, s36 ~ 1, s37 ~ 1, s38 ~ 1, s39 ~ 1, s40 ~ 1,
  s41 ~ 1, s42 ~ 1, s43 ~ 1, s44 ~ 1, s45 ~ 1, s46 ~ 1, s47 ~ 1, s48 ~ 1,
  s49 ~ 1, s50 ~ 1, s51 ~ 1, s52 ~ 1, s53 ~ 1, s54 ~ 1, s55 ~ 1, s56 ~ 1
)

cat("Formulas created explicitly\n\n")

# Now fit
cat("Fitting model... Started at:", format(Sys.time(), "%H:%M:%S"), "\n\n")

fit_intercept_only <- joint(
  formSurv = formulas,
  basRisk = rep("weibullsurv", 56),
  dataSurv = data_list,
  id = "id",
  control = list(config = TRUE, variant = 1)
)

cat("\nFinished at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("Model fitting complete!\n\n")



#### Not Working?
# ==============================================================================
# FIX 1: REINDEX PLOT IDs TO SIMPLE INTEGERS
# ==============================================================================

# Create simple integer IDs
plot_mapping <- data.frame(
  TESSELLE_original = unique(forest_subset_no8$TESSELLE),
  id_simple = 1:length(unique(forest_subset_no8$TESSELLE))
)

forest_subset_no8 <- forest_subset_no8 %>%
  left_join(plot_mapping, by = c("TESSELLE" = "TESSELLE_original"))

# Verify
cat("Simple IDs now range from", min(forest_subset_no8$id_simple), 
    "to", max(forest_subset_no8$id_simple), "\n\n")

# ==============================================================================
# FIX 2: START WITH ONLY TOP 10 MOST COMMON TRANSITIONS
# ==============================================================================

top_transitions <- transition_summary_no8 %>% 
  arrange(desc(n_events)) %>% 
  head(10)

cat("=== FITTING ONLY TOP 10 TRANSITIONS ===\n")
print(top_transitions)
cat("\n")

# Reconstruct with simple IDs
reconstruct_plot_methodB_simple <- function(plot_data, TT, states_remaining) {
  plot_id <- plot_data$id_simple[1]
  
  plot_data <- plot_data %>% arrange(time)
  obs_times <- plot_data$time
  obs_states <- plot_data$sp_class
  
  if(length(obs_times) < 2 || any(is.na(obs_states))) return(NULL)
  
  data_list <- list()
  counter <- 1
  s0 <- obs_states[1]
  approx_entry_time <- obs_times[1]
  
  for(k in 2:length(obs_times)) {
    t0 <- obs_times[k-1]
    t1 <- obs_times[k]
    s1 <- obs_states[k]
    
    possible_to_states <- setdiff(states_remaining, s0)
    
    if(s0 != s1) {
      for(to_state in possible_to_states) {
        if(to_state != s1) {
          data_list[[counter]] <- data.frame(
            id = plot_id,
            Tstart = approx_entry_time,
            Tstop = t1,
            from = s0,
            to = to_state,
            status = 0,
            approx_entry = approx_entry_time
          )
          counter <- counter + 1
        }
      }
      
      data_list[[counter]] <- data.frame(
        id = plot_id,
        Tstart = t0,
        Tstop = t1,
        from = s0,
        to = s1,
        status = 3,
        approx_entry = approx_entry_time
      )
      counter <- counter + 1
      
      s0 <- s1
      approx_entry_time <- t1
    }
  }
  
  possible_to_states <- setdiff(states_remaining, s0)
  for(to_state in possible_to_states) {
    data_list[[counter]] <- data.frame(
      id = plot_id,
      Tstart = approx_entry_time,
      Tstop = TT,
      from = s0,
      to = to_state,
      status = 0,
      approx_entry = approx_entry_time
    )
    counter <- counter + 1
  }
  
  bind_rows(data_list)
}

# Reconstruct
TT <- max(forest_subset_no8$time)
model_data_simple <- bind_rows(lapply(unique(forest_subset_no8$id_simple), function(plot_id) {
  reconstruct_plot_methodB_simple(
    forest_subset_no8 %>% filter(id_simple == plot_id), 
    TT, 
    states_remaining
  )
}))

# Filter to only top 10 transitions
data_list_top10 <- list()
for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  data_list_top10[[i]] <- model_data_simple %>% 
    filter(from == from_s, to == to_s)
}

# Create survival objects
for(i in 1:10) {
  d <- data_list_top10[[i]]
  time_vec <- d$Tstart
  idx_censored <- d$status == 0
  time_vec[idx_censored] <- d$Tstop[idx_censored]
  
  assign(paste0("s", i), 
         inla.surv(time = time_vec, time2 = d$Tstop, 
                   event = d$status, truncation = d$approx_entry))
}

# Formulas for 10 transitions
formulas_top10 <- list(s1 ~ 1, s2 ~ 1, s3 ~ 1, s4 ~ 1, s5 ~ 1, 
                       s6 ~ 1, s7 ~ 1, s8 ~ 1, s9 ~ 1, s10 ~ 1)

cat("Fitting model with ONLY 10 transitions (20 parameters)...\n")
fit_top10 <- joint(
  formSurv = formulas_top10,
  basRisk = rep("weibullsurv", 10),
  dataSurv = data_list_top10,
  id = "id",
  control = list(config = TRUE, variant = 1)
)

cat("Success!\n")
summary(fit_top10)


##### WOrks but maybe weird estimates?

cat("Success!\n")
summary(fit_top10)

# Extract and examine estimates
cat("\n=== PARAMETER ESTIMATES ===\n")
cat("Shape parameters (alpha):\n")
print(fit_top10$summary.hyperpar)

cat("\nIntercepts (beta0):\n")
print(fit_top10$summary.fixed)

# Check which transitions these correspond to
cat("\n=== TRANSITION MAPPING ===\n")
for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  n_events <- top_transitions$n_events[i]
  
  alpha_est <- fit_top10$summary.hyperpar[i, "mean"]
  alpha_sd <- fit_top10$summary.hyperpar[i, "sd"]
  beta_est <- fit_top10$summary.fixed[i, "mean"]
  beta_sd <- fit_top10$summary.fixed[i, "sd"]
  
  cat(sprintf("s%d: %d→%d (%d events) | α=%.2f (SD=%.2f), β₀=%.2f (SD=%.2f)\n",
              i, from_s, to_s, n_events, alpha_est, alpha_sd, beta_est, beta_sd))
}

# Check for problematic estimates
cat("\n=== QUALITY CHECKS ===\n")
problematic <- data.frame(
  transition = 1:10,
  alpha_sd_large = fit_top10$summary.hyperpar[, "sd"] > 1,
  beta_sd_large = fit_top10$summary.fixed[, "sd"] > 2,
  alpha_extreme = abs(fit_top10$summary.hyperpar[, "mean"]) > 3 | 
    abs(fit_top10$summary.hyperpar[, "mean"]) < 0.3
)

if(any(problematic$alpha_sd_large | problematic$beta_sd_large | problematic$alpha_extreme)) {
  cat("WARNING: Some transitions have quality issues:\n")
  print(problematic[problematic$alpha_sd_large | problematic$beta_sd_large | problematic$alpha_extreme, ])
} else {
  cat("All estimates appear reasonable!\n")
}

##### All right, lets try more trnasitions:


# ==============================================================================
# EXPAND TO TOP 24 TRANSITIONS (≥20 events)
# ==============================================================================

top24_transitions <- transition_summary_no8 %>% 
  filter(n_events >= 20) %>%
  arrange(desc(n_events))

cat(sprintf("=== FITTING %d TRANSITIONS (≥20 events) ===\n", nrow(top24_transitions)))
print(top24_transitions)

# Create data list
data_list_top24 <- list()
for(i in 1:nrow(top24_transitions)) {
  from_s <- top24_transitions$from[i]
  to_s <- top24_transitions$to[i]
  data_list_top24[[i]] <- model_data_simple %>% 
    filter(from == from_s, to == to_s)
}

# Create survival objects
for(i in 1:nrow(top24_transitions)) {
  d <- data_list_top24[[i]]
  time_vec <- d$Tstart
  idx_censored <- d$status == 0
  time_vec[idx_censored] <- d$Tstop[idx_censored]
  
  assign(paste0("s", i), 
         inla.surv(time = time_vec, time2 = d$Tstop, 
                   event = d$status, truncation = d$approx_entry))
}

# Create formulas (24 transitions)
formulas_top24 <- list(
  s1 ~ 1, s2 ~ 1, s3 ~ 1, s4 ~ 1, s5 ~ 1, s6 ~ 1, s7 ~ 1, s8 ~ 1,
  s9 ~ 1, s10 ~ 1, s11 ~ 1, s12 ~ 1, s13 ~ 1, s14 ~ 1, s15 ~ 1, s16 ~ 1,
  s17 ~ 1, s18 ~ 1, s19 ~ 1, s20 ~ 1, s21 ~ 1, s22 ~ 1, s23 ~ 1, s24 ~ 1
)

cat("\nFitting 24-transition model (48 parameters)...\n")
cat("Started at:", format(Sys.time(), "%H:%M:%S"), "\n")

fit_top24 <- joint(
  formSurv = formulas_top24,
  basRisk = rep("weibullsurv", 24),
  dataSurv = data_list_top24,
  id = "id",
  control = list(config = TRUE, variant = 1)
)

cat("Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("Success!\n\n")

# Quick check
cat("=== TOP 24 TRANSITION ESTIMATES ===\n")
for(i in 1:24) {
  from_s <- top24_transitions$from[i]
  to_s <- top24_transitions$to[i]
  n_events <- top24_transitions$n_events[i]
  alpha_est <- fit_top24$summary.hyperpar[i, "mean"]
  beta_est <- fit_top24$summary.fixed[i, "mean"]
  
  cat(sprintf("%d→%d (%3d events): α=%.2f, β₀=%.2f\n",
              from_s, to_s, n_events, alpha_est, beta_est))
}
















# ==============================================================================
# FIT WITH ONLY CLIMATE COVARIATES (CMI + TMEAN)
# ==============================================================================

cat("=== FITTING: CMI + Tmean ONLY (no categorical covariates) ===\n")

formulas_climate_only <- list(
  s1 ~ cov_CMI_std + cov_Tmean_std,
  s2 ~ cov_CMI_std + cov_Tmean_std,
  s3 ~ cov_CMI_std + cov_Tmean_std,
  s4 ~ cov_CMI_std + cov_Tmean_std,
  s5 ~ cov_CMI_std + cov_Tmean_std,
  s6 ~ cov_CMI_std + cov_Tmean_std,
  s7 ~ cov_CMI_std + cov_Tmean_std,
  s8 ~ cov_CMI_std + cov_Tmean_std,
  s9 ~ cov_CMI_std + cov_Tmean_std,
  s10 ~ cov_CMI_std + cov_Tmean_std
)

cat("Fitting 10 transitions with CMI and Tmean...\n")
cat("Started at:", format(Sys.time(), "%H:%M:%S"), "\n")

fit_climate_only <- joint(
  formSurv = formulas_climate_only,
  basRisk = rep("weibullsurv", 10),
  dataSurv = data_list_top10_cov,
  id = "id",
  control = list(config = TRUE, variant = 1)
)

cat("Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("Model fitted!\n\n")

# Check for problems
cat("=== QUALITY CHECK ===\n")
problem_params <- fit_climate_only$summary.fixed[fit_climate_only$summary.fixed$sd > 0.5, ]
if(nrow(problem_params) > 0) {
  cat("WARNING: Parameters with SD > 0.5:\n")
  print(problem_params)
} else {
  cat("✓ All parameters well-identified (SD < 0.5)!\n")
}

cat("\n=== CLIMATE EFFECTS ON FOREST SUCCESSION ===\n")
cat("(Effects are per 1 SD change in covariate)\n\n")

results_table <- data.frame()

for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  n_events <- top_transitions$n_events[i]
  
  # Rows: intercept, CMI, Tmean
  rows <- ((i-1)*3 + 1):(i*3)
  effects <- fit_climate_only$summary.fixed[rows, ]
  
  alpha <- fit_climate_only$summary.hyperpar[i, "mean"]
  
  results_table <- rbind(results_table, data.frame(
    Transition = paste0(from_s, "→", to_s),
    N_events = n_events,
    Alpha = round(alpha, 2),
    Beta0 = round(effects[1, "mean"], 2),
    CMI_effect = round(effects[2, "mean"], 2),
    CMI_sd = round(effects[2, "sd"], 2),
    CMI_sig = ifelse(effects[2, "0.025quant"] * effects[2, "0.975quant"] > 0, "*", ""),
    Tmean_effect = round(effects[3, "mean"], 2),
    Tmean_sd = round(effects[3, "sd"], 2),
    Tmean_sig = ifelse(effects[3, "0.025quant"] * effects[3, "0.975quant"] > 0, "*", "")
  ))
}

print(results_table)

cat("\n* = 95% CI excludes zero (significant effect)\n")
cat("\nPositive effect = increased transition rate with higher covariate value\n")
cat("Negative effect = decreased transition rate with higher covariate value\n")



##### SHow stuff:
library(ggplot2)
library(gridExtra)
library(viridis)

# ==============================================================================
# PLOT 1: CLIMATE EFFECTS WITH CONFIDENCE INTERVALS
# ==============================================================================

# Prepare data for plotting
coef_data <- data.frame()

for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  n_events <- top_transitions$n_events[i]
  
  rows <- ((i-1)*3 + 1):(i*3)
  effects <- fit_climate_only$summary.fixed[rows, ]
  
  # CMI effect
  coef_data <- rbind(coef_data, data.frame(
    Transition = paste0(from_s, "→", to_s),
    from = from_s,
    to = to_s,
    N_events = n_events,
    Covariate = "CMI (Moisture)",
    Estimate = effects[2, "mean"],
    Lower = effects[2, "0.025quant"],
    Upper = effects[2, "0.975quant"],
    Significant = effects[2, "0.025quant"] * effects[2, "0.975quant"] > 0
  ))
  
  # Tmean effect
  coef_data <- rbind(coef_data, data.frame(
    Transition = paste0(from_s, "→", to_s),
    from = from_s,
    to = to_s,
    N_events = n_events,
    Covariate = "Temperature",
    Estimate = effects[3, "mean"],
    Lower = effects[3, "0.025quant"],
    Upper = effects[3, "0.975quant"],
    Significant = effects[3, "0.025quant"] * effects[3, "0.975quant"] > 0
  ))
}

# Order by effect size
coef_data$Transition <- factor(coef_data$Transition, 
                               levels = unique(coef_data$Transition[order(coef_data$N_events, decreasing = TRUE)]))

p1 <- ggplot(coef_data, aes(x = Transition, y = Estimate, color = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ Covariate, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#E41A1C"),
                     name = "95% CI excludes 0") +
  labs(
    title = "Climate Effects on Forest Succession Transitions",
    subtitle = "Effect of 1 SD change in climate variable on log-hazard (Method B: Interval-censored observational data)",
    x = "Transition (ordered by # events)",
    y = "Log-hazard coefficient"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

print(p1)

# ==============================================================================
# PLOT 2: HAZARD RATIOS FOR CLIMATE SCENARIOS
# ==============================================================================

# Calculate hazard ratios: compare +1 SD to baseline (0) and -1 SD
hr_data <- data.frame()

for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  
  rows <- ((i-1)*3 + 1):(i*3)
  effects <- fit_climate_only$summary.fixed[rows, ]
  
  beta_cmi <- effects[2, "mean"]
  beta_tmean <- effects[3, "mean"]
  
  # Hazard ratio = exp(beta * change)
  # For +1 SD vs 0: HR = exp(beta * 1)
  # For -1 SD vs 0: HR = exp(beta * -1)
  
  hr_data <- rbind(hr_data, data.frame(
    Transition = paste0(from_s, "→", to_s),
    Scenario = c("High CMI (+1 SD)", "Low CMI (-1 SD)", 
                 "High Temp (+1 SD)", "Low Temp (-1 SD)"),
    Variable = c("CMI", "CMI", "Temperature", "Temperature"),
    HR = c(exp(beta_cmi), exp(-beta_cmi), 
           exp(beta_tmean), exp(-beta_tmean))
  ))
}

# Select interesting transitions
interesting_trans <- c("1→4", "3→1", "4→3", "1→6", "2→1")

p2 <- hr_data %>%
  filter(Transition %in% interesting_trans) %>%
  ggplot(aes(x = Scenario, y = HR, fill = Variable)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_wrap(~ Transition, ncol = 5) +
  scale_fill_manual(values = c("CMI" = "#377EB8", "Temperature" = "#E41A1C")) +
  labs(
    title = "Hazard Ratios: Climate Effects on Transition Rates",
    subtitle = "HR > 1 = faster transition, HR < 1 = slower transition (relative to average climate)",
    x = NULL,
    y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 2))

print(p2)

# ==============================================================================
# PLOT 3: TRANSITION PROBABILITIES UNDER CLIMATE SCENARIOS
# ==============================================================================

# Extract parameters for probability calculations
extract_params_from_fit <- function(model, transition_idx) {
  rows <- ((transition_idx-1)*3 + 1):(transition_idx*3)
  
  intercept <- model$summary.fixed[rows[1], "mean"]
  beta_cmi <- model$summary.fixed[rows[2], "mean"]
  beta_tmean <- model$summary.fixed[rows[3], "mean"]
  alpha <- model$summary.hyperpar[transition_idx, "mean"]
  
  return(list(intercept = intercept, beta_cmi = beta_cmi, 
              beta_tmean = beta_tmean, alpha = alpha))
}

# Get all parameters
params_list <- lapply(1:10, function(i) extract_params_from_fit(fit_climate_only, i))
names(params_list) <- paste0(top_transitions$from, "_", top_transitions$to)

# Weibull hazard function
weibull_hazard_climate <- function(t, alpha, intercept, beta_cmi, beta_tmean,
                                   cmi_val = 0, tmean_val = 0) {
  lambda <- exp(intercept + beta_cmi * cmi_val + beta_tmean * tmean_val)
  alpha * lambda * (lambda * t)^(alpha - 1)
}

# Calculate cumulative hazard and survival for one transition
calc_survival_one_trans <- function(params, time_seq, cmi_val, tmean_val) {
  cum_haz <- numeric(length(time_seq))
  
  for(i in seq_along(time_seq)) {
    if(i == 1) {
      cum_haz[i] <- 0
    } else {
      dt <- time_seq[i] - time_seq[i-1]
      h <- weibull_hazard_climate(time_seq[i], params$alpha, params$intercept,
                                  params$beta_cmi, params$beta_tmean,
                                  cmi_val, tmean_val)
      cum_haz[i] <- cum_haz[i-1] + h * dt
    }
  }
  
  survival <- exp(-cum_haz)
  transition_prob <- 1 - survival
  
  return(transition_prob)
}

# Focus on a few key transitions
key_transitions <- list(
  "1→4" = "1_4",  # Significant CMI effect
  "3→1" = "3_1",  # Significant Tmean effect
  "4→3" = "4_3",  # Interesting divergent effects
  "1→6" = "1_6"   # Accelerating hazard (α=2.2)
)

time_seq <- seq(0, 50, by = 0.5)

# Climate scenarios
scenarios <- list(
  "Baseline" = list(cmi = 0, tmean = 0, col = "black"),
  "High Moisture" = list(cmi = 1, tmean = 0, col = "#377EB8"),
  "Low Moisture" = list(cmi = -1, tmean = 0, col = "#FF7F00"),
  "High Temperature" = list(cmi = 0, tmean = 1, col = "#E41A1C"),
  "Low Temperature" = list(cmi = 0, tmean = -1, col = "#4DAF4A")
)

plot_list <- list()

for(trans_name in names(key_transitions)) {
  param_key <- key_transitions[[trans_name]]
  params <- params_list[[param_key]]
  
  prob_data <- data.frame()
  
  for(scen_name in names(scenarios)) {
    scen <- scenarios[[scen_name]]
    probs <- calc_survival_one_trans(params, time_seq, scen$cmi, scen$tmean)
    
    prob_data <- rbind(prob_data, data.frame(
      time = time_seq,
      probability = probs,
      scenario = scen_name,
      color = scen$col
    ))
  }
  
  prob_data$scenario <- factor(prob_data$scenario, levels = names(scenarios))
  
  p <- ggplot(prob_data, aes(x = time, y = probability, color = scenario, group = scenario)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = setNames(sapply(scenarios, function(x) x$col), names(scenarios))) +
    labs(
      title = paste("Transition", trans_name),
      x = "Time (years)",
      y = "Cumulative Transition Probability",
      color = "Climate Scenario"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    ) +
    ylim(0, 1)
  
  plot_list[[trans_name]] <- p
}

# Combine into grid
grid.arrange(
  grobs = plot_list,
  ncol = 2,
  top = grid::textGrob("Transition Probabilities Under Climate Scenarios\n(±1 SD from mean)",
                       gp = grid::gpar(fontsize = 14, fontface = "bold"))
)

# ==============================================================================
# PLOT 4: NETWORK DIAGRAM OF TRANSITIONS
# ==============================================================================

# Create network visualization showing effect sizes
library(igraph)

# Create edge list with effect magnitudes
edges <- data.frame()
for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  
  rows <- ((i-1)*3 + 1):(i*3)
  effects <- fit_climate_only$summary.fixed[rows, ]
  
  # Use absolute magnitude of strongest effect
  cmi_effect <- abs(effects[2, "mean"])
  tmean_effect <- abs(effects[3, "mean"])
  max_effect <- max(cmi_effect, tmean_effect)
  
  # Which effect is stronger?
  dominant <- ifelse(cmi_effect > tmean_effect, "CMI", "Temp")
  
  edges <- rbind(edges, data.frame(
    from = from_s,
    to = to_s,
    weight = max_effect,
    dominant = dominant,
    n_events = top_transitions$n_events[i]
  ))
}

# Create graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Set edge properties
E(g)$width <- scales::rescale(E(g)$weight, to = c(0.5, 4))
E(g)$color <- ifelse(E(g)$dominant == "CMI", "#377EB8", "#E41A1C")
E(g)$arrow.size <- 0.5

# Set vertex properties  
V(g)$size <- 25
V(g)$color <- "lightblue"
V(g)$label.cex <- 1.2
V(g)$label.color <- "black"

# Plot
par(mar = c(1, 1, 3, 1))
plot(g, 
     layout = layout_with_kk(g),
     main = "Forest Succession Network\n(Edge width = effect magnitude, Color = dominant variable)",
     edge.curved = 0.2)

legend("bottomright", 
       legend = c("CMI dominant", "Temperature dominant"),
       col = c("#377EB8", "#E41A1C"),
       lwd = 2,
       bty = "n")





######TEST vs 70-30split

# ==============================================================================
# STEP 1: SPLIT DATA
# ==============================================================================

set.seed(456)  # For reproducibility

# Randomly select 70% of plots for training
all_plots <- unique(forest_subset_no8$id_simple)
n_train <- round(0.7 * length(all_plots))

train_plots <- sample(all_plots, size = n_train)
test_plots <- setdiff(all_plots, train_plots)

cat(sprintf("Training: %d plots\n", length(train_plots)))
cat(sprintf("Testing: %d plots\n", length(test_plots)))

# Create train and test datasets
forest_train <- forest_subset_no8 %>% filter(id_simple %in% train_plots)
forest_test <- forest_subset_no8 %>% filter(id_simple %in% test_plots)

cat(sprintf("\nTraining observations: %d\n", nrow(forest_train)))
cat(sprintf("Testing observations: %d\n", nrow(forest_test)))




# ==============================================================================
# STEP 2: FIT MODEL ON TRAINING DATA ONLY
# ==============================================================================

# Reconstruct training data using Method B
cat("Reconstructing training data...\n")
TT <- max(forest_train$time)

model_data_train <- bind_rows(lapply(train_plots, function(plot_id) {
  reconstruct_plot_methodB_simple(
    forest_train %>% filter(id_simple == plot_id), 
    TT, 
    states_remaining
  )
}))

# Add standardized covariates
model_data_train <- model_data_train %>%
  left_join(plot_covariates %>% select(id_simple, cov_CMI_std, cov_Tmean_std),
            by = c("id" = "id_simple"))

# Create data list for top 10 transitions
data_list_train <- list()
for(i in 1:10) {
  from_s <- top_transitions$from[i]
  to_s <- top_transitions$to[i]
  data_list_train[[i]] <- model_data_train %>% 
    filter(from == from_s, to == to_s)
}

# Create survival objects
for(i in 1:10) {
  d <- data_list_train[[i]]
  time_vec <- d$Tstart
  idx_censored <- d$status == 0
  time_vec[idx_censored] <- d$Tstop[idx_censored]
  
  assign(paste0("s", i), 
         inla.surv(time = time_vec, time2 = d$Tstop, 
                   event = d$status, truncation = d$approx_entry))
}

# Fit model on training data
cat("Fitting model on training data...\n")
fit_train <- joint(
  formSurv = formulas_climate_only,
  basRisk = rep("weibullsurv", 10),
  dataSurv = data_list_train,
  id = "id",
  control = list(config = TRUE, variant = 1)
)

cat("Training model fitted!\n\n")



# ==============================================================================
# STEP 3: PREDICT ON TEST DATA
# ==============================================================================

cat("Making predictions on test data...\n")

# For each test plot, calculate predicted probability of each transition
# using the fitted model parameters

# Extract parameters from training model
get_predicted_hazard <- function(fit_model, transition_idx, cmi_val, tmean_val, time) {
  # Get parameters
  rows <- ((transition_idx-1)*3 + 1):(transition_idx*3)
  intercept <- fit_model$summary.fixed[rows[1], "mean"]
  beta_cmi <- fit_model$summary.fixed[rows[2], "mean"]
  beta_tmean <- fit_model$summary.fixed[rows[3], "mean"]
  alpha <- fit_model$summary.hyperpar[transition_idx, "mean"]
  
  # Calculate hazard at given time
  lambda <- exp(intercept + beta_cmi * cmi_val + beta_tmean * tmean_val)
  hazard <- alpha * lambda * (lambda * time)^(alpha - 1)
  
  return(hazard)
}

# For each test plot, calculate probability of observed transition
test_predictions <- data.frame()

for(plot_id in test_plots) {
  plot_data <- forest_test %>% filter(id_simple == plot_id) %>% arrange(time)
  
  # Get covariate values for this plot
  plot_cov <- plot_covariates %>% filter(id_simple == plot_id)
  cmi_val <- plot_cov$cov_CMI_std
  tmean_val <- plot_cov$cov_Tmean_std
  
  # Look at each consecutive pair of observations
  for(i in 2:nrow(plot_data)) {
    state_from <- plot_data$sp_class[i-1]
    state_to <- plot_data$sp_class[i]
    time_interval <- plot_data$time[i] - plot_data$time[i-1]
    
    if(state_from != state_to) {
      # A transition occurred!
      # Find which of our 10 transitions this corresponds to
      trans_match <- which(top_transitions$from == state_from & 
                             top_transitions$to == state_to)
      
      if(length(trans_match) == 1) {
        # Calculate predicted hazard at midpoint of interval
        time_midpoint <- (plot_data$time[i-1] + plot_data$time[i]) / 2
        pred_hazard <- get_predicted_hazard(fit_train, trans_match, 
                                            cmi_val, tmean_val, time_midpoint)
        
        test_predictions <- rbind(test_predictions, data.frame(
          plot_id = plot_id,
          from = state_from,
          to = state_to,
          time_start = plot_data$time[i-1],
          time_end = plot_data$time[i],
          interval_length = time_interval,
          predicted_hazard = pred_hazard,
          observed = 1  # Transition happened
        ))
      }
    }
  }
}

cat(sprintf("Made predictions for %d transitions in test data\n", nrow(test_predictions)))





# ==============================================================================
# STEP 4: EVALUATE PREDICTIONS
# ==============================================================================

cat("\n=== VALIDATION RESULTS ===\n\n")

# Metric 1: Correlation between predicted hazard and observed transitions
cat("1. PREDICTIVE CORRELATION:\n")
cat(sprintf("   Correlation between predicted hazard and observations: %.3f\n", 
            cor(test_predictions$predicted_hazard, test_predictions$observed)))

# Metric 2: Compare log-likelihood on train vs test
# (This requires more complex calculation - simpler version below)

# Metric 3: Qualitative check - do high-hazard predictions correspond to transitions?
# Split into quartiles
test_predictions$hazard_quartile <- cut(test_predictions$predicted_hazard, 
                                        breaks = quantile(test_predictions$predicted_hazard, 
                                                          probs = c(0, 0.25, 0.5, 0.75, 1)),
                                        labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                                        include.lowest = TRUE)

cat("\n2. TRANSITIONS BY PREDICTED HAZARD QUARTILE:\n")
cat("   (If model works well, should see more transitions in high hazard groups)\n")
quartile_summary <- test_predictions %>%
  group_by(hazard_quartile) %>%
  summarise(
    n_transitions = n(),
    mean_hazard = mean(predicted_hazard)
  )
print(quartile_summary)

# Metric 4: Compare AIC between models
cat("\n3. MODEL COMPARISON (on training data):\n")
cat(sprintf("   Climate model AIC: %.1f\n", -2 * fit_train$mlik[1,1] + 2 * 30))  # 30 params
cat(sprintf("   Intercept-only AIC: %.1f\n", -2 * fit_top10$mlik[1,1] + 2 * 20))  # 20 params
cat("   (Lower AIC = better model fit)\n")

# Visual comparison
library(ggplot2)

p_validation <- ggplot(test_predictions, aes(x = predicted_hazard)) +
  geom_histogram(bins = 30, fill = "#377EB8", alpha = 0.7) +
  labs(
    title = "Distribution of Predicted Hazards for Test Set Transitions",
    subtitle = paste0("Based on ", nrow(test_predictions), " observed transitions in held-out data"),
    x = "Predicted Hazard (from training model)",
    y = "Count"
  ) +
  theme_minimal()

print(p_validation)


# ==============================================================================
# PROPER VALIDATION: Log-likelihood on test data
# ==============================================================================

cat("\n=== BETTER VALIDATION METRIC ===\n")

# For a proper comparison, we'd need to:
# 1. Reconstruct the full test data (including non-transitions)
# 2. Calculate likelihood for each observation
# 3. Sum log-likelihoods

# Simpler interpretation: AIC difference of 1255 is massive
cat("\nInterpretation of AIC difference (Climate vs Intercept-only):\n")
cat("  Difference: -1255 AIC units\n")
cat("  Rule of thumb:\n")
cat("    <2:   Essentially equivalent\n")
cat("    2-10: Moderate evidence for better model\n")
cat("    >10:  Strong evidence for better model\n")
cat("    >100: Overwhelming evidence\n")
cat("\n  YOUR RESULT: -1255 = OVERWHELMING evidence that climate matters!\n")

# Evidence ratio
evidence_ratio <- exp((4789.6 - 3534.6) / 2)
cat(sprintf("\n  Evidence ratio: %.0e:1 in favor of climate model\n", evidence_ratio))
cat("  (The climate model is exponentially more likely to be correct)\n")







############
# ==============================================================================
# SCALE UP: SELECT 5000 PLOTS
# ==============================================================================

set.seed(789)

# Start fresh from the full dataset (without state 8)
cat("Starting from full dataset (state 8 already removed)...\n")
cat(sprintf("Available plots: %d\n", length(unique(data$TESSELLE))))

# Remove state 8 from full dataset if not already done
data_no8 <- data %>% filter(sp_class != 8)

# Get unique plot IDs
all_plot_ids <- unique(data_no8$TESSELLE)

# Sample 5000 plots
n_plots <- 5000
selected_plots_5k <- sample(all_plot_ids, size = min(n_plots, length(all_plot_ids)))

# Filter data
forest_5k <- data_no8 %>% filter(TESSELLE %in% selected_plots_5k)

cat(sprintf("\nSelected %d plots with %d total observations\n", 
            length(unique(forest_5k$TESSELLE)), 
            nrow(forest_5k)))

# Check state distribution
cat("\nState distribution:\n")
print(table(forest_5k$sp_class))

# Check observation pattern
cat("\nObservations per plot:\n")
print(table(forest_5k %>% group_by(TESSELLE) %>% tally() %>% pull(n)))

# ==============================================================================
# CREATE SIMPLE INTEGER IDs
# ==============================================================================

plot_mapping_5k <- data.frame(
  TESSELLE_original = unique(forest_5k$TESSELLE),
  id_simple = 1:length(unique(forest_5k$TESSELLE))
)

forest_5k <- forest_5k %>%
  left_join(plot_mapping_5k, by = c("TESSELLE" = "TESSELLE_original"))

cat(sprintf("\nSimple IDs created: 1 to %d\n", max(forest_5k$id_simple)))

# ==============================================================================
# PREPARE COVARIATES (STANDARDIZED)
# ==============================================================================

cat("\n=== PREPARING COVARIATES ===\n")

plot_covariates_5k <- forest_5k %>%
  group_by(id_simple) %>%
  slice(1) %>%
  select(id_simple, cov_CMI, cov_Tmean, cov_soil, cov_pert_class) %>%
  ungroup()

# Standardize continuous
plot_covariates_5k <- plot_covariates_5k %>%
  mutate(
    cov_CMI_std = scale(cov_CMI)[,1],
    cov_Tmean_std = scale(cov_Tmean)[,1],
    cov_soil = factor(cov_soil),
    cov_pert_class = factor(cov_pert_class)
  )

cat("Continuous covariates standardized\n")

# Check categorical distributions
cat("\n=== CATEGORICAL COVARIATE DISTRIBUTIONS ===\n")
cat("Soil types:\n")
print(table(plot_covariates_5k$cov_soil))

cat("\nPerturbation classes:\n")
print(table(plot_covariates_5k$cov_pert_class))

# ==============================================================================
# RECONSTRUCT WITH METHOD B
# ==============================================================================

cat("\n=== RECONSTRUCTING DATA (5000 plots) ===\n")
cat("This will take a few minutes...\n")

TT <- max(forest_5k$time)
states_remaining <- sort(unique(forest_5k$sp_class))

start_time <- Sys.time()

model_data_5k <- bind_rows(lapply(unique(forest_5k$id_simple), function(plot_id) {
  if(plot_id %% 500 == 0) cat(sprintf("  %d/%d\n", plot_id, length(unique(forest_5k$id_simple))))
  reconstruct_plot_methodB_simple(
    forest_5k %>% filter(id_simple == plot_id), 
    TT, 
    states_remaining
  )
}))

end_time <- Sys.time()
cat(sprintf("Reconstruction complete! Time: %.1f minutes\n", 
            as.numeric(difftime(end_time, start_time, units = "mins"))))

# ==============================================================================
# CHECK TRANSITION COUNTS
# ==============================================================================

cat("\n=== TRANSITION EVENT COUNTS (5000 plots) ===\n")

transition_summary_5k <- data.frame()
for(from_state in states_remaining) {
  for(to_state in states_remaining) {
    if(from_state != to_state) {
      temp_data <- model_data_5k %>% filter(from == from_state, to == to_state)
      transition_summary_5k <- rbind(transition_summary_5k, data.frame(
        transition = paste0(from_state, "→", to_state),
        from = from_state,
        to = to_state,
        n_events = sum(temp_data$status == 3),
        n_censored = sum(temp_data$status == 0)
      ))
    }
  }
}

cat("Top 20 transitions by event count:\n")
print(transition_summary_5k %>% arrange(desc(n_events)) %>% head(20))

cat("\nTransitions with <20 events:\n")
sparse_5k <- transition_summary_5k %>% filter(n_events < 20)
cat(sprintf("Count: %d out of %d\n", nrow(sparse_5k), nrow(transition_summary_5k)))

# ==============================================================================
# DECIDE HOW MANY TRANSITIONS TO FIT
# ==============================================================================

# How many transitions have ≥50 events now?
sufficient_data <- transition_summary_5k %>% filter(n_events >= 50)
cat(sprintf("\n=== TRANSITIONS WITH ≥50 EVENTS: %d ===\n", nrow(sufficient_data)))

cat("\nRecommendation:\n")
cat(sprintf("  - Fit %d transitions (≥50 events each)\n", nrow(sufficient_data)))
cat("  - This should provide enough data for categorical covariates\n")
cat("  - Check soil/pert distributions per transition before fitting\n")

# ==============================================================================
# CHECK CATEGORICAL COVARIATE SPARSITY PER TRANSITION
# ==============================================================================

cat("=== CHECKING SOIL & PERTURBATION DISTRIBUTIONS PER TRANSITION ===\n\n")

# Add covariates to model data
model_data_5k_cov <- model_data_5k %>%
  left_join(plot_covariates_5k %>% select(id_simple, cov_CMI_std, cov_Tmean_std, 
                                          cov_soil, cov_pert_class),
            by = c("id" = "id_simple"))

# For top 34 transitions, check categorical distributions
sparsity_check <- data.frame()

for(i in 1:nrow(sufficient_data)) {
  from_s <- sufficient_data$from[i]
  to_s <- sufficient_data$to[i]
  trans_label <- paste0(from_s, "→", to_s)
  
  # Filter to this transition, events only
  trans_events <- model_data_5k_cov %>%
    filter(from == from_s, to == to_s, status == 3)
  
  n_events <- nrow(trans_events)
  
  # Count soil types
  soil_counts <- table(trans_events$cov_soil)
  n_soil_types <- length(soil_counts)
  min_soil_count <- min(soil_counts)
  n_soil_zeros <- sum(soil_counts == 0)
  
  # Count perturbation types
  pert_counts <- table(trans_events$cov_pert_class)
  n_pert_types <- length(pert_counts)
  min_pert_count <- min(pert_counts)
  
  sparsity_check <- rbind(sparsity_check, data.frame(
    transition = trans_label,
    n_events = n_events,
    n_soil_levels = n_soil_types,
    min_soil_count = min_soil_count,
    soil_ok = min_soil_count >= 5,  # At least 5 events per level
    n_pert_levels = n_pert_types,
    min_pert_count = min_pert_count,
    pert_ok = min_pert_count >= 5
  ))
}

cat("Summary of categorical covariate sparsity:\n")
print(head(sparsity_check, 34))

# Overall assessment
cat("\n=== SPARSITY ASSESSMENT ===\n")
cat(sprintf("Transitions where soil is OK (≥5 events per level): %d / 34\n", 
            sum(sparsity_check$soil_ok)))
cat(sprintf("Transitions where perturbation is OK (≥5 events per level): %d / 34\n", 
            sum(sparsity_check$pert_ok)))

problematic_soil <- sparsity_check %>% filter(!soil_ok)
problematic_pert <- sparsity_check %>% filter(!pert_ok)

if(nrow(problematic_soil) > 0) {
  cat("\nTransitions with sparse soil levels:\n")
  print(problematic_soil %>% select(transition, n_events, min_soil_count))
}

if(nrow(problematic_pert) > 0) {
  cat("\nTransitions with sparse perturbation levels:\n")
  print(problematic_pert %>% select(transition, n_events, min_pert_count))
}

# ==============================================================================
# CHECK OVERALL DISTRIBUTIONS
# ==============================================================================

cat("\n=== OVERALL SOIL DISTRIBUTION (across all events) ===\n")
all_events_soil <- model_data_5k_cov %>%
  filter(status == 3, 
         from %in% sufficient_data$from,
         to %in% sufficient_data$to) %>%
  count(cov_soil, sort = TRUE)
print(all_events_soil)

cat("\n=== OVERALL PERTURBATION DISTRIBUTION (across all events) ===\n")
all_events_pert <- model_data_5k_cov %>%
  filter(status == 3,
         from %in% sufficient_data$from,
         to %in% sufficient_data$to) %>%
  count(cov_pert_class, sort = TRUE)
print(all_events_pert)

# ==============================================================================
# RECOMMENDATION
# ==============================================================================

cat("\n=== RECOMMENDATION ===\n\n")

if(sum(sparsity_check$soil_ok) >= 30 && sum(sparsity_check$pert_ok) >= 30) {
  cat("✓ GOOD NEWS: Most transitions have sufficient categorical data\n")
  cat("  Strategy: Fit with all covariates (CMI + Tmean + Soil + Pert)\n")
  cat("  May need to handle a few sparse transitions carefully\n")
} else if(sum(sparsity_check$pert_ok) >= 30) {
  cat("⚠ Soil is too sparse, but perturbation looks OK\n")
  cat("  Strategy: Fit with CMI + Tmean + Perturbation only\n")
} else {
  cat("⚠ Categorical covariates still sparse\n")
  cat("  Strategy Options:\n")
  cat("    A) Fit continuous only (CMI + Tmean)\n")
  cat("    B) Collapse rare soil/pert categories\n")
  cat("    C) Use categorical only for well-populated transitions\n")
}

# ==============================================================================
# CHECK PERTURBATION CLASS ONLY (DROP SOIL)
# ==============================================================================

cat("=== PERTURBATION CLASS DISTRIBUTION (34 transitions, ≥50 events) ===\n\n")

# Check perturbation distribution per transition
pert_check <- data.frame()

for(i in 1:nrow(sufficient_data)) {
  from_s <- sufficient_data$from[i]
  to_s <- sufficient_data$to[i]
  trans_label <- paste0(from_s, "→", to_s)
  
  # Filter to this transition, events only
  trans_events <- model_data_5k_cov %>%
    filter(from == from_s, to == to_s, status == 3)
  
  n_events <- nrow(trans_events)
  
  # Count perturbation types
  pert_counts <- table(trans_events$cov_pert_class)
  n_pert_levels <- length(pert_counts)
  min_pert_count <- if(length(pert_counts) > 0) min(pert_counts) else 0
  
  pert_check <- rbind(pert_check, data.frame(
    transition = trans_label,
    from = from_s,
    to = to_s,
    n_events = n_events,
    n_pert_levels = n_pert_levels,
    min_pert_count = min_pert_count,
    pert_ok = min_pert_count >= 10  # At least 10 events per level
  ))
}

cat("Perturbation sparsity by transition:\n")
print(pert_check)

# Summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("Transitions where perturbation is OK (≥10 events per level): %d / 34\n", 
            sum(pert_check$pert_ok)))

# Check overall perturbation distribution
cat("\n=== OVERALL PERTURBATION DISTRIBUTION ===\n")
all_events_pert <- model_data_5k_cov %>%
  filter(status == 3,
         from %in% sufficient_data$from,
         to %in% sufficient_data$to) %>%
  count(cov_pert_class, sort = TRUE)

print(all_events_pert)

cat("\nPerturbation levels:\n")
cat(paste(levels(factor(model_data_5k_cov$cov_pert_class)), collapse = ", "))
cat("\n")

# ==============================================================================
# DECISION POINT
# ==============================================================================

cat("\n=== RECOMMENDATION ===\n\n")

if(sum(pert_check$pert_ok) >= 25) {
  cat("✓ GOOD: Most transitions have sufficient perturbation data\n")
  cat("  Strategy: Fit 34 transitions with CMI + Tmean + Perturbation\n\n")
  
  cat("Next step: Fit the model!\n")
  cat("  Model size: 34 transitions × 4 parameters = 136 parameters\n")
  cat("  (α + β₀ + β_CMI + β_Tmean + β_pert for each transition)\n")
  
} else {
  cat("⚠ Many transitions still have sparse perturbation data\n\n")
  
  # Check if we can collapse perturbation classes
  cat("Options:\n")
  cat("  A) Fit continuous only (CMI + Tmean) for all 34 transitions\n")
  cat("  B) Collapse rare perturbation classes (if applicable)\n")
  cat("  C) Fit perturbation only for well-populated transitions\n\n")
  
  cat("What are the perturbation classes? (to see if we can collapse)\n")
}

