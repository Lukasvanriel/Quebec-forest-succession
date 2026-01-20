#!/usr/bin/env Rscript
# ==============================================================================
# FOREST SUCCESSION DATA PREPARATION - FULL DATASET
# Method B Reconstruction with Parallelization
# ==============================================================================

# Load required libraries
library(tidyverse)
library(data.table)
library(parallel)
library(here)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths
DATA_PATH <- here("Data", "BTE", "bte_msm_ready.csv")  # INPUT: Your raw data
OUTPUT_DIR <- here("Data-Output", "INLA-realData")        # OUTPUT: Where to save results
PROGRESS_FILE <- file.path(OUTPUT_DIR, "progress_log.txt")
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

# Processing parameters
N_CORES <- 8            # Number of cores to use (adjust for your cluster)
CHECKPOINT_INTERVAL <- 100  # Save checkpoint every N plots #10000
SEED <- 42

# Create directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PROGRESS LOGGING FUNCTION
# ==============================================================================

log_progress <- function(message, progress_file = PROGRESS_FILE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_message <- sprintf("[%s] %s\n", timestamp, message)
  cat(log_message)
  cat(log_message, file = progress_file, append = TRUE)
}

# Initialize progress file
cat("", file = PROGRESS_FILE)  # Clear file
log_progress("=== FOREST SUCCESSION DATA PREPARATION - STARTED ===")
log_progress(sprintf("Using %d cores", N_CORES))

# ==============================================================================
# LOAD DATA AND SELECT COLUMNS
# ==============================================================================

log_progress("Loading raw data...")
data_raw <- fread(DATA_PATH)
# Remove first column if it's an index
if(colnames(data_raw)[1] == "V1") data_raw[, V1 := NULL]
log_progress(sprintf("Loaded %d observations with %d columns", 
                     nrow(data_raw), ncol(data_raw)))

# ==============================================================================
# SELECT ONLY NECESSARY COLUMNS (MEMORY OPTIMIZATION)
# ==============================================================================

log_progress("Selecting relevant columns to reduce memory footprint...")

required_cols <- c(
  # ID and geography
  "TESSELLE", "LATIT", "LONGI",
  
  # Ecological regions (for reference)
  "SDOM_ECO", "SREG_ECO",
  
  # State and time
  "sp_class", "time", "dom_sp",
  
  # Climate covariates
  "cov_CMI", "cov_Tmean",
  
  # Categorical covariates
  "cov_soil", "cov_pert_class", "cov_pert_sev", "cov_time_pert"
)

# Check which columns exist
missing_cols <- setdiff(required_cols, colnames(data_raw))
if(length(missing_cols) > 0) {
  log_progress(sprintf("WARNING: Missing columns: %s", 
                       paste(missing_cols, collapse = ", ")))
  log_progress("Continuing with available columns only...")
}

existing_cols <- intersect(required_cols, colnames(data_raw))
data_raw <- data_raw[, existing_cols]

log_progress(sprintf("Reduced to %d columns (from ~35)", ncol(data_raw)))
log_progress(sprintf("Estimated memory: %.1f GB", 
                     as.numeric(object.size(data_raw)) / 1024^3))


# ==============================================================================
# DATA PREPARATION
# ==============================================================================

log_progress("Preparing data...")

# Remove state 8
data_clean <- data_raw %>% filter(sp_class != 8)
log_progress(sprintf("After removing state 8: %d observations", nrow(data_clean)))

# Create simple integer IDs
unique_plots <- unique(data_clean$TESSELLE)
n_plots <- length(unique_plots)
log_progress(sprintf("Total plots to process: %d", n_plots))

plot_mapping <- data.table(
  TESSELLE_original = unique_plots,
  id_simple = 1:n_plots
)

data_clean <- as.data.table(data_clean)
data_clean <- merge(data_clean, plot_mapping, 
                    by.x = "TESSELLE", by.y = "TESSELLE_original", 
                    all.x = TRUE)

# Prepare covariates at plot level
log_progress("Preparing covariates...")
plot_covariates <- data_clean[, .(
  cov_CMI = first(cov_CMI),
  cov_Tmean = first(cov_Tmean),
  cov_soil = first(cov_soil),
  cov_pert_class = first(cov_pert_class)
), by = id_simple]

# Standardize continuous covariates
plot_covariates[, cov_CMI_std := scale(cov_CMI)[,1]]
plot_covariates[, cov_Tmean_std := scale(cov_Tmean)[,1]]

log_progress("Covariates prepared and standardized")

# Get time and state info
TT <- max(data_clean$time)
states_remaining <- sort(unique(data_clean$sp_class))
log_progress(sprintf("States: %s", paste(states_remaining, collapse = ", ")))
log_progress(sprintf("Max time: %d years", TT))

# ==============================================================================
# METHOD B RECONSTRUCTION FUNCTION (OPTIMIZED)
# ==============================================================================

reconstruct_plot_methodB_fast <- function(plot_data, TT, states_remaining) {
  # Get plot info
  plot_id <- plot_data$id_simple[1]
  
  # Order by time
  setorder(plot_data, time)
  obs_times <- plot_data$time
  obs_states <- plot_data$sp_class
  
  # Safety checks
  if(length(obs_times) < 2 || any(is.na(obs_states))) return(NULL)
  
  # Pre-allocate list
  max_rows <- length(obs_times) * length(states_remaining) * 2
  result_list <- vector("list", max_rows)
  counter <- 1
  
  s0 <- obs_states[1]
  approx_entry_time <- obs_times[1]
  
  # Loop through observation periods
  for(k in 2:length(obs_times)) {
    t0 <- obs_times[k-1]
    t1 <- obs_times[k]
    s1 <- obs_states[k]
    
    possible_to_states <- setdiff(states_remaining, s0)
    
    if(s0 != s1) {
      # Transition occurred
      
      # Compacted competing risk rows
      for(to_state in possible_to_states) {
        if(to_state != s1) {
          result_list[[counter]] <- list(
            id = plot_id,
            Tstart = approx_entry_time,
            Tstop = t1,
            from = s0,
            to = to_state,
            status = 0L,
            approx_entry = approx_entry_time
          )
          counter <- counter + 1
        }
      }
      
      # Actual transition
      result_list[[counter]] <- list(
        id = plot_id,
        Tstart = t0,
        Tstop = t1,
        from = s0,
        to = s1,
        status = 3L,
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
    result_list[[counter]] <- list(
      id = plot_id,
      Tstart = approx_entry_time,
      Tstop = TT,
      from = s0,
      to = to_state,
      status = 0L,
      approx_entry = approx_entry_time
    )
    counter <- counter + 1
  }
  
  # Convert to data.table
  rbindlist(result_list[1:(counter-1)])
}

# ==============================================================================
# PARALLEL PROCESSING WITH PROGRESS TRACKING
# ==============================================================================

log_progress("Starting parallel reconstruction...")
log_progress(sprintf("Splitting %d plots into %d chunks", n_plots, N_CORES))

# Split plots into chunks for parallel processing
plot_ids <- unique(data_clean$id_simple)
chunks <- split(plot_ids, ceiling(seq_along(plot_ids) / ceiling(length(plot_ids) / N_CORES)))

log_progress(sprintf("Created %d chunks", length(chunks)))

# Function to process a chunk with progress tracking
process_chunk <- function(chunk_ids, chunk_num, total_chunks) {
  # Create chunk-specific progress file
  chunk_progress_file <- file.path(OUTPUT_DIR, sprintf("chunk_%03d_progress.txt", chunk_num))
  
  write(sprintf("Chunk %d/%d: Starting with %d plots", 
                chunk_num, total_chunks, length(chunk_ids)), 
        chunk_progress_file)
  
  results <- vector("list", length(chunk_ids))
  
  for(i in seq_along(chunk_ids)) {
    plot_id <- chunk_ids[i]
    
    # Get plot data
    plot_data <- data_clean[id_simple == plot_id]
    
    # Reconstruct
    results[[i]] <- reconstruct_plot_methodB_fast(plot_data, TT, states_remaining)
    
    # Progress update every 100 plots
    if(i %% 100 == 0) {
      write(sprintf("Chunk %d/%d: Processed %d/%d plots (%.1f%%)", 
                    chunk_num, total_chunks, i, length(chunk_ids), 
                    100 * i / length(chunk_ids)), 
            chunk_progress_file, append = TRUE)
    }
  }
  
  # Combine results for this chunk
  chunk_result <- rbindlist(results, use.names = TRUE, fill = TRUE)
  
  write(sprintf("Chunk %d/%d: COMPLETE - %d rows generated", 
                chunk_num, total_chunks, nrow(chunk_result)), 
        chunk_progress_file, append = TRUE)
  
  return(chunk_result)
}

# ==============================================================================
# RUN PARALLEL PROCESSING
# ==============================================================================

log_progress("Launching parallel workers...")

start_time <- Sys.time()

# Determine if we can use mclapply (Unix-like systems)
use_mclapply <- .Platform$OS.type == "unix"

if(use_mclapply) {
  log_progress(sprintf("Using mclapply (fork-based parallelization) with %d cores", N_CORES))
  
  # Process chunks and save to disk immediately
  chunk_files <- mclapply(seq_along(chunks), function(i) {
    
    # Run calculation
    chunk_result <- process_chunk(chunks[[i]], i, length(chunks))
    
    # Save to checkpoint directory
    checkpoint_file <- file.path(CHECKPOINT_DIR, sprintf("chunk_%03d_result.rds", i))
    saveRDS(chunk_result, checkpoint_file, compress = TRUE)
    
    # Clean up worker memory
    rm(chunk_result)
    gc()
    
    return(checkpoint_file)
  }, mc.cores = N_CORES)
  
  # Check for errors
  if(any(sapply(chunk_files, function(x) inherits(x, "try-error")))) {
    log_progress("ERROR: One or more parallel workers failed")
    quit(status = 1)
  }
  
} else {
  log_progress(sprintf("Using makeCluster (socket-based) with %d cores", N_CORES))
  
  # Set up cluster
  cl <- makeCluster(N_CORES)
  
  # Export necessary objects
  clusterExport(cl, c("data_clean", "TT", "states_remaining", 
                      "reconstruct_plot_methodB_fast", "OUTPUT_DIR",
                      "CHECKPOINT_DIR", "process_chunk", "chunks"))
  
  # Load packages on workers
  clusterEvalQ(cl, {
    library(data.table)
  })
  
  # Process chunks and save to disk immediately
  chunk_files <- parLapply(cl, seq_along(chunks), function(i) {
    
    # Run calculation
    chunk_result <- process_chunk(chunks[[i]], i, length(chunks))
    
    # Save to checkpoint directory
    checkpoint_file <- file.path(CHECKPOINT_DIR, sprintf("chunk_%03d_result.rds", i))
    saveRDS(chunk_result, checkpoint_file, compress = TRUE)
    
    # Clean up worker memory
    rm(chunk_result)
    gc()
    
    return(checkpoint_file)
  })
  
  stopCluster(cl)
}

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

log_progress(sprintf("Parallel processing complete! Time: %.1f minutes", elapsed))
log_progress(sprintf("All chunks saved to: %s", CHECKPOINT_DIR))

# ==============================================================================
# COMBINE RESULTS SEQUENTIALLY
# ==============================================================================

log_progress("Combining results from disk one-by-one...")
log_progress("(This avoids loading all chunks into RAM simultaneously)")

model_data <- rbindlist(lapply(unlist(chunk_files), function(f) {
  log_progress(sprintf("  Loading: %s", basename(f)))
  chunk_data <- readRDS(f)
  return(chunk_data)
}))

log_progress(sprintf("Combined data: %d rows", nrow(model_data)))

# ==============================================================================
# ADD COVARIATES
# ==============================================================================

log_progress("Adding covariates to reconstructed data...")
model_data <- merge(model_data, 
                    plot_covariates[, .(id_simple, cov_CMI_std, cov_Tmean_std, 
                                        cov_soil, cov_pert_class)],
                    by.x = "id", by.y = "id_simple",
                    all.x = TRUE)

log_progress(sprintf("Final data dimensions: %d rows × %d columns", 
                     nrow(model_data), ncol(model_data)))

# ==============================================================================
# CALCULATE TRANSITION SUMMARY
# ==============================================================================

log_progress("Calculating transition summary statistics...")

transition_summary <- model_data[status == 3, .N, by = .(from, to)]
setorder(transition_summary, -N)
transition_summary[, transition := paste0(from, "→", to)]

log_progress("Top 20 transitions by event count:")
for(i in 1:min(20, nrow(transition_summary))) {
  log_progress(sprintf("  %s: %d events", 
                       transition_summary$transition[i], 
                       transition_summary$N[i]))
}

# Save transition summary
fwrite(transition_summary, file.path(OUTPUT_DIR, "transition_summary.csv"))
log_progress("Saved transition summary to CSV")

# ==============================================================================
# SAVE FINAL DATA
# ==============================================================================

log_progress("Saving reconstructed data as RDS...")

output_file <- file.path(OUTPUT_DIR, "forest_reconstruction_complete.rds")
saveRDS(model_data, output_file, compress = TRUE)

file_size_mb <- file.size(output_file) / 1024^2
log_progress(sprintf("Data saved! File size: %.1f MB", file_size_mb))

# Also save plot covariates separately
saveRDS(plot_covariates, file.path(OUTPUT_DIR, "plot_covariates.rds"))
log_progress("Saved plot covariates separately")

# Save metadata
metadata <- list(
  n_plots = n_plots,
  n_rows = nrow(model_data),
  n_transitions = nrow(transition_summary),
  states = states_remaining,
  max_time = TT,
  processing_time_mins = elapsed,
  timestamp = Sys.time()
)
saveRDS(metadata, file.path(OUTPUT_DIR, "metadata.rds"))
log_progress("Saved metadata")

# ==============================================================================
# CLEANUP CHECKPOINTS (OPTIONAL)
# ==============================================================================

# Optionally remove checkpoint files after successful completion
KEEP_CHECKPOINTS <- FALSE  # Set to TRUE to keep checkpoints

if(!KEEP_CHECKPOINTS) {
  log_progress("Removing checkpoint files...")
  checkpoint_files <- list.files(CHECKPOINT_DIR, 
                                 pattern = "^chunk_.*_result\\.rds$", 
                                 full.names = TRUE)
  file.remove(checkpoint_files)
  log_progress(sprintf("Removed %d checkpoint files", length(checkpoint_files)))
} else {
  log_progress(sprintf("Checkpoint files preserved in: %s", CHECKPOINT_DIR))
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

log_progress("=== RECONSTRUCTION COMPLETE ===")
log_progress(sprintf("Total plots processed: %d", n_plots))
log_progress(sprintf("Total rows generated: %d", nrow(model_data)))
log_progress(sprintf("Total processing time: %.1f minutes (%.1f hours)", 
                     total_time, total_time / 60))
log_progress(sprintf("Average time per 1000 plots: %.2f minutes", 
                     total_time / (n_plots / 1000)))
log_progress(sprintf("Output file: %s", output_file))
log_progress("===================================")

# Print to console for job completion
cat("\n\n")
cat("##############################################\n")
cat("# RECONSTRUCTION COMPLETE\n")
cat(sprintf("# Processed: %d plots\n", n_plots))
cat(sprintf("# Generated: %d rows\n", nrow(model_data)))
cat(sprintf("# Time: %.1f hours\n", total_time / 60))
cat(sprintf("# Output: %s\n", output_file))
cat("##############################################\n")
