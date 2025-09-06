# --- FULLY ROBUST SIMULATION SCRIPT (PARALLELIZED FOR CODESPACES) ---

# --- 1. SETUP ---
# Libraries are installed by the devcontainer, but we still need to load them
library(mice)
library(dplyr)
library(ggplot2)
library(tidyr)
library(future)
library(furrr)

# Load your custom imputation function
# source("mice.impute.gkp_pmm.R") # Good practice to save your function in a file

# --- 2. DEFINE SIMULATION FUNCTION ---
# It's best practice to wrap one full simulation run in a single function
run_one_simulation <- function(sim_id, size, mech, rate, patt_name, patt_matrix, methods) {
  
  # Bootstrap data
  data_complete_sized <- data_complete_base[sample(1:nrow(data_complete_base), size, replace = TRUE), ]
  
  # Ampute data
  data_amputed_obj <- ampute(data_complete_sized, prop = rate, patterns = patt_matrix, mech = mech)
  data_missing <- data_amputed_obj$amp
  missing_pattern <- is.na(data_missing) & !is.na(data_complete_sized)
  
  # This inner loop over methods is fast, so we don't parallelize it.
  results_for_this_run <- list()
  for (method in methods) {
    method_name <- method; impute_method <- "gkp_pmm"; k_pre_clusters_val <- 3
    if (method == "gkp_pmm_no_clust") k_pre_clusters_val <- 0
    else if (method != "gkp_pmm") impute_method <- method

    start_time <- Sys.time()
    imputed_obj <- mice(data_missing, m = 1, maxit = 5, method = impute_method, 
                        k_pre_clusters = k_pre_clusters_val, printFlag = FALSE)
    end_time <- Sys.time()
    data_imputed <- complete(imputed_obj, 1)

    for (col_name in colnames(data_complete_sized)) {
      if (any(missing_pattern[, col_name])) {
        true_vals <- data_complete_sized[missing_pattern[, col_name], col_name]
        imputed_vals <- data_imputed[missing_pattern[, col_name], col_name]
        metric_name <- if(is.numeric(true_vals)) "NRMSE" else "Misclassification"
        error_val <- if(is.numeric(true_vals)) {
          sqrt(mean((true_vals - imputed_vals)^2, na.rm=TRUE)) / sd(data_complete_sized[, col_name], na.rm=TRUE)
        } else { mean(true_vals != imputed_vals, na.rm=TRUE) }
        
        results_for_this_run[[length(results_for_this_run) + 1]] <- data.frame(
          sample_size = size, mechanism = mech, rate = rate, pattern = patt_name,
          sim_run = sim_id, method = method_name, variable = col_name,
          metric = metric_name, error = error_val,
          time_sec = as.numeric(difftime(end_time, start_time, units = "secs"))
        )
      }
    }
  }
  return(bind_rows(results_for_this_run))
}

# --- 3. EXECUTE THE PARALLEL SIMULATION ---

# Define all parameters (same as before)
# ... [SAMPLE_SIZES, MECHANISMS, MISSING_RATES, patterns_list, N_SIM, METHODS_TO_COMPARE] ...
data_complete_base <- na.omit(nhanes[, c("age", "bmi", "hyp", "chl")])
SAMPLE_SIZES <- c(100, 500, 1000, 10000)
MECHANISMS <- c("MCAR", "MAR") # Start smaller
MISSING_RATES <- c(0.20, 0.40)
patterns_list <- list("All_Vars_Missing" = NULL)
N_SIM <- 50 # We can increase N_SIM now since it runs in parallel
METHODS_TO_COMPARE <- c("gkp_pmm", "gkp_pmm_no_clust", "pmm", "rf")


# Create a grid of all simulation conditions
conditions <- expand.grid(
  sim_id = 1:N_SIM,
  size = SAMPLE_SIZES,
  mech = MECHANISMS,
  rate = MISSING_RATES,
  patt_name = names(patterns_list),
  stringsAsFactors = FALSE
)

# Set up the parallel backend. This will use all available cores.
plan(multisession) 

cat("Starting parallel simulation...\n")
# Use future_map_dfr to run the simulation in parallel and combine results
all_results_robust <- future_map_dfr(
  1:nrow(conditions),
  ~ run_one_simulation(
      sim_id = conditions$sim_id[.x],
      size = conditions$size[.x],
      mech = conditions$mech[.x],
      rate = conditions$rate[.x],
      patt_name = conditions$patt_name[.x],
      patt_matrix = patterns_list[[conditions$patt_name[.x]]],
      methods = METHODS_TO_COMPARE
    ),
  .progress = TRUE # Shows a nice progress bar!
)
cat("Parallel simulation complete!\n")

# Save the results
saveRDS(all_results_robust, "simulation_results.rds")


# --- 4. ANALYZE AND VISUALIZE ENHANCED RESULTS ---

# --- 4a. Computational Cost ---
time_summary <- all_results %>%
  distinct(rate, pattern, method, sim_run, time_sec) %>%
  group_by(rate, pattern, method) %>%
  summarise(mean_time = mean(time_sec), .groups = 'drop')

cat("\n--- Computational Cost Summary ---\n")
print(time_summary)

ggplot(time_summary, aes(x = method, y = mean_time, fill = method)) +
  geom_bar(stat = "identity") +
  facet_grid(rate ~ pattern) +
  labs(
    title = "Computational Cost Across Conditions",
    x = "Imputation Method", y = "Mean Time (seconds)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- 4b. Imputation Quality ---
accuracy_summary <- all_results %>%
  group_by(rate, pattern, method, variable, metric) %>%
  summarise(mean_error = mean(error), .groups = 'drop') %>%
  arrange(rate, pattern, variable, mean_error)

cat("\n--- Imputation Quality Summary ---\n")
print(accuracy_summary, n=100) # Print more rows

ggplot(all_results, aes(x = reorder(method, error), y = error, fill = method)) +
  geom_boxplot() +
  # Create a grid of plots: rows are variables, columns are patterns and rates
  facet_grid(variable ~ pattern + rate, scales = "free_y") +
  labs(
    title = "Imputation Quality Across Varying Rates and Patterns",
    subtitle = "Lower error is better. Rows = Variables, Columns = Conditions.",
    x = "Imputation Method",
    y = "Error (NRMSE or Misclassification)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") # Hide legend to save space