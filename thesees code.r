# =============================================================================
# Setup Script for GKP-PMM and Two-Phase Sampling Analysis
# =============================================================================
#
# This script installs and loads all required packages for the Gower-kPrototypes
# Predictive Mean Matching (GKP-PMM) imputation method and the two-phase sampling
# methodology with Gower's pre-selection and LASSO.
#
# Author: [Your Name]
# Date: September 6, 2025
# =============================================================================

# Install required packages (uncomment if needed)
# install.packages(c("mice", "cluster", "clustMixType", "MASS", "nnet",
#                    "glmnet", "survey", "dplyr", "ggplot2", "tidyr",
#                    "microbenchmark"))

# Load libraries
library(mice)          # For multiple imputation
library(cluster)       # For Gower's distance
library(clustMixType)  # For k-prototypes clustering
library(MASS)          # For ordered logistic regression
library(nnet)          # For multinomial regression
library(glmnet)        # For LASSO regression
library(survey)        # For survey analysis
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(tidyr)         # For data reshaping
library(microbenchmark)# For performance benchmarking

cat("All required packages loaded successfully.\n")
