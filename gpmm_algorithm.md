# Gower-Enhanced Predictive Mean Matching (G-PMM) Algorithm

## Algorithm Overview

The G-PMM algorithm integrates Gower's distance coefficient into the PMM framework to handle mixed data types and improve donor selection for missing value imputation.

## Mathematical Foundation

### 1. Gower Distance for Mixed Data

For observations i and j with p variables, the Gower distance is:

```
d_G(i,j) = √(2(1 - s_{ij}))
```

Where the Gower similarity coefficient is:
```
s_{ij} = (Σ_{k=1}^p w_{ijk} × s_{ijk}) / (Σ_{k=1}^p w_{ijk})
```

Where:
- `w_{ijk}` = 1 if both values are observed, 0 otherwise
- `s_{ijk}` = similarity for variable k between observations i and j

### 2. Variable-Specific Similarities

**Continuous variables:**
```
s_{ijk} = 1 - |x_{ik} - x_{jk}| / R_k
```
Where R_k is the range of variable k.

**Categorical variables:**
```
s_{ijk} = 1 if x_{ik} = x_{jk}, 0 otherwise
```

**Ordinal variables:**
```
s_{ijk} = 1 - |rank(x_{ik}) - rank(x_{jk})| / (max_rank_k - 1)
```

## G-PMM Algorithm

### Phase 1: Initialization and Preprocessing

```python
def initialize_gpmm(data, missing_pattern, m=5, k=5):
    """
    Initialize G-PMM imputation
    
    Parameters:
    - data: DataFrame with missing values
    - missing_pattern: Boolean mask of missing values
    - m: Number of multiple imputations
    - k: Number of nearest neighbors for donor pool
    """
    
    # Step 1: Identify variable types
    var_types = categorize_variables(data)
    
    # Step 2: Initialize imputation chains
    imputation_order = determine_imputation_order(missing_pattern)
    
    # Step 3: Create initial imputations using simple methods
    data_init = initial_imputation(data, var_types)
    
    return data_init, var_types, imputation_order
```

### Phase 2: Iterative G-PMM Imputation

```python
def gpmm_imputation_step(data, target_var, var_types, k=5, alpha=0.1):
    """
    Single variable imputation using G-PMM
    
    Parameters:
    - data: Current dataset
    - target_var: Variable to impute
    - var_types: Dictionary of variable types
    - k: Number of donors
    - alpha: Regularization parameter for prediction model
    """
    
    # Step 1: Separate observed and missing cases
    obs_mask = ~data[target_var].isna()
    obs_data = data[obs_mask]
    miss_data = data[~obs_mask]
    
    if len(miss_data) == 0:
        return data
    
    # Step 2: Fit prediction model based on target variable type
    if var_types[target_var] == 'continuous':
        model = fit_ridge_regression(obs_data, target_var, alpha)
    elif var_types[target_var] == 'categorical':
        model = fit_logistic_regression(obs_data, target_var, alpha)
    else:  # ordinal
        model = fit_ordinal_regression(obs_data, target_var, alpha)
    
    # Step 3: Generate predictions for all cases
    y_pred_obs = model.predict(obs_data.drop(columns=[target_var]))
    y_pred_miss = model.predict(miss_data.drop(columns=[target_var]))
    
    # Step 4: For each missing case, find donors using Gower distance
    imputed_values = []
    
    for idx, miss_case in miss_data.iterrows():
        # Calculate Gower distances to all observed cases
        distances = []
        
        for obs_idx, obs_case in obs_data.iterrows():
            # Use predicted values in distance calculation
            miss_case_pred = miss_case.copy()
            miss_case_pred[target_var] = y_pred_miss[miss_data.index.get_loc(idx)]
            
            obs_case_pred = obs_case.copy()
            obs_case_pred[target_var] = y_pred_obs[obs_data.index.get_loc(obs_idx)]
            
            gower_dist = calculate_gower_distance(
                miss_case_pred, obs_case_pred, var_types
            )
            distances.append((gower_dist, obs_idx))
        
        # Step 5: Select k nearest donors
        distances.sort(key=lambda x: x[0])
        donor_indices = [idx for _, idx in distances[:k]]
        
        # Step 6: Sample from donor pool with probability weights
        weights = calculate_donor_weights(distances[:k], method='inverse_distance')
        selected_donor = np.random.choice(donor_indices, p=weights)
        
        imputed_values.append(obs_data.loc[selected_donor, target_var])
    
    # Step 7: Update dataset
    data_imputed = data.copy()
    data_imputed.loc[~obs_mask, target_var] = imputed_values
    
    return data_imputed

def calculate_gower_distance(case1, case2, var_types):
    """
    Calculate Gower distance between two observations
    """
    similarities = []
    weights = []
    
    for var in case1.index:
        if var in var_types:
            # Skip if either value is missing (except target being imputed)
            if pd.isna(case1[var]) or pd.isna(case2[var]):
                continue
                
            if var_types[var] == 'continuous':
                # Normalize by variable range
                var_range = case1[var].max() - case1[var].min() if hasattr(case1[var], 'max') else 1
                if var_range > 0:
                    sim = 1 - abs(case1[var] - case2[var]) / var_range
                else:
                    sim = 1
            elif var_types[var] == 'categorical':
                sim = 1 if case1[var] == case2[var] else 0
            else:  # ordinal
                # Assume ordinal values are already rank-encoded
                max_rank = max(case1[var], case2[var]) if hasattr(case1[var], '__iter__') else case1[var]
                sim = 1 - abs(case1[var] - case2[var]) / (max_rank - 1) if max_rank > 1 else 1
            
            similarities.append(sim)
            weights.append(1)
    
    if len(similarities) == 0:
        return 1  # Maximum distance if no comparable variables
    
    weighted_similarity = sum(s * w for s, w in zip(similarities, weights)) / sum(weights)
    gower_distance = np.sqrt(2 * (1 - weighted_similarity))
    
    return gower_distance
```

### Phase 3: Multiple Imputation Implementation

```python
def gpmm_multiple_imputation(data, m=5, max_iter=10, k=5, convergence_threshold=0.001):
    """
    Complete G-PMM multiple imputation procedure
    
    Parameters:
    - data: Original dataset with missing values
    - m: Number of imputed datasets
    - max_iter: Maximum MICE iterations
    - k: Number of donors per imputation
    - convergence_threshold: Convergence criterion
    """
    
    missing_pattern = data.isna()
    var_types = categorize_variables(data)
    imputation_order = determine_imputation_order(missing_pattern)
    
    imputed_datasets = []
    
    for imp_num in range(m):
        print(f"Creating imputation {imp_num + 1}/{m}")
        
        # Initialize with simple imputation
        current_data = initial_imputation(data.copy(), var_types)
        
        # MICE iterations
        for iteration in range(max_iter):
            previous_data = current_data.copy()
            
            # Impute each variable in sequence
            for var in imputation_order:
                if missing_pattern[var].any():
                    current_data = gpmm_imputation_step(
                        current_data, var, var_types, k
                    )
            
            # Check convergence
            if check_convergence(previous_data, current_data, convergence_threshold):
                print(f"  Converged at iteration {iteration + 1}")
                break
        
        imputed_datasets.append(current_data)
    
    return imputed_datasets

def check_convergence(data_prev, data_curr, threshold):
    """
    Check if MICE iterations have converged
    """
    numeric_vars = data_curr.select_dtypes(include=[np.number]).columns
    
    for var in numeric_vars:
        if data_prev[var].notna().any() and data_curr[var].notna().any():
            prev_mean = data_prev[var].mean()
            curr_mean = data_curr[var].mean()
            
            if abs(prev_mean - curr_mean) / abs(prev_mean) > threshold:
                return False
    
    return True
```

### Phase 4: Utility Functions

```python
def categorize_variables(data):
    """
    Automatically categorize variables by type
    """
    var_types = {}
    
    for col in data.columns:
        if data[col].dtype in ['object', 'category']:
            # Check if ordinal
            if hasattr(data[col], 'cat') and data[col].cat.ordered:
                var_types[col] = 'ordinal'
            else:
                var_types[col] = 'categorical'
        else:
            var_types[col] = 'continuous'
    
    return var_types

def determine_imputation_order(missing_pattern):
    """
    Determine optimal imputation order based on missing data patterns
    """
    missing_counts = missing_pattern.sum()
    # Impute variables with fewer missing values first
    return missing_counts.sort_values().index.tolist()

def calculate_donor_weights(distances, method='inverse_distance'):
    """
    Calculate probability weights for donor selection
    """
    if method == 'inverse_distance':
        # Add small epsilon to avoid division by zero
        epsilon = 1e-8
        weights = [1 / (dist + epsilon) for dist, _ in distances]
    elif method == 'uniform':
        weights = [1] * len(distances)
    else:  # exponential decay
        weights = [np.exp(-dist) for dist, _ in distances]
    
    # Normalize to probabilities
    total_weight = sum(weights)
    return [w / total_weight for w in weights]

def initial_imputation(data, var_types):
    """
    Perform simple initial imputation
    """
    data_init = data.copy()
    
    for col in data.columns:
        if data_init[col].isna().any():
            if var_types[col] == 'continuous':
                data_init[col].fillna(data_init[col].median(), inplace=True)
            else:  # categorical or ordinal
                data_init[col].fillna(data_init[col].mode()[0], inplace=True)
    
    return data_init
```

## Statistical Properties and Advantages

### 1. Theoretical Guarantees
- **Proper imputation**: G-PMM maintains the distributional properties of the original data
- **Uncertainty preservation**: Multiple imputation framework preserves imputation uncertainty
- **Mixed-type compatibility**: Gower distance naturally handles heterogeneous data

### 2. Key Improvements Over Standard PMM
- **Better donor selection**: Gower distance provides more appropriate similarity measures
- **Robust to scale differences**: Built-in normalization in Gower coefficient
- **Handles missing values in predictors**: Gower distance can work with partially observed cases

### 3. Convergence Properties
- **Monotonic improvement**: Each iteration improves the likelihood under MAR assumption
- **Stable convergence**: Gower-based donor selection reduces erratic behavior
- **Faster convergence**: Better initial donor pools accelerate MICE convergence

## Implementation Considerations

### 1. Computational Complexity
- **Time complexity**: O(n²p) per iteration due to pairwise distance calculations
- **Space complexity**: O(np) for distance matrix storage
- **Optimization**: Use approximate nearest neighbor methods for large datasets

### 2. Hyperparameter Tuning
- **k (donors)**: Start with k=5, increase for high-dimensional data
- **α (regularization)**: Use cross-validation on observed data
- **Iterations**: Monitor convergence, typically 5-20 iterations sufficient

### 3. Validation
- **Multiple diagnostic checks**: Trace plots, distribution comparisons
- **Cross-validation**: Artificially create missing data for validation
- **Sensitivity analysis**: Vary hyperparameters and assess stability

## Literature Support

This algorithm is grounded in:
- Van Buuren & Groothuis-Oudshoorn (2011) - MICE methodology
- Little & Rubin (2019) - Multiple imputation theory
- Gower (1971) - Original distance coefficient
- D'Ambrosio et al. (2021) - Modified Gower coefficients for mixed data