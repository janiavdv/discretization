library(stats)

dicretize <- function(X, y, threshold = 0.01, continuous_cols) {
  # X: data frame
  # y: target variable
  # threshold: % threshold for NLL improvement
  # continuous_cols: list of continuous columns in X

  quantiles <- seq(0, 1, length.out = 11) # 10 quantiles
  
  # initialize the discretization with the quantiles
  for (col in continuous_cols) {
      X[[col]] <- cut(X[[col]], quantile(X[[col]], probs = quantiles), 
                      include.lowest = TRUE, labels = FALSE)
  }
  
  # fit logistic regression and calculate NLL
  model <- glm(y ~ ., data = cbind(X, y = y), family = binomial)
  NLL <- -logLik(model)
  
  # shuffle columns order
  continuous_cols <- sample(continuous_cols)
  
  for (col in continuous_cols) {
    while (TRUE) {
      # find the next best cut point for col
      X_cut = next_best_cut(X, y, col, sort(unique(X[[col]])))
      
      model_cut <- glm(y ~ ., data = cbind(X_cut, y = y), family = binomial)
      NLL_cut <- -logLik(model)
      
      if (NLL_cut < ((1 - threshold) * NLL)) {
        X <- X_cut
        NLL <- NLL_cut
      } else {
        break # Model is no longer improving with new cuts
      }
    }
  }
  return(X)
}

next_best_cut <- function(X, y, col, bins) {
  # X: data frame
  # y: target variable
  # col: col to make cuts on next
  # bins: current bins we're allowed to make a cut on
  
  n_bins <- length(bins)
  if (n_bins <= 1) {
    return(X)  # No cuts possible (recursive base case)
  }
  
  best_obj_value <- Inf
  prev_beta <- NULL # for warm-starting the model optimization
  
  for (i in 1:(n_bins-1)) {
    # Make a cut point by merging bins
    merged_bins <- bins
    merged_bins[i+1] <- merged_bins[i] # Merge bins i and i+1
    X_temp <- X
    X_temp[[col]] <- factor(X_temp[[col]], levels = bins) # create factor
    X_temp[[col]] <- factor(X_temp[[col]], levels = merged_bins) # change levels of factor 
    
    # Fit risk model on the new data
    if (is.null(prev_beta)) {
      mod <- riskscores::risk_mod(X = X_temp, y = y)
    } else {
      mod <- riskscores::risk_mod(X = X_temp, y = y, beta = prev_beta)
    }
    
    obj_value <- obj_fn(mod$X, mod$y, mod$gamma, mod$beta, mod$weights, mod$lambda0)
    
    # Check if new cut improved obj function
    if (obj_value < best_obj_value) {
      
      best_obj_value <- obj_value
      prev_beta <- mod$beta
      X <- X_temp
      
      # Recursively try cuts on either side of the new cut point
      # ISSUE: doesn't it introduce bias to do post-order traversal?
      L <- X
      L_bins <- bins[bins <= merged_bins[i]]
      L <- next_best_cut(L, y, col, L_bins)
      
      R <- X
      R_bins <- bins[bins > merged_bins[i]]
      R <- next_best_cut(R, y, col, R_bins)
      
      # Combine results from left and right into X
      # if it belongs to L, use L[[col]], else R[[col]]
      X[[col]] <- factor(ifelse(X[[col]] %in% L_bins, L[[col]], R[[col]]), levels = bins)
    }
  }
  return(X)
}