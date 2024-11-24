library(stats)

dicretize <- function(X, y, threshold = 0.01, continuous_cols, quantiles = seq(0, 1, length.out = 11)) {
  # X: data frame
  # y: target variable
  # threshold: % threshold for NLL improvement
  # continuous_cols: list of continuous columns in X
  # quantiles: possible cuts for discretization
  
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
    bins <- c(0, 1)
    
    while (TRUE) {
      # find the next best cut point for col
      X_split <- next_best_split(X, y, col, bins, quantiles)
      
      model_split <- glm(y ~ ., data = cbind(X_split, y = y), family = binomial)
      NLL_split <- -logLik(model_split)
      
      if (NLL_split < ((1 - threshold) * NLL)) {
        X <- X_split
        NLL <- NLL_split
        bins <- unique(as.numeric(levels(X[[col]])))
      } else {
        break # Model is no longer improving with new cuts
      }
    }
  }
  return(X)
}

next_best_split <- function(X, y, col, bins, quantiles) {
  # X: data frame
  # y: target variable
  # col: col to make cuts on next
  # bins: current bins we're allowed to make a cut on
  # quantiles: possible split boundaries
  
  best_obj_value <- Inf
  prev_beta <- NULL # for warm-starting the model optimization
  
  for (q in quantiles) {
    # Make a cut point at q by splitting the bin
    new_bins <- sort(unique(c(bins, q)))
    
    X_temp <- X
    X_temp[[col]] <- cut(X[[col]], breaks = quantile(X[[col]], probs = new_bins, na.rm = TRUE),
                         include.lowest = TRUE, labels = FALSE)
    
    # Fit risk model and calculate objective function
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
    }
  }
  return(X)
}