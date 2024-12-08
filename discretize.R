library(stats)
source("../riskscores/R/risk_mod.R")
source("../riskscores/R/utils.R")

discretize <- function(X, y, threshold = 0.01, continuous_cols, n_quantiles = NULL) {
  # X: data frame
  # y: target variable
  # threshold: % threshold for NLL improvement
  # continuous_cols: list of continuous columns in X
  # n_quantiles: list of positive integers of quantiles for each column
  
  n_cols = length(continuous_cols)
  
  # If quantiles is NULL, use 10 quantiles for every column
  if (is.null(n_quantiles)) {
    n_quantiles <- rep(10, n_cols)
  }
  
  if (length(n_quantiles) != n_cols) {
    stop ("continuous_cols and n_quantiles must be the same length!")
  }
  
  if (any(n_quantiles <= 0) || any(n_quantiles %% 1 != 0)) {
    stop("n_quantiles must contain only positive integers")
  }

  quantiles <- list()
  for (i in 1:n_cols) {
    # Generate quantiles for each column
    quantiles[[i]] <- seq(0, 1, length.out = n_quantiles[i] + 1)
    
    # Initialize the discretization with the quantiles
    col <- continuous_cols[i]
    breaks <- unique(quantile(X[[col]], probs = quantiles[[i]], na.rm = TRUE))
    X[[col]] <- cut(X[[col]], breaks = breaks, 
                    include.lowest = TRUE, labels = FALSE)
  }
  
  # fit logistic regression and calculate NLL
  model <- glm(y ~ ., data = cbind(X, y = y), family = binomial)
  NLL <- -logLik(model)
  
  # Shuffle columns and quantiles based on the same indices
  shuffle <- sample(n_cols)
  continuous_cols <- continuous_cols[shuffle]
  quantiles <- quantiles[shuffle]
  
  for (i in 1:n_cols) {
    bins <- c(0, 1)
    col <- continuous_cols[i]
    
    while (TRUE) {
      # find the next best cut point for col
      X_split <- next_best_split(X, y, col, bins, quantiles[[i]], threshold)
      
      model_split <- glm(y ~ ., data = cbind(X_split, y = y), family = binomial)
      NLL_split <- -logLik(model_split)
      
      if (NLL_split >= ((1 + threshold) * NLL)) {
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

next_best_split <- function(X, y, col, bins, col_quantiles, threshold) {
  # X: data frame
  # y: target variable
  # col: col to make cuts on next
  # bins: current bins we're allowed to make a cut on
  # col_quantiles: possible split boundaries for this column
  # threshold: % threshold for objective function improvement
  
  best_obj_value <- Inf
  prev_beta <- NULL # for warm-starting the model optimization
  
  for (q in col_quantiles) {
    # Make a cut point at q by splitting the bin
    new_bins <- sort(unique(c(bins, q)))
    
    X_temp <- X
    breaks <- unique(quantile(X[[col]], breaks = new_bins, na.rm = TRUE))
    if (length(breaks) < 2) next
    X_temp[[col]] <- cut(X[[col]], breaks = breaks,
                         include.lowest = TRUE, labels = FALSE)
    rm_df <- X_temp
    rm_df$y <- y
    
    # model matrix and save
    rm_mat <- model.matrix(y ~ ., data = rm_df)
    rm_mat <- cbind(rm_mat, y = rm_df$y)
    X_rm <- as.matrix(rm_mat[,-ncol(rm_mat)])
    y_rm <- as.matrix(rm_mat[,ncol(rm_mat)])
    
    # Fit risk model and calculate objective function
    if (is.null(prev_beta)) {
      mod <- risk_mod(X = X_rm, y = y_rm)
    } else {
      mod <- risk_mod(X = X_rm, y = y_rm, beta = prev_beta)
    }
    
    obj_value <- obj_fcn(mod$X, mod$y, mod$gamma, mod$beta, mod$weights, mod$lambda0)
    # Check if new cut improved obj function by threshold %
    if (obj_value >= ((1 + threshold) * best_obj_value)) {
      best_obj_value <- obj_value
      prev_beta <- mod$beta
      X <- X_temp
    }
  }
  return(X)
}