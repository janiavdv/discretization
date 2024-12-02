library(stats)
source("../riskscores/R/risk_mod.R")
source("../riskscores/R/utils.R")

discretize <- function(X, y, threshold = 0.01, continuous_cols, quantiles = seq(0, 1, length.out = 11)) {
  # X: data frame
  # y: target variable
  # threshold: % threshold for NLL improvement
  # continuous_cols: list of continuous columns in X
  # quantiles: possible cuts for discretization
  
  # initialize the discretization with the quantiles
  for (col in continuous_cols) {
    breaks <- unique(quantile(X[[col]], probs = quantiles, na.rm = TRUE))
    X[[col]] <- cut(X[[col]], breaks = breaks, 
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
      X_split <- next_best_split(X, y, col, bins, quantiles, threshold)
      
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

next_best_split <- function(X, y, col, bins, quantiles, threshold) {
  # X: data frame
  # y: target variable
  # col: col to make cuts on next
  # bins: current bins we're allowed to make a cut on
  # quantiles: possible split boundaries
  # threshold: % threshold for objective function improvement
  
  best_obj_value <- Inf
  prev_beta <- NULL # for warm-starting the model optimization
  
  for (q in quantiles) {
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