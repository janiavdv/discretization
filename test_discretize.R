source('discretize.R')
library(dplyr)

raw_data <- read.csv('data/Peru_TB_data.csv')

floor_median <- function(x) {
  return(floor(median(x, na.rm=TRUE)))
}

preprocess <- function(data) {
  # define binary outcome (y)
  data$adherence_outcome <- ifelse(data$PCTadherence_sensi < 0.95, 1, 0)
  
  # drop subsequent regimens
  data <- data %>%
    filter(regimen == "first") %>%
    select(-regimen)
  
  # drop X, patient id, and other variables we don't want to include in our model
  data <- data %>% select(-c(X, # index
                             prop_misseddoses_int, # not in dictionary
                             PCTadherence,
                             PCTadherence_sensi,
                             PTID2, # study ID
                             hunger_freq, # "not including in the model"
                             health_ctr, # 71 health centers ,
                             daily_cont, # not in dictionary
                             post_tb_pt_work)) #  not in dictionary
  
  # family support - simplify by taking median response to get almost integer value
  fam_vars <- c("fam_affection", "fam_parent_getalong", "fam_others_getalong",
                "fam_confide", "fam_parent_emosupport", "fam_others_emosupport",
                "fam_support_others", "fam_satisfied")
  data$family_median <- apply(data[, fam_vars], 1, floor_median)
  data <- data %>% select(-all_of(c(fam_vars, "fam_support")))
  
  # evaluation of health services
  health_serv_vars <- c("aten_wait", "aten_respect", "aten_explain", "aten_space",
                        "aten_concern","aten_satis_hours")
  data$health_svc_median <- apply(data[, health_serv_vars],
                                   1, floor_median)
  data <- data %>% select(-all_of(c(health_serv_vars, "healthsvc_satis")))
  
  # motivation
  motivation_vars <- c("motiv_no_transmit", "motiv_fam_worry", "motiv_study_work",
                       "motiv_activities")
  data$motiv_median <- apply(data[, motivation_vars], 1, floor_median)
  data <- data %>% select(-all_of(c(motivation_vars,  "motiv_summary")))
  
  # tb disinformation
  knowledge_vars <- c("conoc_cure", "conoc_missed_doses", "conoc_default")
  data$knowledge_median <- apply(data[, knowledge_vars], 1, floor_median)
  data <- data %>% select(-all_of(c(knowledge_vars,  "tb_knowledge")))
  
  # pills and adr_freq need to multiplied by 4 - does not match data dictionary
  data$pills <- 4*data$pills
  data$adr_freq <- 4*data$adr_freq
  
  # change to categorical
  data$current_sx_none <- as.factor(data$current_sx_none)
  data$tto_anterior_tb <- as.factor(data$tto_anterior_tb)
  data$monitor1 <- as.factor(data$monitor1)
  
  data <- na.omit(data)
  
  # some categorical variable have levels with few observations
  # drop ram and regular_drug since only 5 observations in class 1
  data <- data %>% select(-c(ram, regular_drug))

  # education levels could be combined but we drop since not
  # included in data documentation and unclear best way to relevel
  data <- data %>% select(-c(edu_level_mom, edu_level_dad))
  
  non_numeric_columns <- names(data)[!sapply(data, is.numeric)]
  data[non_numeric_columns] <- lapply(data[non_numeric_columns], as.factor)
  
  return(data)
}

df <- preprocess(raw_data)
X <- df %>% select(-adherence_outcome) %>% select(-audit_tot) # TODO: figure out what's wrong with this col
y <- df %>% pull(adherence_outcome) # get target column as vector

cont_cols <- c("self_eff", "tx_mos", "stig_tot",
              "phq9_tot", "age_BLchart", "ace_score")

disc_df <- discretize(X, y, continuous_cols = cont_cols)
