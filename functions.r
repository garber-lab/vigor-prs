library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(broom)
library(here)


# Function to recode specified SNPs so that risk alleles match weights
recode_snps <- function(df, snps_to_recode) {
  df %>%
    mutate(across(all_of(snps_to_recode), 
                  ~ ifelse(. == 0, 2,
                           ifelse(. == 2, 0,
                                  ifelse(. == 1, 1, NA)))))
}

# Function to perform matrix multiplication to generate the main PRS
calculate_main_prs <- function(risk_genos = risk_genos_recoded,
                               weights = weights,
                               iid_col = "IID",
                               weights_col = "Weights_risk_direction_only") {
  weight_mat <- weights %>%
    column_to_rownames(var = "SNP") %>%
    select(weights = !!sym(weights_col)) %>%
    as.matrix()
  
  geno_mat <- risk_genos %>%
    select(rownames(weight_mat)) %>%
    as.matrix()
  
  rownames(geno_mat) <- risk_genos %>% pull(!!sym(iid_col))
  
  res <- geno_mat %*% weight_mat
  res <- as.data.frame(res) %>%
    rownames_to_column(., var="tmp_id") %>%
    transmute(!!sym(iid_col) := tmp_id,
              score := weights)
  return(res)
}

# Function to calculate Class II risk score
calculate_class2_risk_score <- function(risk_genos) {
  high_risk_hap_ln_or <- 2.0541237
  generic_risk_ln_or <- 0.5721089
  
  risk_genos %>%
    mutate(mhc_class2_only_prs = case_when(
      RS145954018 >= 1 ~ high_risk_hap_ln_or,
      RS114448410 == 2 ~ 2 * generic_risk_ln_or,
      RS114448410 == 1 ~ generic_risk_ln_or,
      RS114448410 == 0 & RS145954018 == 0 ~ 0,
      TRUE ~ NA_real_
    )) %>%
    select(IID, mhc_class2_only_prs)
}

# write a function that scales the PRS based on controls
library(dplyr)

scale_prs <- function(data = analysis_df,
                      pheno_col = "vitiligo",
                      pheno_cond = 0,
                      scale_columns = c("mhc_class2_only_prs", "CONFIRMED_prs")) {
  for (i in scale_columns) {
    mean_val <- mean(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    std_dev_val <- sd(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    
    new_col_name <- paste0(i, "_scaled")
    percentile_col_name <- paste0(i, "_percentile")
    
    # Calculate percentile
    control_prs <- na.omit(data[data[[pheno_col]] == pheno_cond, i])
    ecdf_func <- ecdf(control_prs)
    percentiles <- ifelse(data$case_control == "control", ecdf_func(data[[i]]), NA)
    
    data <- data %>% 
      mutate(!!sym(new_col_name) := (data[[i]] - mean_val) / std_dev_val,
             !!sym(percentile_col_name) := ecdf_func(!!sym(i)))
  }
  return(data)
}
