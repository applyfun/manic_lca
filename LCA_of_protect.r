### Latent class analysis of responses to manic or irritable symptoms
### 05/10/2020
### Ryan Arathimos

### load packages
library(data.table)
library(ukbkings)
library(tidyverse)
library(ggplot2)
library(stringi)
library(poLCA)
library(parallel)
library(foreach)
library(doParallel)
library(janitor)
library(broom)
library(gridExtra)
library(scales)

set.seed(123456)

### SET UP DIRECTORIES

dirs <- read.table("~/brc_scratch/scripts/paths.txt")

dir_prs <- as.character(dirs[1, ])
new_2019_qc_dir <- as.character(dirs[3, ])
qc_dir <- as.character(dirs[3, ])
project_dir <- as.character(dirs[9, ])

data_dir <- as.character(paste0(dirs[10, ], "data/manic_lca"))
output_dir <- as.character(paste0(dirs[10, ], "output/manic_lca"))
scripts_dir <- as.character(paste0(dirs[10, ], "scripts/manic_lca"))

### pick a color palette for plots

palette_choice <- "Set1" # options from rcolorbrewer palette eg. "Dark2"

### to run LCA - 2 classes represents separation between manic/not-manic individuals in full dataset, >2 needed

n_series <- 7 # maximum number of classes to attempt

df <-
  readRDS(paste0(output_dir, "/", "manic_data_derived_protect_mhq.rds"))

dim(df)

### subset responders of manic symptoms question - drop non-responders
main_manic_df <-
  df[which(
    df$Ever.extreme.irritability == 1 |
      df$Ever.manic.or.excitable == 1 |
      df$Ever.extreme.irritability == 0 |
      df$Ever.manic.or.excitable == 0
  ), ]

dim(main_manic_df)

manic_derived_vars <-
  c(
    "More.active.manic.episode.derived",
    "More.confident.manic.episode.derived",
    "Easily.distracted.manic.episode.derived",
    "More.creative.manic.episode.derived",
    "Less.sleep.manic.episode.derived",
    "Thoughts.racing.manic.episode.derived",
    "More.restless.manic.episode.derived",
    "More.talkative.manic.episode.derived"
  )

covars <- c("Ever.extreme.irritability", "Ever.manic.or.excitable")

manic_derived_vars_and_covariates <- c(manic_derived_vars, covars)

### create an ID list for later merge

main_manic_df_ids <- main_manic_df$TecID

#############################################################
################# DEFINE ANALYTICAL SAMPLES #################

### define full sample subset - no IDs no covariates

manic_df_full_sample <- main_manic_df[c(manic_derived_vars)]

manic_df_full_sample[] <- lapply(manic_df_full_sample, factor)

dim(manic_df_full_sample)

### define full sample with covariates subset - no IDs

manic_df_full_sample_covariate_sensitivity <-
  main_manic_df[c(manic_derived_vars_and_covariates)]

manic_df_full_sample_covariate_sensitivity <-
  manic_df_full_sample_covariate_sensitivity[which(
    main_manic_df$Ever.extreme.irritability == 1 |
      main_manic_df$Ever.manic.or.excitable == 1 |
      main_manic_df$Ever.extreme.irritability == 0 |
      main_manic_df$Ever.manic.or.excitable == 0
  ), ]

manic_df_full_sample_covariate_sensitivity[] <-
  lapply(manic_df_full_sample_covariate_sensitivity, factor)

dim(manic_df_full_sample_covariate_sensitivity)

### subset cases only - if answered yes to either manic episode of irritability - no IDs no covariates

manic_df_cases_with_ids <-
  main_manic_df[which(main_manic_df$Ever.extreme.irritability == 1 |
    main_manic_df$Ever.manic.or.excitable == 1), ]

manic_df_cases <- manic_df_cases_with_ids[c(manic_derived_vars)]

manic_df_cases[] <- lapply(manic_df_cases, factor)

saveRDS(manic_df_cases_with_ids, file = paste0(output_dir, "/", "lca_cases_only_analytical_dataset_protect.rds"))

### define cases only with covariates - no IDs

manic_df_cases_covariate_sensitivity <-
  main_manic_df[which(
    main_manic_df$Ever.extreme.irritability == 1 | main_manic_df$Ever.manic.or.excitable == 1
    & (!is.na(main_manic_df$Ever.extreme.irritability) & !is.na(main_manic_df$Ever.manic.or.excitable)
      )
  ), ]

manic_df_cases_covariate_sensitivity <-
  manic_df_cases_covariate_sensitivity[c(manic_derived_vars_and_covariates)]

manic_df_cases_covariate_sensitivity[] <-
  lapply(manic_df_cases_covariate_sensitivity, factor)

dim(manic_df_cases_covariate_sensitivity)

### sensitivity subsets (adjustment for covariates) will be different sample size to main analyses
### since includes only individuals that responded to both questions (manic episode and irritability)

#### checks

dim(manic_df_cases)
str(manic_df_cases)
tabyl(manic_df_cases[, 3])

dim(manic_df_cases_covariate_sensitivity)

tabyl(manic_df_cases_covariate_sensitivity[, 10])
dim(manic_df_cases_covariate_sensitivity)
dim(manic_df_full_sample_covariate_sensitivity)

###########################################################################
######################## DESCRIBE ANALYTICAL SUBSETS ######################

sensitivity_subsets <- c(
  "manic_df_cases_covariate_sensitivity"
)

base_subsets <- c("manic_df_cases", "manic_df_full_sample")

all_subsets <- c(base_subsets, sensitivity_subsets)

for (subsetn in all_subsets) {

  tb <- list()

  for (var in manic_derived_vars) {

    tb1 <- tabyl(get(subsetn)[[var]])

    names(tb1)[1] <- var

    tb[[var]] <- tb1

  }

  ### collapse list and rename variables

  tb2 <- bind_rows(tb, .id = c("n"))[, 1:4]

  names(tb2) <- c("Response", "Value", "N", "percent")

  tb2$`N (percent)` <- paste0(tb2$N, " (", round(tb2$percent, digits = 2), ")")

  assign(paste0("stats_", subsetn), tb2)

  write.table(
    tb2,
    file = paste0(output_dir, "/", "stats_", subsetn, "_protect.txt"),
    row.names = F
  )

}


### descriptive stats of covariates in sensitivity analysis subsets

for (subsetn in sensitivity_subsets) {

  tb <- list()

  for (var in covars) {

    tb1 <- tabyl(get(subsetn)[[var]])

    names(tb1)[1] <- var

    tb[[var]] <- tb1

  }

  ### collapse list and rename variables

  tb2 <- bind_rows(tb, .id = c("n"))

  names(tb2) <- c("Response", "Value", "N", "percent", "valid_percent")

  tb2$`N (percent)` <- paste0(tb2$N, " (", round(tb2$percent, digits = 2), ")")

  assign(paste0("covar_stats_", subsetn), tb)

  write.table(
    tb2,
    file = paste0(output_dir, "/", "covar_stats_", subsetn, "_protect.txt"),
    row.names = F
  )

}


##########################################################################
############################## ANALYSIS ##################################

### create formulas for all models
measurevar <- "1"

groupvars <- manic_derived_vars

### create formula for base model
f <- as.formula(paste0("cbind(", paste(groupvars, collapse = ","), ")~", measurevar))

print(f)

### initiate cluster for parallel process

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

all_lca_cases <- foreach(j = 2:(n_series), .packages = "poLCA") %dopar% {
  res <-
    poLCA(
      formula = f,
      data = manic_df_cases,
      nclass = j,
      maxiter = 10000,
      nrep = 20
    )
  return(res)
}

### save the main results at this point

saveRDS(all_lca_cases, file = paste0(output_dir, "/", "lca_cases_only_main_results_protect.rds"))

stopImplicitCluster()

gc()


