### Latent class analysis of responses to manic or irritable symptoms
### 05/10/2020
### Ryan Arathimos

### load packages
library(data.table)
library(ukbkings)
library(tidyverse)
library(ggplot2)
library(stringi)
library(FactoMineR)
library(factoextra)
library(poLCA)
library(parallel)
library(foreach)
library(doParallel)
library(janitor)
library(broom)
library(gridExtra)

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

### sensitivty analysis with full sample - both adjusted for mania/irritability and not [x2 models]
### cases only, with and without adjustment for mania and irritability [x2 models]
### to run LCA - 2 classes represents separation between manic/not-manic individuals in full dataset, >2 needed

n_series <- 7 # maximum number of classes to attempt

df <-
  readRDS(paste0(output_dir, "/", "manic_data_derived_subset_mhq.rds"))

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
main_manic_df_ids <- main_manic_df$eid

#############################################################
################# DEFINE ANALYTICAL SAMPLES #################

### define full sample subset - no IDs no covariates

manic_df_full_sample <- main_manic_df[c(manic_derived_vars)]

manic_df_full_sample[] <- lapply(manic_df_full_sample, factor)

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

### subset cases only - if answered yes to either manic episode of irritability - no IDs no covariates

manic_df_cases_with_ids <-
  main_manic_df[which(main_manic_df$Ever.extreme.irritability == 1 |
    main_manic_df$Ever.manic.or.excitable == 1), ]

manic_df_cases <- manic_df_cases_with_ids[c(manic_derived_vars)]

manic_df_cases[] <- lapply(manic_df_cases, factor)

saveRDS(manic_df_cases_with_ids, file = paste0(output_dir, "/", "lca_cases_only_analytical_dataset.rds"))

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
  "manic_df_cases_covariate_sensitivity",
  "manic_df_full_sample_covariate_sensitivity"
)

base_subsets <- c("manic_df_cases", "manic_df_full_sample")

all_subsets <- c(base_subsets, sensitivity_subsets)

for (subsetn in all_subsets) {
  tball <- list()

  for (var in manic_derived_vars) {
    tb1 <- tabyl(get(subsetn)[[var]])
    names(tb1)[1] <- var

    tball[[var]] <- tb1
  }

  ### collapse list and rename variables
  tball2 <- bind_rows(tball, .id = c("n"))[, 1:4]
  names(tball2) <- c("Response", "Value", "N", "percent")
  tball2$`N (percent)` <-
    paste0(tball2$N, " (", round(tball2$percent, digits = 2), ")")

  assign(paste0("stats_", subsetn), tball2)

  write.table(
    tball2,
    file = paste0(output_dir, "/", "stats_", subsetn, ".txt"),
    row.names = F
  )
}


### descriptive stats of covariates in sensitivity analysis subsets
for (subsetn in sensitivity_subsets) {
  tball <- list()

  for (var in covars) {
    tb1 <- tabyl(get(subsetn)[[var]])
    names(tb1)[1] <- var

    tball[[var]] <- tb1
  }

  ### collapse list and rename variables
  tball2 <- bind_rows(tball, .id = c("n"))
  names(tball2) <- c("Response", "Value", "N", "percent", "valid_percent")
  tball2$`N (percent)` <-
    paste0(tball2$N, " (", round(tball2$percent, digits = 2), ")")

  assign(paste0("covar_stats_", subsetn), tball)

  write.table(
    tball2,
    file = paste0(output_dir, "/", "covar_stats_", subsetn, ".txt"),
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

###  create formula for mddel sensitivity with adjustments for yes/no mani episode and yes/no irritability
#### for models with covariates, replace the ~1 with the desired function of predictors iv1,iv2,iv3... eg cbind(dv1,dv2,dv3)~iv1+iv2*iv3.

f_adjusted <- as.formula(paste0("cbind(", paste(groupvars, collapse = ","), ")~", paste0(paste(covars, collapse = " + "))))

# initiate cluster for parallel process
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

saveRDS(all_lca_cases, file = paste0(output_dir, "/", "lca_cases_only_main_results_tmp.rds"))

all_lca_cases_covariate_sensitivity <-
  foreach(j = 2:(n_series), .packages = "poLCA") %dopar% {
    res <-
      poLCA(
        formula = f_adjusted,
        data = manic_df_cases_covariate_sensitivity,
        nclass = j,
        maxiter = 10000,
        nrep = 20
      )
    return(res)
  }

saveRDS(all_lca_cases_covariate_sensitivity, file = paste0(output_dir, "/", "lca_cases_covariate_sensitivity_results_tmp.rds"))

all_lca_full_sample <-
  foreach(j = 2:(n_series), .packages = "poLCA") %dopar% {
    res <-
      poLCA(
        formula = f,
        data = manic_df_full_sample,
        nclass = j,
        maxiter = 10000,
        nrep = 20
      )
    return(res)
  }

saveRDS(all_lca_full_sample, file = paste0(output_dir, "/", "lca_full_sample_sensitivity_results_tmp.rds"))

# all_lca_full_sample_covariate_sensitivity <-
#  foreach(j = 2:(n_series), .packages = 'poLCA') %dopar% {
#    res <-
#      poLCA(
#        formula = f_adjusted,
#        data = manic_df_full_sample_covariate_sensitivity,
#        nclass = j,
#        maxiter = 10000,
#        nrep = 20
#      )
#    return(res)
#  }

stopImplicitCluster()

gc()

### list the analyses conducted separately to loop through for descriptions
list_of_analyses_filenames <- c(
  "lca_cases_only_main", "lca_cases_covariate_sensitivity",
  "lca_full_sample_sensitivity"
)


##################################################################################
######################### PLOT & DESCRIBE RESULTS ################################

#lca <- readRDS(file = paste0(output_dir, "/", "lca_cases_only_main_results_tmp.rds"))

### relative entropy function
### https://link.springer.com/article/10.1007%2FBF01246098 - relative entropy

relative_entropy <- function(lcares) {
  lcares2 <- lcares
  entropy <- function(p) sum(-p * log(p))
  error_prior <- entropy(lcares2[["P"]]) # Class proportions
  error_post <- mean(apply(lcares2[["posterior"]], 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  R2_entropy
}

### function to get lca results in to shape for plotting
prepare_lca_res_for_plotting <- function(lcares) {
  res_optimum <- lcares

  cond_df <- as.data.frame(res_optimum$probs)
  cond_dfse <- as.data.frame(res_optimum$probs.se)

  names1 <- c(paste0(names(res_optimum$probs), ".1"))

  cond_df1 <- as.data.frame(t(cond_df[, names1]))
  cond_dfse1 <- as.data.frame(t(cond_dfse[, names1]))
  cond_dfse1$response <- row.names(cond_dfse1)
  cond_df1$response <- row.names(cond_df1)

  cond_dfse1$stat <- "SE"
  cond_df1$stat <- "probability"

  colnames(cond_dfse1) <- colnames(cond_df1)

  cond_df2 <- melt(cond_df1, idvars = c("response"))
  cond_dfse2 <- melt(cond_dfse1, idvars = c("response"))

  cond_df3 <- merge(cond_df2, cond_dfse2, by = c("response", "variable"))

  names(cond_df3)[names(cond_df3) == "value.x"] <- "Probability"

  names(cond_df3)[names(cond_df3) == "value.y"] <- "SE"

  names(cond_df3)[names(cond_df3) == "variable"] <- "Classes"

  cond_df3$Response <- paste0(sapply(strsplit(cond_df3$response, "[.]"), "[[", 1), " ", sapply(strsplit(cond_df3$response, "[.]"), "[[", 2))

  cond_df3$percentage <- round(res_optimum$P * 100, digits = 1)

  cond_df3$Classes <- paste0(cond_df3$Classes, cond_df3$percentage, "%")

  return(cond_df3)
}

for (analysis_i in list_of_analyses_filenames) {

  all_lca <- readRDS(file = paste0(output_dir, "/",analysis_i ,"_results_tmp.rds"))

  ### storage vectors for fit stats
  bic_vector <- numeric()
  aic_vector <- numeric()
  numiter_vector <- numeric()
  start_attempts_vector <- numeric
  relative_entropy_vector <- numeric()
  poLCA_entropy_vector <- numeric()
  number_of_classes <- numeric()
  N <- numeric()

  ### calculate fit stats
  for (i in 1:length(all_lca)) {
    bic_vector <- c(bic_vector, all_lca[[i]]$bic)
    aic_vector <- c(aic_vector, all_lca[[i]]$aic)
    numiter_vector <- c(numiter_vector, all_lca[[i]]$numiter)
    # start_attempts_vector <- c(start_attempts_vector, NROW(all_lca[[i]]$attempts))
    poLCA_entropy_vector <- c(poLCA_entropy_vector, poLCA.entropy(all_lca[[i]]))
    relative_entropy_vector <- c(relative_entropy_vector, relative_entropy(all_lca[[i]]))
    number_of_classes <- c(number_of_classes, NROW(unique(all_lca[[i]]$predclass)))
    N <- c(N, (all_lca[[i]]$N))
    
    currentclassn <- NROW(unique(all_lca[[i]]$predclass))

    ### apply preparation functionreally inefficient function, could be recoded
    lca_res_plots <- prepare_lca_res_for_plotting(all_lca[[i]])

    ### Manually relevel responses for plot x axis
    lca_res_plots$Response <- factor(lca_res_plots$Response, levels = c(
      "More active", "More talkative", "More confident",
      "More creative", "Less sleep", "Thoughts racing",
      "More restless", "Easily distracted"
    ))

    gg1 <- ggplot(lca_res_plots, aes(Response, Probability, group = Classes, color = Classes)) +
      geom_line() +
      geom_errorbar(aes(ymin = (Probability - (1.96 * SE)), ymax = (Probability + (1.96 * SE))),
        width = .4,
        position = position_dodge(.001)
      ) +
      labs(colour = "Classes (% of sample)") +
      scale_color_brewer(palette = palette_choice) +
      theme_classic()


    png(filename = paste0(output_dir, "/", analysis_i, "_line_plot_", currentclassn, "_object", i, ".png"), width = 24, height = 12, units = "cm", res = 300)
    print(gg1)
    dev.off()
  }

  fit_stats <- rbind(
    as.character(number_of_classes), round(as.numeric(bic_vector), digits = 1), round(as.numeric(aic_vector), digits = 1),
    as.character(numiter_vector), round(as.numeric(poLCA_entropy_vector), digits = 2),
    round(as.numeric(relative_entropy_vector), digits = 3), as.numeric(N)
  )


  colnames(fit_stats) <- c(paste0("model_", 2:n_series, "_classes"))

  row.names(fit_stats) <- c("Number of classes", "BIC", "AIC", "N iterations", "Entropy", "Relative entropy")


  gg1 <- ggplot(as.data.frame(t(fit_stats)), aes(
    x = as.numeric(as.character(`Number of classes`)), y =
      as.numeric(as.character(`BIC`))
  )) +
    geom_point(color = "cornflowerblue", alpha = 0.8) +
    geom_path(color = "cornflowerblue") +
    theme_classic() +
    xlab("Number of classes") +
    ylab("BIC")

  gg2 <- ggplot(as.data.frame(t(fit_stats)), aes(
    x = as.numeric(as.character(`Number of classes`)), y =
      as.numeric(as.character(`AIC`))
  )) +
    geom_point(color = "cornflowerblue", alpha = 0.8) +
    geom_path(color = "cornflowerblue") +
    theme_classic() +
    xlab("Number of classes") +
    ylab("AIC")

  png(filename = paste0(output_dir, "/", analysis_i, "_bic_aic_skreeplot_.png"), width = 24, height = 10, units = "cm", res = 300)
  grid.arrange(gg1, gg2, nrow = 1)
  dev.off()

  write.table(fit_stats, paste0(output_dir, "/", analysis_i, "_fit_stats.txt"), col.names=F)
  
}
