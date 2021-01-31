### Post LCA estimation analyses in PROTECT
### Associations of latent classes with other variables

### Ryan Arathimos
### 23/11/2020

library(data.table)
library(ukbkings)
library(tidyverse)
library(ggplot2)
library(janitor)
library(broom)
library(gridExtra)
library(RColorBrewer)
library(caret)
library(nnet)
library(readstata13)
library(RNOmni)

set.seed(123456)

### SET UP DIRECTORIES

dirs <- read.table("~/brc_scratch/scripts/paths.txt")

prs_dir <- as.character(dirs[1, ])
new_2019_qc_dir <- as.character(dirs[3, ])
qc_dir <- as.character(dirs[3, ])
project_dir <- as.character(dirs[9, ])
dep_phen_dir <- as.character(dirs[6, ])
psychosis_phen_dir <- as.character(dirs[7, ])
prs_protect_dir <- as.character(dirs[11, ])
protect_qc_dir <- as.character(dirs[12, ])

data_dir <- as.character(paste0(dirs[10, ], "data/manic_lca"))
output_dir <- as.character(paste0(dirs[10, ], "output/manic_lca"))
scripts_dir <- as.character(paste0(dirs[10, ], "scripts/manic_lca"))

### pick a color palette for plots - options from rcolorbrewer palette eg. "Dark2"

palette_choice <- "Set1"

### set up functions for multinomial logistic regressions

source(paste0(scripts_dir, "/function_categorical_multinomial.r"))

source(paste0(scripts_dir, "/function_continuous_multinomial.r"))

source(paste0(scripts_dir, "/function_binary_multinomial.r"))

`%notin%` <- Negate(`%in%`)

### optimum model - i.e. number of classes which optimally fit the data

optimum_classes <- 5
optimum_model_in_res <- optimum_classes - 1 # the object number in the list is numbered (N classes -1)

### load phenotype data

sociodem <- readRDS(paste0(output_dir, "/", "sociodem_data_derived_protect.rds"))

### load lca results

lca <- readRDS(file = paste0(output_dir, "/", "lca_cases_only_main_results_protect.rds"))

### load the subset of the dataset with the ids

ids <- readRDS(file = paste0(output_dir, "/", "lca_cases_only_analytical_dataset_protect.rds"))


################################################################################
######## TABULATE NUMBERS OF IRRITABLE AND MANIC IN EACH CLASS WITH % ##########

if (NROW(lca[[optimum_model_in_res]]$posterior) != NROW(ids)) {
  print("Number of rows in the ID file do not match the number of samples in the LCA results!~~~")
}

posterior_probs <- as.data.frame(lca[[optimum_model_in_res]]$posterior)

posterior_probs$predclass <- as.numeric(as.character(lca[[optimum_model_in_res]]$predclass))

### manually relevel class numbering to match discovery - this will change with each LCA rerun

posterior_probs[["predclass"]] <- recode(posterior_probs[["predclass"]], `1` = 2, `2` = 3, `3` = 4, `4` = 1, `5` = 5)

### also manually rename probabilities to match new class numberings

names(posterior_probs) <- c("V2","V3","V4","V1","V5","predclass")

posterior_probs <- posterior_probs[,c("V1","V2","V3","V4","V5","predclass")]

### bind phenotypes and LCA results

posterior_probs <- cbind(posterior_probs, ids)

tabyl(posterior_probs$Ever.manic.or.excitable)
tabyl(posterior_probs$Ever.extreme.irritability)

### merge sociodemographic phenotypes

dim(posterior_probs)

posterior_probs <- merge(posterior_probs, sociodem, by = c("ResultsID", "TecID"), all.x = T)

dim(posterior_probs)

### create new composite variable for tabulations

posterior_probs$Ever.manic.and.irritable <- NA
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable == 1] <- 1
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.extreme.irritability == 1] <- 2
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable == 1 & posterior_probs$Ever.extreme.irritability == 1] <- 3

tabyl(posterior_probs$Ever.manic.and.irritable)

###

tb <- list()

for (i in 1:optimum_classes) {
  tb1 <- tabyl(posterior_probs$Ever.manic.and.irritable[posterior_probs$predclass == i])
  names(tb1)[1] <- paste0("class_", i)
  tb1$categories <- c("Ever manic", "Ever irritable", "Ever manic and irritable")
  tb[[i]] <- tb1
}

tb2 <- bind_rows(tb, .id = c("class"))
tb2 <- tb2[c("categories", "n", "percent", "class")]
tb2$`N (percent)` <-
  paste0(tb2$n, " (", round(tb2$percent*100, digits = 2), ")")
tb2$Response <- as.factor(as.character(gsub("\\.", " ", as.character(tb2$categories)))) # remove periods if any exist



### plot numbers/percentages of irritable, manic or manic AND irritable in each class

png(
  filename = paste0(output_dir, "/", "percent_manic_or_irritable_", optimum_classes, "class_protect_barplot.png"),
  width = 20, height = 13, units = "cm", res = 300
)

ggplot(tb2, aes(class, percent)) +
  geom_bar(aes(fill = Response), position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Proportion of class") +
  theme_classic() +
  geom_text(
    aes(
      class,
      y = percent + 0.03,
      group = Response,
      label = format(
        percent,
        nsmall = 0,
        digits = 1,
        scientific = FALSE, size = 4
      )
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
  )

dev.off()

png(
  filename = paste0(output_dir, "/", "numbers_manic_or_irritable_", optimum_classes, "class_protect_barplot.png"),
  width = 20, height = 13, units = "cm", res = 300
)

ggplot(tb2, aes(class, n)) +
  geom_bar(aes(fill = Response), position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Number of individuals in class") +
  theme_classic()


dev.off()

### make a table

print((tb2[c("class", "Response", "N (percent)")]))

write.table(
  tb2[c("class", "Response", "N (percent)")],
  file = paste0(
    output_dir,
    "/",
    "tabulations_manic_or_irritable_",
    optimum_classes,
    "class_protect.csv"
  ),
  sep = ",",
  row.names = F,
  quote = F
)

### could also create a composite variable based on contribution by summed class probabilities


########################################################################
######## ASSOCIATIONS WITH DURATION OF MANIC SYMPTOMATOLOGY ############

### dummy code the predicted class variable - comparisons between each class and all others

class(posterior_probs$predclass)
posterior_probs$predclass <- as.factor(as.character(posterior_probs$predclass))

res_dummy <- as.data.frame(model.matrix(~ predclass - 1, data = posterior_probs))

### checks
head(res_dummy)
str(res_dummy)
tabyl(res_dummy[, c("predclass2")])

posterior_probs <- cbind(posterior_probs, res_dummy)
tabyl(posterior_probs$Brief.duration.mania.or.irritability)

### reshape data to long so that participants have as many rows as classes and probability for each row
posterior_long <- reshape(posterior_probs,
  idvar = "ResultsID", direction = "long",
  varying = 3:7, v.names = (f <- c("probabilities")), timevar = "predclass"
)

### check reshape to long for a random participant
posterior_long$probabilities[posterior_long$TecID == sample(posterior_long$TecID, 1)] # should have 5 probabilities sum to 1
if (round(sum(posterior_long$probabilities[posterior_long$TecID == sample(posterior_long$TecID, 1)]), digits = 3) != 1) print("Error with reshape!")

posterior_long$predclass[posterior_long$TecID == sample(posterior_long$TecID, 1)] # should have values for X number of classes

duration_naive_res <- list()
duration_weighted_res <- list()

posterior_probs$predclass <- as.factor(as.character(posterior_probs$predclass))
posterior_long$predclass <- as.factor(as.character(posterior_long$predclass))

for (reference_i in c("main", "relevel")) {

  ### if loop is 'relevel' then change the reference class to "3"

  if (reference_i == "relevel") {
    posterior_probs$predclass <-
      relevel(as.factor(as.character(posterior_probs$predclass)), ref = "3")

    posterior_long$predclass <-
      relevel(as.factor(as.character(posterior_long$predclass)), ref = "3")
  } else {
    print("This is the main analysis (reference level '1')")
  }


  ### create formula for multinomial regression
  for (duration_i in c(
    "Brief.duration.mania.or.irritability",
    "Moderate.duration.mania.or.irritability",
    "Extended.duration.mania.or.irritability"
  )) {
    formula_duration <- paste0("predclass ~ ", duration_i)

    ### multinomial logistic regression naive to probabilities
    mres <- multinom(formula_duration, data = posterior_probs)
    tidy_mres <- tidy(mres, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)
    tidy_mres$model <- "naive"
    tidy_mres$exposure <- duration_i

    duration_naive_res[[duration_i]] <- tidy_mres

    ### multinomial logistic regression weighted for probabilities
    mres_weighted <- multinom(formula_duration, data = posterior_long, weights = probabilities)

    # extract and exponentiate coefficients
    tidy_mres_weighted <- tidy(mres_weighted, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)
    tidy_mres_weighted$model <- "weighted"
    tidy_mres_weighted$exposure <- duration_i

    duration_weighted_res[[duration_i]] <- tidy_mres_weighted
  }

  ### bind and drop intercept from results
  duration_naive_res_df <- bind_rows(duration_naive_res)
  duration_naive_res_df <- duration_naive_res_df[which(duration_naive_res_df$term != "(Intercept)"), ]

  duration_weighted_res_df <- bind_rows(duration_weighted_res)
  duration_weighted_res_df <- duration_weighted_res_df[which(duration_weighted_res_df$term != "(Intercept)"), ]

  ### create a 'comparison' column of statistical test for plot
  duration_naive_res_df$comparison <- paste0(levels(posterior_probs[["predclass"]])[1], "v", duration_naive_res_df$`y.level`)
  duration_weighted_res_df$comparison <- paste0(levels(posterior_long[["predclass"]])[1], "v", duration_naive_res_df$`y.level`)

  ### remove periods if any exist from names
  duration_naive_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(duration_naive_res_df$exposure)))
  duration_weighted_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(duration_weighted_res_df$exposure)))
 
  write.csv(duration_weighted_res_df, paste0(output_dir,"/","multinomial_duration_weighted_res_",reference_i,"_protect.csv"), col.names=T, row.names=F, quote = F)
  
  ### plot multinomial results
  png(
    filename = paste0(output_dir, "/", "multinomial_naive_duration_", reference_i, "_", optimum_classes, "class_protect_scatterplot.png"),
    width = 20, height = 8, units = "cm", res = 300
  )

  print(ggplot(duration_naive_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Estimate (RR)") +
    theme_classic() +
    facet_wrap(~exposure) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white")
    ))

  dev.off()

  ###

  png(
    filename = paste0(output_dir, "/", "multinomial_weighted_duration_", reference_i, "_", optimum_classes, "class_protect_scatterplot.png"),
    width = 20, height = 8, units = "cm", res = 300
  )

  print(ggplot(duration_weighted_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Estimate (RR)") +
    theme_classic() +
    facet_wrap(~exposure) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white")
    ))

  dev.off()

  print("Plotted mutinomial results!")

  ### duration summarise by most probable class

  tb_duration <- list()

  for (i in 1:optimum_classes) {
    tb1 <- tabyl(posterior_probs$Duration.mania.or.irritability[posterior_probs$predclass == i])

    names(tb1)[1] <- paste0("class_", i)

    tb1$categories <- c("1 - Brief", "2 - Moderate", "3 - Extended", "Don't know NA")

    tb_duration[[i]] <- tb1
  }

  tb_duration2 <- bind_rows(tb_duration, .id = c("class"))
  tb_duration2 <- tb_duration2[c("categories", "n", "percent", "class")]
  tb_duration2$`N (percent)` <-
    paste0(tb_duration2$n, " (", round(tb_duration2$percent*100, digits = 2), ")")
  tb_duration2$Duration <- as.factor(gsub("\\.", " ", as.character(tb_duration2$categories))) # remove periods if any exist

  ### tabulate number of missing rows (hence individuals with NA for duration) and subset
  missing_duration <- NROW(posterior_probs$Duration.mania.or.irritability) - NROW(tb_duration2[which(posterior_probs$Duration.mania.or.irritability != "NA"), ])
  tb_duration2 <- tb_duration2[which(tb_duration2$categories != "NA"), ]

  write.table(tb_duration2[c("categories", "n", "percent", "class", "N (percent)")],
              file = paste0(output_dir, "/",  "duration_protect", optimum_classes, "class_tabulations.csv"),
              row.names = FALSE, quote = FALSE, sep = ","
  )
  
  ### plot raw comparisons by most probable class
  png(
    filename = paste0(output_dir, "/", "percent_duration_", optimum_classes, "class_protect_barplot.png"),
    width = 20, height = 13, units = "cm", res = 300
  )

  print(ggplot(tb_duration2, aes(class, percent)) +
    geom_bar(aes(fill = Duration), position = "dodge", stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Proportion of class") +
    theme_classic() +
    xlab(paste0("Class (N= ", NROW(posterior_probs[which(posterior_probs[["Duration.mania.or.irritability"]] != "NA"), ]), ")")) +
    geom_text(aes(class,
      y = percent + 0.03,
      group = Duration,
      label = format(
        percent,
        nsmall = 0,
        digits = 1,
        scientific = FALSE, size = 4
      )
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
    ))

  dev.off()
}


##############################################################################
######## ASSOCIATIONS WITH DISRUPTIVENESS OF MANIC SYMPTOMATOLOGY ############

# assumes there is a var 'predclass' and var 'probabilities' in df for regressions
source(paste0(scripts_dir, "/function_binary_multinomial.r"))

### rename levels for plot before applying function
posterior_probs$Problematic.mania.or.irritability[posterior_probs$Problematic.mania.or.irritability == 0] <- "Not disruptive"
posterior_probs$Problematic.mania.or.irritability[posterior_probs$Problematic.mania.or.irritability == 1] <- "Disruptive"

tabyl(posterior_probs$Problematic.mania.or.irritability)

analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Problematic.mania.or.irritability",
  traitname = "disruptiveness_protect"
)

####################################################################################
############### ASSOCIATIONS WITH SOCIODEMOGRAPHIC FACTORS #########################

### create subset of data for sociodemographic analyses

sociodemographic_vars <- c(
  "Sex", "Education.level",
  "Smoking", "Alcohol.intake.frequency"
)


posterior_probs_sociodem <- na.omit(posterior_probs[, c(
  "TecID", "ResultsID",
  sociodemographic_vars, "predclass"
)])

dim(posterior_probs_sociodem)

posterior_long_sociodem <- na.omit(posterior_long[, c("TecID", sociodemographic_vars, "predclass", "probabilities")])

dim(posterior_long_sociodem)

if ((NROW(posterior_long_sociodem) / optimum_classes) != NROW(posterior_probs_sociodem)) print("Error: Long and wide datasets have different numbers of participants!")

### tabulate/describe sample of sociodemographics

tabulated_sociodem <- list()

tabulated_sociodem[["categorical"]] <- bind_rows(
  tabyl(posterior_probs_sociodem$Sex),
  tabyl(posterior_probs_sociodem$Education.level),
  tabyl(posterior_probs_sociodem$Alcohol.intake.frequency),
  tabyl(posterior_probs_sociodem$Smoking)
)

### apply multinomial analysis function to sex

Sex_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Sex",
  traitname = "sex_protect"
)


### apply multinomial analysis function to education level (categorical)

tabyl(posterior_probs_sociodem$Education.level)
class(posterior_probs_sociodem$Education.level)

Education.level_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Education.level",
  traitname = "education_protect"
)

### apply multinomial analysis function to smoking (categorical)

tabyl(posterior_probs_sociodem$Smoking)
class(posterior_probs_sociodem$Smoking)

Smoking_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Smoking",
  traitname = "smoking_protect"
)

### apply multinomial analysis function to alcohol consumption (categorical)

tabyl(posterior_probs_sociodem$Alcohol.intake.frequency)
class(posterior_probs_sociodem$Alcohol.intake.frequency)

Alcohol.intake.frequency_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Alcohol.intake.frequency",
  traitname = "alcohol_frequency_protect"
)


##############################################################################
######## ASSOCIATIONS WITH PRS OF PSYCHIATRIC DISORDERS ######################

### PRS to be used - with vector of names (in same order)

prs_filenames <- c("BIPO01", "SCHI02", "DEPR07", "ADHD05", "AUTI05", "ANXI02")

prs_names <- c(
  "Bipolar disorder PRS", "Schizophrenia PRS", "Depression PRS",
  "ADHD PRS", "ASD PRS", "Anxiety PRS"
)

### load PRS from files keep only those calculated at pvalue=1 threshold

pvalue_threshold <- "1" # can be "0.05" or other

prs_raw <- list()

for (i_prs in 1:length(prs_filenames)) {
  prs_dat <- setDF(fread(paste0(prs_protect_dir, prs_filenames[i_prs], ".all.score"), header = T))

  prs_dat <- prs_dat[c("IID", paste0(pvalue_threshold))]

  names(prs_dat) <- c("IID", c(paste0(prs_filenames[i_prs], "_", pvalue_threshold)))

  prs_raw[[prs_filenames[i_prs]]] <- prs_dat
}

### ensure same IDs across PRS files (no withdrawn IDs in some and not others)

id_freqs <- data.frame(table(unlist(lapply(prs_raw, "[[", 1))))

exclude_ids <- as.character(id_freqs$Var1[which(as.numeric(as.character(id_freqs$Freq)) < NROW(prs_raw))])

### remove any IDs that do not fit

prs_raw <- lapply(prs_raw, function(x) x[which(as.character(x[, 1]) %notin% exclude_ids), ])

prs_df <- bind_cols(prs_raw)

### load PCs and covs

covs <-
  setDF(fread(
    paste0(protect_qc_dir, "pcs/europeans/top20_pcs_reformatted_europeans.csv")
  ))

cov_names <- paste0("pc", rep(1:6))

### merge PRS and covariates

for_pcs_regression <- merge(prs_df, covs, by = c("IID"))

### regress out PCs

prs.res <- list()

for (i_trait in c(paste0(prs_filenames, "_", pvalue_threshold))) {
  prs.lm <- lm(paste0(i_trait, " ~ ", paste0(cov_names, collapse = " + ")),
    data = for_pcs_regression
  )

  prs.res[[i_trait]] <- resid(prs.lm)
}

prs.res_df <- as.data.frame(prs.res)

hist(prs.res_df[[paste0("SCHI02_", pvalue_threshold)]], breaks = 500)

### standardise PRS residuals

prs.res_df[names(prs.res_df)] <-
  lapply(prs.res_df[names(prs.res_df)], function(x) {
    c(scale(x))
  })

### add ID column back in

prs.res_df$IID <- for_pcs_regression$IID

### check conforms to normal distribution with SD=1

hist(prs.res_df[[paste0("SCHI02_", pvalue_threshold)]], breaks = 500)

sd(prs.res_df[[paste0("SCHI02_", pvalue_threshold)]])

### merge with posterior probs and class memembership

posterior_probs_prs <- merge(posterior_probs, prs.res_df, by.x = "TecID", by.y = "IID")

### find which samples weren't found in PRS

samples_without_prs <- NROW(posterior_probs) - NROW(posterior_probs_prs)

print(paste0("There are ", samples_without_prs, " number of samples in latent class analysis without PRS"))

ids_without_prs <- posterior_probs$TecID[posterior_probs$TecID %notin% posterior_probs_prs$TecID]

print(NROW(ids_without_prs))

### reshape data to long so that participants have as many rows as classes and a probability for each row/class

posterior_prs_long <- reshape(posterior_probs_prs,
  idvar = "TecID", direction = "long",
  varying = 3:7, v.names = (f <- c("probabilities")), timevar = "predclass"
)

### check reshape to long correct by picking a random participant ID

posterior_prs_long$probabilities[posterior_prs_long$TecID == sample(posterior_prs_long$TecID, 1)] # should have 5 probabilities sum to 1

if (round(sum(posterior_prs_long$probabilities[posterior_prs_long$TecID == sample(posterior_prs_long$TecID, 1)]), digits = 6) != 1) print("Error with reshape!")

### check that values for X number of classes present (ie 5 if optimum model is 5 class)

posterior_prs_long$predclass[posterior_prs_long$TecID == sample(posterior_prs_long$TecID, 1)]

### save stata dataset for cross-checking results using mlogit in STATA
# save.dta13(posterior_prs_long[, c("TecID", "probabilities", "BIPO01_1", "predclass")],
#  file = "posterior_prs_long.dta"
# )

### multinom regressions of PRS on classes - create blank lists for results

prs_naive_res <- list()

prs_weighted_res <- list()

### loop over either having class 1 or class 3 as reference class

posterior_probs_prs$predclass <- as.factor(as.character(posterior_probs_prs$predclass))
levels(posterior_probs_prs$predclass)

posterior_prs_long$predclass <- as.factor(as.character(posterior_prs_long$predclass))
levels(posterior_prs_long$predclass)

posterior_probs_prs$predclass <-
  relevel(as.factor(as.character(posterior_probs_prs$predclass)), ref = "1")

posterior_long$predclass <-
  relevel(as.factor(as.character(posterior_long$predclass)), ref = "1")

for (reference_i in c("main", "relevel")) {

  ### if loop is 'relevel' then change the reference class to "3"

  if (reference_i == "relevel") {
    posterior_probs_prs$predclass <-
      relevel(as.factor(as.character(posterior_probs_prs$predclass)), ref = "3")

    posterior_prs_long$predclass <-
      relevel(as.factor(as.character(posterior_prs_long$predclass)), ref = "3")
  } else {
    posterior_probs_prs$predclass <-
      relevel(as.factor(as.character(posterior_probs_prs$predclass)), ref = "1")

    posterior_prs_long$predclass <-
      relevel(as.factor(as.character(posterior_prs_long$predclass)), ref = "1")

    print("This is the main analysis (reference level '1')")
  }

  print(levels(posterior_probs_prs$predclass))
  print(levels(posterior_prs_long$predclass))

  for (prs_num in 1:NROW(prs_filenames)) {

    ### pull out filename and associated trait name (same order vectors)

    prs_trait_i <- prs_filenames[prs_num]

    prs_trait_name <- prs_names[prs_num]

    ### create regression formula

    formula_prs <- formula(paste0("predclass ~ ", prs_trait_i, "_", pvalue_threshold))

    ### multinomial logistic regression naive to probabilities

    mres <- multinom(formula_prs, data = posterior_probs_prs)
    tidy_mres <-
      tidy(mres,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )

    tidy_mres$model <- "naive"

    tidy_mres$exposure <- prs_trait_name

    prs_naive_res[[prs_trait_i]] <- tidy_mres

    ### multinomial logistic regression weighted for probabilities (bias adjusted)

    mres_weighted <-
      multinom(formula_prs, data = posterior_prs_long, weights = probabilities)

    ### extract and exponentiate coefficients, add columns for labelling purposes

    tidy_mres_weighted <-
      tidy(mres_weighted,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )

    tidy_mres_weighted$model <- "weighted"

    tidy_mres_weighted$exposure <- prs_trait_name

    prs_weighted_res[[prs_trait_i]] <- tidy_mres_weighted
  }

  ### bind and drop intercept term from results
  prs_naive_res_df <- bind_rows(prs_naive_res)

  prs_naive_res_df <-
    prs_naive_res_df[which(prs_naive_res_df$term != "(Intercept)"), ]

  prs_weighted_res_df <- bind_rows(prs_weighted_res)

  prs_weighted_res_df <-
    prs_weighted_res_df[which(prs_weighted_res_df$term != "(Intercept)"), ]

  ### create a 'comparison' column of statistical test for plot labels
  prs_naive_res_df$comparison <-
    paste0(levels(posterior_prs_long$predclass)[1], "v", prs_naive_res_df$`y.level`)

  prs_weighted_res_df$comparison <-
    paste0(levels(posterior_prs_long$predclass)[1], "v", prs_weighted_res_df$`y.level`)

  print(prs_weighted_res_df)

  ### remove periods if any exist from names
  prs_naive_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_naive_res_df$exposure)))

  prs_weighted_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_weighted_res_df$exposure)))

  write.csv(prs_weighted_res_df, paste0(output_dir,"/","multinomial_prs_weighted_res_",reference_i,"_protect.csv"), col.names=T, row.names=F, quote = F)
  
  ### plot multinomial results
  png(
    filename = paste0(output_dir, "/", "multinomial_naive_prs_", reference_i, "_", optimum_classes, "class_protect_scatterplot.png"),
    width = 20, height = 10, units = "cm", res = 300
  )

  print(ggplot(prs_naive_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Risk ratio (per SD increase in PRS)") +
    theme_classic() +
    facet_wrap(~exposure) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(hjust = 0)
    ))

  dev.off()

  ###

  png(
    filename = paste0(output_dir, "/", "multinomial_weighted_prs_", reference_i, "_", optimum_classes, "class_protect_scatterplot.png"),
    width = 20, height = 10, units = "cm", res = 300
  )

  print(ggplot(prs_weighted_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Risk ratio (per SD increase in PRS)") +
    theme_classic() +
    facet_wrap(~exposure) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(hjust = 0)
    ))

  dev.off()
}

gc()

### plot density plots of PRS by most likely class

posterior_probs_prs$predclass <- relevel(posterior_probs_prs$predclass, ref = "1")

posterior_prs_long$predclass <- relevel(posterior_prs_long$predclass, ref = "1")

for (prs_num in 1:NROW(prs_filenames)) {

  ### pull out filename and associated trait name (same order vectors)

  prs_trait_i <- prs_filenames[prs_num]

  prs_trait_name <- prs_names[prs_num]

  prs_trait_i_p <- paste0(prs_trait_i, "_", pvalue_threshold)

  gg1 <- ggplot(posterior_probs_prs, aes(x = get(prs_trait_i_p), fill = predclass)) +
    geom_density(alpha = 0.3) +
    xlab(paste0("PRS N=", NROW(posterior_probs_prs))) +
    ggtitle(paste0(prs_trait_name)) +
    labs(fill = "Class") +
    theme_classic()

  png(
    filename = paste0(output_dir, "/", "prs_", prs_trait_i_p, "_", optimum_classes, "class_protect_densityplot.png"),
    width = 20, height = 15, units = "cm", res = 300
  )

  print(gg1)

  dev.off()
}

##############################################################################################
################## ASSOCIATIONS WITH DIAGNOSIS OF PSYCHIATRIC DISORDERS ######################

diagnoses_to_check <- c(
  "Self.reported.ADHD", "Self.reported.GAD.and.others", "Self.reported.ASD",
  "Self.reported.Depression", "Self.reported.Mania.or.bipolar.disorder", "Self.reported.PsychosisAny"
)

for (i in diagnoses_to_check) {
  print(i)
  print(tabyl(ids[[i]]))
}

posterior_probs$Self.reported.Schizophrenia.or.psychosis <- factor(
  posterior_probs$Self.reported.PsychosisAny,
  levels = sort(unique(posterior_probs$Self.reported.PsychosisAny)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.Schizophrenia.or.psychosis <- factor(
  posterior_long$Self.reported.PsychosisAny,
  levels = sort(unique(posterior_long$Self.reported.PsychosisAny)),
  labels = c("No", "Yes")
)

Self.reported.Psychosis.any_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.Schizophrenia.or.psychosis",
  traitname = "schizophrenia_psychosis_protect"
)

posterior_probs$Self.reported.Mania.or.bipolar.disorder <- factor(
  posterior_probs$Self.reported.Mania.or.bipolar.disorder,
  levels = sort(unique(posterior_probs$Self.reported.Mania.or.bipolar.disorder)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.Mania.or.bipolar.disorder <- factor(
  posterior_long$Self.reported.Mania.or.bipolar.disorder,
  levels = sort(unique(posterior_long$Self.reported.Mania.or.bipolar.disorder)),
  labels = c("No", "Yes")
)

Self.reported.Mania.or.bipolar.disorder_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.Mania.or.bipolar.disorder",
  traitname = "mania_bipolar_protect"
)

posterior_probs$Self.reported.Depression <- factor(
  posterior_probs$Self.reported.Depression,
  levels = sort(unique(posterior_probs$Self.reported.Depression)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.Depression <- factor(
  posterior_long$Self.reported.Depression,
  levels = sort(unique(posterior_long$Self.reported.Depression)),
  labels = c("No", "Yes")
)

Self.reported.Depression_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.Depression",
  traitname = "depression_protect"
)

posterior_probs$Self.reported.GAD.and.others <- factor(
  posterior_probs$Self.reported.GAD.and.others,
  levels = sort(unique(posterior_probs$Self.reported.GAD.and.others)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.GAD.and.others <- factor(
  posterior_long$Self.reported.GAD.and.others,
  levels = sort(unique(posterior_long$Self.reported.GAD.and.others)),
  labels = c("No", "Yes")
)

Self.reported.GAD.and.others_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.GAD.and.others",
  traitname = "gad_protect"
)
