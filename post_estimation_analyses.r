### Post LCA estimation analyses
### Associations of latent classes with other variables

### Ryan Arathimos
### 13/10/2020

library(data.table)
library(ukbkings)
library(ggplot2)
library(janitor)
library(broom)
library(gridExtra)
library(RColorBrewer)
library(caret)
library(nnet)
library(png)
library(cowplot)
library(readstata13)
library(RNOmni)
library(tidyverse)
library(ggpubr)
library(tidyselect)
library(dplyr)

set.seed(123456)

### SET UP DIRECTORIES

dirs <- read.table("~/brc_scratch/scripts/paths.txt")

prs_dir <- as.character(dirs[1, ])
new_2019_qc_dir <- as.character(dirs[3, ])
qc_dir <- as.character(dirs[3, ])
project_dir <- as.character(dirs[9, ])
dep_phen_dir <- as.character(dirs[6, ])
psychosis_phen_dir <- as.character(dirs[7, ])

data_dir <- as.character(paste0(dirs[10, ], "data/manic_lca"))
output_dir <- as.character(paste0(dirs[10, ], "output/manic_lca"))
scripts_dir <- as.character(paste0(dirs[10, ], "scripts/manic_lca"))

### pick a color palette for plots - options from rcolorbrewer palette eg. "Dark2"

palette_choice <- "Set1"

palette_choice_ma <- c(
  "MA v AR" = "#E41A1C", "MA v EA" = "#377EB8",
  "MA v IR" = "#984EA3", "MA v FC" = "#4DAF4A")

palette_choice_ir <- c("IR v AR" = "#E41A1C", "IR v EA" = "#377EB8",
"IR v MA" = "#FF7F00", "IR v FC" = "#4DAF4A"
)

### set up functions for multinomial logistic regressions

source(paste0(scripts_dir, "/function_categorical_multinomial.r"))

source(paste0(scripts_dir, "/function_continuous_multinomial.r"))

source(paste0(scripts_dir, "/function_binary_multinomial.r"))

`%notin%` <- Negate(`%in%`)

### optimum model - i.e. number of classes which optimally fit the data

optimum_classes <- 5
optimum_model_in_res <- optimum_classes - 1 # the object number in the list is numbered (N classes -1)

### names classes of optimum model - manually rename if ordering or number changes

class_names <- data.frame(
  class = as.character(c(1, 2, 3, 4, 5)),
  Class = as.character(c("class 1", "class 2", "class 3", "class 4", "class 5")),
  class_name = as.character(c("Inactive restless", "Extensively affected", "Minimally affected", "Focused creative", "Active restless")),
  class_abbreviation = as.character(c("IR", "EA", "MA", "FC", "AR"))
)

class_names$predclass <- class_names$class

### load phenotype data

df <- readRDS(paste0(output_dir, "/", "manic_data_derived_subset_mhq.rds"))

dim(df)

### load ICD diagnoses phenotype data from HES

icd_hes <- readRDS(paste0(output_dir, "/", "ICD_diagnoses_hesin_cleaned_data.rds"))

dim(icd_hes)

### load lca results

lca <- readRDS(file = paste0(output_dir, "/", "lca_cases_only_main_results_tmp.rds"))

### load the subset of the dataset with the ids

ids <- readRDS(file = paste0(output_dir, "/", "lca_cases_only_analytical_dataset.rds"))


################################################################################
######## TABULATE NUMBERS OF IRRITABLE AND MANIC IN EACH CLASS WITH % ##########

if (NROW(lca[[optimum_model_in_res]]$posterior) != NROW(ids)) {
  print("Number of rows in the ID file do not match the number of samples in the LCA results!~~~")
}

posterior_probs <- as.data.frame(lca[[optimum_model_in_res]]$posterior)

posterior_probs$predclass <- (lca[[optimum_model_in_res]]$predclass)

# write to tidy file for UKB returns

posterior_probs_clean <- posterior_probs

names(posterior_probs_clean) <- c("Class_1_probability","Class_2_probability","Class_3_probability","Class_4_probability","Class_5_probability","Most_likely_class")

posterior_probs_clean <- cbind(posterior_probs_clean, ids)

posterior_probs_clean <- posterior_probs_clean[,c("eid","Class_1_probability","Class_2_probability","Class_3_probability","Class_4_probability","Class_5_probability","Most_likely_class")]

write.csv(posterior_probs_clean, file=paste0(output_dir,"/LCA_classes_with_posterior_probabilities.csv"))

# bind with phenos

posterior_probs <- cbind(posterior_probs, ids)

tabyl(posterior_probs$Ever.manic.or.excitable)
tabyl(posterior_probs$Ever.extreme.irritability)

# add class names

posterior_probs$predclass <- factor(posterior_probs$predclass)

posterior_probs <- left_join(posterior_probs, class_names, by = "predclass")

### merge in ICD data from HES

dim(posterior_probs)

posterior_probs <- merge(posterior_probs, icd_hes, by = "eid", all.x = TRUE)

dim(posterior_probs)

### create new composite variable for tabulations

posterior_probs$Ever.manic.and.irritable <- NA
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable == 1] <- 1
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.extreme.irritability == 1] <- 2
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable == 1 & posterior_probs$Ever.extreme.irritability == 1] <- 3

tabyl(posterior_probs$Ever.manic.and.irritable)

###

tbx <- list()

for (i in 1:optimum_classes) {
  tb1 <- tabyl(posterior_probs$Ever.manic.and.irritable[posterior_probs$predclass == i])
  names(tb1)[1] <- as.character(posterior_probs$class_abbreviation[posterior_probs$predclass == i][1])
  tb1$categories <- c("Ever manic", "Ever irritable", "Ever manic and irritable")
  tbx[[i]] <- tb1
}

tb2 <- bind_rows(tbx, .id = c("class"))

tb2 <- tb2[c("categories", "n", "percent", "class")]

tb2$`N (percent)` <-
  paste0(tb2$n, " (", round(tb2$percent * 100, digits = 2), ")")

tb2$Response <- as.factor(as.character(gsub("\\.", " ", as.character(tb2$categories)))) # remove periods if any exist

tb2 <- left_join(tb2, class_names[, c("class", "class_abbreviation")])

tb2$class <- tb2$class_abbreviation

tb2$class_abbreviation <- NULL

### plot numbers/percentages of irritable, manic or manic AND irritable in each class

png(
  filename = paste0(output_dir, "/", "percent_manic_or_irritable_", optimum_classes, "class_barplot.png"),
  width = 20, height = 13, units = "cm", res = 300
)

ggplot(tb2, aes(Response, percent)) +
  geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Proportion of class") +
  xlab(paste0("Response (N= ", NROW(posterior_probs[which(posterior_probs[["predclass"]] != "NA"), ]), ")")) +
  theme_classic() +
  geom_text(
    aes(
      Response,
      y = percent + 0.03,
      group = class,
      label = format(percent, nsmall = 0, digits = 1, scientific = FALSE, size = 4)
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
  ) +
  labs(fill = "Class")


dev.off()


png(
  filename = paste0(output_dir, "/", "numbers_manic_or_irritable_", optimum_classes, "class_barplot.png"),
  width = 20, height = 13, units = "cm", res = 300
)

ggplot(tb2, aes(Response, n)) +
  geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Number of individuals in class") +
  xlab(paste0("Response (N= ", NROW(posterior_probs[which(posterior_probs[["predclass"]] != "NA"), ]), ")")) +
  theme_classic() +
  labs(fill = "Class")


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
    "class.csv"
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

### check
head(res_dummy)
str(res_dummy)
tabyl(res_dummy[, c("predclass2")])

posterior_probs <- cbind(posterior_probs, res_dummy)
tabyl(posterior_probs$Brief.duration.mania.or.irritability)

# summary(glm(data = posterior_probs, predclass1 ~ Brief.duration.mania.or.irritability, family = "binomial"))
# summary(glm(data = posterior_probs, predclass2 ~ Brief.duration.mania.or.irritability, family = "binomial"))
# summary(glm(data = posterior_probs, predclass3 ~ Brief.duration.mania.or.irritability, family = "binomial"))
# summary(glm(data = posterior_probs, predclass4 ~ Brief.duration.mania.or.irritability, family = "binomial"))
# summary(glm(data = posterior_probs, predclass5 ~ Brief.duration.mania.or.irritability, family = "binomial"))
### logistic regression with dummy coded contrasts of classes

### reshape data to long so that participants have as many rows as classes and probability for each row
posterior_long <- reshape(posterior_probs,
  idvar = "eid", direction = "long",
  varying = 2:6, v.names = (f <- c("probabilities")), timevar = "predclass"
)

### check reshape to long for a random participant
posterior_long$probabilities[posterior_long$eid == sample(posterior_long$eid, 1)] # should have 5 probabilities sum to 1
if (round(sum(posterior_long$probabilities[posterior_long$eid == sample(posterior_long$eid, 1)]), digits = 3) != 1) print("Error with reshape!")

posterior_long$predclass[posterior_long$eid == sample(posterior_long$eid, 1)] # should have values for X number of classes

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
    
    palette_choice4 <- palette_choice_ma
  } else {
    print("This is the main analysis (reference level '1')")
    palette_choice4 <- palette_choice_ir
    
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

  # join names
  duration_naive_res_df <- left_join(duration_naive_res_df, class_names, by = c("y.level" = "class"))
  duration_weighted_res_df <- left_join(duration_weighted_res_df, class_names, by = c("y.level" = "class"))

  refclass <- as.character(class_names$class_abbreviation)[which(class_names$class == levels(posterior_probs$predclass)[1])]
  refclasslong <- as.character(class_names$class_abbreviation)[which(class_names$class == levels(posterior_long$predclass)[1])]

  ### create a 'comparison' column of statistical test for plot
  duration_naive_res_df$comparison <- paste0(refclass, " v ", duration_naive_res_df$class_abbreviation)
  duration_weighted_res_df$comparison <- paste0(refclasslong, " v ", duration_naive_res_df$class_abbreviation)

  ### remove periods if any exist from names
  duration_naive_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(duration_naive_res_df$exposure)))
  duration_weighted_res_df$exposure <- as.factor(gsub("\\.", " ", as.character(duration_weighted_res_df$exposure)))

  write.csv(duration_weighted_res_df, paste0(output_dir, "/", "multinomial_duration_weighted_res_", reference_i, ".csv"), row.names = F, quote = F)

  ### plot multinomial results
  png(
    filename = paste0(output_dir, "/", "multinomial_naive_duration_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
    width = 20, height = 8, units = "cm", res = 300
  )

  print(ggplot(duration_naive_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_colour_manual(values = palette_choice4) +
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
    filename = paste0(output_dir, "/", "multinomial_weighted_duration_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
    width = 20, height = 8, units = "cm", res = 300
  )

  print(ggplot(duration_weighted_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_colour_manual(values = palette_choice4) +
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

  ### set "Don't know" to NA

  posterior_probs$Duration.mania.or.irritability[posterior_probs$Duration.mania.or.irritability == 4] <- NA
  posterior_long$Duration.mania.or.irritability[posterior_long$Duration.mania.or.irritability == 4] <- NA

  for (i in 1:optimum_classes) {
    tb1 <- tabyl(posterior_probs$Duration.mania.or.irritability[posterior_probs$predclass == i])

    names(tb1)[1] <- paste0("class_", i)

    tb1$categories <- c("1 - Brief", "2 - Moderate", "3 - Extended", "Missing - NA")

    tb_duration[[i]] <- tb1
  }

  tb_duration2 <- bind_rows(tb_duration, .id = c("class"))

  tb_duration2 <- left_join(tb_duration2, class_names[, c("class", "class_abbreviation")])

  tb_duration2$class <- tb_duration2$class_abbreviation

  tb_duration2$class_abbreviation <- NULL

  tb_duration2 <- tb_duration2[c("categories", "n", "percent", "valid_percent", "class")]

  tb_duration2$`N (percent)` <-
    paste0(tb_duration2$n, " (", round(tb_duration2$percent * 100, digits = 2), ")")

  tb_duration2$`N (valid percent)` <-
    paste0(tb_duration2$n, " (", round(tb_duration2$valid_percent * 100, digits = 2), ")")

  tb_duration2$Duration <- as.factor(gsub("\\.", " ", as.character(tb_duration2$categories))) # remove periods if any exist

  write.table(tb_duration2[c("categories", "n", "percent", "class", "N (percent)")],
    file = paste0(output_dir, "/", "duration_", optimum_classes, "class_tabulations.csv"),
    row.names = FALSE, quote = FALSE, sep = ","
  )

  ### tabulate number of missing rows (hence individuals with NA for duration) and subset

  missing_duration <- NROW(posterior_probs$Duration.mania.or.irritability) - NROW(tb_duration2[which(posterior_probs$Duration.mania.or.irritability != "NA"), ])

  ### drop NAs

  tb_duration2_valid <- tb_duration2[which(tb_duration2$categories != "Missing - NA"), ]

  ### plot bar plots of comparisons by most probable class

  png(
    filename = paste0(output_dir, "/", "valid_percent_duration_", optimum_classes, "class_barplot.png"),
    width = 20, height = 13, units = "cm", res = 300
  )

  print(ggplot(tb_duration2_valid, aes(Duration, valid_percent)) +
    geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Proportion of class") +
    theme_classic() +
    xlab(paste0("Response (N= ", NROW(posterior_probs[which(posterior_probs[["Duration.mania.or.irritability"]] != "NA"), ]), ")")) +
    geom_text(aes(Duration,
      y = valid_percent + 0.03,
      group = class,
      label = format(valid_percent, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
    ) +
    labs(fill = "Class"))

  dev.off()

  ###

  png(
    filename = paste0(output_dir, "/", "percent_duration_NAs_", optimum_classes, "class_barplot.png"),
    width = 20, height = 13, units = "cm", res = 300
  )

  print(ggplot(tb_duration2, aes(Duration, percent)) +
    geom_bar(aes(fill = class), position = "dodge", stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Proportion of class") +
    theme_classic() +
    xlab(paste0("Response (N= ", NROW(posterior_probs[which(posterior_probs[["Duration.mania.or.irritability"]] != "NA"), ]), ")")) +
    geom_text(aes(Duration,
      y = percent + 0.03,
      group = class,
      label = format(percent, nsmall = 0, digits = 2, scientific = FALSE, size = 4)
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
    ) +
    labs(fill = "Class"))

  dev.off()
}


##############################################################################
######## ASSOCIATIONS WITH DISRUPTIVENESS OF MANIC SYMPTOMATOLOGY ############

# assumes there is a var 'predclass' and var 'probabilities' in df for regressions

### rename levels for plot before applying function
posterior_probs$Problematic.mania.or.irritability[posterior_probs$Problematic.mania.or.irritability == 0] <- "Not disruptive"
posterior_probs$Problematic.mania.or.irritability[posterior_probs$Problematic.mania.or.irritability == 1] <- "Disruptive"

tabyl(posterior_probs$Problematic.mania.or.irritability)

analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Problematic.mania.or.irritability",
  traitname = "disruptiveness",
  missingnessplot = TRUE
)

####################################################################################
############### ASSOCIATIONS WITH SOCIODEMOGRAPHIC FACTORS #########################

### create subset of data for sociodemographic analyses

sociodemographic_vars <- c(
  "Sex", "Townsend.deprivation", "Education.level",
  "Smoking", "Alcohol.intake.frequency"
)

sociodemographic_vars_auxillary <- c("Age")

### create subset to make comparions between LCA subset and whole MHQ

posterior_probs_sociodem_description <- (posterior_probs[, c(
  "eid",
  sociodemographic_vars, sociodemographic_vars_auxillary, "predclass"
)])

tabulated_sociodem <- list()

tabulated_sociodem[["continuous"]] <- bind_rows(
  c(
    summary(posterior_probs_sociodem_description$Townsend.deprivation),
    sd(na.omit(posterior_probs_sociodem_description$Townsend.deprivation))
  ),
  c(
    summary(posterior_probs_sociodem_description$Age),
    sd(na.omit(posterior_probs_sociodem_description$Age))
  )
)

tabulated_sociodem[["categorical"]] <- bind_rows(
  tabyl(posterior_probs_sociodem_description$Sex),
  tabyl(posterior_probs_sociodem_description$Education.level),
  tabyl(posterior_probs_sociodem_description$Alcohol.intake.frequency),
  tabyl(posterior_probs_sociodem_description$Smoking)
)

### full MHQ for description

all_mhq_description <- (df[, c(
  "eid",
  sociodemographic_vars, sociodemographic_vars_auxillary
)])

tabulated_sociodem_mhq <- list()

tabulated_sociodem_mhq[["continuous"]] <- bind_rows(
  c(
    summary(all_mhq_description$Townsend.deprivation),
    sd(na.omit(all_mhq_description$Townsend.deprivation))
  ),
  c(
    summary(all_mhq_description$Age),
    sd(na.omit(all_mhq_description$Age))
  )
)

tabulated_sociodem_mhq[["categorical"]] <- bind_rows(
  tabyl(all_mhq_description$Sex),
  tabyl(all_mhq_description$Education.level),
  tabyl(all_mhq_description$Alcohol.intake.frequency),
  tabyl(all_mhq_description$Smoking)
)

### form descriptive stats in to table

tabulated_sociodem_mhq[["continuous"]] <- as.data.frame(tabulated_sociodem_mhq[["continuous"]])

tabulated_sociodem_mhq[["continuous"]][is.na(tabulated_sociodem_mhq[["continuous"]])] <- " "

names(tabulated_sociodem_mhq[["continuous"]]) <- c(names(tabulated_sociodem_mhq[["continuous"]][1:7]), "sd1")

tabulated_sociodem_mhq[["continuous"]]$sd <- as.numeric(paste0(
  tabulated_sociodem_mhq[["continuous"]]["sd1"]$sd1
))

tabulated_sociodem_mhq[["continuous"]] <- tabulated_sociodem_mhq[["continuous"]][, c("Mean", "NA's", "sd")]

tabulated_sociodem_mhq[["continuous"]]$Mean <- paste0(round(tabulated_sociodem_mhq[["continuous"]]$Mean, digits = 1), "(", round(tabulated_sociodem_mhq[["continuous"]]$sd, digits = 1), ")")

#

tabulated_sociodem[["continuous"]] <- as.data.frame(tabulated_sociodem[["continuous"]])

tabulated_sociodem[["continuous"]][is.na(tabulated_sociodem[["continuous"]])] <- " "

names(tabulated_sociodem[["continuous"]]) <- c(names(tabulated_sociodem[["continuous"]][1:7]), "sd1")

tabulated_sociodem[["continuous"]]$sd <- as.numeric(paste0(
  (tabulated_sociodem[["continuous"]]["sd1"]$sd1)
))

tabulated_sociodem[["continuous"]] <- tabulated_sociodem[["continuous"]][, c("Mean", "NA's", "sd")]

tabulated_sociodem[["continuous"]]$Mean <- paste0(round(tabulated_sociodem[["continuous"]]$Mean, digits = 1), "(", round(tabulated_sociodem[["continuous"]]$sd, digits = 1), ")")

write.table(tabulated_sociodem_mhq[["continuous"]], file = paste0(output_dir, "/summaries_table1_mhq_subset.tab"), sep = "\t", quote = F, row.names = F)

write.table(tabulated_sociodem[["continuous"]], file = paste0(output_dir, "/summaries_table1_lca_subset.tab"), sep = "\t", quote = F, row.names = F)

### categorical

tabulated_sociodem_cat <- setDF(tabulated_sociodem[["categorical"]])

tabulated_sociodem_cat$Measure <- tabulated_sociodem_cat$`posterior_probs_sociodem_description$Sex`

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Education.level"]]))] <- tabulated_sociodem_cat[["posterior_probs_sociodem_description$Education.level"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Education.level"]]))]

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Smoking"]]))] <- tabulated_sociodem_cat[["posterior_probs_sociodem_description$Smoking"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Smoking"]]))]

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Alcohol.intake.frequency"]]))] <- tabulated_sociodem_cat[["posterior_probs_sociodem_description$Alcohol.intake.frequency"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["posterior_probs_sociodem_description$Alcohol.intake.frequency"]]))]

tabulated_sociodem_tidy <- tabulated_sociodem_cat[, c("Measure", "n", "percent", "valid_percent")]

tabulated_sociodem_tidy$Measure <- as.factor(gsub("\\.", " ", as.character(tabulated_sociodem_tidy$Measure)))

tabulated_sociodem_tidy$percent <- signif((tabulated_sociodem_tidy$percent * 100), digits = 2)

### mhq

tabulated_sociodem_cat <- setDF(tabulated_sociodem_mhq[["categorical"]])

tabulated_sociodem_cat$Measure <- tabulated_sociodem_cat$`all_mhq_description$Sex`

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Education.level"]]))] <- tabulated_sociodem_cat[["all_mhq_description$Education.level"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Education.level"]]))]

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Smoking"]]))] <- tabulated_sociodem_cat[["all_mhq_description$Smoking"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Smoking"]]))]

tabulated_sociodem_cat$Measure[which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Alcohol.intake.frequency"]]))] <- tabulated_sociodem_cat[["all_mhq_description$Alcohol.intake.frequency"]][which(is.na(tabulated_sociodem_cat$Measure) &
  !is.na(tabulated_sociodem_cat[["all_mhq_description$Alcohol.intake.frequency"]]))]

tabulated_sociodem_tidy_mhq <- tabulated_sociodem_cat[, c("Measure", "n", "percent", "valid_percent")]

tabulated_sociodem_tidy_mhq$Measure <- as.factor(gsub("\\.", " ", as.character(tabulated_sociodem_tidy_mhq$Measure)))

tabulated_sociodem_tidy_mhq$percent <- signif((tabulated_sociodem_tidy_mhq$percent * 100), digits = 2)

### write

write.table(tabulated_sociodem_tidy_mhq, file = paste0(output_dir, "/tabulations_table1_mhq_subset.tab"), sep = "\t", quote = F, row.names = F)

write.table(tabulated_sociodem_tidy, file = paste0(output_dir, "/tabulations_table1_lca_subset.tab"), sep = "\t", quote = F, row.names = F)

### subset for analysis with NAs omitted

posterior_probs_sociodem <- na.omit(posterior_probs[, c(
  "eid",
  sociodemographic_vars, "predclass"
)])

dim(posterior_probs_sociodem)

posterior_long_sociodem <- na.omit(posterior_long[, c("eid", sociodemographic_vars, "predclass", "probabilities")])

dim(posterior_long_sociodem)

if ((NROW(posterior_long_sociodem) / optimum_classes) != NROW(posterior_probs_sociodem)) print("Error: Long and wide datasets have different numbers of participants!")

### first rank normalise continuous townsend deprivation

class(posterior_long_sociodem$Townsend.deprivation)
summary(posterior_long_sociodem$Townsend.deprivation)

posterior_long_sociodem$Townsend.deprivation.rank.normalised <- RankNorm(posterior_long_sociodem$Townsend.deprivation)
posterior_probs_sociodem$Townsend.deprivation.rank.normalised <- RankNorm(posterior_probs_sociodem$Townsend.deprivation)

### check rank normal var

hist(posterior_probs_sociodem$Townsend.deprivation.rank.normalised)
summary(posterior_long_sociodem$Townsend.deprivation.rank.normalised)

sd(posterior_probs_sociodem$Townsend.deprivation.rank.normalised)
mean(posterior_probs_sociodem$Townsend.deprivation.rank.normalised)

### tabulate/describe sample of sociodemographics

tabulated_sociodem <- list()

tabulated_sociodem[["continuous"]] <- c(summary(posterior_probs_sociodem$Townsend.deprivation.rank.normalised), sd(posterior_probs_sociodem$Townsend.deprivation.rank.normalised))

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
  traitname = "sex"
)

### apply multinomial analysis function to townsend deprivation

summary(posterior_probs_sociodem$Townsend.deprivation.rank.normalised)
class(posterior_probs_sociodem$Townsend.deprivation.rank.normalised)

Townsend.deprivation.rank.normalised_res <- analyse_continuous_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Townsend.deprivation.rank.normalised",
  traitname = "TDI"
)

### apply multinomial analysis function to education level (categorical)

tabyl(posterior_probs_sociodem$Education.level)
class(posterior_probs_sociodem$Education.level)

Education.level_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Education.level",
  traitname = "education"
)

### apply multinomial analysis function to smoking (categorical)

tabyl(posterior_probs_sociodem$Smoking)
class(posterior_probs_sociodem$Smoking)

Smoking_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Smoking",
  traitname = "smoking"
)

### apply multinomial analysis function to alcohol consumption (categorical)

tabyl(posterior_probs_sociodem$Alcohol.intake.frequency)
class(posterior_probs_sociodem$Alcohol.intake.frequency)

Alcohol.intake.frequency_res <- analyse_categorical_multinomial(
  datasetwide = posterior_probs_sociodem,
  datasetlong = posterior_long_sociodem,
  varname = "Alcohol.intake.frequency",
  traitname = "alcohol_frequency"
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
  prs_dat <- setDF(fread(paste0(prs_dir, prs_filenames[i_prs], "_header.all.score")))

  prs_raw[[i_prs]] <- prs_dat[c("IID", paste0(prs_filenames[i_prs], "_", pvalue_threshold))]
}

### ensure same IDs across PRS files (no withdrawn IDs in some and not others)

id_freqs <- data.frame(table(unlist(lapply(prs_raw, "[[", 1))))

exclude_ids <- as.character(id_freqs$Var1[which(as.numeric(as.character(id_freqs$Freq)) < NROW(prs_raw))])

### remove any IDs that do not fit

prs_raw <- lapply(prs_raw, function(x) x[which(as.character(x[, 1]) %notin% exclude_ids), ])

#prs_df <- rbindlist( prs_raw, idcol = TRUE )
prs_df <- bind_cols(prs_raw)

prs_df$IID <- prs_df$IID...1

### load PCs and covs

ukb_covs <-
  setDF(fread(
    paste0(new_2019_qc_dir, "ukb18177_glanville_covariates.txt")
  ))

cov_names <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

### merge PRS and covariates

for_pcs_regression <- merge(prs_df, ukb_covs, by = c("IID"))

### regress out PCs (top 6)

prs.res <- list()

for (i_trait in c(paste0(prs_filenames, "_", pvalue_threshold))) {
  prs.lm <- lm(paste0(i_trait, " ~ ", paste0(cov_names, collapse = " + ")),
    data = for_pcs_regression
  )

  prs.res[[i_trait]] <- resid(prs.lm)
}

prs.res_df <- as.data.frame(prs.res)

### standardise PRS residuals

prs.res_df[names(prs.res_df)] <-
  lapply(prs.res_df[names(prs.res_df)], function(x) {
    c(scale(x))
  })

### add ID column back in

prs.res_df$IID <- for_pcs_regression$IID

### check conforms to normal distribution with SD=1

hist(prs.res_df[[paste0("SCHI02_", pvalue_threshold)]], breaks = 100)

sd(prs.res_df[[paste0("SCHI02_", pvalue_threshold)]])

### merge with posterior probs and class memembership

posterior_probs_prs <- merge(posterior_probs, prs.res_df, by.x = "eid", by.y = "IID")

### find which samples weren't found in PRS

samples_without_prs <- NROW(posterior_probs) - NROW(posterior_probs_prs)

print(paste0("There are ", samples_without_prs, " number of samples in latent class analysis without PRS"))

ids_without_prs <- posterior_probs$eid[posterior_probs$eid %notin% posterior_probs_prs$eid]

print(NROW(ids_wthout_prs))

### reshape data to long so that participants have as many rows as classes and a probability for each row/class

posterior_prs_long <- reshape(posterior_probs_prs,
  idvar = "eid", direction = "long",
  varying = 2:6, v.names = (f <- c("probabilities")), timevar = "predclass"
)

### check reshape to long correct by picking a random participant ID

posterior_prs_long$probabilities[posterior_prs_long$eid == sample(posterior_prs_long$eid, 1)] # should have 5 probabilities sum to 1

if (round(sum(posterior_prs_long$probabilities[posterior_prs_long$eid == sample(posterior_prs_long$eid, 1)]), digits = 6) != 1) print("Error with reshape!")

### check that values for X number of classes present (ie 5 if optimum model is 5 class)

posterior_prs_long$predclass[posterior_prs_long$eid == sample(posterior_prs_long$eid, 1)]

### save stata dataset for cross-checking results using mlogit in STATA
# save.dta13(posterior_prs_long[, c("eid", "probabilities", "BIPO01_1", "predclass")],
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

  # join names
  prs_naive_res_df <- left_join(prs_naive_res_df, class_names, by = c("y.level" = "class"))
  prs_weighted_res_df <- left_join(prs_weighted_res_df, class_names, by = c("y.level" = "class"))

  refclass <- as.character(class_names$class_abbreviation)[which(class_names$class == levels(posterior_probs_prs[["predclass"]])[1])]
  refclasslong <- as.character(class_names$class_abbreviation)[which(class_names$class == levels(posterior_prs_long[["predclass"]])[1])]

  ### create a 'comparison' column of statistical test for plot
  prs_naive_res_df$comparison <- paste0(refclass, " v ", prs_naive_res_df$class_abbreviation)
  prs_weighted_res_df$comparison <- paste0(refclasslong, " v ", prs_weighted_res_df$class_abbreviation)

  print(prs_weighted_res_df)

  ### remove periods if any exist from names
  prs_naive_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_naive_res_df$exposure)))

  prs_weighted_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_weighted_res_df$exposure)))

  write.csv(prs_weighted_res_df, paste0(output_dir, "/", "multinomial_prs_weighted_res_", reference_i, ".csv"), col.names = T, row.names = F, quote = F)

  plot_prs_ls <- list()

  ### plot multinomial results
  png(
    filename = paste0(output_dir, "/", "multinomial_naive_prs_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
    width = 20, height = 10, units = "cm", res = 300
  )

  print(ggplot(prs_naive_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    facet_wrap(~exposure)) +
    scale_colour_manual(values = palette_choice4) +
    ylab("Risk ratio (per SD increase in PRS)") +
    theme_classic() +
    theme(
      strip.text.x = element_text(hjust = 0.01),
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(colour = "") +
    guides(colour = FALSE)

  dev.off()

  for (k in unique(prs_naive_res_df$exposure)) {
    prs_naive_res_df1 <- prs_naive_res_df[prs_naive_res_df$exposure == k, ]

    plot_prs_ls[[paste0("naive_", k)]] <- ggplot(prs_naive_res_df1, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      scale_colour_manual(values = palette_choice4) +
      ylab("Risk ratio (per SD increase in PRS)") +
      theme_classic() +
      theme(legend.position = "none") +
      labs(colour = "", title = k) +
      guides(colour = FALSE)
  }

  ###

  png(
    filename = paste0(output_dir, "/", "multinomial_weighted_prs_", reference_i, "_", optimum_classes, "class_scatterplot.png"),
    width = 20, height = 10, units = "cm", res = 300
  )

  print(ggplot(prs_weighted_res_df, aes(comparison, estimate, colour = comparison)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    facet_wrap(~exposure)) +
    scale_colour_manual(values = palette_choice4) +
    theme_classic() +
    ylab("Risk ratio (per SD increase in PRS)") +
    theme(
      strip.text.x = element_text(hjust = 0.01),
      panel.background = element_rect(fill = "white", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "white")
    ) +
    labs(colour = " ") +
    guides(colour = FALSE)

  dev.off()


  for (k in unique(prs_weighted_res_df$exposure)) {
    prs_weighted_res_df1 <- prs_weighted_res_df[prs_weighted_res_df$exposure == k, ]

    plot_prs_ls[[paste0("weighted_", k)]] <- ggplot(prs_weighted_res_df1, aes(comparison, estimate, colour = comparison)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, color = "coral", alpha = 0.4) +
      geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
      xlab("Class comparisons") +
      scale_colour_manual(values = palette_choice4) +
      ylab("Risk ratio (per SD increase in PRS)") +
      theme_classic() +
      theme(legend.position = "none") +
      labs(colour = "", title = k) +
      guides(colour = FALSE)
  }
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
    filename = paste0(output_dir, "/", "prs_", prs_trait_i_p, "_", optimum_classes, "class_densityplot.png"),
    width = 20, height = 15, units = "cm", res = 300
  )

  print(gg1)

  dev.off()
}


png(
  filename = paste0(output_dir, "/", "multinomial_weighted_relevel", "_all_prs.png"),
  width = 25, height = 12, units = "cm", res = 300
)

annotate_figure(grid.arrange(
  (plot_prs_ls[["weighted_ADHD PRS"]] + rremove("ylab") + rremove("xlab")),
  (plot_prs_ls[["weighted_Anxiety PRS"]] + rremove("ylab") + rremove("xlab")),
  (plot_prs_ls[["weighted_ASD PRS"]] + rremove("ylab") + rremove("xlab")),
  (plot_prs_ls[["weighted_Bipolar disorder PRS"]] + rremove("ylab") + rremove("xlab")),
  (plot_prs_ls[["weighted_Depression PRS"]] + rremove("ylab") + rremove("xlab")),
  (plot_prs_ls[["weighted_Schizophrenia PRS"]] + rremove("ylab") + rremove("xlab")),
  ncol = 3
),
bottom = text_grob("Class comparison", color = "black", size = 12),
left = text_grob("Estimate - RR", color = "black", rot = 90, size = 12)
)

dev.off()

##############################################################################################
################## ASSOCIATIONS WITH DIAGNOSIS OF PSYCHIATRIC DISORDERS ######################

diagnoses_to_check <- c(
  "Self.reported.ADHD", "Self.reported.GAD.and.others", "Self.reported.ASD",
  "Self.reported.Depression", "Self.reported.Mania.or.bipolar.disorder", "Self.reported.PsychosisAny"
)

diagnosis_personality <- c("Self.reported.Personality.disorder", "Neuroticism.score")

for (i in diagnoses_to_check) {
  print(i)
  print(tabyl(df[[i]]))
}

### convert to factors with correct leveling for Yes/No response

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

Self.reported.Schizophrenia.or.psychosis_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.Schizophrenia.or.psychosis",
  traitname = "schizophrenia_psychosis"
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
  traitname = "mania_bipolar"
)

posterior_probs$Self.reported.ADHD <- factor(
  posterior_probs$Self.reported.ADHD,
  levels = sort(unique(posterior_probs$Self.reported.ADHD)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.ADHD <- factor(
  posterior_long$Self.reported.ADHD,
  levels = sort(unique(posterior_long$Self.reported.ADHD)),
  labels = c("No", "Yes")
)

Self.reported.ADHD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.ADHD",
  traitname = "adhd"
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
  traitname = "depression"
)

posterior_probs$Self.reported.ASD <- factor(
  posterior_probs$Self.reported.ASD,
  levels = sort(unique(posterior_probs$Self.reported.ASD)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.ASD <- factor(
  posterior_long$Self.reported.ASD,
  levels = sort(unique(posterior_long$Self.reported.ASD)),
  labels = c("No", "Yes")
)

Self.reported.ASD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.ASD",
  traitname = "asd"
)

posterior_probs$Self.reported.GAD <- factor(
  posterior_probs$Self.reported.GAD.and.others,
  levels = sort(unique(posterior_probs$Self.reported.GAD.and.others)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.GAD <- factor(
  posterior_long$Self.reported.GAD.and.others,
  levels = sort(unique(posterior_long$Self.reported.GAD.and.others)),
  labels = c("No", "Yes")
)

Self.reported.GAD.and.others_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.GAD",
  traitname = "gad"
)

### grid plots of self-reported diagnoses

png(
  filename = paste0(output_dir, "/", "multinomial_relevel_log", "_all_SR_diagnoses.png"),
  width = 25, height = 12, units = "cm", res = 300
)

annotate_figure(grid.arrange(
  (Self.reported.ADHD_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.GAD.and.others_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.ASD_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Mania.or.bipolar.disorder_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Depression_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Schizophrenia.or.psychosis_res[[1]]$weightedlog_relevel + rremove("ylab") + rremove("xlab")),
  ncol = 3
),
bottom = text_grob("Class comparison", color = "black", size = 11),
left = text_grob("Estimate - log(RR)", color = "black", rot = 90, size = 11)
)

dev.off()

###

png(
  filename = paste0(output_dir, "/", "multinomial_main_log", "_all_SR_diagnoses.png"),
  width = 25, height = 12, units = "cm", res = 300
)

annotate_figure(grid.arrange(
  (Self.reported.ADHD_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.GAD.and.others_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.ASD_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Mania.or.bipolar.disorder_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Depression_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Schizophrenia.or.psychosis_res[[1]]$weightedlog_main + rremove("ylab") + rremove("xlab")),
  ncol = 3
),
bottom = text_grob("Class comparison", color = "black", size = 11),
left = text_grob("Estimate- log(RR)", color = "black", rot = 90, size = 11)
)

dev.off()

###

png(
  filename = paste0(output_dir, "/", "multinomial_relevel", "_all_SR_diagnoses.png"),
  width = 25, height = 12, units = "cm", res = 300
)

annotate_figure(grid.arrange(
  (Self.reported.ADHD_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.GAD.and.others_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.ASD_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Mania.or.bipolar.disorder_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Depression_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  (Self.reported.Schizophrenia.or.psychosis_res[[1]]$weighted_relevel + rremove("ylab") + rremove("xlab")),
  ncol = 3
),
bottom = text_grob("Class comparison", color = "black", size = 11),
left = text_grob("Estimate - RR", color = "black", rot = 90, size = 11)
)

# hjust = 1, x = 1, face = "italic",

dev.off()

###

png(
  filename = paste0(output_dir, "/", "multinomial_main", "_all_SR_diagnoses.png"),
  width = 25, height = 12, units = "cm", res = 300
)

annotate_figure(grid.arrange(
  (Self.reported.ADHD_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.GAD.and.others_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.ASD_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Mania.or.bipolar.disorder_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Depression_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  (Self.reported.Schizophrenia.or.psychosis_res[[1]]$weighted_main + rremove("ylab") + rremove("xlab")),
  ncol = 3
),
bottom = text_grob("Class comparison", color = "black", size = 11),
left = text_grob("Estimate - RR", color = "black", rot = 90, size = 11)
)

dev.off()

### prepare data for upset plot creation (external script)

diagnoses <- c(
  "Self.reported.ADHD", "Self.reported.GAD.and.others", "Self.reported.ASD",
  "Self.reported.Depression", "Self.reported.Mania.or.bipolar.disorder", "Self.reported.Schizophrenia.or.psychosis"
)

working_df <- posterior_probs[, c(diagnoses, "predclass")]

for (self_disorders in diagnoses) {
  working_df[[self_disorders]] <- as.character(working_df[[self_disorders]])
  working_df[[self_disorders]][working_df[[self_disorders]] == "Yes"] <- 1
  working_df[[self_disorders]][working_df[[self_disorders]] == "No"] <- 0
}

working_df$predclass <- as.character(working_df$predclass)

working_df[] <- lapply(working_df, as.numeric)

names(working_df) <- gsub("Self.reported.", "", names(working_df))
names(working_df) <- gsub("\\.", " ", names(working_df))

working_df$predclass <- as.character(working_df$predclass)

### save data to plot upsets as job on cluster

saveRDS(working_df, file = paste0(output_dir, "/working_df_diagnoses_for_plots_tmp.rds"))


### Personality traits

posterior_probs$Self.reported.Personality.disorder <- factor(
  posterior_probs$Self.reported.Personality.disorder,
  levels = sort(unique(posterior_probs$Self.reported.Personality.disorder)),
  labels = c("No", "Yes")
)

posterior_long$Self.reported.Personality.disorder <- factor(
  posterior_long$Self.reported.Personality.disorder,
  levels = sort(unique(posterior_long$Self.reported.Personality.disorder)),
  labels = c("No", "Yes")
)

Self.reported.Personality.disorder_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Self.reported.Personality.disorder",
  traitname = "personality_disorder"
)

### rescale neuroticism score to mean zero and sd of 1

posterior_probs$Neuroticism.score <- scale(posterior_probs$Neuroticism.score)

posterior_long$Neuroticism.score <- scale(posterior_long$Neuroticism.score)

Neuroticism.score_res <- analyse_continuous_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Neuroticism.score",
  traitname = "neuroticism"
)

### ICD10 diagnosis

### tabulate overlap with self-reported diagnosis of counterpart disorder

table(posterior_probs$Depression.ICD, posterior_probs$Self.reported.Depression)

table(posterior_probs$Schizophrenia.or.psychotic.disorder.ICD, posterior_probs$Self.reported.Schizophrenia)

table(posterior_probs$Bipolar.or.mood.disorder.ICD, posterior_probs$Self.reported.Mania.or.bipolar.disorder)

posterior_probs$Dementia.ICD <- factor(
  posterior_probs$Dementia.ICD,
  levels = sort(unique(posterior_probs$Dementia.ICD)),
  labels = c("No", "Yes")
)

posterior_long$Dementia.ICD <- factor(
  posterior_long$Dementia.ICD,
  levels = sort(unique(posterior_long$Dementia.ICD)),
  labels = c("No", "Yes")
)

Dementia.ICD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Dementia.ICD",
  traitname = "dementia_ICD"
)

posterior_probs$Depression.ICD <- factor(
  posterior_probs$Depression.ICD,
  levels = sort(unique(posterior_probs$Depression.ICD)),
  labels = c("No", "Yes")
)

posterior_long$Depression.ICD <- factor(
  posterior_long$Depression.ICD,
  levels = sort(unique(posterior_long$Depression.ICD)),
  labels = c("No", "Yes")
)

Depression.ICD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Depression.ICD",
  traitname = "depression_ICD"
)

posterior_probs$Schizophrenia.or.psychotic.disorder.ICD <- factor(
  posterior_probs$Schizophrenia.or.psychotic.disorder.ICD,
  levels = sort(unique(posterior_probs$Schizophrenia.or.psychotic.disorder.ICD)),
  labels = c("No", "Yes")
)

posterior_long$Schizophrenia.or.psychotic.disorder.ICD <- factor(
  posterior_long$Schizophrenia.or.psychotic.disorder.ICD,
  levels = sort(unique(posterior_long$Schizophrenia.or.psychotic.disorder.ICD)),
  labels = c("No", "Yes")
)

Schizophrenia.or.psychotic.disorder.ICD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Schizophrenia.or.psychotic.disorder.ICD",
  traitname = "schizophrenia_ICD"
)

posterior_probs$Mania.or.bipolar.disorder.ICD <- factor(
  posterior_probs$Bipolar.or.mood.disorder.ICD,
  levels = sort(unique(posterior_probs$Bipolar.or.mood.disorder.ICD)),
  labels = c("No", "Yes")
)

posterior_long$Mania.or.bipolar.disorder.ICD <- factor(
  posterior_long$Bipolar.or.mood.disorder.ICD,
  levels = sort(unique(posterior_long$Bipolar.or.mood.disorder.ICD)),
  labels = c("No", "Yes")
)

Mania.or.bipolar.disorder.ICD_res <- analyse_binary_multinomial(
  datasetwide = posterior_probs,
  datasetlong = posterior_long,
  varname = "Mania.or.bipolar.disorder.ICD",
  traitname = "bipolar_ICD"
)

### grid plots of ICD diagnoses

png(
  filename = paste0(output_dir, "/", "multinomial_relevel", "_all_ICD.png"),
  width = 30, height = 18, units = "cm", res = 300
)

grid.arrange(
  Depression.ICD_res[[1]]$weighted_relevel,
  Schizophrenia.or.psychotic.disorder.ICD_res[[1]]$weighted_relevel,
  Mania.or.bipolar.disorder.ICD_res[[1]]$weighted_relevel,
  Dementia.ICD_res[[1]]$weighted_relevel,
  
  ncol = 2
)

dev.off()

###

png(
  filename = paste0(output_dir, "/", "multinomial_main", "_all_ICD.png"),
  width = 30, height = 18, units = "cm", res = 300
)

grid.arrange(
  Depression.ICD_res[[1]]$weighted_main,
  Schizophrenia.or.psychotic.disorder.ICD_res[[1]]$weighted_main,
  Mania.or.bipolar.disorder.ICD_res[[1]]$weighted_main,
  Dementia.ICD_res[[1]]$weighted_main,
  ncol = 2
)

dev.off()


## create Figure 4

comb2pngs <- function(imgs, bottom_text = NULL) {
  img1 <- grid::rasterGrob(as.raster(readPNG(imgs[1])),
                           interpolate = FALSE
  )
  img2 <- grid::rasterGrob(as.raster(readPNG(imgs[2])),
                           interpolate = FALSE
  )
  # grid.arrange(img1, img2, ncol = 1, bottom = bottom_text)
  plot_grid(img1, img2, labels = c("A", "B"), ncol = 1, label_size = 15, hjust = 0, scale = 1.01)
}

png1_dest <- paste0(output_dir, "/", "multinomial_relevel_log", "_all_SR_diagnoses.png")
png2_dest <- paste0(output_dir, "/", "multinomial_weighted_relevel", "_all_prs.png")

png(
  filename = paste0(output_dir, "/", "figure_4_panel_plot_diagnoses_PRS.png"),
  width = 24, height = 28, units = "cm", res = 300
)

comb2pngs(c(png1_dest, png2_dest))

dev.off()

