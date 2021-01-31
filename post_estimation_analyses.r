### Post LCA estimation analyses
### Associations of latent classes with other variables
### Ryan Arathimos
### 13/10/2020

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

set.seed(123456)

### SET UP DIRECTORIES
data_dir <- "~/brc_scratch/data/genotype_lca"
output_dir <- "~/brc_scratch/output/genotype_lca"
scripts_dir <- "~/brc_scratch/scripts/genotype_lca"

project_dir <- "/scratch/datasets/ukbiobank/ukb18177_glanville"
qc_dir <-  "/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/"
dir_prs <-  "/mnt/lustre/groups/ukbiobank/sumstats/PRS/ukb18177_glanville/PRS_noheader/"

dirs <- read.table('~/brc_scratch/data/directories_list.txt')
dir_prs <- as.character(dirs[1,])
new_2019_qc_dir <- as.character(dirs[3,])
new_2019_qc_dir <- "/mnt/lustre/groups/ukbiobank/ukb18177_glanville/2019_new_qc/"

palette_choice <- "Set1" #rcolorbrewer palette "Dark2"
`%notin%` <- Negate(`%in%`)

### optimum model - i.e. number of classes which optimally fit the data

optimum_classes <- 5
optimum_model_in_res <- optimum_classes - 1 #the object number in the list is numbered (N classes -1)

### names classes of optimum model

class_names <- factor(c("Focused creative","Active restless","Unaffected","Invariable affected","Sedentary restless") )      
# undiversified, invariable, functional, perception invariable undifferentiated

### load phenotype data

df <- readRDS(paste0(output_dir, "/", "manic_data_derived_subset_mhq.rds"))

dim(df)

### load lca results

lca <- readRDS(file = paste0(output_dir,"/manic_lca/","lca_cases_only_main_results_tmp.rds"))

### load the subset of the dataset with the ids

ids <- readRDS(file = paste0(output_dir,"/manic_lca/","lca_cases_only_analytical_dataset.rds"))


################################################################################
######## TABULATE NUMBERS OF IRRITABLE AND MANIC IN EACH CLASS WITH % ##########

if (NROW(lca[[optimum_model_in_res]]$posterior) != NROW(ids)) 
  print("Number of rows in the ID file do not match the number of samples in the LCA results!~~~")

posterior_probs <- as.data.frame(lca[[optimum_model_in_res]]$posterior)

posterior_probs$predclass <- (lca[[optimum_model_in_res]]$predclass)

posterior_probs <- cbind(posterior_probs,ids)

tabyl(posterior_probs$Ever.manic.or.excitable)
tabyl(posterior_probs$Ever.extreme.irritability)

### create new composite variable for tabulations

posterior_probs$Ever.manic.and.irritable <- NA
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable==1 ] <- 1
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.extreme.irritability==1 ] <- 2
posterior_probs$Ever.manic.and.irritable[posterior_probs$Ever.manic.or.excitable==1 & posterior_probs$Ever.extreme.irritability==1] <- 3

tabyl(posterior_probs$Ever.manic.and.irritable)

### 

tball <- list()

for (i in 1:optimum_classes) {
  
  tb1 <- tabyl(posterior_probs$Ever.manic.and.irritable[posterior_probs$predclass==i])
  names(tb1)[1] <- paste0("class_",i)
  tb1$categories <- c("Ever manic","Ever irritable","Ever manic and irritable")
  tball[[i]] <- tb1
  
}

tball2 <- bind_rows(tball, .id=c("class"))
tball2 <- tball2[c ("categories", "n", "percent","class")]
tball2$`N (percent)` <-
  paste0(tball2$n, " (", round(tball2$percent, digits = 2), ")")
tball2$Response <- as.factor(gsub("\\."," ", as.character(tball2$categories))) # remove periods if any exist



### plot numbers/percentages of irritable, manic or manic AND irritable in each class

png(filename = paste0(output_dir,"/manic_lca/","percent_manic_or_irritable_", optimum_classes,"class_barplot.png"), 
    width = 20, height = 13, units = "cm", res = 300)

ggplot(tball2, aes(class, percent)) +
  geom_bar(aes(fill = Response), position = "dodge", stat = "identity", alpha=0.9) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Percent of class") + theme_bw() +
  geom_text(
    aes(
      class,
      y = percent + 0.03,
      group = Response,
      label = format(
        percent,
        nsmall = 0,
        digits = 1,
        scientific = FALSE, size=4
      )
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
  )

dev.off()

png(filename = paste0(output_dir,"/manic_lca/","numbers_manic_or_irritable_", optimum_classes,"class_barplot.png"), 
    width = 20, height = 13, units = "cm", res = 300)

  ggplot(tball2, aes(class, n)) +   
    geom_bar(aes(fill = Response), position = "dodge", stat="identity", alpha=0.9)    +
    scale_fill_brewer(palette=palette_choice) +
    ylab("Number of individuals in class") + theme_bw() 
    
  
dev.off()

### make a table

print((tball2[c("class","Response","N (percent)")]))

write.table(
  tball2[c("class", "Response", "N (percent)")],
  file = paste0(
    output_dir,
    "/manic_lca/",
    "tabulations_manic_or_irritable_",
    optimum_classes,
    "class.csv"
  ),
  sep = "," ,
  row.names = F,
  quote = F
)

### could also create a composite variable based on contribution by summed class probabilities


########################################################################
######## ASSOCIATIONS WITH DURATION OF MANIC SYMPTOMATOLOGY ############

### dummy code the predicted class variable - comparisons between each class and all others

class(posterior_probs$predclass)
posterior_probs$predclass <- as.factor(posterior_probs$predclass)

res_dummy <- as.data.frame(model.matrix(~predclass -1 , data = posterior_probs))

### check
head(res_dummy)
str(res_dummy)
tabyl(res_dummy[,c("predclass2")])

posterior_probs <- cbind(posterior_probs, res_dummy)
tabyl(posterior_probs$Brief.duration.mania.or.irritability)

summary(glm(data=posterior_probs, predclass1 ~ Brief.duration.mania.or.irritability , family="binomial"))
summary(glm(data=posterior_probs, predclass2 ~ Brief.duration.mania.or.irritability , family="binomial"))
summary(glm(data=posterior_probs, predclass3 ~ Brief.duration.mania.or.irritability , family="binomial"))
summary(glm(data=posterior_probs, predclass4 ~ Brief.duration.mania.or.irritability , family="binomial"))
summary(glm(data=posterior_probs, predclass5 ~ Brief.duration.mania.or.irritability , family="binomial"))
### logistic regression with dummy coded contrasts of classes

### reshape data to long so that participants have as many rows as classes and probability for each row
posterior_long <- reshape(posterior_probs, idvar="eid", direction="long", 
                          varying=1:5,  v.names=(f=c("probabilities")), timevar="predclass")

### check reshape to long for a random participant
posterior_long$probabilities[posterior_long$eid=="5798307"] # should have 5 probabilities sum to 1
if(round(sum(posterior_long$probabilities[posterior_long$eid=="5798307"]), digits=3) != 1) print("Error with reshape!")

posterior_long$predclass[posterior_long$eid=="5798307"] #should have values for X number of classes

duration_naive_res <- list()
duration_weighted_res <- list()

### create formula for multinomial regression
for(duration_i in c("Brief.duration.mania.or.irritability",
                    "Moderate.duration.mania.or.irritability",
                    "Extended.duration.mania.or.irritability")){
    
formula_duration <- paste0("predclass ~ ",duration_i)

### multinomial logistic regression naive to probabilities
mres <- multinom(formula_duration , data = posterior_probs)
tidy_mres <- tidy(mres,conf.int = TRUE, conf.level = 0.95, exponentiate=TRUE)
tidy_mres$model <- "naive"
tidy_mres$exposure <- duration_i

duration_naive_res[[duration_i]] <- tidy_mres

### multinomial logistic regression weighted for probabilities
mres_weighted <- multinom(formula_duration , data = posterior_long, weights = probabilities)

#extract and exponentiate coefficients
tidy_mres_weighted <- tidy(mres_weighted,conf.int = TRUE, conf.level = 0.95, exponentiate=TRUE)
tidy_mres_weighted$model <- "weighted"
tidy_mres_weighted$exposure <- duration_i

duration_weighted_res[[duration_i]] <- tidy_mres_weighted

}

### bind and drop intercept from results
duration_naive_res_df <- bind_rows(duration_naive_res)
duration_naive_res_df <- duration_naive_res_df[which(duration_naive_res_df$term!="(Intercept)"),]

duration_weighted_res_df <- bind_rows(duration_weighted_res)
duration_weighted_res_df <- duration_weighted_res_df[which(duration_weighted_res_df$term!="(Intercept)"),]

### create a 'comparison' column of statistical test for plot
duration_naive_res_df$comparison <- paste0("1v",duration_naive_res_df$`y.level`)
duration_weighted_res_df$comparison <- paste0("1v",duration_naive_res_df$`y.level`)

### remove periods if any exist from names
duration_naive_res_df$exposure <- as.factor(gsub("\\."," ", as.character(duration_naive_res_df$exposure))) 
duration_weighted_res_df$exposure <- as.factor(gsub("\\."," ", as.character(duration_weighted_res_df$exposure))) 

### plot multinomial results
png(filename = paste0(output_dir,"/manic_lca/","multinomial_naive_duration_", optimum_classes,"class_scatterplot.png"), 
    width = 20, height = 8, units = "cm", res = 300)

ggplot(duration_naive_res_df, aes(comparison, estimate, colour=comparison)) +
  geom_point(alpha=0.9) +
  geom_hline(yintercept=1, color="coral", alpha=0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  xlab("Class comparisons") +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Estimate (RR)") + theme_bw() +
  facet_wrap(~exposure) +
  theme(panel.background=element_rect(fill='white', colour='black'),
         strip.background=element_rect(fill='white', colour='white'))

dev.off()

###

png(filename = paste0(output_dir,"/manic_lca/","multinomial_weighted_duration_", optimum_classes,"class_scatterplot.png"), 
    width = 20, height = 8, units = "cm", res = 300)

ggplot(duration_weighted_res_df, aes(comparison, estimate, colour=comparison)) +
  geom_point(alpha=0.9) +
  geom_hline(yintercept=1, color="coral", alpha=0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  xlab("Class comparisons") +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Estimate (RR)") + theme_bw() +
  facet_wrap(~exposure) +
  theme(panel.background=element_rect(fill='white', colour='black'),
        strip.background=element_rect(fill='white', colour='white'))

dev.off()


### duration summarise by most probable class

tball_duration <- list()

for (i in 1:optimum_classes) {
  
  tb1 <- tabyl(posterior_probs$Duration.mania.or.irritability[posterior_probs$predclass==i])
  
  names(tb1)[1] <- paste0("class_",i)
  
  tb1$categories <- c("1 - Brief","2 - Moderate","3 - Extended", "Don't know","NA")
  
  tball_duration[[i]] <- tb1
  
}

tball_duration2 <- bind_rows(tball_duration, .id=c("class"))
tball_duration2 <- tball_duration2[c("categories", "n", "percent","class")]
tball_duration2$`N (percent)` <-
  paste0(tball_duration2$n, " (", round(tball_duration2$percent, digits = 2), ")")
tball_duration2$Duration <- as.factor(gsub("\\."," ", as.character(tball_duration2$categories))) # remove periods if any exist

### tabulate number of missing rows (hence individuals with NA for duration) and subset
missing_duration <- NROW(posterior_probs$Duration.mania.or.irritability) - NROW(tball_duration2[which(posterior_probs$Duration.mania.or.irritability!="NA"),])
tball_duration2 <- tball_duration2[which(tball_duration2$categories!="NA"),]

### plot raw comparisons by most probable class 
png(filename = paste0(output_dir,"/manic_lca/","percent_duration_", optimum_classes,"class_barplot.png"), 
    width = 20, height = 13, units = "cm", res = 300)

ggplot(tball_duration2, aes(class, percent)) +
  geom_bar(aes(fill = Duration), position = "dodge", stat = "identity", alpha=0.9) +
  scale_fill_brewer(palette = palette_choice) +
  ylab("Percent of class") + theme_bw() +
  geom_text(aes(class,
      y = percent + 0.03,
      group = Duration,
      label = format(
        percent,
        nsmall = 0,
        digits = 1,
        scientific = FALSE, size=4
      )
    ),
    color = "cornflowerblue",
    position = position_dodge(.9),
    hjust = .5
  )

dev.off()

##############################################################################
######## ASSOCIATIONS WITH DISRUPTIVENESS OF MANIC SYMPTOMATOLOGY ############


##############################################################################
######## ASSOCIATIONS WITH PRS OF PSYCHIATRIC DISORDERS ######################

### PRS to be used

prs_filenames <- c("BIPO01", "SCHI02", "DEPR07", "ADHD05", "AUTI05", "ANXI02")

### load PRS from files keep only those calculated at pvalue=1 threshold

pvalue_threshold <- "1" #can be "0.05" or other

prs_raw <- list()

for (i_prs in 1:length(prs_filenames)) { 
  
  prs_dat <- setDF(fread(paste0(dir_prs , prs_filenames[i_prs] , "_header.all.score")))
  
  prs_raw[[i_prs]] <- prs_dat[c("IID",paste0(prs_filenames[i_prs],"_", pvalue_threshold))]
  
}

### ensure same IDs across PRS files (no withdrawn IDs in some and not others)

id_freqs <- data.frame(table(unlist(lapply(prs_raw, '[[', 1))))

exclude_ids <- as.character(id_freqs$Var1[which(as.numeric(as.character(id_freqs$Freq)) < NROW(prs_raw))])

### remove any IDs that do not fit

prs_raw <- lapply(prs_raw, function(x) x[which(as.character(x[,1]) %notin% exclude_ids),]) 

prs_df <- bind_cols(prs_raw)

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

for (i_trait in c(paste0(prs_filenames, "_", pvalue_threshold)) ) {
  
  prs.lm <- lm(paste0(i_trait, " ~ ", paste0(cov_names, collapse = " + ")),
               data = for_pcs_regression)
  
  prs.res[[i_trait]] <- resid(prs.lm)
  
}

prs.res_df <- as.data.frame(prs.res)

### standardise PRS residuals

prs.res_df[names(prs.res_df)] <-
  lapply(prs.res_df[names(prs.res_df)], function(x)
    c(scale(x)))

### add ID column back in

prs.res_df$IID <- for_pcs_regression$IID

### check conforms to normal distribution with SD=1

hist(prs.res_df[[paste0("SCHI02_",pvalue_threshold)]], breaks = 100)

sd(prs.res_df[[paste0("SCHI02_",pvalue_threshold)]])

### merge with posterior probs and class memembership

posterior_probs_prs <- merge(posterior_probs, prs.res_df, by.x="eid", by.y="IID")

### find which samples weren't found in PRS

samples_without_prs <- NROW(posterior_probs) - NROW(posterior_probs_prs)

print(paste0("There are ", samples_without_prs," number of samples in latent class analysis without PRS"))

ids_without_prs <- posterior_probs$eid[posterior_probs$eid %notin%  posterior_probs_prs$eid]

### reshape data to long so that participants have as many rows as classes and a probability for each row/class

posterior_prs_long <- reshape(posterior_probs_prs, idvar="eid", direction="long", 
                          varying=2:6,  v.names=(f=c("probabilities")), timevar="predclass")

### check reshape to long correct by picking a random participant ID

random_eid <- sample(posterior_prs_long$eid, 1)

posterior_prs_long$probabilities[posterior_prs_long$eid==random_eid] # should have 5 probabilities sum to 1

if(round(sum(posterior_prs_long$probabilities[posterior_prs_long$eid==random_eid]), digits=6) != 1) print("Error with reshape!")

### check that values for X number of classes present (ie 5 if optimum model is 5 class)

posterior_prs_long$predclass[posterior_prs_long$eid=="5798307"] 

### save stata dataset for cross-checking results using mlogit in STATA

save.dta13(posterior_prs_long[,c("eid","probabilities","BIPO01_1","predclass")], 
           file="posterior_prs_long.dta")

### multinom regressions of PRS on classes - create blank lists for results

prs_naive_res <- list()

prs_weighted_res <- list()

### loop over either having class 1 or class 3 as reference class 

for(reference_i in c("main", "relevel")) {
  
  ### if loop is 'relevel' then change the reference class to "3"
  
  if(reference_i=="relevel") {
    
    posterior_probs_prs$predclass <- 
      relevel(as.factor(posterior_probs_prs$predclass), ref=3)
    
    posterior_prs_long$predclass <- 
      relevel(as.factor(posterior_prs_long$predclass), ref=3)
    
  } else { print("This is the main analysis (reference level '1')")}
  
  
  for (prs_trait_i in prs_filenames) {
    
    formula_prs <- paste0("predclass ~ ", prs_trait_i, "_",pvalue_threshold)
    
    ### multinomial logistic regression naive to probabilities
    
    mres <- multinom(formula_prs , data = posterior_probs_prs)
    tidy_mres <-  
      tidy(mres,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )
    
    tidy_mres$model <- "naive"
    
    tidy_mres$exposure <- prs_trait_i
    
    prs_naive_res[[prs_trait_i]] <- tidy_mres
    
    ### multinomial logistic regression weighted for probabilities (bias adjusted)
    
    mres_weighted <-
      multinom(formula_prs , data = posterior_prs_long, weights = probabilities)
    
    ### extract and exponentiate coefficients, add columns for labelling purposes
    
    tidy_mres_weighted <-
      tidy(mres_weighted,
        conf.int = TRUE,
        conf.level = 0.95,
        exponentiate = TRUE
      )
    
    tidy_mres_weighted$model <- "weighted"
    
    tidy_mres_weighted$exposure <- prs_trait_i
    
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
    paste0(levels(posterior_prs_long$predclass)[1],"v", prs_naive_res_df$`y.level`)
  
  prs_weighted_res_df$comparison <-
    paste0(levels(posterior_prs_long$predclass)[1],"v", prs_naive_res_df$`y.level`)
  
  ### remove periods if any exist from names
  prs_naive_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_naive_res_df$exposure)))
  
  prs_weighted_res_df$exposure <-
    as.factor(gsub("\\.", " ", as.character(prs_weighted_res_df$exposure)))
  
  ### plot multinomial results
  png(filename = paste0(output_dir,"/manic_lca/","multinomial_naive_prs_",reference_i,"_", optimum_classes,"class_scatterplot.png"), 
      width = 20, height = 10, units = "cm", res = 300)
  
  print(ggplot(prs_naive_res_df, aes(comparison, estimate, colour=comparison)) +
    geom_point(alpha=0.9) +
    geom_hline(yintercept=1, color="coral", alpha=0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Risk ratio (per SD increase in PRS)") + theme_bw() +
    facet_wrap(~exposure) +
    theme(panel.background=element_rect(fill='white', colour='black'),
          strip.background=element_rect(fill='white', colour='white')))
  
  dev.off()
  
  ###
  
  png(filename = paste0(output_dir,"/manic_lca/","multinomial_weighted_prs_", reference_i,"_", optimum_classes,"class_scatterplot.png"), 
      width = 20, height = 10, units = "cm", res = 300)
  
  print(ggplot(prs_weighted_res_df, aes(comparison, estimate, colour=comparison)) +
    geom_point(alpha=0.9) +
    geom_hline(yintercept=1, color="coral", alpha=0.4) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    xlab("Class comparisons") +
    scale_fill_brewer(palette = palette_choice) +
    ylab("Risk ratio (per SD increase in PRS)") + theme_bw() +
    facet_wrap(~exposure) +
    theme(panel.background=element_rect(fill='white', colour='black'),
          strip.background=element_rect(fill='white', colour='white')))
  
  dev.off()
  
  
}

gc()
  
### plot histograms of PRS by most likely class
  
for (prs_trait_i in prs_filenames) {
    
    prs_trait_i_p <- paste0(prs_trait_i, "_", pvalue_threshold)
    
    gg1 <- ggplot(posterior_probs_prs, aes(x = get(prs_trait_i_p), fill = predclass)) +
      geom_density(alpha = 0.3) +
      xlab(paste0("PRS ", prs_trait_i_p)) +
      theme_bw()
  
    png(filename = paste0(output_dir,"/manic_lca/","prs_", prs_trait_i_p , "_", optimum_classes,"class_densityplot.png"),
    width = 20, height = 15, units = "cm", res = 300)
    
      print(gg1)
    
    dev.off()
    
}
  
 
