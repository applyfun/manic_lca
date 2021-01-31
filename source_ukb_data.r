### Source MHQ data through ukbkings
### 30/09/2020
### R ARATHIMOS

### SET UP DIRECTORIES

dirs <- read.table("~/brc_scratch/scripts/paths.txt")

dir_prs <- as.character(dirs[1, ])
new_2019_qc_dir <- as.character(dirs[3, ])
qc_dir <- as.character(dirs[3, ])
project_dir <- as.character(dirs[9, ])

data_dir <- as.character(paste0(dirs[10, ], "data/manic_lca"))
output_dir <- as.character(paste0(dirs[10, ], "output/manic_lca"))
scripts_dir <- as.character(paste0(dirs[10, ], "scripts/manic_lca"))

#### load packages
library(data.table)
library(ukbkings)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringi)
library(janitor)
library(poLCA)
library(RNOmni)

######################### EXTRACT DATA USING UKBKINGS PACKAGE

f <- bio_field(project_dir)

setwd(data_dir)

### create a vector of field IDs to extract to write to file for extraction from UKB

vars204 <-
  c(paste0("204", "00"), paste0("2040", 1:9), paste0("204", 10:99))
vars205 <-
  c(paste0("205", "00"), paste0("2050", 1:9), paste0("205", 10:99))

### include sociodemographic fields also

fields_to_extract <- c("31", "189", "21022", "1558", "1239", "1249", "26410", "6138", 
                       "20544", "20499", "41270","20127", vars204, vars205)

### write a file of field IDs to extract one per line, no header

write.table(
  fields_to_extract,
  file = "fieldIDs_to_extract.txt",
  quote = F,
  row.names = F,
  col.names = F
)

### extract field ids using ukbkings

bio_phen(
  project_dir = project_dir,
  field = "fieldIDs_to_extract.txt",
  out = "mhq_phenotypes_subset"
) # e.g. "data/ukb" writes "data/ukb.rds")

### HES data

hesin <- bio_hesin(project_dir = project_dir, record = "hesin", hesin_dir = "raw/")

diag <- bio_hesin(project_dir = project_dir, record = "diag", hesin_dir = "raw/")

psych <- bio_hesin(project_dir = project_dir, record = "psych", hesin_dir = "raw/")

### load ICD10 codes
icd10_coding <- fread(paste0(data_dir, "/coding19.tsv"))

### subset to psych codes

icd10_F_coding <- icd10_coding[which(sapply(strsplit(icd10_coding$coding, ""), "[[", 1) == "F")]

dim(icd10_F_coding)

# Schizophrenia, schizotypal and delusional disorders F20-F29
Schizophrenia.or.psychotic.disorder.ICD <- icd10_F_coding[155:185, ]

# Mood [affective] disorders including bipolar F30, F31, F34, F38, F39
Bipolar.or.mood.disorder.ICD <- icd10_F_coding[c(186:202, 218:227), ]

# Mood [affective] disorders depressive type - F32-F33
Depression.ICD <- icd10_F_coding[c(203:217), ]

write.table(Schizophrenia.or.psychotic.disorder.ICD, file = paste0(output_dir, "/Schizophrenia.or.psychotic.disorder.ICD10_codes_reference_derived.tsv"), row.names = F, sep = "\t")

write.table(Bipolar.or.mood.disorder.ICD, file = paste0(output_dir, "/Bipolar.or.mood.disorder.ICD10_codes_reference_derived.tsv"), row.names = F, sep = "\t")

write.table(Depression.ICD, file = paste0(output_dir, "/Depression.ICD10_codes_reference_derived.tsv"), row.names = F, sep = "\t")

### load extracted data

MHQ <- readRDS(paste0(data_dir, "/mhq_phenotypes_subset.rds"))

### explore data

str(MHQ)
names(MHQ)

print(paste0("The number of rows/individuals with an MHQ completion date is ", NROW(MHQ[which(MHQ$`20400-0.0` != ""), ])))

### load var codings

mhq <- fread(paste0(data_dir, "/variable_codes_ukb_mhq.csv"))

names(mhq)
mhq[["UK Biobank Code"]]

### match name in coding names file to the names of MHQ dataframe

names(MHQ)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 20548 can be split in to multiple variables since multiple choices were allowed

# Coding	- Meaning
# 1 -	I was more talkative than usual
# 2 - 	I was more restless than usual
# 3 -	My thoughts were racing
# 5 - 	I needed less sleep than usual
# 6 -	I was more creative or had more ideas than usual
# 7 - 	I was easily distracted
# 8 -	I was more confident than usual
# 9 -	I was more active than usual
#-818	- Prefer not to answer


tabyl(MHQ$`20548-0.1`) # example response - instance 1
tabyl(MHQ$`20548-0.2`) # all -818 prefer not to answer are in instance 1
tabyl(MHQ$`20548-0.3`)
tabyl(MHQ$`20548-0.8`)

MHQ$More.talkative.manic.episode.derived[MHQ$`20548-0.1` == 1 |
  MHQ$`20548-0.2` == 1 |
  MHQ$`20548-0.3` == 1 |
  MHQ$`20548-0.4` == 1 |
  MHQ$`20548-0.5` == 1 |
  MHQ$`20548-0.6` == 1 |
  MHQ$`20548-0.7` == 1 |
  MHQ$`20548-0.8` == 1] <- 1


MHQ$More.restless.manic.episode.derived[MHQ$`20548-0.1` == 2 |
  MHQ$`20548-0.2` == 2 |
  MHQ$`20548-0.3` == 2 |
  MHQ$`20548-0.4` == 2 |
  MHQ$`20548-0.5` == 2 |
  MHQ$`20548-0.6` == 2 |
  MHQ$`20548-0.7` == 2 |
  MHQ$`20548-0.8` == 2] <- 1

MHQ$Thoughts.racing.manic.episode.derived[MHQ$`20548-0.1` == 3 |
  MHQ$`20548-0.2` == 3 |
  MHQ$`20548-0.3` == 3 |
  MHQ$`20548-0.4` == 3 |
  MHQ$`20548-0.5` == 3 |
  MHQ$`20548-0.6` == 3 |
  MHQ$`20548-0.7` == 3 |
  MHQ$`20548-0.8` == 3] <- 1

MHQ$Less.sleep.manic.episode.derived[MHQ$`20548-0.1` == 5 |
  MHQ$`20548-0.2` == 5 |
  MHQ$`20548-0.3` == 5 |
  MHQ$`20548-0.4` == 5 |
  MHQ$`20548-0.5` == 5 |
  MHQ$`20548-0.6` == 5 |
  MHQ$`20548-0.7` == 5 |
  MHQ$`20548-0.8` == 5] <- 1

MHQ$More.creative.manic.episode.derived[MHQ$`20548-0.1` == 6 |
  MHQ$`20548-0.2` == 6 |
  MHQ$`20548-0.3` == 6 |
  MHQ$`20548-0.4` == 6 |
  MHQ$`20548-0.5` == 6 |
  MHQ$`20548-0.6` == 6 |
  MHQ$`20548-0.7` == 6 |
  MHQ$`20548-0.8` == 6] <- 1

MHQ$Easily.distracted.manic.episode.derived[MHQ$`20548-0.1` == 7 |
  MHQ$`20548-0.2` == 7 |
  MHQ$`20548-0.3` == 7 |
  MHQ$`20548-0.4` == 7 |
  MHQ$`20548-0.5` == 7 |
  MHQ$`20548-0.6` == 7 |
  MHQ$`20548-0.7` == 7 |
  MHQ$`20548-0.8` == 7] <- 1

MHQ$More.confident.manic.episode.derived[MHQ$`20548-0.1` == 8 |
  MHQ$`20548-0.2` == 8 |
  MHQ$`20548-0.3` == 8 |
  MHQ$`20548-0.4` == 8 |
  MHQ$`20548-0.5` == 8 |
  MHQ$`20548-0.6` == 8 |
  MHQ$`20548-0.7` == 8 |
  MHQ$`20548-0.8` == 8] <- 1

MHQ$More.active.manic.episode.derived[MHQ$`20548-0.1` == 9 |
  MHQ$`20548-0.2` == 9 |
  MHQ$`20548-0.3` == 9 |
  MHQ$`20548-0.4` == 9 |
  MHQ$`20548-0.5` == 9 |
  MHQ$`20548-0.6` == 9 |
  MHQ$`20548-0.7` == 9 |
  MHQ$`20548-0.8` == 9] <- 1

### question was asked when Field 20501 was Yes or Field 20502 was Yes.

manic_derived_vars <- c(
  "More.active.manic.episode.derived", "More.confident.manic.episode.derived", "Easily.distracted.manic.episode.derived",
  "More.creative.manic.episode.derived", "Less.sleep.manic.episode.derived", "Thoughts.racing.manic.episode.derived",
  "More.restless.manic.episode.derived", "More.talkative.manic.episode.derived"
)

### add negative responses (controls) back in to individual manic symptoms responses
for (i in manic_derived_vars) {
  MHQ[[i]][which(MHQ$`20501-0.0` == 0 & MHQ$`20502-0.0` == 0)] <- 0
  MHQ[[i]][which(is.na(MHQ[[i]]) & (MHQ$`20502-0.0` == 1 | MHQ$`20501-0.0` == 1))] <- 0

  print(i)
  print(tabyl(MHQ[[i]]))
}

### make NA particpants that prefered not to answer symptoms subquestion

tabyl(MHQ$More.active.manic.episode.derived)

MHQ$More.talkative.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$More.restless.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$Thoughts.racing.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$Less.sleep.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$More.creative.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$Easily.distracted.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$More.confident.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA
MHQ$More.active.manic.episode.derived[MHQ$`20548-0.1` == -818] <- NA

tabyl(MHQ$More.active.manic.episode.derived)

###

### subset to individuals with response to question on manic or irritable episode 

MHQ <- MHQ[which(MHQ$`20501-0.0` == 1 | MHQ$`20502-0.0` == 1 | MHQ$`20501-0.0` == 0 | MHQ$`20502-0.0` == 0), ]
dim(MHQ)

### subset to individuals with a response to subquestion symptoms also

MHQ <- MHQ[which(!is.na(MHQ$More.active.manic.episode.derived)), ]
dim(MHQ)

### create an ID list for later merge
MHQ_ids <- MHQ$eid

### kep only responses to manic symtoms, no ID
manic_check_MHQ <- MHQ[c(manic_derived_vars)]

#### checks
dim(manic_check_MHQ)
str(manic_check_MHQ)
tabyl(manic_check_MHQ[, 3])

### convert to factors instead of numeric
manic_check_MHQ[] <- lapply(manic_check_MHQ, factor)
str(manic_check_MHQ)

### create binary variables of ever manic and ever irritable
MHQ$Ever.manic.or.excitable <- MHQ$`20501-0.0`
tabyl(MHQ$Ever.manic.or.excitable)

MHQ$Ever.manic.or.excitable[which(MHQ$Ever.manic.or.excitable < 0)] <- NA
tabyl(MHQ$Ever.manic.or.excitable)

MHQ$Ever.extreme.irritability <- MHQ$`20502-0.0`
tabyl(MHQ$Ever.extreme.irritability)

MHQ$Ever.extreme.irritability[which(MHQ$Ever.extreme.irritability < 0)] <- NA
tabyl(MHQ$Ever.extreme.irritability)

### create variables for irritability duration
MHQ$Brief.duration.mania.or.irritability <- NA
MHQ$Brief.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 1)] <- 1
tabyl(MHQ$Brief.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report moderate or extended duration
MHQ$Brief.duration.mania.or.irritability[which(MHQ$`20501-0.0` == 0 & MHQ$`20502-0.0` == 0)] <- 0
MHQ$Brief.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 2)] <- 0
MHQ$Brief.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 3)] <- 0
tabyl(MHQ$Brief.duration.mania.or.irritability)

### now moderate duration irritability
MHQ$Moderate.duration.mania.or.irritability <- NA
MHQ$Moderate.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 2)] <- 1
tabyl(MHQ$Moderate.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report brief or extended duration
MHQ$Moderate.duration.mania.or.irritability[which(MHQ$`20501-0.0` == 0 & MHQ$`20502-0.0` == 0)] <- 0
MHQ$Moderate.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 1)] <- 0
MHQ$Moderate.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 3)] <- 0
tabyl(MHQ$Moderate.duration.mania.or.irritability)

### now extended duration
MHQ$Extended.duration.mania.or.irritability <- NA
MHQ$Extended.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 3)] <- 1
tabyl(MHQ$Extended.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report moderate or brief duration
MHQ$Extended.duration.mania.or.irritability[which(MHQ$`20501-0.0` == 0 & MHQ$`20502-0.0` == 0)] <- 0
MHQ$Extended.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 1)] <- 0
MHQ$Extended.duration.mania.or.irritability[which(MHQ$`20492-0.0` == 2)] <- 0
tabyl(MHQ$Extended.duration.mania.or.irritability)

### now a categorical variable
MHQ$Duration.mania.or.irritability <- MHQ$`20492-0.0`
MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == -818)] <- NA
MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == -121)] <- NA
tabyl(MHQ$Duration.mania.or.irritability)

### 20493 is problems caused by manic or irritable episode
### How much of a problem have these "high" or "irritable" periods caused you?
MHQ$Problematic.mania.or.irritability <- MHQ$`20493-0.0`
MHQ$Problematic.mania.or.irritability[which(MHQ$Problematic.mania.or.irritability < 0)] <- NA
tabyl(MHQ$Problematic.mania.or.irritability)

### define fields "31","189","21022","1558","1239","1249" sociodemographics

MHQ$Sex <- MHQ$`31-0.0`

tabyl(MHQ$Sex)

MHQ$Sex[MHQ$Sex == 0] <- "Female"
MHQ$Sex[MHQ$Sex == 1] <- "Male"

###

MHQ$Townsend.deprivation <- MHQ$`189-0.0`

summary(MHQ$Townsend.deprivation)
hist(MHQ$Townsend.deprivation)

MHQ$Age <- MHQ$`21022-0.0`

summary(MHQ$Age)
hist(MHQ$Age)

### define alcohol intake frequency

MHQ$Alcohol.intake.frequency <- NA

MHQ$Alcohol.intake.frequency[MHQ$`1558-0.0` == 1] <- "Daily"

MHQ$Alcohol.intake.frequency[MHQ$`1558-0.0` == 2 | MHQ$`1558-0.0` == 3] <- "Weekly"

MHQ$Alcohol.intake.frequency[MHQ$`1558-0.0` == 4 | MHQ$`1558-0.0` == 5] <- "Occasionally"

MHQ$Alcohol.intake.frequency[MHQ$`1558-0.0` == 6] <- "Never"

tabyl(MHQ$`1558-0.0`)
tabyl(MHQ$Alcohol.intake.frequency)

head(MHQ$Alcohol.intake.frequency)

### define smoking

MHQ$Smoking <- NA

MHQ$Smoking[which(MHQ$`1239-0.0` == 1 | MHQ$`1239-0.0` == 2)] <- "Current"

MHQ$Smoking[which(MHQ$`1249-0.0` == 1 | MHQ$`1249-0.0` == 2)] <- "Past"

MHQ$Smoking[which(MHQ$`1249-0.0` == 3 | MHQ$`1249-0.0` == 4)] <- "Never"

### check

tabyl(MHQ$Smoking)
head(MHQ$Smoking)

### define education levels - field #6138

names(MHQ)

# "6138-0.0" "6138-0.1" eg instances
# create 4 levels of education:
# none - 4
# O-levels or CSE or equivalent - 3
# A-levels/AS-levels, National Vocational Qualification, HND/HNC or equivalent or other professional qualification - 2
# University degree - 1

MHQ$Education.level <- NA

MHQ$Education.level[MHQ$`6138-0.0` == 1 |
  MHQ$`6138-0.1` == 1 |
  MHQ$`6138-0.2` == 1 |
  MHQ$`6138-0.3` == 1 |
  MHQ$`6138-0.4` == 1 |
  MHQ$`6138-0.5` == 1] <- "University.degree"

MHQ$Education.level[MHQ$`6138-0.0` == 2 |
  MHQ$`6138-0.1` == 2 |
  MHQ$`6138-0.2` == 2 |
  MHQ$`6138-0.3` == 2 |
  MHQ$`6138-0.4` == 2 |
  MHQ$`6138-0.5` == 2] <- "A.levels.NVQ.HNC.or.HND"

MHQ$Education.level[MHQ$`6138-0.0` == 3 |
  MHQ$`6138-0.1` == 3 |
  MHQ$`6138-0.2` == 3 |
  MHQ$`6138-0.3` == 3 |
  MHQ$`6138-0.4` == 3 |
  MHQ$`6138-0.5` == 3] <- "O.levels.or.CSE"

MHQ$Education.level[MHQ$`6138-0.0` == 4 |
  MHQ$`6138-0.1` == 4 |
  MHQ$`6138-0.2` == 4 |
  MHQ$`6138-0.3` == 4 |
  MHQ$`6138-0.4` == 4 |
  MHQ$`6138-0.5` == 4] <- "O.levels.or.CSE"

MHQ$Education.level[MHQ$`6138-0.0` == 5 |
  MHQ$`6138-0.1` == 5 |
  MHQ$`6138-0.2` == 5 |
  MHQ$`6138-0.3` == 5 |
  MHQ$`6138-0.4` == 5 |
  MHQ$`6138-0.5` == 5] <- "A.levels.NVQ.HNC.or.HND"

MHQ$Education.level[MHQ$`6138-0.0` == 6 |
  MHQ$`6138-0.1` == 6 |
  MHQ$`6138-0.2` == 6 |
  MHQ$`6138-0.3` == 6 |
  MHQ$`6138-0.4` == 6 |
  MHQ$`6138-0.5` == 6] <- "A.levels.NVQ.HNC.or.HND"

MHQ$Education.level[MHQ$`6138-0.0` == -7 |
  MHQ$`6138-0.1` == -7 |
  MHQ$`6138-0.2` == -7 |
  MHQ$`6138-0.3` == -7 |
  MHQ$`6138-0.4` == -7 |
  MHQ$`6138-0.5` == -7] <- "None"

tabyl(MHQ$Education.level)

### MDI
# MHQ$`26410-0.0`


################################### DIAGNOSES ###########################################

##### Adapted from script by Kylie Glanville 09/11/2020

tabyl(MHQ$`20499-0.0`)

MHQ$Ever.sought.or.received.professional.help.for.mental.distress <- MHQ$`20499-0.0`

tabyl(MHQ$Ever.sought.or.received.professional.help.for.mental.distress)

### rename all 20544 array indices

names(MHQ)

names(MHQ) <- gsub("20544\\-0", "Mental.health.problems.ever.diagnosed.by.a.professional", names(MHQ))

### Physician Diagnosed Disorders

# To the question: "Have you been diagnosed with one or more of the following", there are 16 possible responses (i.e. array=16,  see MHQ pMHQ for details). Therefore, an individual can report up to 16 diagnoses. So for each individual, there are 16 columns of data for this question, such that each column contains an integer representing a possible response (e.g. Social Phobia =1) or NA.
# The r code is creating a column for social anxiety or social phobia by first looking at the column Ever.sought.or.received.professional.help.for.mental.distress and adding an NA if the response was NA, then looking at each array column for Mental.health.problems.ever.diagnosed.by.a.professional and looking for the integer associated with a phenotype (e.g. 1 = social anxiety). If the column does not contain NA and does contain the relevant integer, participant is coded as 1 for case, else 0.

MHQ$Self.reported.SocPhobia <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 1) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 1), 1, 0)
))

MHQ$Self.reported.Schizophrenia <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 2) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 2), 1, 0)
))

MHQ$Self.reported.PsychosisOther <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 3) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 3), 1, 0)
))

# The next line of code is creating a column for Self.reported.PsychosisAny based on the last two columns created, i.e. if participant has scz OR psychosis(other), they will be a case in this column.

MHQ$Self.reported.PsychosisAny <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Self.reported.Schizophrenia) & Self.reported.Schizophrenia == 1) | (!is.na(Self.reported.PsychosisOther) & Self.reported.PsychosisOther == 1), 1, 0)
))


MHQ$Self.reported.Personality.disorder <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 4) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 4), 1, 0)
))

MHQ$Self.reported.OtherPhobia <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 5) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 5), 1, 0)
))

MHQ$Self.reported.PanicAttacks <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 6) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 6), 1, 0)
))

MHQ$Self.reported.OCD <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 7) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 7), 1, 0)
))

MHQ$Self.reported.Mania.or.bipolar.disorder <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 10) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 10), 1, 0)
))

MHQ$Self.reported.Depression <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 11) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 11), 1, 0)
))

MHQ$Self.reported.Mood <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse(!is.na(Self.reported.Mania.or.bipolar.disorder) & Self.reported.Mania.or.bipolar.disorder == 1 |
    !is.na(Self.reported.Depression) & Self.reported.Depression == 1, 1, 0)
))

MHQ$Self.reported.BulimiaNervosa <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 12) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 12), 1, 0)
))

MHQ$Self.reported.BingeEating <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 13) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 13), 1, 0)
))

MHQ$Self.reported.ASD <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 14) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 14), 1, 0)
))

MHQ$Self.reported.GAD.and.others <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 15) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 15), 1, 0)
))

MHQ$Self.reported.AnorexiaNervosa <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 16) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 16), 1, 0)
))

MHQ$Self.reported.EatingDisorderAny <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Self.reported.AnorexiaNervosa) & Self.reported.AnorexiaNervosa == 1) |
    (!is.na(Self.reported.BulimiaNervosa) & Self.reported.BulimiaNervosa == 1) |
    (!is.na(Self.reported.BingeEating) & Self.reported.BingeEating == 1), 1, 0)
))

MHQ$Self.reported.Agoraphobia <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 17) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 17), 1, 0)
))

MHQ$Self.reported.AnxietyAny <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Self.reported.SocPhobia) & Self.reported.SocPhobia == 1) |
    (!is.na(Self.reported.GAD.and.others) & Self.reported.GAD.and.others == 1) |
    (!is.na(Self.reported.PanicAttacks) & Self.reported.PanicAttacks == 1) |
    (!is.na(Self.reported.Agoraphobia) & Self.reported.Agoraphobia == 1) |
    (!is.na(Self.reported.OtherPhobia) & Self.reported.OtherPhobia == 1) |
    (!is.na(Self.reported.OCD) & Self.reported.OCD == 1), 1, 0)
))

MHQ$Self.reported.ADHD <- with(MHQ, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.1) & Mental.health.problems.ever.diagnosed.by.a.professional.1 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.2) & Mental.health.problems.ever.diagnosed.by.a.professional.2 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.3) & Mental.health.problems.ever.diagnosed.by.a.professional.3 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.4) & Mental.health.problems.ever.diagnosed.by.a.professional.4 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.5) & Mental.health.problems.ever.diagnosed.by.a.professional.5 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.6) & Mental.health.problems.ever.diagnosed.by.a.professional.6 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.7) & Mental.health.problems.ever.diagnosed.by.a.professional.7 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.8) & Mental.health.problems.ever.diagnosed.by.a.professional.8 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.9) & Mental.health.problems.ever.diagnosed.by.a.professional.9 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.10) & Mental.health.problems.ever.diagnosed.by.a.professional.10 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.11) & Mental.health.problems.ever.diagnosed.by.a.professional.11 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.12) & Mental.health.problems.ever.diagnosed.by.a.professional.12 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.13) & Mental.health.problems.ever.diagnosed.by.a.professional.13 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.14) & Mental.health.problems.ever.diagnosed.by.a.professional.14 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.15) & Mental.health.problems.ever.diagnosed.by.a.professional.15 == 18) |
    (!is.na(Mental.health.problems.ever.diagnosed.by.a.professional.16) & Mental.health.problems.ever.diagnosed.by.a.professional.16 == 18), 1, 0)
))

diagnoses_vars <- c(
  "Self.reported.ADHD", "Self.reported.AnxietyAny", "Self.reported.Agoraphobia", "Self.reported.EatingDisorderAny", "Self.reported.AnorexiaNervosa", "Self.reported.GAD.and.others", "Self.reported.ASD", "Self.reported.BingeEating",
  "Self.reported.BulimiaNervosa", "Self.reported.Mood", "Self.reported.Depression", "Self.reported.Mania.or.bipolar.disorder", "Self.reported.OCD", "Self.reported.OtherPhobia", "Self.reported.Personality.disorder", "Self.reported.PsychosisAny",
  "Self.reported.Schizophrenia", "Self.reported.SocPhobia", "Self.reported.PsychosisOther"
)

### describe self-reported diagnoses

for (i in diagnoses_vars) {
  print(i)
  print(tabyl(MHQ[[i]]))
}

### sort neuroticism score
MHQ$Neuroticism.score <- MHQ$`20127-0.0`

######################################

Schizophrenia.or.psychotic.disorder.ICD_df <- fread(paste0(output_dir, "/Schizophrenia.or.psychotic.disorder.ICD10_codes_reference_derived.tsv"))

Depression.ICD_df <- fread(paste0(output_dir, "/Depression.ICD10_codes_reference_derived.tsv"))

Bipolar.or.mood.disorder.ICD_df <- fread(paste0(output_dir, "/Bipolar.or.mood.disorder.ICD10_codes_reference_derived.tsv"))

### pull out IDs with specific ICD10 codes

scz_diag <- diag[diag$diag_icd10 %in% Schizophrenia.or.psychotic.disorder.ICD_df$coding, ]
NROW(unique(scz_diag$eid))

bip_diag <- diag[diag$diag_icd10 %in% Bipolar.or.mood.disorder.ICD_df$coding, ]
NROW(unique(bip_diag$eid))

dep_diag <- diag[diag$diag_icd10 %in% Depression.ICD_df$coding, ]
NROW(unique(dep_diag$eid))

### tabulate instances of each ICD10 code

instances_scz <- tabyl(scz_diag$diag_icd10)

instances_dep <- tabyl(dep_diag$diag_icd10)

instances_bip <- tabyl(bip_diag$diag_icd10)

### pull out unique IDs and write to file

icd10_cases_eid <- as.data.frame(unique(diag$eid))
names(icd10_cases_eid) <- c("eid")

icd10_cases_eid$Schizophrenia.or.psychotic.disorder.ICD <- "No"
icd10_cases_eid$Schizophrenia.or.psychotic.disorder.ICD[which(icd10_cases_eid$eid %in% unique(scz_diag$eid))] <- "Yes"

icd10_cases_eid$Depression.ICD <- "No"
icd10_cases_eid$Depression.ICD[which(icd10_cases_eid$eid %in% unique(dep_diag$eid))] <- "Yes"

icd10_cases_eid$Bipolar.or.mood.disorder.ICD <- "No"
icd10_cases_eid$Bipolar.or.mood.disorder.ICD[which(icd10_cases_eid$eid %in% unique(bip_diag$eid))] <- "Yes"

### check

tabyl(icd10_cases_eid$Bipolar.or.mood.disorder.ICD)

### save ICD hesin data in order to merge in to main MHQ dataset at multinomial analysis stage

icd_hesin_vars <- c("eid", "Bipolar.or.mood.disorder.ICD", "Depression.ICD", "Schizophrenia.or.psychotic.disorder.ICD")

saveRDS(icd10_cases_eid[, c(icd_hesin_vars)], file = paste0(output_dir, "/", "ICD_diagnoses_hesin_cleaned_data.rds"))

### write file

subcols <- c(
  "eid", "Sex", "Age", "Smoking", "Alcohol.intake.frequency", "Townsend.deprivation",
  "Extended.duration.mania.or.irritability", "Moderate.duration.mania.or.irritability",
  "Brief.duration.mania.or.irritability", "Problematic.mania.or.irritability",
  "Duration.mania.or.irritability", c(manic_derived_vars), "Ever.extreme.irritability",
  "Ever.manic.or.excitable", "Education.level", "Neuroticism.score", c(diagnoses_vars)
)

saveRDS(MHQ[, c(subcols)], file = paste0(output_dir, "/", "manic_data_derived_subset_mhq.rds"))

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
