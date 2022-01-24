### Source PROTECT data for LCA
### 19/11/2020
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
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringi)
library(janitor)
library(RNOmni)

######################### LOAD DATASETS

MHQ <- setDF(fread(paste0(data_dir, "/Data/MHQ_PROT_DA_052_fullsample.csv")))

print(paste0("There are ", NROW(unique(MHQ$ResultsID)), " individuals in the mhq dataset..."))

dim(MHQ)

### sex and education phenotypes - keep first response by oldest date

SOCIODEM <- setDF(fread(paste0(data_dir, "/Data/DEMOGS_PROT_DA_052_fullsample.csv")))

SOCIODEM$Submitted_SOCIODEM <- as.Date(SOCIODEM$Submitted, "%d/%m/%Y")

print(paste0("There are ", NROW(unique(SOCIODEM$ResultsID)), " individuals in the SOCIODEM dataset..."))

dim(SOCIODEM)

SOCIODEM <- as.data.frame(SOCIODEM %>%
  group_by(ResultsID) %>%
  arrange(Submitted_SOCIODEM) %>%
  slice(1))

dim(SOCIODEM)

### smoking and alcohol phenotypes - keep first response by oldest date

LQ <- setDF(fread(paste0(data_dir, "/Data/LQ_Sleep_PROT_DA_052_fullsample.csv")))

print(paste0("There are ", NROW(unique(LQ$ResultsID)), " individuals in the LQ dataset..."))

LQ$Submitted_LQ <- as.Date(LQ$Submitted, "%d/%m/%Y")

dim(LQ)

LQ <- as.data.frame(LQ %>%
  group_by(ResultsID) %>%
  arrange(Submitted_LQ) %>%
  slice(1))

dim(LQ)

### explore data

names(MHQ)

### repeat question: X5.1 and X4.1

tabyl((MHQ[["X5.1.Have.you.ever.had.a.period.of.time.when.you.were.feeling.so.good...high....excited...or..hyper..that"]]))

tabyl((MHQ[["X4.1.Have.you.ever.had.a.period.of.time.when.you.were.feeling.so.good...high....excited...or..hyper..that"]]))

###

tabyl((MHQ[["X4.2.Have.you.ever.had.a.period.of.time.when.you.were.so.irritable.that.you.found.yourself.shouting.at.pe"]]))

tabyl(MHQ[["X5.2.Please.try.to.remember.a.period.when.you.were.in.a..high..or..irritable..state.and.select.all.of.the"]])

### subset data

main_fields <- c(
  "X4.2.Have.you.ever.had.a.period.of.time.when.you.were.so.irritable.that.you.found.yourself.shouting.at.pe",
  "X4.1.Have.you.ever.had.a.period.of.time.when.you.were.feeling.so.good...high....excited...or..hyper..that",
  "X5.2.Please.try.to.remember.a.period.when.you.were.in.a..high..or..irritable..state.and.select.all.of.the",
  "X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted.",
  "X5.4.How.much.of.a.problem.have.these..high..or..irritable..periods.caused.you.",
  "X1.3.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional.",
  "X1.4.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional.",
  "X1.2.In.your.life..did.you.seek.or.receive.help.from.a.professional..medical.doctor..psychologist..social"
)

MHQ <- MHQ[c("TecID","ResultsID", "Submitted", c(main_fields))]

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Response to questions - values of each response:

# Yes 01, Do not know NA, Prefer not to answer DA, No 00

# Do not know NA, Prefer not to answer DA, Less than 24 hours 01, At least a day, but less than a week 02, A week or more 03,

# Do not know NA, Prefer not to answer DA, No problems 00, Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life 01,

# Symptoms responses: 
# Prefer not to answer DA, None of the above 00, I was more active than usual 01,
# I was more talkative than usual 02, I needed less sleep than usual 03,
# I was more creative or had more ideas than usual 04, I was more restless than usual 05,
# I was more confident than usual 06, My thoughts were racing 07, I was easily distracted 08,


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### define new variables and set missing/did not respond to NA

MHQ$Ever.manic.or.excitable <- MHQ[["X4.1.Have.you.ever.had.a.period.of.time.when.you.were.feeling.so.good...high....excited...or..hyper..that"]]

### values "DA", ".", "" <blank> and "NA" <character> are considered NA/missing responses in all fields

MHQ$Ever.manic.or.excitable[which(MHQ$Ever.manic.or.excitable == "DA")] <- NA
MHQ$Ever.manic.or.excitable[which(MHQ$Ever.manic.or.excitable == "")] <- NA
MHQ$Ever.manic.or.excitable[which(MHQ$Ever.manic.or.excitable == ".")] <- NA
MHQ$Ever.manic.or.excitable[which(MHQ$Ever.manic.or.excitable == "NA")] <- NA

tabyl(MHQ$Ever.manic.or.excitable)
class(MHQ$Ever.manic.or.excitable)

MHQ$Ever.manic.or.excitable <- as.numeric(MHQ$Ever.manic.or.excitable)

tabyl(MHQ$Ever.manic.or.excitable)

###

MHQ$Ever.extreme.irritability <- ((MHQ[["X4.2.Have.you.ever.had.a.period.of.time.when.you.were.so.irritable.that.you.found.yourself.shouting.at.pe"]]))

MHQ$Ever.extreme.irritability[which(MHQ$Ever.extreme.irritability == "DA")] <- NA
MHQ$Ever.extreme.irritability[which(MHQ$Ever.extreme.irritability == "")] <- NA
MHQ$Ever.extreme.irritability[which(MHQ$Ever.extreme.irritability == ".")] <- NA
MHQ$Ever.extreme.irritability[which(MHQ$Ever.extreme.irritability == "NA")] <- NA

tabyl(MHQ$Ever.extreme.irritability)
class(MHQ$Ever.extreme.irritability)

MHQ$Ever.extreme.irritability <- as.numeric(MHQ$Ever.extreme.irritability)

table(MHQ$Ever.extreme.irritability, MHQ$Ever.manic.or.excitable)

### responses to symptoms

MHQ$bulksymptoms <- MHQ[["X5.2.Please.try.to.remember.a.period.when.you.were.in.a..high..or..irritable..state.and.select.all.of.the"]]

### create a temp "bulksymptoms" vairbale which to split on delimiter

tabyl(MHQ$bulksymptoms)

MHQ$bulksymptoms[which(MHQ$bulksymptoms == ".")] <- NA
MHQ$bulksymptoms[which(MHQ$bulksymptoms == "DA")] <- NA
MHQ$bulksymptoms[which(MHQ$bulksymptoms == "NA")] <- NA
MHQ$bulksymptoms[which(MHQ$bulksymptoms == "")] <- NA

tabyl(MHQ$bulksymptoms)

### split the strings containing symptom codes on ","

split_list <- strsplit(MHQ$bulksymptoms, ",")

str(split_list)

### individuals that have selected "DA" (i.e. 'Don't Know' response) along with a symptom code are considered valid responses (DA is ignored)

n.obs <- sapply(split_list, length)

seq.max <- seq_len(max(n.obs))

mat <- as.data.frame(t(sapply(split_list, "[", i = seq.max)), stringsAsFactors = F)

### bind individual symptom codes in matrix - which are now in individual columns back in to MHQ

MHQ <- cbind(MHQ, mat)

class(MHQ$`V1`)

### create individual variable response by pulling out symptom codes from the split matrix

MHQ$More.active.manic.episode.derived[MHQ$`V1` == "1" |
  MHQ$`V2` == "1" |
  MHQ$`V3` == "1" |
  MHQ$`V4` == "1" |
  MHQ$`V5` == "1" |
  MHQ$`V6` == "1" |
  MHQ$`V7` == "1" |
  MHQ$`V8` == "1"] <- 1


MHQ$More.talkative.manic.episode.derived[MHQ$`V1` == "2" |
  MHQ$`V2` == "2" |
  MHQ$`V3` == "2" |
  MHQ$`V4` == "2" |
  MHQ$`V5` == "2" |
  MHQ$`V6` == "2" |
  MHQ$`V7` == "2" |
  MHQ$`V8` == "2"] <- 1

MHQ$Less.sleep.manic.episode.derived[MHQ$`V1` == "3" |
  MHQ$`V2` == "3" |
  MHQ$`V3` == "3" |
  MHQ$`V4` == "3" |
  MHQ$`V5` == "3" |
  MHQ$`V6` == "3" |
  MHQ$`V7` == "3" |
  MHQ$`V8` == "3"] <- 1

MHQ$More.creative.manic.episode.derived[MHQ$`V1` == "4" |
  MHQ$`V2` == "4" |
  MHQ$`V3` == "4" |
  MHQ$`V4` == "4" |
  MHQ$`V5` == "4" |
  MHQ$`V6` == "4" |
  MHQ$`V7` == "4" |
  MHQ$`V8` == "4"] <- 1

MHQ$More.restless.manic.episode.derived[MHQ$`V1` == "5" |
  MHQ$`V2` == "5" |
  MHQ$`V3` == "5" |
  MHQ$`V4` == "5" |
  MHQ$`V5` == "5" |
  MHQ$`V6` == "5" |
  MHQ$`V7` == "5" |
  MHQ$`V8` == "5"] <- 1

MHQ$More.confident.manic.episode.derived[MHQ$`V1` == "6" |
  MHQ$`V2` == "6" |
  MHQ$`V3` == "6" |
  MHQ$`V4` == "6" |
  MHQ$`V5` == "6" |
  MHQ$`V6` == "6" |
  MHQ$`V7` == "6" |
  MHQ$`V8` == "6"] <- 1

MHQ$Thoughts.racing.manic.episode.derived[MHQ$`V1` == "7" |
  MHQ$`V2` == "7" |
  MHQ$`V3` == "7" |
  MHQ$`V4` == "7" |
  MHQ$`V5` == "7" |
  MHQ$`V6` == "7" |
  MHQ$`V7` == "7" |
  MHQ$`V8` == "7"] <- 1

MHQ$Easily.distracted.manic.episode.derived[MHQ$`V1` == "8" |
  MHQ$`V2` == "8" |
  MHQ$`V3` == "8" |
  MHQ$`V4` == "8" |
  MHQ$`V5` == "8" |
  MHQ$`V6` == "8" |
  MHQ$`V7` == "8" |
  MHQ$`V8` == "8"] <- 1

## repeat for responses with a zero in fron of the symptom code

MHQ$More.active.manic.episode.derived[MHQ$`V1` == "01" |
  MHQ$`V2` == "01" |
  MHQ$`V3` == "01" |
  MHQ$`V4` == "01" |
  MHQ$`V5` == "01" |
  MHQ$`V6` == "01" |
  MHQ$`V7` == "01" |
  MHQ$`V8` == "01"] <- 1


MHQ$More.talkative.manic.episode.derived[MHQ$`V1` == "02" |
  MHQ$`V2` == "02" |
  MHQ$`V3` == "02" |
  MHQ$`V4` == "02" |
  MHQ$`V5` == "02" |
  MHQ$`V6` == "02" |
  MHQ$`V7` == "02" |
  MHQ$`V8` == "02"] <- 1

MHQ$Less.sleep.manic.episode.derived[MHQ$`V1` == "03" |
  MHQ$`V2` == "03" |
  MHQ$`V3` == "03" |
  MHQ$`V4` == "03" |
  MHQ$`V5` == "03" |
  MHQ$`V6` == "03" |
  MHQ$`V7` == "03" |
  MHQ$`V8` == "03"] <- 1

MHQ$More.creative.manic.episode.derived[MHQ$`V1` == "04" |
  MHQ$`V2` == "04" |
  MHQ$`V3` == "04" |
  MHQ$`V4` == "04" |
  MHQ$`V5` == "04" |
  MHQ$`V6` == "04" |
  MHQ$`V7` == "04" |
  MHQ$`V8` == "04"] <- 1

MHQ$More.restless.manic.episode.derived[MHQ$`V1` == "05" |
  MHQ$`V2` == "05" |
  MHQ$`V3` == "05" |
  MHQ$`V4` == "05" |
  MHQ$`V5` == "05" |
  MHQ$`V6` == "05" |
  MHQ$`V7` == "05" |
  MHQ$`V8` == "05"] <- 1

MHQ$More.confident.manic.episode.derived[MHQ$`V1` == "06" |
  MHQ$`V2` == "06" |
  MHQ$`V3` == "06" |
  MHQ$`V4` == "06" |
  MHQ$`V5` == "06" |
  MHQ$`V6` == "06" |
  MHQ$`V7` == "06" |
  MHQ$`V8` == "06"] <- 1

MHQ$Thoughts.racing.manic.episode.derived[MHQ$`V1` == "07" |
  MHQ$`V2` == "07" |
  MHQ$`V3` == "07" |
  MHQ$`V4` == "07" |
  MHQ$`V5` == "07" |
  MHQ$`V6` == "07" |
  MHQ$`V7` == "07" |
  MHQ$`V8` == "07"] <- 1

MHQ$Easily.distracted.manic.episode.derived[MHQ$`V1` == "08" |
  MHQ$`V2` == "08" |
  MHQ$`V3` == "08" |
  MHQ$`V4` == "08" |
  MHQ$`V5` == "08" |
  MHQ$`V6` == "08" |
  MHQ$`V7` == "08" |
  MHQ$`V8` == "08"] <- 1


###

manic_derived_vars <- c(
  "More.active.manic.episode.derived", "More.confident.manic.episode.derived", "Easily.distracted.manic.episode.derived",
  "More.creative.manic.episode.derived", "Less.sleep.manic.episode.derived", "Thoughts.racing.manic.episode.derived",
  "More.restless.manic.episode.derived", "More.talkative.manic.episode.derived"
)

### add negative responses (controls) back in to individual manic symptoms responses

for (i in manic_derived_vars) {

  MHQ[[i]][which(MHQ$Ever.manic.or.excitable == "0" & MHQ$Ever.extreme.irritability == "0")] <- 0

  MHQ[[i]][which(is.na(MHQ[[i]]) & (MHQ$Ever.manic.or.excitable == "1" | MHQ$Ever.extreme.irritability == "1"))] <- 0

  MHQ[[i]][which(MHQ$bulksymptoms == "DA")] <- NA

  print(i)
  print(tabyl(MHQ[[i]]))

}

### duration of episode question responses

tabyl(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]])
class(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]])

### create variables for duration

MHQ$Brief.duration.mania.or.irritability <- NA

MHQ$Brief.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "1")] <- 1
tabyl(MHQ$Brief.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report moderate or extended duration

MHQ$Brief.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "." |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "DA" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "NA" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "" |
  is.na(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]]))] <- NA

MHQ$Brief.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "2")] <- 0
MHQ$Brief.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "3")] <- 0
tabyl(MHQ$Brief.duration.mania.or.irritability)

### now moderate duration irritability

MHQ$Moderate.duration.mania.or.irritability <- NA

MHQ$Moderate.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "2")] <- 2
tabyl(MHQ$Moderate.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report moderate or extended duration

MHQ$Moderate.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "." |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "DA" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "NA" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "" |
  is.na(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]]))] <- NA

MHQ$Moderate.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "1")] <- 0

MHQ$Moderate.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "3")] <- 0
tabyl(MHQ$Moderate.duration.mania.or.irritability)

### now extended duration

MHQ$Extended.duration.mania.or.irritability <- NA

MHQ$Extended.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "3")] <- 3
tabyl(MHQ$Extended.duration.mania.or.irritability)

### controls are individuals who reported no manic episode and no extreme irritability plus those who report moderate or extended duration

MHQ$Extended.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "." |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "DA" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "" |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "NA" |
  is.na(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]]))] <- NA

MHQ$Extended.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "2")] <- 0

MHQ$Extended.duration.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "1")] <- 0

tabyl(MHQ$Extended.duration.mania.or.irritability)

### now a categorical variable

MHQ$Duration.mania.or.irritability <- MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]]

MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == "")] <- NA
MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == ".")] <- NA
MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == "DA")] <- NA
MHQ$Duration.mania.or.irritability[which(MHQ$Duration.mania.or.irritability == "NA")] <- NA

tabyl(MHQ$Duration.mania.or.irritability)
class(MHQ$Duration.mania.or.irritability)

MHQ$Duration.mania.or.irritability <- as.numeric(MHQ$Duration.mania.or.irritability)

### How much of a problem have these "high" or "irritable" periods caused you?

tabyl(MHQ[["X5.4.How.much.of.a.problem.have.these..high..or..irritable..periods.caused.you."]])

MHQ$Problematic.mania.or.irritability <- MHQ[["X5.4.How.much.of.a.problem.have.these..high..or..irritable..periods.caused.you."]]

MHQ$Problematic.mania.or.irritability[which(MHQ$Problematic.mania.or.irritability == "DA")] <- NA

MHQ$Problematic.mania.or.irritability[which(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "." |
  MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]] == "" |
    is.na(MHQ[["X5.3.What.is.the.longest.time.that.these..high..or..irritable..periods.have.lasted."]]))] <- NA

tabyl(MHQ$Problematic.mania.or.irritability)

##############################################################
####################### DIAGNOSES ############################

# Diagnosis by a professional X1.4

# Do not know NA, Prefer not to answer DA, None of the above 00, Anorexia nervosa 01,
# Bulimia nervosa 02, Psychological over-eating or binge-eating 03, Schizophrenia 04,
# Any other type of psychosis or psychotic illness 05, A personality disorder 06,
# Autism, Asperger's or autistic spectrum disorder 07,
# Attention deficit or attention deficit and hyperactivity disorder (ADD/ADHD) 08,

# Diagnosis by a professional X1.3

# Do not know NA, Prefer not to answer DA, Depression 01, Mania, hypomania, bipolar or manic-depression 02,
# Anxiety, nerves or generalized anxiety disorder 03, Social anxiety or social phobia 04, Agoraphobia 05,
# Panic attacks 06, Obsessive compulsive disorder (OCD) 07, None of the above 00,

tabyl(MHQ[["X1.3.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional."]])

### second field with diagnoses

tabyl(MHQ[["X1.4.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional."]])

### ever sought help

tabyl(MHQ[["X1.2.In.your.life..did.you.seek.or.receive.help.from.a.professional..medical.doctor..psychologist..social"]])

MHQ$Ever.sought.or.received.professional.help.for.mental.distress <- as.character(MHQ[["X1.2.In.your.life..did.you.seek.or.receive.help.from.a.professional..medical.doctor..psychologist..social"]])

### coding for this var: Yes 1, No 2

MHQ$Ever.sought.or.received.professional.help.for.mental.distress[which(MHQ$Ever.sought.or.received.professional.help.for.mental.distress == "DA" |
  MHQ$Ever.sought.or.received.professional.help.for.mental.distress == "" |
  MHQ$Ever.sought.or.received.professional.help.for.mental.distress == "NA" |
  MHQ$Ever.sought.or.received.professional.help.for.mental.distress == ".")] <- NA

MHQ$Ever.sought.or.received.professional.help.for.mental.distress[which(MHQ$Ever.sought.or.received.professional.help.for.mental.distress == "2")] <- 0

MHQ$Ever.sought.or.received.professional.help.for.mental.distress <- as.numeric(MHQ$Ever.sought.or.received.professional.help.for.mental.distress)

tabyl(MHQ$Ever.sought.or.received.professional.help.for.mental.distress)

### split diagnoses fields in to individuals codes

MHQ$bulkdiagnoses1.3 <- MHQ[["X1.3.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional."]]

MHQ$bulkdiagnoses1.4 <- MHQ[["X1.4.Have.you.been.diagnosed.with.one.or.more.of.the.following.mental.health.problems.by.a.professional."]]

### split the strings containing symptom codes on ","

split_list <- strsplit(MHQ$bulkdiagnoses1.4, ",")

n.obs <- sapply(split_list, length)

seq.max <- seq_len(max(n.obs))

Bulk_diagnoses1 <- as.data.frame(t(sapply(split_list, "[", i = seq.max)), stringsAsFactors = F)

MHQ_help_tmp <- MHQ[, c("ResultsID", "Ever.sought.or.received.professional.help.for.mental.distress")]

Bulk_diagnoses1 <- cbind(Bulk_diagnoses1, MHQ_help_tmp)

### Physician Diagnosed Disorders

head(Bulk_diagnoses1)
str(Bulk_diagnoses1)

Bulk_diagnoses1$Self.reported.Schizophrenia <- with(Bulk_diagnoses1, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "4") |
    (!is.na(V2) & V2 == "4") |
    (!is.na(V3) & V3 == "4") |
    (!is.na(V4) & V4 == "4") |
    (!is.na(V5) & V5 == "4") |
    (!is.na(V6) & V6 == "4") |
    (!is.na(V7) & V7 == "4") |
    (!is.na(V8) & V8 == "4") |
    (!is.na(V1) & V1 == "04") |
    (!is.na(V2) & V2 == "04") |
    (!is.na(V3) & V3 == "04") |
    (!is.na(V4) & V4 == "04") |
    (!is.na(V5) & V5 == "04") |
    (!is.na(V6) & V6 == "04") |
    (!is.na(V7) & V7 == "04") |
    (!is.na(V8) & V8 == "04"), 1, 0)
))

Bulk_diagnoses1$Self.reported.PsychosisOther <- with(Bulk_diagnoses1, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "5") |
    (!is.na(V2) & V2 == "5") |
    (!is.na(V3) & V3 == "5") |
    (!is.na(V4) & V4 == "5") |
    (!is.na(V5) & V5 == "5") |
    (!is.na(V6) & V6 == "5") |
    (!is.na(V7) & V7 == "5") |
    (!is.na(V8) & V8 == "5") |
    (!is.na(V1) & V1 == "05") |
    (!is.na(V2) & V2 == "05") |
    (!is.na(V3) & V3 == "05") |
    (!is.na(V4) & V4 == "05") |
    (!is.na(V5) & V5 == "05") |
    (!is.na(V6) & V6 == "05") |
    (!is.na(V7) & V7 == "05") |
    (!is.na(V8) & V8 == "05"), 1, 0)
))


# The next line of code is creating a column for Self.reported.PsychosisAny based on the last two columns created, i.e. if participant has scz OR psychosis(other), they will be a case in this column.

Bulk_diagnoses1$Self.reported.PsychosisAny <- with(Bulk_diagnoses1, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(Self.reported.Schizophrenia) & Self.reported.Schizophrenia == 1) | (!is.na(Self.reported.PsychosisOther) & Self.reported.PsychosisOther == 1), 1, 0)
))


Bulk_diagnoses1$Self.reported.ASD <- with(Bulk_diagnoses1, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "7") |
    (!is.na(V2) & V2 == "7") |
    (!is.na(V3) & V3 == "7") |
    (!is.na(V4) & V4 == "7") |
    (!is.na(V5) & V5 == "7") |
    (!is.na(V6) & V6 == "7") |
    (!is.na(V7) & V7 == "7") |
    (!is.na(V8) & V8 == "7") |
    (!is.na(V1) & V1 == "07") |
    (!is.na(V2) & V2 == "07") |
    (!is.na(V3) & V3 == "07") |
    (!is.na(V4) & V4 == "07") |
    (!is.na(V5) & V5 == "07") |
    (!is.na(V6) & V6 == "07") |
    (!is.na(V7) & V7 == "07") |
    (!is.na(V8) & V8 == "07"), 1, 0)
))


Bulk_diagnoses1$Self.reported.ADHD <- with(Bulk_diagnoses1, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "8") |
    (!is.na(V2) & V2 == "8") |
    (!is.na(V3) & V3 == "8") |
    (!is.na(V4) & V4 == "8") |
    (!is.na(V5) & V5 == "8") |
    (!is.na(V6) & V6 == "8") |
    (!is.na(V7) & V7 == "8") |
    (!is.na(V8) & V8 == "8") |
    (!is.na(V1) & V1 == "08") |
    (!is.na(V2) & V2 == "08") |
    (!is.na(V3) & V3 == "08") |
    (!is.na(V4) & V4 == "08") |
    (!is.na(V5) & V5 == "08") |
    (!is.na(V6) & V6 == "08") |
    (!is.na(V7) & V7 == "08") |
    (!is.na(V8) & V8 == "08"), 1, 0)
))


diagnoses_vars1 <- c(
  "Self.reported.ADHD", "Self.reported.ASD", "Self.reported.PsychosisAny",
  "Self.reported.Schizophrenia", "Self.reported.PsychosisOther"
)

MHQ <- cbind(MHQ, Bulk_diagnoses1[, c(diagnoses_vars1)])

### 2 - second field with diagnoses

split_list <- strsplit(MHQ$bulkdiagnoses1.3, ",")

n.obs <- sapply(split_list, length)

seq.max <- seq_len(max(n.obs))

Bulk_diagnoses2 <- as.data.frame(t(sapply(split_list, "[", i = seq.max)), stringsAsFactors = F)

Bulk_diagnoses2 <- cbind(Bulk_diagnoses2, MHQ_help_tmp)

str(Bulk_diagnoses2)

Bulk_diagnoses2$Self.reported.Depression <- with(Bulk_diagnoses2, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "1") |
    (!is.na(V2) & V2 == "1") |
    (!is.na(V3) & V3 == "1") |
    (!is.na(V4) & V4 == "1") |
    (!is.na(V5) & V5 == "1") |
    (!is.na(V6) & V6 == "1") |
    (!is.na(V7) & V7 == "1") |
    (!is.na(V8) & V8 == "1") |
    (!is.na(V1) & V1 == "01") |
    (!is.na(V2) & V2 == "01") |
    (!is.na(V3) & V3 == "01") |
    (!is.na(V4) & V4 == "01") |
    (!is.na(V5) & V5 == "01") |
    (!is.na(V6) & V6 == "01") |
    (!is.na(V7) & V7 == "01") |
    (!is.na(V8) & V8 == "01"), 1, 0)
))


Bulk_diagnoses2$Self.reported.Mania.or.bipolar.disorder <- with(Bulk_diagnoses2, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "2") |
    (!is.na(V2) & V2 == "2") |
    (!is.na(V3) & V3 == "2") |
    (!is.na(V4) & V4 == "2") |
    (!is.na(V5) & V5 == "2") |
    (!is.na(V6) & V6 == "2") |
    (!is.na(V7) & V7 == "2") |
    (!is.na(V8) & V8 == "2") |
    (!is.na(V1) & V1 == "02") |
    (!is.na(V2) & V2 == "02") |
    (!is.na(V3) & V3 == "02") |
    (!is.na(V4) & V4 == "02") |
    (!is.na(V5) & V5 == "02") |
    (!is.na(V6) & V6 == "02") |
    (!is.na(V7) & V7 == "02") |
    (!is.na(V8) & V8 == "02"), 1, 0)
))

Bulk_diagnoses2$Self.reported.GAD.and.others <- with(Bulk_diagnoses2, ifelse(is.na(Ever.sought.or.received.professional.help.for.mental.distress), NA,
  ifelse((!is.na(V1) & V1 == "3") |
    (!is.na(V2) & V2 == "3") |
    (!is.na(V3) & V3 == "3") |
    (!is.na(V4) & V4 == "3") |
    (!is.na(V5) & V5 == "3") |
    (!is.na(V6) & V6 == "3") |
    (!is.na(V7) & V7 == "3") |
    (!is.na(V8) & V8 == "3") |
    (!is.na(V1) & V1 == "03") |
    (!is.na(V2) & V2 == "03") |
    (!is.na(V3) & V3 == "03") |
    (!is.na(V4) & V4 == "03") |
    (!is.na(V5) & V5 == "03") |
    (!is.na(V6) & V6 == "03") |
    (!is.na(V7) & V7 == "03") |
    (!is.na(V8) & V8 == "03"), 1, 0)
))


diagnoses_vars2 <- c(
  "Self.reported.GAD.and.others", "Self.reported.Depression", "Self.reported.Mania.or.bipolar.disorder"
)

MHQ <- cbind(MHQ, Bulk_diagnoses2[, c(diagnoses_vars2)])

### describe self-reported diagnoses

for (i in c(diagnoses_vars1, diagnoses_vars2)) {
  print(i)
  print(tabyl(MHQ[[i]]))
}

### some IDs are not unique, duplicated entries

NROW(unique(MHQ$ResultsID))

n_occur <- data.frame(table(MHQ$ResultsID))

duplicated_MHQ <- MHQ[MHQ$ResultsID %in% n_occur$Var1[n_occur$Freq > 1], ]

### keep only earliest response

dim(MHQ)
class(MHQ$Submitted)

MHQ$Submitted <- as.Date(MHQ$Submitted, "%d/%m/%Y")

MHQ <- as.data.frame(MHQ %>%
  group_by(ResultsID) %>%
  arrange(Submitted) %>%
  slice(1))

dim(MHQ)

### subset to individuals with response to question on manic or irritable episode

MHQ <- MHQ[which(MHQ$Ever.extreme.irritability == "1" | MHQ$Ever.manic.or.excitable == "1" | MHQ$Ever.extreme.irritability == "0" | MHQ$Ever.manic.or.excitable == "0"), ]

MHQ <- MHQ[which(!is.na(MHQ$Ever.extreme.irritability) | !is.na(MHQ$Ever.manic.or.excitable)), ]

dim(MHQ)

### subset to individuals with a response to subquestion symptoms also

MHQ <- MHQ[which(!is.na(MHQ$More.active.manic.episode.derived)), ]
dim(MHQ)

tabyl(MHQ$Ever.extreme.irritability)
tabyl(MHQ$Ever.manic.or.excitable)
tabyl(MHQ$More.active.manic.episode.derived)

### subset

MHQ <- as.data.frame(MHQ[c(
  "TecID","ResultsID", "Submitted", c(manic_derived_vars), c(diagnoses_vars1, diagnoses_vars2), "Ever.manic.or.excitable", "Ever.extreme.irritability", "Problematic.mania.or.irritability",
  "Brief.duration.mania.or.irritability", "Moderate.duration.mania.or.irritability", "Extended.duration.mania.or.irritability", "Duration.mania.or.irritability"
)])

### save MHQ

saveRDS(MHQ, file = paste0(output_dir, "/", "manic_data_derived_protect_mhq.rds"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### clean sociodemographic phenotypes

names(SOCIODEM)

### sex

tabyl(SOCIODEM[["X1.2.Are.you.male.or.female."]])

SOCIODEM$Sex <- as.character(SOCIODEM[["X1.2.Are.you.male.or.female."]])

SOCIODEM$Sex[SOCIODEM$Sex == "F"] <- "Female"
SOCIODEM$Sex[SOCIODEM$Sex == "M"] <- "Male"

tabyl(SOCIODEM$Sex)

### education levels

tabyl(SOCIODEM$X1.5.What.is.the.highest.level.of.education.you.have.completed.)

### coding for this var is:
### Secondary Education (GCSE/O-Levels) 1, Post-Secondary Education (College, A-Levels, NVQ3 or below, or similar) 2,
### Vocational Qualification (Diploma, Certificate, BTEC, NVQ 4 and above, or similar) 3,
### Undergraduate Degree (BA, BSc etc.) 4, Post-graduate Degree (MA, MSc etc.) 5, Doctorate (PhD) 6,

SOCIODEM$Education.level <- NA

SOCIODEM[["SOCIODEM$X1.5.What.is.the.highest.level.of.education.you.have.completed."]] <- as.numeric(as.character(SOCIODEM$X1.5.What.is.the.highest.level.of.education.you.have.completed.))

SOCIODEM$Education.level[SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 6 |
  SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 5 |
  SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 4] <- "University.degree"

SOCIODEM$Education.level[SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 2] <- "College.A.Levels.NVQ3.or.below.or.similar"

SOCIODEM$Education.level[SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 1] <- "Secondary.education.GCSE.or.O.Levels"

SOCIODEM$Education.level[SOCIODEM[["X1.5.What.is.the.highest.level.of.education.you.have.completed."]] == 3] <- "Vocational.or.BTEC.or.similar"

tabyl(SOCIODEM$Education.level)

## age - one individual with age 116

tabyl(SOCIODEM$Age)

### smoking

names(LQ)

LQ$Smoking <- NA

tabyl(LQ$X1.1.Do.you.smoke.)

LQ[["X1.1.Do.you.smoke."]] <- as.numeric(as.character(LQ[["X1.1.Do.you.smoke."]]))

tabyl(LQ$X1.1.Do.you.smoke.)

tabyl(LQ$X1.3.Have.you.ever.smoked.)

LQ[["X1.3.Have.you.ever.smoked."]] <- as.numeric(as.character(LQ[["X1.3.Have.you.ever.smoked."]]))

LQ$Smoking[LQ[["X1.1.Do.you.smoke."]] == 1] <- "Current"

LQ$Smoking[LQ[["X1.3.Have.you.ever.smoked."]] == 0 & LQ[["X1.1.Do.you.smoke."]] == 0] <- "Never"

LQ$Smoking[LQ[["X1.3.Have.you.ever.smoked."]] == 1 & LQ[["X1.1.Do.you.smoke."]] == 0] <- "Past"

tabyl(LQ$Smoking)

### alcohol consumption

tabyl(LQ$X2.1.How.often.do.you.normally.have.a.drink.of.something.with.alcohol.in.)

LQ[["X2.1.How.often.do.you.normally.have.a.drink"]] <- as.numeric(as.character(LQ[["X2.1.How.often.do.you.normally.have.a.drink.of.something.with.alcohol.in."]]))

### coding for this var is: at least weekly 3, less than once a week 2, less than once a month 1, never 0

LQ$Alcohol.intake.frequency <- NA

LQ$Alcohol.intake.frequency[LQ[["X2.1.How.often.do.you.normally.have.a.drink"]] == 3] <- "Weekly"

LQ$Alcohol.intake.frequency[LQ[["X2.1.How.often.do.you.normally.have.a.drink"]] == 1 | LQ[["X2.1.How.often.do.you.normally.have.a.drink"]] == 2] <- "Occasionally"

LQ$Alcohol.intake.frequency[LQ[["X2.1.How.often.do.you.normally.have.a.drink"]] == 0] <- "Never"

tabyl(LQ$Alcohol.intake.frequency)

### merge

SOCIODEM <- merge(SOCIODEM, LQ, by = c("ResultsID","TecID"), all.x = T, all.y = T)

dim(SOCIODEM)

SOCIODEM <- SOCIODEM[c("ResultsID", "TecID","Alcohol.intake.frequency", "Sex", "Education.level", "Smoking","Age")]

saveRDS(SOCIODEM, file = paste0(output_dir, "/", "sociodem_data_derived_protect.rds"))

##
