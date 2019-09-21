########################################################
# 
# Project title: "Heterogeneous Effects of Poverty on Cognition"
# Authors: Farbmacher, Koegel, Spindler
#
# This script replicates the results of the analysis in the RAND American Life Panel
# 
#
# Steps in the script:
# 1) Define functions that are used in the analysis
# 2) Prepare the data set for the analysis (the data come from Caravalho et al.'s (2016)
#    Flanker task experiment in the RAND American Life Panel)
# 3) Run the OLS regressions for the Flanker task
# 4) Save all of the results into the folder defined by "setwd" below
#
#
# Notes:
# - We conducted the analysis using
#   - R version 3.4.4
#   - RStudio version 1.1.447
#   - Windows 7
#
########################################################

# Remove all objects from the workspace
rm( list=ls() )


########################################################
# >>> The following two commands are the only commands
# in this script that need to be changed by the user <<<

# Set the directory where you want to save the results
# Example: setwd("C:/replication")
setwd("your_path")

# Define the path where Carvalho et al.'s (2016) data folder "20140481_dataset_&_programs" from
# the AER website is located
# Example: path.data <- "C:/20140481_data/"
path.data <- "your_path/"
########################################################


# Create a folder inside the main folder, which was determined by setwd above
for(subfolder in c("ALP_results")) {
  dir.create(subfolder)
}  

# Load packages
library(lmtest)
library(sandwich)
library(multiwayvcov)
library(foreign)
library(tidyr)
library(stringr)

############################
# 1) Define functions that are needed for the analyses that follow
############################
# Define a function for obtaining HC1 clustered standard errors (Stata's default clustered standard errors),
# along with the number of observations and number of individuals used in the analysis
se_clustered <- function(lm0, df0) {
  regression <- coeftest(lm0, vcov.=cluster.vcov(lm0, df0[,"prim_key"]))
  obs <- nobs(lm0)
  clusters <- length(unique(df0[,"prim_key"]))
  result <- list(regression=regression, obs=obs, clusters=clusters)
  return(result)
}

############################
# 2) Create the data set for our analysis, based on Caravalho et al.'s (2016) American Life Panel data from the AER website
############################

# Prepare the dataset Controls.dta
controls <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Controls.dta"))
# Keep only the individuals from the American Life Panel sample
controls <- controls[which(controls$sample=="ALP"),]
# Keep only the variables which are needed for the subgroup analysis
controls <- controls[,c("prim_key", "age")]

# Prepare the dataset Baseline.dta
baseline <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Baseline.dta"))
# Keep only the individuals from the American Life Panel sample
baseline <- baseline[which(baseline$sample=="ALP"),]
# For each individual, keep only the observation which refers to the payday used for the experiment
baseline <- baseline[which(baseline$experimental_payment==1),]
# Keep only the variables which are needed for the subgroup analysis
baseline <- baseline[, c("prim_key", "total_payamount")]

# Prepare the dataset Followup.dta
followup.orig <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Followup.dta"))
# Keep only the individuals from the American Life Panel sample
followup <- followup.orig[which(followup.orig$sample=="ALP"),]
# Keep only the variables which are needed for the analysis
# Get the names of the 20 flanker_correct variables
flanker.correct.trials <- grep("^flanker_correct",names(followup),value=TRUE)
# Get the names of the 20 flanker_lntime variables
flanker.lntime.trials <- grep("^flanker_lntime", names(followup), value=TRUE)
followup <- followup[,c("prim_key", "treatment", flanker.correct.trials, flanker.lntime.trials)]
# Recode the treatment dummy
followup$treatment <- ifelse(followup$treatment=="Before Payday", 1, 0)
# Reshape the data set into long format
followup.long <- gather(followup, key="var_name", value="value", -c("treatment", "prim_key"))
followup.long$trial <- str_extract(followup.long$var_name, "[0-9]{1,2}$")
followup.long$trial <- as.numeric(followup.long$trial)
followup.long$var_name <- str_replace(followup.long$var_name, "[0-9]{1,2}$", "")
followup.long <- spread(followup.long, key="var_name", value="value")
followup.long <- followup.long[order(followup.long$prim_key, followup.long$trial),]


# Create the raw dataset for our analysis by merging the followup, baseline, and controls data frames
full.alp <- merge(followup.long, controls, by="prim_key", all.x=T)
full.alp <- merge(full.alp, baseline, by="prim_key", all.x=T)

############################
# 3) Run the OLS regressions for the Flanker task
############################

###
# 3a) Full sample
###
# Outcome: correct
df.correct <- full.alp[,c("flanker_correct", "treatment", "trial", "prim_key")]
df.correct <- na.omit(df.correct)
lm.correct <- lm(flanker_correct~treatment+factor(trial), data=df.correct)
lm.correct.clustered <- se_clustered(lm.correct, df.correct)

# Outcome: log response time
df.lntime <- full.alp[,c("flanker_lntime", "treatment", "trial", "prim_key")]
df.lntime <- na.omit(df.lntime)
lm.lntime <- lm(flanker_lntime~treatment+factor(trial), data=df.lntime)
lm.lntime.clustered <- se_clustered(lm.lntime, df.lntime)

###
# 3b) Subgroup: All individuals whose current income is below $750 and whose age is either below 30 or above 70 years
###
# Outcome: correct
df.correct.subgroup <- subset(full.alp, subset= total_payamount<750 & (age<30 | age>70),
                              select=c("flanker_correct", "treatment", "trial", "prim_key", "total_payamount", "age"))
df.correct.subgroup <- na.omit(df.correct.subgroup)
lm.correct.subgroup <- lm(flanker_correct~treatment+factor(trial), data=df.correct.subgroup)
lm.correct.subgroup.clustered <- se_clustered(lm.correct.subgroup, df.correct.subgroup)

# Outcome: lntime
df.lntime.subgroup <- subset(full.alp, subset= total_payamount<750 & (age<30 | age>70),
                              select=c("flanker_lntime", "treatment", "trial", "prim_key", "total_payamount", "age"))
df.lntime.subgroup <- na.omit(df.lntime.subgroup)
lm.lntime.subgroup <- lm(flanker_lntime~treatment+factor(trial), data=df.lntime.subgroup)
lm.lntime.subgroup.clustered <- se_clustered(lm.lntime.subgroup, df.lntime.subgroup)


############################
# 4) Save all of the results
############################
# Full sample OLS results
cat("lm.correct.clustered \n", file="ALP_results/Table_3_flanker_analysis.txt")
capture.output(round(lm.correct.clustered[[1]], digits=3), file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n N=", lm.correct.clustered[[2]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n Number of individuals=", lm.correct.clustered[[3]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)

cat("\n\nlm.lntime.clustered \n", file="ALP_results/Table_3_flanker_analysis.txt", append=T)
capture.output(round(lm.lntime.clustered[[1]], digits=3), file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n N=", lm.lntime.clustered[[2]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n Number of individuals=", lm.lntime.clustered[[3]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)

# Subgroup OLS results
cat("\n\nlm.correct.subgroup.clustered \n", file="ALP_results/Table_3_flanker_analysis.txt", append=T)
capture.output(round(lm.correct.subgroup.clustered[[1]], digits=3), file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n N=", lm.correct.subgroup.clustered[[2]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n Number of individuals=", lm.correct.subgroup.clustered[[3]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)

cat("\n\nlm.lntime.subgroup.clustered \n", file="ALP_results/Table_3_flanker_analysis.txt", append=T)
capture.output(round(lm.lntime.subgroup.clustered[[1]], digits=3), file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n N=", lm.lntime.subgroup.clustered[[2]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)
cat("\n Number of individuals=", lm.lntime.subgroup.clustered[[3]], file="ALP_results/Table_3_flanker_analysis.txt", append=T)


######################################
#### Script end
######################################

