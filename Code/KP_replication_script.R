########################################################
# 
# Project title: "Heterogeneous Effects of Poverty on Cognition"
# Authors: Farbmacher, Koegel, Spindler
#
# This script replicates the results of the analysis in the GfK KnowledgePanel
# 
#
# Steps in the script:
# 1) Define functions that are used in the analysis
# 2) Prepare the data set for the analysis (the data come from Caravalho et al.'s (2016)
#    experiment in the GfK Knowledge Panel)
# 3) Estimate all of the models based on the GfK KnowledgePanel in the paper
# 4) Create the matrices for the heatmaps and typical individuals and estimate the corresponding effects
# 5) Create the descriptive statistics and causal forest results
# 6) Save all of the results into the folder defined by "setwd" below
# 7) Create the graphs in the paper and save them into the folder defined by "setwd" below
#
# Notes:
# - We conducted the analysis using
#   - R version 3.4.4
#   - RStudio version 1.1.447
#   - Windows 7
#   - grf package version 0.9.6
#
# - To install the grf package version 0.9.6, proceed as follows:
#   a) Load the R package devtools
#   b) Run the command: install_version("grf, version="0.9.6")
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


# Create the required folders inside the main folder, which was determined by setwd above
for(subfolder in c("KP_results", "KP_figures")) {
  dir.create(subfolder)
}  

# Load packages
library(sandwich)
library(lmtest)
library(data.table)
library(plyr)
library(dplyr)
library(clubSandwich)
library(mltools)
library(gtools)
library(grf)
library(ggplot2)
library(gridExtra)
library(pryr)
library(parallel)
library(stringr)
library(cowplot)
library(colorspace)
library(fBasics)
library(gtable)
library(grid)
library(quantreg)
library(foreign)

# Asser that grf package version 0.9.6 is installed
stopifnot(packageVersion("grf")=="0.9.6")

############################
# 1) Define functions that are needed for the analyses that follow
############################

###
# 1a) Define functions for finding the most frequently occuring category for each dummy variable and dummy variable class
###

# Marital status
f.marstat<- function(y) {
  y <- y %>%
    mutate(
      mar.stat.max=pmax(married,divorced,widowed,never_married),
      married=ifelse(married==mar.stat.max, 1, 0),
      divorced=ifelse(divorced==mar.stat.max, 1, 0),
      widowed=ifelse(widowed==mar.stat.max, 1, 0),
      never_married=ifelse(never_married==mar.stat.max, 1, 0)
    )
  y
}

# Race
f.race<-function(y) {
  y <- y %>%
    mutate(
      race.max=pmax(white,black,other_race,hispanic),
      white=ifelse(white==race.max,1,0),
      black=ifelse(black==race.max,1,0),
      other_race=ifelse(other_race==race.max,1,0),
      hispanic=ifelse(hispanic==race.max,1,0)
    )
  y
}

# Employment status
f.empstat <- function(y) {
  y <- y %>%
    mutate(
      emp.stat.max=pmax(working2,unemployed2,disabled,retired,other_empstatus),
      working2=ifelse(working2==emp.stat.max,1,0),
      unemployed2=ifelse(unemployed2==emp.stat.max,1,0),
      disabled=ifelse(disabled==emp.stat.max,1,0),
      retired=ifelse(retired==emp.stat.max,1,0),
      other_empstatus=ifelse(other_empstatus==emp.stat.max,1,0)
    )
  y
}

# Education
f.educ <- function(y) {
  y <- y %>%
    mutate(
      educ.max=pmax(college_graduate,some_college,high_school,less_high_school),
      college_graduate=ifelse(college_graduate==educ.max,1,0),
      some_college=ifelse(some_college==educ.max,1,0),
      high_school=ifelse(high_school==educ.max,1,0),
      less_high_school=ifelse(less_high_school==educ.max,1,0)
    )
  y
}

# Income
f.inc <- function(y) {
  y <- y %>%
    mutate(
      inc.max=pmax(income_less_5k, income_btw_5k_10k, income_btw_10k_15k, income_btw_15k_20k, income_btw_20k_25k, income_btw_25k_30k, income_btw_30k_35k, income_btw_35k_40k),
      income_less_5k=ifelse(income_less_5k==inc.max,1,0),
      income_btw_5k_10k=ifelse(income_btw_5k_10k==inc.max,1,0),
      income_btw_10k_15k=ifelse(income_btw_10k_15k==inc.max,1,0),
      income_btw_15k_20k=ifelse(income_btw_15k_20k==inc.max,1,0),
      income_btw_20k_25k=ifelse(income_btw_20k_25k==inc.max,1,0),
      income_btw_25k_30k=ifelse(income_btw_25k_30k==inc.max,1,0),
      income_btw_30k_35k=ifelse(income_btw_30k_35k==inc.max,1,0),
      income_btw_35k_40k=ifelse(income_btw_35k_40k==inc.max,1,0)
    )
  y
}

# Other dummies
f.other <- function(y) {
  y$male <- round(y$male)
  y$metro <- round(y$metro)
  y$hhld_head <- round(y$hhld_head)
  y$children <- round(y$children)
  y
}

# Financial strain variables
f.strain <- function(y) {
  y$livecheckbycheck <- round(y$livecheckbycheck)
  y$caloric_crunch <- round(y$caloric_crunch)
  y$liquidity_constrained <- round(y$liquidity_constrained)
  y$hardship <- round(y$hardship)
  y
}

###
# 1b) Define functions which call the functions from 1a) above to determine the most frequently occuring categories
###

# All dummy class variables (marital status, race, employment status, education, income)
f.class <- function(x0) {
  x0 <- f.marstat(x0)
  x0 <- f.race(x0)
  x0 <- f.empstat(x0)
  x0 <- f.educ(x0)
  x0 <- f.inc(x0)
  x0
}

# All single dummy variables (e.g. metro, male, hardship)
f.dum <- function(x0) {
  x0 <- f.strain(x0)
  x0 <- f.other(x0)
  x0
}


###
# 1c) Define functions that create dummies indicating if there are tied categories (equally frequently occuring characteristics)
# in the dummy classes
###

# Marital status
ties.marstat <- function(y) {
  marstat.names <- c("married", "divorced", "widowed", "never_married")
  marstat.names.mean <- paste0(marstat.names,".mean")
  # If more than one marital status dummy has the same mean as the maximum martital status dummy mean, then set tie.marstat=1
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,marstat.names.mean]==y[i0,"mar.stat.max"])
    if(length(tie.cats)>1) {
      y[i0, "tie.marstat"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.marstat"] <- 0
    }
  }
  # Assert that only one of the marital status dummies is set to one (among the rows that do not have missing marital status
  # information, as captured by !is.na(y$married), and that do not have ties)
  stopifnot(rowSums(y[which(!is.na(y$married) & y$tie.marstat==0),marstat.names])==1)
  y
}

# Race
ties.race <- function(y) {
  race.names <- c("white", "black", "hispanic", "other_race")
  race.names.mean <- paste0(race.names,".mean")
  # If more than one race status dummy has the same mean as the maximum race status dummy mean, then set tie.race=1
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,race.names.mean]==y[i0,"race.max"])
    if(length(tie.cats)>1) {
      y[i0, "tie.race"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.race"] <- 0
    }
  }
  # Assert that only one of the race dummies is set to one (among the rows that do not have missing race information, as captured
  # by !is.na(y$white), and that do not have ties)
  stopifnot(rowSums(y[which(!is.na(y$white) & y$tie.race==0),race.names])==1)
  y
}

# Employment status
ties.empstat <- function(y) {
  empstat.names <- c("working2", "unemployed2", "disabled", "retired", "other_empstatus")
  empstat.names.mean <- paste0(empstat.names,".mean")
  # If more than one employment status dummy has the same mean as the maximum employment status dummy mean, then set tie.empstat=1
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,empstat.names.mean]==y[i0,"emp.stat.max"])
    if(length(tie.cats)>1) {
      y[i0, "tie.empstat"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.empstat"] <- 0
    }
  }
  # Assert that only one of the employment status dummies is set to one (among the rows that do not have missing
  # employment status information, as captured by !is.na(y$working2), and that do not have ties)
  stopifnot(rowSums(y[which(!is.na(y$working2) & y$tie.empstat==0),empstat.names])==1)
  y
}


###
# 1d) Define functions that solve tied categories in ordinal dummy variable classes (and household size)
# and that create dummies indicating if there are tied categories in the classes 
# (the characteristics/classes are education, income, household size)
###

# Functions for the matrix based on which we perform the predictions below

# Education
fsolt.educ <- function(y) {
  educ.sorted <- c("less_high_school", "high_school", "some_college", "college_graduate")
  educ.sorted.mean <- paste0(educ.sorted,".mean")
  for(i0 in 1:nrow(y)) {
    # Set the first dummy category whose mean value equals the maximum mean value over all dummies within that class to 1
    # (Because the dummies are ordered in increasing order, the lowest category in case of ties is the first. Thus, the lowest
    # category is set to 1 in case of ties. which.max(x) gives the index of the first maximum, i.e., the first 1)
    y[i0,educ.sorted][which.max(y[i0,educ.sorted.mean]==y[i0,"educ.max"])] <- 1
    # Set all other dummies to zero
    y[i0,educ.sorted][-which.max(y[i0,educ.sorted.mean]==y[i0,"educ.max"])] <- 0
    # Create a dummy indicating if there are ties in the education categories
    if(length(which(y[i0,educ.sorted.mean]==y[i0,"educ.max"]))>1) {
      y[i0,"tie.educ"] <- 1
    }
    if(length(which(y[i0,educ.sorted.mean]==y[i0,"educ.max"]))==1) {
      y[i0,"tie.educ"] <- 0
    }
  }
  stopifnot(rowSums(y[!is.na(y$less_high_school),educ.sorted])==1)
  y
}

# Income
fsolt.inc <- function(y) {
  inc.sorted <- c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k", "income_btw_20k_25k",
                  "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")
  inc.sorted.mean <- paste0(inc.sorted,".mean")
  for(i0 in 1:nrow(y)) {
    # Set the first dummy category whose mean value equals the maximum mean value over all dummies within that class to 1
    # (Because the dummies are ordered in increasing order, the lowest category in case of ties is the first. Thus, the lowest
    # category is set to 1 in case of ties. which.max(x) gives the index of the first maximum, i.e., the first 1)
    y[i0,inc.sorted][which.max(y[i0,inc.sorted.mean]==y[i0,"inc.max"])] <- 1
    # Set all other dummies to zero
    y[i0,inc.sorted][-which.max(y[i0,inc.sorted.mean]==y[i0,"inc.max"])] <- 0
    # Create a dummy indicating if there are ties in the income categories
    if(length(which(y[i0,inc.sorted.mean]==y[i0,"inc.max"]))>1) {
      y[i0,"tie.inc"] <- 1
    }
    if(length(which(y[i0,inc.sorted.mean]==y[i0,"inc.max"]))==1) {
      y[i0,"tie.inc"] <- 0
    }
  }
  stopifnot(rowSums(y[!is.na(y$income_less_5k),inc.sorted])==1)
  y
}

# Household size
# If there are tie categories in the household size in the sense that the household size is between two integer values
# (which yields a non-integer household size, e.g., median(c(3,3,4,4,))=3.5), then set the household size to the integer
# part of the non-integer household size (i.e., choose the "lower category")
fsolt.hhld_size <- function(y) {
  y["tie.hhldsize"] <- ifelse(y["hhld_size"]%%1!=0, 1, 0)
  y[,"hhld_size"] <- floor(y[,"hhld_size"])
  y
}


# Define functions that set the dummies for the categories in which there are ties to NA (in the extended matrix defined below, based
# on which we solve the ties in the non-ordinal characteristics marital status, race, employment status; we use these functions so
# it becomes clear when estimating effects that there are still unresolved ties)

# Marital status
fsolt.marstat.ext <- function(y) {
  marstat.names <- c("married", "divorced", "widowed", "never_married")
  marstat.names.mean <- paste0(marstat.names,".mean")
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,marstat.names.mean]==y[i0,"mar.stat.max"])
    if(length(tie.cats)>1) {
      y[i0,marstat.names] <- NA
      y[i0, "tie.marstat"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.marstat"] <- 0
    }
  }
  # Assert that only one of the marital status dummies is set to 1 if the martital status is not missing
  stopifnot(rowSums(y[!is.na(y$married),marstat.names])==1)
  y
}

# Race
fsolt.race.ext <- function(y) {
  race.names <- c("white", "black", "hispanic", "other_race")
  race.names.mean <- paste0(race.names,".mean")
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,race.names.mean]==y[i0,"race.max"])
    if(length(tie.cats)>1) {
      y[i0,race.names] <- NA
      y[i0, "tie.race"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.race"] <- 0
    }
  }
  # Assert that only one of the marital status dummies is set to 1 if the martital status is not missing
  stopifnot(rowSums(y[!is.na(y$white),race.names])==1)
  y
}

# Employment status
fsolt.empstat.ext <- function(y) {
  empstat.names <- c("working2", "unemployed2", "disabled", "retired", "other_empstatus")
  empstat.names.mean <- paste0(empstat.names,".mean")
  for(i0 in 1:nrow(y)) {
    tie.cats <- which(y[i0,empstat.names.mean]==y[i0,"emp.stat.max"])
    if(length(tie.cats)>1) {
      y[i0,empstat.names] <- NA
      y[i0, "tie.empstat"] <- 1
    }
    if(length(tie.cats)==1) {
      y[i0, "tie.empstat"] <- 0
    }
  }
  # Assert that only one of the marital status dummies is set to 1 if the martital status is not missing
  stopifnot(rowSums(y[!is.na(y$working2),empstat.names])==1)
  y
}


###
# 1e) Define a function for computing the means and medians for all variables in a given window centered at a given
# age-current income combination
# window.age0: width of the window in the age direction
# window.tpa0: width of the window in the current income direction
# The function f.typical calls this function
###

f.roll <- function(window.age0, window.tpa0, fun0) {
  # Define the window size (in the age and current income direction)
  win.age0 <- (window.age0-1)/2
  win.tpa0 <- (window.tpa0-1)/2
  # Select all covariates that we use in the analysis. vars.full is defined in the function f.pred.custom below
  df.tmp0 <- full.kp.mod[,vars.full]
  
  # Define a function that selects all observations in a given window and performs calculations based on them
  f.select <- function(age0, tpa0, win.age0, win.tpa0, fun0){
    df.tmp1 <- df.tmp0[df.tmp0$age>=age0-win.age0 & df.tmp0$age<=age0+win.age0 &
                         df.tmp0$total_payamount>=tpa0-win.tpa0 & df.tmp0$total_payamount<=tpa0+win.tpa0,]
    # Define a variable which gives the number of observations in a given window
    df.tmp1 <- df.tmp1 %>% mutate(n=n())
    # Calculate some other measure based on the observations in a given window (we calculate median and mean using this
    # function below in f.typical), and add to the data frame the age and current income where the window is centered at
    df.tmp1 <- df.tmp1 %>% summarise_all(fun0) %>% mutate(age=age0, total_payamount=tpa0)
  }
  
  # Calculate the measure of interest (for us, median and mean below) for all current income-age combinations given by seq.tpa, seq.age
  seq.tpa <<- seq(400,500, by=25)
  seq.age <<- c(seq(18,22, by=1),seq(73,77,by=1))
  # For a given current income ltpa0, loop over all ages lage0. Do the looping for all current incomes.
  list0<- lapply(seq.tpa, function(ltpa0) {
    tmp <-lapply(seq.age, function(lage0) f.select(lage0, ltpa0, win.age0, win.tpa0, fun0))
    # For a given current income and all ages, create the data frame tmp. Put all the data frames for all current incomes
    # into the list list0
    do.call(rbind,tmp) })
  # Create a data frame which contains the measure of interest for all current income-age combinations
  df0 <- do.call(rbind, list0)
}


###
# 1f) Define a function which creates a data frame for the typical individual estimations
# agewin1, tpawin1 refer to the age-current income window that we perferably use in our analysis
# agewin2, tpawin2 refer to the extended age-current income window that we use to solve ties in the non-ordinal characteristics
###

f.typical <- function(agewin1, tpawin1, agewin2, tpawin2) {
  # Create two data frames which contain all mean and median values, respectively, for each age-current income combination, using
  # the function f.roll defined immediately above
  x.mean.spec <- f.roll(agewin1, tpawin1, mean)
  x.median.spec <- f.roll(agewin1, tpawin1, median)
  # Assert that the order of the age-current income combinations in the mean and median data frame is identical
  stopifnot(identical(x.mean.spec[c("age","total_payamount")],x.median.spec[c("age","total_payamount")]))
  
  # For the extended age-current income window: Create a data frame which contains all mean values
  # for each age-current income combination
  x.mean.spec.ext <- f.roll(agewin2, tpawin2, mean)

  # Combine the mean and median matrices
  # (for the dummy variables, use the mean matrix values and for the continuous variables use the median matrix values)
  # For the main window size (agewin1, tpawin1)
  x.spec <- x.mean.spec
  x.spec$payamount_fraction <- x.median.spec$payamount_fraction
  x.spec$hhld_size <- x.median.spec$hhld_size
  
  # For the extended window size (agewin2, tpawin2)
  x.spec.ext <- x.mean.spec.ext
  
  # Keep the original mean value variables for the dummies
  # For the main window size (agewin1, tpawin1)
  keep <- c("married", "divorced", "widowed", "never_married", "white", "black", "hispanic", "other_race", "working2", "unemployed2",
            "disabled", "retired", "other_empstatus", "college_graduate", "some_college", "high_school", "less_high_school", "income_less_5k",
            "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k", "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k",
            "income_btw_35k_40k", "livecheckbycheck", "caloric_crunch", "liquidity_constrained", "hardship",
            "male", "children", "metro", "hhld_head")
  keep.mean <- paste0(keep, ".mean")
  x.spec[,keep.mean] <- x.spec[,keep]
  
  # For the extended window size (agewin2, tpawin2)
  keep.ext <- c("married", "divorced", "widowed", "never_married", "white", "black", "hispanic", "other_race", "working2", "unemployed2",
            "disabled", "retired", "other_empstatus", "livecheckbycheck", "caloric_crunch", "liquidity_constrained", "hardship",
            "male", "children", "metro", "hhld_head")
  keep.mean.ext <- paste0(keep.ext, ".mean")
  x.spec.ext[,keep.mean.ext] <- x.spec.ext[,keep.ext]
  
  # Set the variables' values according to the most frequently occuring values within each age-current income window
  # For the main window size (agewin1, tpawin1)
  # Financial strain dummies, and male, metro, children, household head
  x.spec <- f.dum(x.spec)
  # Marital status, race, employment status, education, income
  x.spec <- f.class(x.spec)
  
  # For the extended window size (agewin2, tpawin2)
  # Financial strain dummies, and male, metro, children, household head
  x.spec.ext <- f.dum(x.spec.ext)
  # Marital status, race, employment status, education, income
  x.spec.ext <- f.class(x.spec.ext)
  
  # Solve the ties for the ordinal characteristics by selecting the lowest tie category for each class in a given window
  # (Do this for the classes/characteristics education, income, and also household size; the functions additionally
  # create dummies indicating if there are tie categories)
  x.spec <- fsolt.educ(x.spec)
  x.spec <- fsolt.inc(x.spec)
  x.spec <- fsolt.hhld_size(x.spec)
  
  # Create dummies indicating if there are ties in the non-ordinal dummy classes (marital status, race, employment status)
  x.spec <- ties.marstat(x.spec)
  x.spec <- ties.race(x.spec)
  x.spec <- ties.empstat(x.spec)
  
  # In the extended window size data frame (agewin2, tpawin2): Set dummies referring to a tied dummy class to missing;
  # and create a dummy indicating if there are ties in a respective class
  x.spec.ext <- fsolt.marstat.ext(x.spec.ext)
  x.spec.ext <- fsolt.race.ext(x.spec.ext)
  x.spec.ext <- fsolt.empstat.ext(x.spec.ext)
  
  # Add the tie dummy indicators to the x.spec data frame created above
  # Main window size data frame
  x.spec$tie.n <- rowSums(x.spec[, grep("^tie.", names(x.spec))])
  # Extended window size data frame
  x.spec.ext$tie.n <- rowSums(x.spec.ext[, grep("^tie.", names(x.spec.ext))])
  
  # Replace the ties in the marital status, race and employment status characteristics in the main window size data frame
  # by the extended window size data frame values
  stopifnot(x.spec[,c("age", "total_payamount")]==x.spec.ext[,c("age","total_payamount")])
  # Marital status
  x.spec[which(x.spec$tie.marstat==1),c("married", "divorced", "widowed", "never_married")] <-
    x.spec.ext[which(x.spec$tie.marstat==1),c("married", "divorced", "widowed", "never_married")]
  # Employment status
  x.spec[which(x.spec$tie.empstat==1),c("working2", "unemployed2", "disabled", "retired", "other_empstatus")] <-
    x.spec.ext[which(x.spec$tie.empstat==1),c("working2", "unemployed2", "disabled", "retired", "other_empstatus")]
  # Race
  x.spec[which(x.spec$tie.race==1),c("white", "black", "hispanic", "other_race")] <-
    x.spec.ext[which(x.spec$tie.race==1),c("white", "black", "hispanic", "other_race")]
  
  # Solve the ties in the non-ordinal single dummy variables by replacing the tie values by the corresponding extended data frame values
  # (If there are also ties in the extended matrix, then set the respective tie values to zero. This is implicitly done, because
  # R's rounding function rounds 0.5 to 0.)
  # Liquidity constrained
   x.spec[which(x.spec$liquidity_constrained.mean==0.5), "liquidity_constrained"] <-
    x.spec.ext[which(x.spec$liquidity_constrained.mean==0.5), "liquidity_constrained"]
  # Caloric crunch
  x.spec[which(x.spec$caloric_crunch.mean==0.5), "caloric_crunch"] <-
    x.spec.ext[which(x.spec$caloric_crunch.mean==0.5), "caloric_crunch"]
  # Live check by check
  x.spec[which(x.spec$livecheckbycheck.mean==0.5), "livecheckbycheck"] <-
    x.spec.ext[which(x.spec$livecheckbycheck.mean==0.5), "livecheckbycheck"]
  # Financial hardship
  x.spec[which(x.spec$hardship.mean==0.5), "hardship"] <-
    x.spec.ext[which(x.spec$hardship.mean==0.5), "hardship"]
  # Male
  x.spec[which(x.spec$male.mean==0.5), "male"] <-
    x.spec.ext[which(x.spec$male.mean==0.5), "male"]
  # Children
  x.spec[which(x.spec$children.mean==0.5), "children"] <-
    x.spec.ext[which(x.spec$children.mean==0.5), "children"]
  # Metro
  x.spec[which(x.spec$metro.mean==0.5), "metro"] <-
    x.spec.ext[which(x.spec$metro.mean==0.5), "metro"]
  # Household head
  x.spec[which(x.spec$hhld_head.mean==0.5), "hhld_head"] <-
    x.spec.ext[which(x.spec$hhld_head.mean==0.5), "hhld_head"]
  
  # Create a matrix for each row in the x.spec matrix where you vary all of the variables sequentially
  # Reminder: the data frame x.spec contains all age-current income combinations along with the corresponding
  # most frequent background variables
  cats<- c("baseline",
           "livecheckbycheck1",
           "livecheckbycheck0",
           "caloric_crunch1",
           "caloric_crunch0",
           "liquidity_constrained1",
           "liquidity_constrained0",
           "hardship1",
           "hardship0",
           "male1",
           "male0",
           "married1",
           "divorced1",
           "widowed1",
           "never_married1",
           "white1",
           "black1",
           "hispanic1",
           "other_race1",
           "working21",
           "unemployed21",
           "disabled1",
           "retired1",
           "other_empstatus1",
           "college_graduate1",
           "some_college1",
           "high_school1",
           "less_high_school1",
           "income_less_5k1",
           "income_btw_5k_10k1",
           "income_btw_10k_15k1",
           "income_btw_15k_20k1",
           "income_btw_20k_25k1",
           "income_btw_25k_30k1",
           "income_btw_30k_35k1",
           "income_btw_35k_40k1",
           "hhld_head1",
           "hhld_head0",
           "children1",
           "children0",
           "metro1",
           "metro0",
           "hhld_size1",
           "hhld_size2",
           "hhld_size3",
           "hhld_size4",
           "hhld_size5",
           "payamount_fraction0.25",
           "payamount_fraction0.5",
           "payamount_fraction0.75",
           "payamount_fraction1"
  )
  # Split the cats vector immediately above into the name part and the value part (the end of each string)
  cats.names <- str_replace(str_replace(cats, "[0-9]$", ""), "0\\.[0-9]*$" ,"")
  cats.vals <- as.numeric(unlist(lapply(1:length(cats.names), function(i0) str_replace(cats[i0], cats.names[i0], ""))))
  
  # Create a separate data frame for each age-current income combination, which contains age, current income and the
  # cats vector named change.var
  x.list <- lapply(unique(x.spec$total_payamount), function(x0) {
    tmp <- lapply(unique(x.spec$age), function(y0) data.frame(change.var=factor(cats, levels=rev(cats)), age=y0, total_payamount=x0))
  })
  x.list <- unlist(x.list, recursive=F)
  
  # Name each data frame in the list according to the age-current income combination that it includes
  names(x.list) <- paste0(rep(paste0("tpa",as.character(unique(x.spec$total_payamount))), each=length(unique(x.spec$age))),"age",
                          as.character(unique(x.spec$age)))
  names.x.list <- names(x.list)
  
  # Add all the covariates from the main window size data frame x.spec (which contains all the covariates with their
  # values set according to our variables setting procedure)
  x.list[names(x.list)] <-
    lapply(1:length(x.list), function(k0) merge(x.list[[k0]],x.spec[,c(vars,"n",grep("\\.max$", names(x.spec), val=T), keep.mean, grep("^tie.", names(x.spec), val=T))], by=c("age", "total_payamount"), all.x=T))
  
  # Switch the covariates sequentially on and off below the typical individual baseline covariate combination (in the first
  # row of each data frame)
  for(i0 in names(x.list)) {
    for(k0 in 2:length(cats.names)) {
      # Set the value of a given row to NA if the row would change a variable value to the value to which the variable is
      # already set for the typical individual ("Create the empty rows in the typical individual plots")
      x.list[[i0]][which(x.list[[i0]]$change.var==cats[k0] & x.list[[i0]][,cats.names[k0]]==cats.vals[k0]),cats.names[k0]] <- NA
      # Vary the variables sequentially (switch all of them on and off relative to the baseline covariates specification)
      # (Only switch a variable on if it is not already switched on in the baseline covariates specification)
      x.list[[i0]][which(x.list[[i0]]$change.var==cats[k0] & x.list[[i0]][,cats.names[k0]]!=cats.vals[k0] &
                           !is.na(x.list[[i0]][,cats.names[k0]])),cats.names[k0]] <- cats.vals[k0]
    }
    # Switch the covariate within a dummy class off that is switched on in the baseline specification if another covariate in the
    # same category has been switch on (Example: if married is switched on in the baseline covariate specification, then switch
    # married off if divorced has been switch on; this procedure is only necessary for the dummies that all refer to a
    # common characteristic, such as marital status)
    # Steps: For each characteristic with multiple dummies: 1) select the dummy which is switched on in the baseline covariate specification,
    # 2) Select the rows where another dummy from the same characteristic is switched on (for these rows, the row sum of the relevant
    # dummies equals two, 3) replace the dummy that is switched on in the baseline covariate specification in these rows by zero
    
    # If all of the variables in the rowSums expression are not missing in the first row of a given data frame i0, then continue
    if(rowSums(is.na(x.list[[i0]][1,c("married", "divorced", "widowed", "never_married", "white", "black", "other_race", "hispanic",
                               "working2", "unemployed2", "disabled", "retired", "other_empstatus",
                               "college_graduate", "some_college", "high_school", "less_high_school",
                               "income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                               "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")]))==0) {
    # Marital status
    marstat.baseline <- names(x.list[[i0]][1,c("married", "divorced", "widowed", "never_married")])[x.list[[i0]][1,c("married", "divorced", "widowed", "never_married")]==1]
    x.list[[i0]][which(rowSums(x.list[[i0]][,c("married", "divorced", "widowed", "never_married")])==2),marstat.baseline] <- 0
    # Race
    race.baseline <- names(x.list[[i0]][1,c("white", "black", "other_race", "hispanic")])[x.list[[i0]][1,c("white", "black", "other_race", "hispanic")]==1]
    x.list[[i0]][which(rowSums(x.list[[i0]][,c("white", "black", "other_race", "hispanic")])==2),race.baseline] <- 0
    # Employment status
    empstat.baseline <- names(x.list[[i0]][1,c("working2", "unemployed2", "disabled", "retired", "other_empstatus")])[x.list[[i0]][1,c("working2", "unemployed2", "disabled", "retired", "other_empstatus")]==1]
    x.list[[i0]][which(rowSums(x.list[[i0]][,c("working2", "unemployed2", "disabled", "retired", "other_empstatus")])==2),empstat.baseline] <- 0
    # Education
    educ.baseline <- names(x.list[[i0]][1,c("college_graduate", "some_college", "high_school", "less_high_school")])[x.list[[i0]][1,c("college_graduate", "some_college", "high_school", "less_high_school")]==1]
    x.list[[i0]][which(rowSums(x.list[[i0]][,c("college_graduate", "some_college", "high_school", "less_high_school")])==2),educ.baseline] <- 0
    # Annual household income
    inc.baseline <- names(x.list[[i0]][1,c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                                          "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")])[x.list[[i0]][1,c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                                                                                                                                                    "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")]==1]
    x.list[[i0]][which(rowSums(x.list[[i0]][,c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                                             "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")])==2),inc.baseline] <- 0
    
    
    # Assert that only one of the dummies for each dummy class is set equal to 1 (do not make this assertion for the row in which the
    # dummy that is switched on in the baseline specification would have been set to 1; this row contains a missing, as defined above)
    stopifnot(rowSums(x.list[[i0]][!is.na(x.list[[i0]][,marstat.baseline]),c("married", "divorced", "widowed", "never_married")])==1)
    stopifnot(rowSums(x.list[[i0]][!is.na(x.list[[i0]][,race.baseline]),c("white", "black", "other_race", "hispanic")])==1)
    stopifnot(rowSums(x.list[[i0]][!is.na(x.list[[i0]][,empstat.baseline]),c("working2", "unemployed2", "disabled", "retired", "other_empstatus")])==1)
    stopifnot(rowSums(x.list[[i0]][!is.na(x.list[[i0]][,educ.baseline]),c("college_graduate", "some_college", "high_school", "less_high_school")])==1)
    stopifnot(rowSums(x.list[[i0]][!is.na(x.list[[i0]][,inc.baseline]),c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                                                                        "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")])==1)
    
    
    # Assert that all of the dummies for a given dummy class, which are not set to 1 in the baseline specification, are set equal to 0
    # if the row is the row where the dummy which is set to 1 in the baseline specification would have been set to 1
    stopifnot(rowSums(x.list[[i0]][is.na(x.list[[i0]][,marstat.baseline]),c("married", "divorced", "widowed", "never_married")], na.rm=T)==0)
    stopifnot(rowSums(x.list[[i0]][is.na(x.list[[i0]][,race.baseline]),c("white", "black", "other_race", "hispanic")], na.rm=T)==0)
    stopifnot(rowSums(x.list[[i0]][is.na(x.list[[i0]][,empstat.baseline]),c("working2", "unemployed2", "disabled", "retired", "other_empstatus")], na.rm=T)==0)
    stopifnot(rowSums(x.list[[i0]][is.na(x.list[[i0]][,educ.baseline]),c("college_graduate", "some_college", "high_school", "less_high_school")], na.rm=T)==0)
    stopifnot(rowSums(x.list[[i0]][is.na(x.list[[i0]][,inc.baseline]),c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                                                                       "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")], na.rm=T)==0)
    }
  }
  x.list
}


###
# 1g) Define a function for estimating our causal forest models
# The function returns a list, which contains three causal forest objects: one object for each of our three outcomes
# The function sets all parameters apart from the ones specified in the function's arguments to their default values
###
f.est.full <- function(num.trees0, sample.fraction0, min.node.size0, seed0, mtry0) {
  # Select the 37 covariates we use in our analysis
  vars<-c(
    "livecheckbycheck",
    "caloric_crunch",
    "liquidity_constrained",
    "hardship",
    "total_payamount",
    "payamount_fraction",
    "age",
    "male",
    "married",
    "divorced",
    "widowed",
    "never_married",
    "white",
    "black",
    "hispanic",
    "other_race",
    "working2",
    "unemployed2",
    "disabled",
    "retired",
    "other_empstatus",
    "hhld_size",
    "college_graduate",
    "some_college",
    "high_school",
    "less_high_school",
    "income_less_5k",
    "income_btw_5k_10k",
    "income_btw_10k_15k",
    "income_btw_15k_20k",
    "income_btw_20k_25k",
    "income_btw_25k_30k",
    "income_btw_30k_35k",
    "income_btw_35k_40k",
    "hhld_head",
    "children",
    "metro")

  # Put the names of the variables that we use in the analysis in a vector in the global environment,
  # because other functions also use this vector
  vars.full <<- vars
  
  # Create the matrix/vectors required for the estimation
  # Matrix containing the regressors
  covariates <<- as.matrix(subset(full.kp.mod, select= vars))
  # Vector containing the treatment indicator
  treatment <<- full.kp.mod$treatment
  # Vectors containing the outcomes
  outcome.ncorrect <- full.kp.mod$ncorrect
  outcome.rt <- full.kp.mod$rt
  outcome.ncor.rt <- full.kp.mod$ncor.rt

  
  # Estimation using the ncorrect outcome (number of correct answers)
  cforest.ncorrect <<- causal_forest(covariates,outcome.ncorrect,treatment,
                                   num.trees=num.trees0, sample.fraction = sample.fraction0, honesty=TRUE,
                                   min.node.size = min.node.size0, seed=seed0, mtry=mtry0)
  
  # Estimation using the rt outcome (total response time)
  cforest.rt <<- causal_forest(covariates,outcome.rt,treatment,
                             num.trees=num.trees0, sample.fraction = sample.fraction0, honesty=TRUE,
                             min.node.size = min.node.size0, seed=seed0, mtry=mtry0)
  
  # Estimation using the ncor.rt outcome (number of correct answers per second)
  cforest.ncor.rt <<- causal_forest(covariates,outcome.ncor.rt,treatment,
                                  num.trees=num.trees0, sample.fraction = sample.fraction0, honesty=TRUE,
                                  min.node.size = min.node.size0, seed=seed0, mtry=mtry0)
  
  # Put the causal forest estimation objects in a list
  cf.list.internal <- list(cforest.ncorrect, cforest.rt, cforest.ncor.rt)
  cf.list.internal
}


###
# 1h) Define a function which gives a convenient matrix representation of the variable importance measure for
# our three causal forest estimations
# The function uses the default parameters for the variable importance function
###
fvi <- function(cforest.ncorrect0, cforest.rt0, cforest.ncor.rt0) {
  # Calculate the variable importance measure for the outcome ncorrect
  vi.ncorrect.tmp <- variable_importance(cforest.ncorrect0)
  # Create an object that includes the covariate names in addition to the variable importance measure
  vi.ncorrect <- cbind(vi.ncorrect.tmp[order(vi.ncorrect.tmp)],attr(cforest.ncorrect0[["X.orig"]],"dimnames")[[2]][order(vi.ncorrect.tmp)])
  
  # Calculate the variable importance measure for the outcome rt
  vi.rt.tmp <- variable_importance(cforest.rt0)
  # Create an object that includes the covariate names in addition to the variable importance measure
  vi.rt <- cbind(vi.rt.tmp[order(vi.rt.tmp)],attr(cforest.rt0[["X.orig"]],"dimnames")[[2]][order(vi.rt.tmp)])
  
  # Calculate the variable importance measure for the outcome ncor.rt
  vi.ncor.rt.tmp <-variable_importance(cforest.ncor.rt0)
  # Create an object that includes the covariate names in addition to the variable importance measure
  vi.ncor.rt <- cbind(vi.ncor.rt.tmp[order(vi.ncor.rt.tmp)],attr(cforest.ncor.rt0[["X.orig"]],"dimnames")[[2]][order(vi.ncor.rt.tmp)])
  
  vi.internal <- matrix(c(vi.ncorrect, vi.rt, vi.ncor.rt),ncol=6, dimnames = list(NULL,c("vi.ncorrect.val","vi.ncorrect.var","vi.rt.val","vi.rt.var","vi.ncor.rt.val","vi.ncor.rt.var")))
  vi.internal
}

###
# 1g) Add the causal forest based augmented inverse propensity weighted estimator, which includes the
# subsetting feature. The subsetting features allows to estimate average treatment effects for subgroups of individuals
# (This function has been available in the grf github repository since August 2018, and is now implemented in the
# grf package version 0.10.1. The grf package version 0.9.6 includes the same function, except that it does not include
# the subsetting feature.)
###
average_treatment_effect_subset <- function (forest, target.sample = c("all", "treated", "control", 
                                    "overlap"), method = c("AIPW", "TMLE"), subset = NULL) 
{
  target.sample <- match.arg(target.sample)
  method <- match.arg(method)
  cluster.se <- length(forest$clusters) > 0
  if (!("causal_forest" %in% class(forest))) {
    stop("Average effect estimation only implemented for causal_forest")
  }
  if (is.null(subset)) {
    subset <- 1:length(forest$Y.hat)
  }
  if (class(subset) == "logical" & length(subset) == length(forest$Y.hat)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:length(forest$Y.hat))) {
    stop(paste("If specified, subset must be a vector contained in 1:n,", 
               "or a boolean vector of length n."))
  }
  subset.W.orig <- forest$W.orig[subset]
  subset.W.hat <- forest$W.hat[subset]
  subset.Y.orig <- forest$Y.orig[subset]
  subset.Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- predict(forest)$predictions[subset]
  if (length(forest$clusters) == 0) {
    subset.clusters <- numeric(0)
  }
  else {
    subset.clusters <- forest$clusters[subset]
  }
  if (target.sample == "overlap") {
    W.residual <- subset.W.orig - subset.W.hat
    Y.residual <- subset.Y.orig - subset.Y.hat
    tau.ols <- lm(Y.residual ~ W.residual)
    tau.est <- coef(summary(tau.ols))[2, 1]
    if (cluster.se) {
      tau.se <- sqrt(sandwich::vcovCL(tau.ols, cluster = subset.clusters)[2, 
                                                                          2])
    }
    else {
      tau.se <- sqrt(sandwich::vcovHC(tau.ols)[2, 2])
    }
    return(c(estimate = tau.est, std.err = tau.se))
  }
  if (!all(subset.W.orig %in% c(0, 1))) {
    stop(paste("Average treatment effect estimation only implemented for binary treatment.", 
               "See `average_partial_effect` for continuous W."))
  }
  if (min(subset.W.hat) <= 0.01 && max(subset.W.hat) >= 0.99) {
    rng = range(subset.W.hat)
    warning(paste0("Estimated treatment propensities take values between ", 
                   round(rng[1], 3), " and ", round(rng[2], 3), " and in particular get very close to 0 and 1. ", 
                   "In this case, using `target.sample=overlap`, or filtering data as in ", 
                   "Crump, Hotz, Imbens, and Mitnik (Biometrika, 2009) may be helpful."))
  }
  else if (min(subset.W.hat) <= 0.01 && target.sample != "treated") {
    warning(paste0("Estimated treatment propensities go as low as ", 
                   round(min(subset.W.hat), 3), " which means that treatment ", 
                   "effects for some controls may not be well identified. ", 
                   "In this case, using `target.sample=treated` may be helpful."))
  }
  else if (max(subset.W.hat) >= 0.99 && target.sample != "control") {
    warning(paste0("Estimated treatment propensities go as high as ", 
                   round(max(subset.W.hat), 3), " which means that treatment ", 
                   "effects for some treated units may not be well identified. ", 
                   "In this case, using `target.sample=control` may be helpful."))
  }
  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)
  if (target.sample == "all") {
    tau.avg.raw <- mean(tau.hat.pointwise)
  }
  else if (target.sample == "treated") {
    tau.avg.raw <- mean(tau.hat.pointwise[treated.idx])
  }
  else if (target.sample == "control") {
    tau.avg.raw <- mean(tau.hat.pointwise[control.idx])
  }
  else {
    stop("Invalid target sample.")
  }
  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise
  if (method == "TMLE") {
    loaded <- requireNamespace("sandwich", quietly = TRUE)
    if (!loaded) {
      warning("To use TMLE, please install the package `sandwich`. Using AIPW instead.")
      method = "AIPW"
    }
  }
  if (method == "AIPW") {
    if (target.sample == "all") {
      gamma.control.raw <- 1/(1 - subset.W.hat[control.idx])
      gamma.treated.raw <- 1/subset.W.hat[treated.idx]
    }
    else if (target.sample == "treated") {
      gamma.control.raw <- subset.W.hat[control.idx]/(1 - 
                                                        subset.W.hat[control.idx])
      gamma.treated.raw <- rep(1, length(treated.idx))
    }
    else if (target.sample == "control") {
      gamma.control.raw <- rep(1, length(control.idx))
      gamma.treated.raw <- (1 - subset.W.hat[treated.idx])/subset.W.hat[treated.idx]
    }
    else {
      stop("Invalid target sample.")
    }
    gamma <- rep(0, length(subset.W.orig))
    gamma[control.idx] <- gamma.control.raw/sum(gamma.control.raw) * 
      length(subset.W.orig)
    gamma[treated.idx] <- gamma.treated.raw/sum(gamma.treated.raw) * 
      length(subset.W.orig)
    dr.correction.all <- subset.W.orig * gamma * (subset.Y.orig - 
                                                    Y.hat.1) - (1 - subset.W.orig) * gamma * (subset.Y.orig - 
                                                                                                Y.hat.0)
    dr.correction <- mean(dr.correction.all)
    if (cluster.se) {
      correction.clust <- Matrix::sparse.model.matrix(~factor(subset.clusters) + 
                                                        0, transpose = TRUE) %*% dr.correction.all
      sigma2.hat <- sum(correction.clust^2)/length(dr.correction.all)/(length(dr.correction.all) - 
                                                                         1)
    }
    else {
      sigma2.hat <- mean(dr.correction.all^2)/(length(dr.correction.all) - 
                                                 1)
    }
  }
  else if (method == "TMLE") {
    if (target.sample == "all") {
      eps.tmle.robust.0 <- lm(B ~ A + 0, data = data.frame(A = 1/(1 - 
                                                                    subset.W.hat[subset.W.orig == 0]), B = subset.Y.orig[subset.W.orig == 
                                                                                                                           0] - Y.hat.0[subset.W.orig == 0]))
      eps.tmle.robust.1 <- lm(B ~ A + 0, data = data.frame(A = 1/subset.W.hat[subset.W.orig == 
                                                                                1], B = subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 
                                                                                                                                      1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0, 
                                     newdata = data.frame(A = mean(1/(1 - subset.W.hat))))
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1, 
                                     newdata = data.frame(A = mean(1/subset.W.hat)))
      dr.correction <- delta.tmle.robust.1 - delta.tmle.robust.0
      if (cluster.se) {
        sigma2.hat <- sandwich::vcovCL(eps.tmle.robust.0, 
                                       cluster = subset.clusters[subset.W.orig == 
                                                                   0]) * mean(1/(1 - subset.W.hat))^2 + sandwich::vcovCL(eps.tmle.robust.1, 
                                                                                                                         cluster = subset.clusters[subset.W.orig == 
                                                                                                                                                     1]) * mean(1/subset.W.hat)^2
      }
      else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * 
          mean(1/(1 - subset.W.hat))^2 + sandwich::vcovHC(eps.tmle.robust.1) * 
          mean(1/subset.W.hat)^2
      }
    }
    else if (target.sample == "treated") {
      eps.tmle.robust.0 <- lm(B ~ A + 0, data = data.frame(A = subset.W.hat[subset.W.orig == 
                                                                              0]/(1 - subset.W.hat[subset.W.orig == 0]), B = subset.Y.orig[subset.W.orig == 
                                                                                                                                             0] - Y.hat.0[subset.W.orig == 0]))
      new.center <- mean(subset.W.hat[subset.W.orig == 
                                        1]/(1 - subset.W.hat[subset.W.orig == 1]))
      delta.tmle.robust.0 <- predict(eps.tmle.robust.0, 
                                     newdata = data.frame(A = new.center))
      dr.correction <- -delta.tmle.robust.0
      if (cluster.se) {
        s.0 <- sandwich::vcovCL(eps.tmle.robust.0, cluster = subset.clusters[subset.W.orig == 
                                                                               0]) * new.center^2
        delta.1 <- Matrix::sparse.model.matrix(~factor(subset.clusters[subset.W.orig == 
                                                                         1]) + 0, transpose = TRUE) %*% (subset.Y.orig[subset.W.orig == 
                                                                                                                         1] - Y.hat.1[subset.W.orig == 1])
        s.1 <- sum(delta.1^2)/sum(subset.W.orig == 1)/(sum(subset.W.orig == 
                                                             1) - 1)
        sigma2.hat <- s.0 + s.1
      }
      else {
        sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * 
          new.center^2 + var(subset.Y.orig[subset.W.orig == 
                                             1] - Y.hat.1[subset.W.orig == 1])/sum(subset.W.orig == 
                                                                                     1)
      }
    }
    else if (target.sample == "control") {
      eps.tmle.robust.1 <- lm(B ~ A + 0, data = data.frame(A = (1 - 
                                                                  subset.W.hat[subset.W.orig == 1])/subset.W.hat[subset.W.orig == 
                                                                                                                   1], B = subset.Y.orig[subset.W.orig == 1] - Y.hat.1[subset.W.orig == 
                                                                                                                                                                         1]))
      new.center <- mean((1 - subset.W.hat[subset.W.orig == 
                                             0])/subset.W.hat[subset.W.orig == 0])
      delta.tmle.robust.1 <- predict(eps.tmle.robust.1, 
                                     newdata = data.frame(A = new.center))
      dr.correction <- delta.tmle.robust.1
      if (cluster.se) {
        delta.0 <- Matrix::sparse.model.matrix(~factor(subset.clusters[subset.W.orig == 
                                                                         0]) + 0, transpose = TRUE) %*% (subset.Y.orig[subset.W.orig == 
                                                                                                                         0] - Y.hat.0[subset.W.orig == 0])
        s.0 <- sum(delta.0^2)/sum(subset.W.orig == 0)/(sum(subset.W.orig == 
                                                             0) - 1)
        s.1 <- sandwich::vcovCL(eps.tmle.robust.1, cluster = subset.clusters[subset.W.orig == 
                                                                               1]) * new.center^2
        sigma2.hat <- s.0 + s.1
      }
      else {
        sigma2.hat <- var(subset.Y.orig[subset.W.orig == 
                                          0] - Y.hat.0[subset.W.orig == 0])/sum(subset.W.orig == 
                                                                                  0) + sandwich::vcovHC(eps.tmle.robust.1) * 
          new.center^2
      }
    }
    else {
      stop("Invalid target sample.")
    }
  }
  else {
    stop("Invalid method.")
  }
  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(sigma2.hat)
  return(c(estimate = tau.avg, std.err = tau.se))
}


###
# 1h) Define a function for predicting the conditional average treatment effects (along with standard errors and
# 90% confidence intervals) in a supplied matrix y
# The function returns the supplied matrix y with the estimates and other related variables as columns attached
###
f.pred.custom <- function(y,cforest.ncorrect0, cforest.rt0, cforest.ncor.rt0) {
  # Select only the covariates which are used in the causal forest estimations; these are the only covariates needed for the effect predictions
  dframe <- subset(y, select=vars.full)
  mat <- as.matrix(dframe)
  
  # Predict the effects and compute the related quantities for the outcome ncorrect
  prediction.ncorrect <- predict(cforest.ncorrect0, newdata=mat, estimate.variance=TRUE)
  dframe$pred.ncorrect <- prediction.ncorrect$predictions
  dframe$var.ncorrect <- prediction.ncorrect$variance.estimates
  dframe$sigma.ncorrect <- sqrt(dframe$var.ncorrect)
  dframe$ub.ncorrect <- dframe$pred.ncorrect + 1.645* dframe$sigma.ncorrect
  dframe$lb.ncorrect <- dframe$pred.ncorrect - 1.645* dframe$sigma.ncorrect
  dframe$tstat.ncorrect <- dframe$pred.ncorrect/dframe$sigma.ncorrect
  
  # Predict the effects and compute the related quantities for the outcome ncorrect
  prediction.rt <- predict(cforest.rt0, newdata=mat, estimate.variance=TRUE)
  dframe$pred.rt <- prediction.rt$predictions
  dframe$var.rt <- prediction.rt$variance.estimates
  dframe$sigma.rt <- sqrt(dframe$var.rt)
  dframe$ub.rt <- dframe$pred.rt + 1.645* dframe$sigma.rt
  dframe$lb.rt <- dframe$pred.rt - 1.645* dframe$sigma.rt
  dframe$tstat.rt <- dframe$pred.rt/dframe$sigma.rt
  
  # Predict the effects and compute the related quantities for the outcome ncorrect
  prediction.ncor.rt <- predict(cforest.ncor.rt0, newdata=mat, estimate.variance=TRUE)
  dframe$pred.ncor.rt <- prediction.ncor.rt$predictions
  dframe$var.ncor.rt <- prediction.ncor.rt$variance.estimates
  dframe$sigma.ncor.rt <- sqrt(dframe$var.ncor.rt)
  dframe$ub.ncor.rt <- dframe$pred.ncor.rt + 1.645* dframe$sigma.ncor.rt
  dframe$lb.ncor.rt <- dframe$pred.ncor.rt - 1.645* dframe$sigma.ncor.rt
  dframe$tstat.ncor.rt <- dframe$pred.ncor.rt/dframe$sigma.ncor.rt
  
  # Select only the variables newly created that refer to the predictions and the related variables in the temporary data frame dframe
  dframe <- dframe[,c("pred.ncorrect", "var.ncorrect", "sigma.ncorrect", "ub.ncorrect", "lb.ncorrect", "tstat.ncorrect",
                      "pred.rt", "var.rt", "sigma.rt", "ub.rt", "lb.rt", "tstat.rt",
                      "pred.ncor.rt", "var.ncor.rt", "sigma.ncor.rt", "ub.ncor.rt", "lb.ncor.rt", "tstat.ncor.rt")]
  
  # Attach the dframe, which includes only the selected variables anymore, to the y data frame that was given as input for the function
  # Set all predictions and related variables to missing if any of the predictor variables are missing
  y <- cbind(y, dframe)
  y[rowSums(is.na(y[vars.full]))>0,c("pred.ncorrect", "var.ncorrect", "sigma.ncorrect", "ub.ncorrect", "lb.ncorrect", "tstat.ncorrect",
                                "pred.rt", "var.rt", "sigma.rt", "ub.rt", "lb.rt", "tstat.rt",
                                "pred.ncor.rt", "var.ncor.rt", "sigma.ncor.rt", "ub.ncor.rt", "lb.ncor.rt", "tstat.ncor.rt")] <- NA
  y
}



############################
# 2) Prepare the data set for the analysis
############################

###
# 2a) Create the raw data set for our analysis, based on Caravalho et al.'s (2016) GfK KnowledgePanel data from the AER website
###

# Prepare the dataset Controls.dta
controls <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Controls.dta"))
# Keep only the individuals from the KnowledgePanel sample
controls <- controls[which(controls$sample=="KP"),]
# Keep only the variables which are needed for the analysis
controls <- controls[,c("prim_key", "age", "male", "married", "divorced", "widowed", "never_married",
                       "white", "black", "other_ethnicity", "hispanic", "mixed_ethnicity",
                       "working", "self_employed", "unemployed", "temp_layoff", "disabled", "retired", "other_empstatus",
                       "hhld_size", "college_graduate", "some_college", "high_school", "less_high_school",
                       "income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                       "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k",
                       "hhld_head", "hhld_children_0_1", "hhld_children_2_5", "hhld_children_6_12", "hhld_children_13_17", "metro")]

# Prepare the dataset Baseline.dta
baseline <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Baseline.dta"))
# Keep only the individuals from the KnowledgePanel sample
baseline <- baseline[which(baseline$sample=="KP"),]
# For each individual, keep only the observation which refers to the payday used for the experiment
baseline <- baseline[which(baseline$experimental_payment==1),]
# Keep only the variables which are needed for the analysis
baseline <- baseline[, c("prim_key", "payamount_fraction", "total_payamount", "livecheckbycheck", "hardship",
                         "caloric_crunch", "liquidity_constrained")]

# Prepare the dataset Followup.dta
followup.orig <- read.dta(paste0(path.data, "20140481_dataset_&_programs/data/secondary/Followup.dta"))
# Keep only the individuals from the KnowledgePanel sample
followup <- followup.orig[which(followup.orig$sample=="KP"),]
# Keep only the variables which are needed for the analysis
# Get the names of the 48 stroop_correct variables
stroop.correct.trials <- grep("^stroop_correct",names(followup),value=TRUE)
# Get the names of the 48 stroop_lntime variables
stroop.lntime.trials <- grep("^stroop_lntime", names(followup), value=TRUE)
followup <- followup[,c("prim_key", "treatment", stroop.correct.trials, stroop.lntime.trials,
                       "total_expenditures", "cash", "balance_checking_savings")]


# Create the raw dataset for our analysis by merging the followup, baseline, and controls data frames
full.kp <- merge(followup, controls, by="prim_key")
full.kp <- merge(full.kp, baseline, by="prim_key")

# Assert that the full.kp data frame includes exactly those 2,723 individuals from the KnowledgePanel who
# Carvalho et al. (2016) use for their analysis (the ordering of the two variables tested for identity is the same)
stopifnot(identical(full.kp$prim_key,followup.orig[which(followup.orig$sample=="KP"),"prim_key"]))

###
# 2a) Create the outcome variables for our analysis
###

# Number of correct answers (over the entire Stroop task)
# Compute the number of correct answers, setting ncorrect to missing for individuals who have not done all 48 trials
full.kp[,"ncorrect"]<-rowSums(full.kp[,stroop.correct.trials])

# Total response time in seconds (over the entire Stroop task)
# Create names for the response time variables in ms
ms.stroop.lntime <- paste("ms",stroop.lntime.trials,sep=".")
# Create response time variables in milliseconds instead of the original log milliseconds
full.kp[,ms.stroop.lntime] <- lapply(stroop.lntime.trials, function(x) exp(full.kp[,x]))
# Create the sum of all response times in ms
full.kp[,"rt.ms"] <- rowSums(full.kp[,ms.stroop.lntime])
# Create the sum of all response times in seconds
full.kp[,"rt"] <- rowSums(full.kp[,ms.stroop.lntime])/1000

# Number of correct answers per second
full.kp$ncor.rt <- full.kp$ncorrect/full.kp$rt


###
# 2b) Select all the variables that we use for creating our analysis sample
###

full.kp.mod <- select(full.kp,
                    prim_key,
                    ncorrect,
                    rt,
                    ncor.rt,
                    treatment,
                    total_payamount,
                    livecheckbycheck,
                    caloric_crunch,
                    liquidity_constrained,
                    hardship,
                    payamount_fraction,
                    age,
                    male,
                    married,
                    divorced,
                    widowed,
                    never_married,
                    white,
                    black,
                    other_ethnicity,
                    hispanic,
                    mixed_ethnicity,
                    working,
                    self_employed,
                    unemployed,
                    temp_layoff,
                    disabled,
                    retired,
                    other_empstatus,
                    hhld_size,
                    college_graduate,
                    some_college,
                    high_school,
                    less_high_school,
                    income_less_5k,
                    income_btw_5k_10k,
                    income_btw_10k_15k,
                    income_btw_15k_20k,
                    income_btw_20k_25k,
                    income_btw_25k_30k,
                    income_btw_30k_35k,
                    income_btw_35k_40k,
                    hhld_head,
                    hhld_children_0_1,
                    hhld_children_2_5,
                    hhld_children_6_12,
                    hhld_children_13_17,
                    metro)


###
# 2c) Create the analysis variables, drop the individuals above the 99th current income percentile, and remove all individuals with missing values in our analysis variables
###

# Drop the 26 individuals with a current income above the 99th percentile of the current income distribution
# This operation also drops the 167 individuals with missing current income (they would have been dropped below for the analysis, too)
full.kp.mod <- subset(full.kp.mod, subset= total_payamount <= quantile(full.kp.mod$total_payamount, 0.99, na.rm=T))

# Recode the treatment variable to be a numerical dummy
full.kp.mod$treatment.orig<- full.kp.mod$treatment
full.kp.mod$treatment<- ifelse(full.kp.mod$treatment=="Before Payday",1,0)

# Number of household children
full.kp.mod$hhld_children<-full.kp.mod$hhld_children_0_1 + full.kp.mod$hhld_children_2_5 + full.kp.mod$hhld_children_6_12 + full.kp.mod$hhld_children_13_17
full.kp.mod$hhld_children<-as.numeric(full.kp.mod$hhld_children)

# Indicator whether a child is present in the household
full.kp.mod$children <- ifelse(full.kp.mod$hhld_children>0, 1, 0)

# Working (working or self employed)
full.kp.mod$working2<-ifelse(full.kp.mod$working==1 | full.kp.mod$self_employed==1, 1,0)

# Unemployed (unemployed or temporary laid off)
full.kp.mod$unemployed2<-ifelse(full.kp.mod$unemployed==1 | full.kp.mod$temp_layoff==1, 1, 0)

# Other race (other ethnicity or mixed ethnicity)
full.kp.mod$other_race<-ifelse(full.kp.mod$other_ethnicity==1 | full.kp.mod$mixed_ethnicity==1, 1, 0)

# Make the factor variables numeric (such that they can be used in the causal forest function)
full.kp.mod$metro.orig<- full.kp.mod$metro
full.kp.mod$metro<- ifelse(full.kp.mod$metro=="Metro", 1, 0)
full.kp.mod$hhld_head.orig<- full.kp.mod$hhld_head
full.kp.mod$hhld_head<- ifelse(full.kp.mod$hhld_head=="Yes", 1, 0)

# Remove all observations with a missing value in any variable
full.kp.mod<-na.omit(full.kp.mod)

###
# 2d) Create the three data sets for the financial circumstances analyses
###

# Select all variables which are needed for the financial resources analyses
full.kp.fres <- select(full.kp,
                    prim_key,
                    treatment,
                    total_expenditures,
                    cash,
                    balance_checking_savings)

# Recode the treatment variable to be a numerical dummy
full.kp.fres$treatment.orig <- full.kp.fres$treatment
full.kp.fres$treatment <- ifelse(full.kp.fres$treatment=="Before Payday",1,0)

# Select all individuals in the analysis sample full.kp.mod who do not have a missing value in the respective outcome

# Cash
full.kp.cash <- full.kp.fres[which(full.kp.fres$prim_key %in% full.kp.mod$prim_key), c("cash", "treatment", "prim_key")]
full.kp.cash <- na.omit(full.kp.cash)

# Total expenditures
full.kp.expenditures <- full.kp.fres[which(full.kp.fres$prim_key %in% full.kp.mod$prim_key), c("total_expenditures", "treatment", "prim_key")]
full.kp.expenditures <- na.omit(full.kp.expenditures)

# Balance and checking savings
full.kp.balance <- full.kp.fres[which(full.kp.fres$prim_key %in% full.kp.mod$prim_key), c("balance_checking_savings", "treatment", "prim_key")]
full.kp.balance <- na.omit(full.kp.balance)


############################
# 3) Estimate all of the models based on the GfK KnowledgePanel in the paper
############################

###
# 3a) Financial circumstances analysis
###

# OLS regressions
# Cash
lm.cash <- lm(cash~treatment, data=full.kp.cash)
lm.cash.rob <- coeftest(lm.cash, vcov=vcovHC(lm.cash, "HC1"))

# Total expenditures
lm.expenditures <- lm(total_expenditures~treatment, data=full.kp.expenditures)
lm.expenditures.rob <- coeftest(lm.expenditures, vcov=vcovHC(lm.expenditures, "HC1"))

# Checking and savings account balances
lm.balance <- lm(balance_checking_savings~treatment, data=full.kp.balance)
lm.balance.rob <- coeftest(lm.balance, vcov=vcovHC(lm.balance, "HC1"))

# Quantile regressions
# Use the defaults for the bootstrap standard errors
# (Solutions may be nonunique, but they are plausible)
# Cash
set.seed(10101)
# Cash
q.cash <- rq(cash~treatment, data=full.kp.cash, tau=.5)
q.cash.boot <- summary(q.cash, se="boot", R=1000)

# Total expenditures
q.expenditures <- rq(total_expenditures~treatment, data=full.kp.expenditures, tau=.5)
q.expenditures.boot <- summary(q.expenditures, se="boot", R=1000)

# Balance checking and savings account
q.balance <- rq(balance_checking_savings~treatment, data=full.kp.balance, tau=.5)
q.balance.boot <- summary(q.balance, se="boot", R=1000)

# Conduct the Wilcoxon rank sum test
# Cash
wilcox.cash <- wilcox.test(cash~treatment, data=full.kp.cash)

# Total expenditures
wilcox.expenditures <- wilcox.test(total_expenditures~treatment, data=full.kp.expenditures)

# Balance checking and savings account
wilcox.balance <- wilcox.test(balance_checking_savings~treatment, data=full.kp.balance)

###
# 3b) OLS estimates for the average effect of the financial circumstances before payday on cognition
###

# Full sample average effect estimates
# Correct answers per second
lm.ncor.rt <- lm(ncor.rt~treatment,data=full.kp.mod)
lm.ncor.rt.rob <- coeftest(lm.ncor.rt,vcov=vcovHC(lm.ncor.rt,"HC1"))

# Number of correct answers
lm.ncorrect <- lm(ncorrect~treatment,data=full.kp.mod)
lm.ncorrect.rob <- coeftest(lm.ncorrect,vcov=vcovHC(lm.ncorrect,"HC1"))

# Total response time
lm.rt <- lm(rt~treatment,data=full.kp.mod)
lm.rt.rob <- coeftest(lm.rt,vcov=vcovHC(lm.rt,"HC1"))

# Average effect estimates for the subgroups which Carvalho et al. (2016) analyze
# Create a list which contains as elements all of the different subsamples to be used (except for the less than $20,000 income subsample, which we define below)
# (The subgroup payamount_fraction=1 refers to the "one payment" subsample of Carvalho et al. (2016))
subgroups <- lapply(c("payamount_fraction", "hardship", "livecheckbycheck", "caloric_crunch",
                      "liquidity_constrained"),
                    function(x0) full.kp.mod[which(full.kp.mod[,x0]==1),])
# Name the elements of the list accordingly
names(subgroups) <- c("payamount_fraction", "hardship", "livecheckbycheck", "caloric_crunch",
                      "liquidity_constrained")
# Add the less than $20,000 income subsample
subgroups[["income_less20k"]] <- full.kp.mod[which(rowSums(full.kp.mod[,c("income_less_5k", "income_btw_5k_10k",
                                                                     "income_btw_10k_15k", "income_btw_15k_20k")])==1),]

# Conduct the subsample analysis for the outcome correct answers per second
lm.subgroups.ncor.rt <- lapply(subgroups, function(x0) lm(ncor.rt~treatment, data=x0))
lm.subgroups.ncor.rt.rob <- lapply(lm.subgroups.ncor.rt, function(x0) round(coeftest(x0, vcov=vcovHC(x0, "HC1")), digits=3))

# Conduct the subsample analysis for the outcome number of correct answers
lm.subgroups.ncorrect <- lapply(subgroups, function(x0) lm(ncorrect~treatment, data=x0))
lm.subgroups.ncorrect.rob <- lapply(lm.subgroups.ncorrect, function(x0) round(coeftest(x0, vcov=vcovHC(x0, "HC1")), digits=3))

# Conduct the subsample analysis for the outcome total response time
lm.subgroups.rt <- lapply(subgroups, function(x0) lm(rt~treatment, data=x0))
lm.subgroups.rt.rob <- lapply(lm.subgroups.rt, function(x0) round(coeftest(x0, vcov=vcovHC(x0, "HC1")), digits=3))


###
# 3c) Estimate the causal forest models
###

# Set the model parameters (all other parameters not specified here are set to the causal_forest function's default values)
# Number of trees
num.trees1 <- c(10000)
# Number of variables tried for each split (the value 27 actually corresponds to the default function value when having 37 covariates, as we do)
mtry1 <- c(27)
# Fraction of the data used to build each tree (corresponds also actually to the function's default value)
sample.fraction1 <- c(0.5)
# Minimum node size target
min.node.size1 <- c(2)
# Seed
seed1 <- c(10101)

# Estimate the causal forest models for our three outcomes
models.cf <- f.est.full(num.trees0 = num.trees1, sample.fraction0 = sample.fraction1, min.node.size0 = min.node.size1, seed0 = seed1, mtry0 = mtry1)
names(models.cf) <- c("ncorrect", "rt", "ncor.rt")

# Put all of the covariates used in the causal forest estimations in the vector vars
# (the vector vars.full has been created by the function f.est.full)
vars <- vars.full

###
# 3d) Conduct the balance checks
###

# Conduct pairwise t-tests for all of our analysis covariates (Welch t-test)
ttests.list <- list()
ttests.list[vars] <- lapply(vars, function(x0) t.test(as.formula(paste0(x0, "~treatment")), data=full.kp.mod))

# Create a data frame containing the results from the t-tests
# Extract the information for which covariate the t-tests have been conducted
ttests.comparisons <- lapply(vars, function(x0) ttests.list[[x0]][["data.name"]])
ttests.comparisons <- unlist(ttests.comparisons)

# Extract the mean values for the after payday group
ttests.mean0 <- lapply(vars, function(x0) ttests.list[[x0]][["estimate"]][[1]])
ttests.mean0 <- unlist(ttests.mean0)

# Extract the mean values for the before payday group
ttests.mean1 <- lapply(vars, function(x0) ttests.list[[x0]][["estimate"]][[2]])
ttests.mean1 <- unlist(ttests.mean1)

# Extract the p-values for the mean comparisons
ttests.pvalue <- lapply(vars, function(x0) ttests.list[[x0]][["p.value"]])
ttests.pvalue <- unlist(ttests.pvalue)

# Create the actual data frame
ttests <- data.frame(comparison=ttests.comparisons, mean.untreated=round(ttests.mean0, digits=3), mean.treated=round(ttests.mean1, digits=3), pvalue=round(ttests.pvalue, digits=3))

# Test if all covariates jointly predict the treatment
# Specify the full model
# Drop the covariates never_married, other_race, other_empstatus, less_high_school, income_btw_35k_40k from the regression due to multicolinearity
vars.reg <- vars[which(!(vars %in% c("never_married", "other_race", "other_empstatus", "less_high_school", "income_btw_35k_40k")))]
model.full <- lm(as.formula(paste0("treatment~", paste0(vars.reg, collapse="+"))), data=full.kp.mod)

# Specify the model which only includes a constant (this is the model we want to test model.full against)
model.constant <- lm(treatment~1, data=full.kp.mod)

# Conduct the joint F-test
ftest <- waldtest(model.full, model.constant, vcov=vcovHC(model.full, "HC1"), test="F")

###
# 3e) Estimate the average treatment effect for the subgroup determined based on the causal forest analysis insights
# (Use the augmented inverse probability weighted estimator; subgroup: all individuals who have a current income below $750
# and who are either below 30 or above 70 years of age)
###

# Number of correct answers per second
# Define the subgroup
ncor.rt.subgroup <- which(models.cf[["ncor.rt"]][["X.orig"]][,"total_payamount"]<750 &
                             (models.cf[["ncor.rt"]][["X.orig"]][,"age"]<30 | models.cf[["ncor.rt"]][["X.orig"]][,"age"]>70))
# Estimate the effect
ate.ncor.rt <- average_treatment_effect_subset(models.cf[["ncor.rt"]], target.sample = "all", method="AIPW",
                                               subset=ncor.rt.subgroup)
# Number of observations
ate.ncor.rt.nobs <- length(ncor.rt.subgroup)

# Number of correct answers
# Define the subgroup
ncorrect.subgroup <- which(models.cf[["ncorrect"]][["X.orig"]][,"total_payamount"]<750 &
                            (models.cf[["ncorrect"]][["X.orig"]][,"age"]<30 | models.cf[["ncorrect"]][["X.orig"]][,"age"]>70))
# Estimate the effect
ate.ncorrect <- average_treatment_effect_subset(models.cf[["ncorrect"]], target.sample = "all", method="AIPW",
                                               subset=ncorrect.subgroup)
# Number of observations
ate.ncorrect.nobs <- length(ncorrect.subgroup)

# Total response time
# Define the subgroup
rt.subgroup <- which(models.cf[["rt"]][["X.orig"]][,"total_payamount"]<750 &
                             (models.cf[["rt"]][["X.orig"]][,"age"]<30 | models.cf[["rt"]][["X.orig"]][,"age"]>70))
# Estimate the effect
ate.rt <- average_treatment_effect_subset(models.cf[["rt"]], target.sample = "all", method="AIPW",
                                               subset=rt.subgroup)
# Number of observations
ate.rt.nobs <- length(rt.subgroup)

# Assert that the defined subgroups are identical
stopifnot(identical(models.cf[["ncor.rt"]][["X.orig"]], models.cf[["ncorrect"]][["X.orig"]]))
stopifnot(identical(models.cf[["ncor.rt"]][["X.orig"]], models.cf[["rt"]][["X.orig"]]))
stopifnot(identical(ncor.rt.subgroup, ncorrect.subgroup))
stopifnot(identical(ncor.rt.subgroup, rt.subgroup))


############################
# 4) Create the matrices for the heatmaps and the typical individuals and estimate the corresponding effects
############################

###
# 4a) Create the matrix based on which we create the heatmaps
###

# Create a data frame which contains all of the age and current income combinations of interest
# Age in one year steps; current income in $25 steps
X.full <- data.frame(total_payamount=rep(seq(0, 1500, by=25), each=(93-18+1)), age=rep(seq(18,93, by=1), times=((1500-0)/25)+1))

# Create a variable that we use for adding the most commonly occuring characteristics in the full
# sample to the X.full data frame
X.full$merge <- 1

# Create data frames which contain the mean and median values for the full sample
x.mean.full <- full.kp.mod %>% select(vars.full) %>% summarise_all(mean) %>% ungroup()
x.median.full <- full.kp.mod %>% select(vars.full) %>% summarise_all(median) %>% ungroup()

# Set the dummies in the x.mean.full data frame according to the most frequently occuring characteristics
# in the full sample
# Marital status, race, employment status, education, income
x.mean.full <- f.class(x.mean.full)
# All dummies which do not refer to a class of dummies (e.g., metro, hardship, children)
x.mean.full <- f.dum(x.mean.full)
# Set the share of payday pay amount relative to current income and the household size
# to their median values (do so by replacing the respective mean values by their median values in
# the data frame x.mean.full)
x.mean.full$payamount_fraction <- x.median.full$payamount_fraction
x.mean.full$hhld_size <- x.median.full$hhld_size

# Create the variable for meging the X.full data frame with the x.mean.full data frame
x.mean.full$merge <- 1

# Merge the X.full and x.mean.full data frames
X.full <- merge(X.full, x.mean.full, by="merge", all.x=T)
X.full$age <- X.full$age.x
X.full$total_payamount <- X.full$total_payamount.x

# Keep only the covariates in the X.full data frame that we use in the causal forest
# analysis
X.full <- X.full[vars.full]

# Assert that all of the dummies are either set to 0 or 1
stopifnot(X.full[vars.full][,!(names(X.full[vars.full]) %in% c("age", "total_payamount", "payamount_fraction", "hhld_size"))] ==0 |
            X.full[vars.full][,!(names(X.full[vars.full]) %in% c("age", "total_payamount", "payamount_fraction", "hhld_size"))] ==1)

# Assert that in each dummy class there is exactly one dummy set to 1 and all other dummies
# are set to zero
stopifnot(rowSums(X.full[,c("married", "divorced", "widowed", "never_married")])==1)
stopifnot(rowSums(X.full[,c("white", "black", "other_race", "hispanic")])==1)
stopifnot(rowSums(X.full[,c("working2", "unemployed2", "disabled", "retired", "other_empstatus")])==1)
stopifnot(rowSums(X.full[,c("college_graduate", "some_college", "high_school", "less_high_school")])==1)
stopifnot(rowSums(X.full[,c("income_less_5k", "income_btw_5k_10k", "income_btw_10k_15k", "income_btw_15k_20k",
                            "income_btw_20k_25k", "income_btw_25k_30k", "income_btw_30k_35k", "income_btw_35k_40k")])==1)

# Assert that the household size is an integer
stopifnot(X.full[,"hhld_size"]%%1==0)

# Assert that the payamount fraction is in the range >0 and 1
stopifnot(X.full[,"payamount_fraction"]>0 & X.full[,"payamount_fraction"]<=1)


###
# 4b) Predict the effects that are plotted in the heatmaps, using the X.full data frame
###

X.full.pred <- f.pred.custom(X.full, models.cf[["ncorrect"]], models.cf[["rt"]],
                             models.cf[["ncor.rt"]])


###
# 4c) Create the data frames for the typical individuals
###

# Use as the standard window size 5 years of age and $250 of current income
# For solving ties in the non-ordinal characteristics, extend the window to 7 years of age and
# $450 of current income; every element in the list X.typical.list corresponds to a data frame for
# a given age-current income combination
X.typical.list <- f.typical(5, 250, 7, 450)

###
# 4d) Estimate the effects for the typical individuals, using the X.typical.list list
# To speed up the estimation process, use parallel computing
###

# Detect the number of cores on the machine and set up the cluster
ncors <- detectCores()
cl <- makeCluster(ncors)

# Export to the cluster the data and functions that are needed for the estimation
clusterExport(cl, c(
  "models.cf",
  "X.typical.list",
  "vars.full",
  "f.pred.custom"
))

clusterEvalQ(cl, c(
  library(plyr),
  library(dplyr),
  library(grf),
  library(pryr)
))

# Estimate the effects
X.typical.list.pred <- parLapply(cl, 1:length(X.typical.list), function(z0) {
  f.pred.custom(X.typical.list[[z0]], models.cf[["ncorrect"]], models.cf[["rt"]],
                models.cf[["ncor.rt"]])
})

# Stop the cluster
stopCluster(cl)

# Name the elements in the list X.typical.list.pred accordingly
for(m0 in 1:length(X.typical.list.pred)) {
  names(X.typical.list.pred) <- names(X.typical.list)
}


############################
# 5) Create the descriptive statistics and causal forest results
############################

###
# 5a) Descriptive statistics
###

# Create a matrix which contains all of the summary statistics for the analysis variables
descriptives.list <- lapply(c("mean", "sd"), function(x) summarise_all(subset(full.kp.mod, select=c(vars.full, "ncor.rt", "ncorrect", "rt", "treatment")), eval(parse(text=x))) %>% mutate(stat=x))
descriptives.wide <- do.call(rbind,c(descriptives.list))
descriptives <- transpose(descriptives.wide)

# Name the rows and columns of the matrix containing the descriptive statistics
row.names(descriptives) <- names(descriptives.wide)
colnames(descriptives) <- descriptives[nrow(descriptives),]

# Remove the last row of the descriptive statistics matrix (the last row contains only the names
# of the descriptive statistics that were calculated)
descriptives <- descriptives[-nrow(descriptives),]

# Round the descriptive statistics to three digits
descriptives[colnames(descriptives)] <- lapply(descriptives, function(x) round(as.numeric(x),digits=3))


###
# 5b) Variable importance measure for the three causal forests
###

# Create a matrix containing the variable importance measures for our three causal forests
vi <- fvi(cforest.ncorrect0=models.cf[["ncorrect"]], cforest.rt0=models.cf[["rt"]],
          cforest.ncor.rt0=models.cf[["ncor.rt"]])
vi[,"vi.ncorrect.val"] <- round(as.numeric(vi[,"vi.ncorrect.val"]), digits=4)
vi[,"vi.rt.val"] <- round(as.numeric(vi[,"vi.rt.val"]), digits=4)
vi[,"vi.ncor.rt.val"] <- round(as.numeric(vi[,"vi.ncor.rt.val"]), digits=4)


###
# 5c) Create a table, which contains the vicinity results
# To do so, loop over all age-current income combinations and extract the estimated effects, standard errors, and
# t-statistics from the list X.typical.list.pred
###

# Set the age and current income ranges for the typical individuals
age.younger.vec <- seq(18,22,1)
age.older.vec <- seq(73,77,1)
inc.vec <- seq(400,500,25)

# Younger individuals
# Correct answers per second
vicinity.ncor.rt <- data.frame(inc=rep(inc.vec, each=3))
vicinity.ncor.rt[,paste0("age",age.younger.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.younger.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.ncor.rt[row.pred, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.ncor.rt"]
    vicinity.ncor.rt[row.sigma, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.ncor.rt"]
    vicinity.ncor.rt[row.t, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.ncor.rt"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.ncor.rt"]
  }
}

# Number of correct answers
vicinity.ncorrect <- data.frame(inc=rep(inc.vec, each=3))
vicinity.ncorrect[,paste0("age",age.younger.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.younger.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.ncorrect[row.pred, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.ncorrect"]
    vicinity.ncorrect[row.sigma, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.ncorrect"]
    vicinity.ncorrect[row.t, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.ncorrect"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.ncorrect"]
  }
}

# Total response time
vicinity.rt <- data.frame(inc=rep(inc.vec, each=3))
vicinity.rt[,paste0("age",age.younger.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.younger.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.rt[row.pred, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.rt"]
    vicinity.rt[row.sigma, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.rt"]
    vicinity.rt[row.t, paste0("age",age.younger.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"pred.rt"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.younger.vec[j0])]][1,"sigma.rt"]
  }
}

# Older individuals
# Correct answers per second
vicinity.ncor.rt.older <- data.frame(inc=rep(inc.vec, each=3))
vicinity.ncor.rt.older[,paste0("age",age.older.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.older.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.ncor.rt.older[row.pred, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.ncor.rt"]
    vicinity.ncor.rt.older[row.sigma, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.ncor.rt"]
    vicinity.ncor.rt.older[row.t, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.ncor.rt"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.ncor.rt"]
  }
}

# Number of correct answers
vicinity.ncorrect.older <- data.frame(inc=rep(inc.vec, each=3))
vicinity.ncorrect.older[,paste0("age",age.older.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.older.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.ncorrect.older[row.pred, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.ncorrect"]
    vicinity.ncorrect.older[row.sigma, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.ncorrect"]
    vicinity.ncorrect.older[row.t, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.ncorrect"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.ncorrect"]
  }
}

# Total response time
vicinity.rt.older <- data.frame(inc=rep(inc.vec, each=3))
vicinity.rt.older[,paste0("age",age.older.vec)] <- NA

for(i0 in 1:length(inc.vec)) {
  for(j0 in 1:length(age.older.vec)) {
    row.pred <- 1 + (i0-1)*3
    row.sigma <- 2 + (i0-1)*3
    row.t <- 3 + (i0-1)*3
    vicinity.rt.older[row.pred, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.rt"]
    vicinity.rt.older[row.sigma, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.rt"]
    vicinity.rt.older[row.t, paste0("age",age.older.vec[j0])] <- X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"pred.rt"]/X.typical.list.pred[[paste0("tpa", inc.vec[[i0]], "age", age.older.vec[j0])]][1,"sigma.rt"]
  }
}


############################
# 6) Save all of the results
############################


###
# 6a) Table displaying the results in Appendix Figure 1: Variable importance measures matrix
###
capture.output(vi, file="KP_results/Table_Figure_A1_vi.txt")

###
# 6b) Appendix Tables 2 and 3: Descriptive statistics
###
capture.output(descriptives, file="KP_results/Table_A2_A3_descriptive_statistics.txt")

###
# 6c) Appendix Table 1: Financial circumstances analysis
###
options("scipen"=999)
# OLS regressions
cat("\n lm.cash.rob \n N=", nobs(lm.cash), file="KP_results/Table_A1_financial_circumstances_ols.txt")
capture.output(round(lm.cash.rob, digits=2), file="KP_results/Table_A1_financial_circumstances_ols.txt", append=T)
cat("\n lm.expenditures.rob \n N=", nobs(lm.expenditures), file="KP_results/Table_A1_financial_circumstances_ols.txt", append=T)
capture.output(round(lm.expenditures.rob, digits=2), file="KP_results/Table_A1_financial_circumstances_ols.txt", append=T)
cat("\n lm.balance.rob \n N=", nobs(lm.balance), file="KP_results/Table_A1_financial_circumstances_ols.txt", append=T)
capture.output(print(lm.balance.rob, digits=6), file="KP_results/Table_A1_financial_circumstances_ols.txt", append=T)

# Quantile regressions
cat("\n q.cash.boot \n", file="KP_results/Table_A1_financial_circumstances_quantile.txt")
capture.output(print(q.cash.boot, digits=2), file="KP_results/Table_A1_financial_circumstances_quantile.txt", append=T)
cat("\n q.expenditures.boot \n", file="KP_results/Table_A1_financial_circumstances_quantile.txt", append=T)
capture.output(print(q.expenditures.boot, digits=2), file="KP_results/Table_A1_financial_circumstances_quantile.txt", append=T)
cat("\n q.balance.boot \n", file="KP_results/Table_A1_financial_circumstances_quantile.txt", append=T)
capture.output(print(q.balance.boot, digits=2), file="KP_results/Table_A1_financial_circumstances_quantile.txt", append=T)

# Wilcoxon test
cat("\n wilcox.cash \n", file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt")
capture.output(print(wilcox.cash, digits=2), file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt", append=T)
cat("\n wilcox.expenditures \n", file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt", append=T)
capture.output(print(wilcox.expenditures, digits=2), file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt", append=T)
cat("\n wilcox.balance \n", file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt", append=T)
capture.output(print(wilcox.balance, digits=2), file="KP_results/Table_A1_financial_circumstances_wilcoxon.txt", append=T)

###
# 6d) OLS Cognition analysis
###

# Table 1: Regressions based on the full sample
cat("\n lm.ncor.rt.rob \n N=", nobs(lm.ncor.rt), file="KP_results/Table_1_cognition_ols.txt")
capture.output(round(lm.ncor.rt.rob, digits=3), file="KP_results/Table_1_cognition_ols.txt", append=T)
cat("\n lm.ncorrect.rob \n N=", nobs(lm.ncorrect), file="KP_results/Table_1_cognition_ols.txt", append=T)
capture.output(round(lm.ncorrect.rob, digits=3), file="KP_results/Table_1_cognition_ols.txt", append=T)
cat("\n lm.rt.rob \n N=", nobs(lm.rt), file="KP_results/Table_1_cognition_ols.txt", append=T)
capture.output(round(lm.rt.rob, digits=3), file="KP_results/Table_1_cognition_ols.txt", append=T)

# Appendix Table 5: Regressions based on the subsamples analyzed by Caravalho et al. (2016)
cat("Number of observations \n", file="KP_results/Table_A5_cognition_subgroups_ols.txt")
capture.output(lapply(lm.subgroups.ncor.rt, function(x0) nobs(x0)), file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
cat("\n lm.subgroups.ncor.rt.rob \n", file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
capture.output(print(lm.subgroups.ncor.rt.rob, digits=1), file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
cat("\n lm.subgroups.ncorrect.rob \n", file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
capture.output(print(lm.subgroups.ncorrect.rob, digits=2), file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
cat("\n lm.subgroups.rt.rob \n", file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)
capture.output(print(lm.subgroups.rt.rob, digits=3), file="KP_results/Table_A5_cognition_subgroups_ols.txt", append=T)

###
# 6e) Table 2: Causal forest average effect estimates for the subgroup determined by the insights from the causal forest analysis
###
cat(paste0("\n ate.ncor.rt \n N=", ate.ncor.rt.nobs, "\n"), file="KP_results/Table_2_cognition_subgroup_cf.txt")
capture.output(round(ate.ncor.rt, digits=3), file="KP_results/Table_2_cognition_subgroup_cf.txt", append=T)
cat(paste0("\n ate.ncorrect \n N=", ate.ncorrect.nobs, "\n"), file="KP_results/Table_2_cognition_subgroup_cf.txt", append=T)
capture.output(round(ate.ncorrect, digits=3), file="KP_results/Table_2_cognition_subgroup_cf.txt", append=T)
cat(paste0("\n ate.rt \n N=", ate.rt.nobs, "\n"), file="KP_results/Table_2_cognition_subgroup_cf.txt", append=T)
capture.output(round(ate.rt, digits=3), file="KP_results/Table_2_cognition_subgroup_cf.txt", append=T)

###
# 6f) Appendix Table 4: Balance checks
###
cat("\n Pairwise mean comparisons \n", file=paste0("KP_results/Table_A4_balance_checks.txt"))
capture.output(ttests, file=paste0("KP_results/Table_A4_balance_checks.txt"), append=T)
cat("\n Joint test \n", file=paste0("KP_results/Table_A4_balance_checks.txt"), append=T)
capture.output(round(ftest, digits=3), file=paste0("KP_results/Table_A4_balance_checks.txt"), append=T)

###
# 6g) Vincinity estimates (including the typical younger and older individual discussed in the text)
###
# Appendix Table 6: Younger individuals
capture.output(
  round(vicinity.ncor.rt, digits=4),
  file="KP_results/Table_A6_vicinity_ncor_rt_younger.txt")
capture.output(
  round(vicinity.ncorrect, digits=3),
  file="KP_results/Table_A6_vicinity_ncorrect_younger.txt")
capture.output(
  round(vicinity.rt, digits=3),
  file="KP_results/Table_A6_vicinity_rt_younger.txt")

# Appendix Table 7: Older individuals
capture.output(
  round(vicinity.ncor.rt.older, digits=4),
  file="KP_results/Table_A7_vicinity_ncor_rt_older.txt")
capture.output(
  round(vicinity.ncorrect.older, digits=3),
  file="KP_results/Table_A7_vicinity_ncorrect_older.txt")
capture.output(
  round(vicinity.rt.older, digits=3),
  file="KP_results/Table_A7_vicinity_rt_older.txt")


############################
# 7) Create the graphs in the paper
############################

# Use the standard font in R for creating the graphs (in the working paper, we use the font Times New Roman)

###
# 7a) Create the variable importance plots
###

# Create a vector containing the labels which we want to use for the variable importance plots
# The order of the elements in var.names corresponds to the order of the variables in vars.full
var.names <- c("Live from check to check", "Caloric crunch", "Liquidity constrained", "Financial hardship",
               "Current income", "Share payday pay amount","Age", "Male", "Married", "Divorced",
               "Widowed", "Never married", "White", "Black", "Hispanic", "Other race", "Working", "Unemployed",
               "Disabled", "Retired", "Other employment status", "Household size", "College graduate", "Some college",
               "High school", "Less than high school", "An. hh inc. less than $5k", "An. hh inc. btw. $5k and $10k",
               "An. hh inc. btw. $10k and $15k", "An. hh inc. btw. $15k and $20k", "An. hh inc. btw. $20k and $25k",
               "An. hh inc. btw. $25k and $30k", "An. hh inc. btw. $30k and $35k", "An. hh inc. btw. $35k and $40k",
               "Household head", "Children", "Metropolitan area")

# Create a data frame which contains the variable importance measure for each of the three causal forest models

# Numer of correct answers per second
vi.ncor.rt <- variable_importance(models.cf[["ncor.rt"]])
vi.ncor.rt <- data.frame(vi.ncor.rt)
vi.ncor.rt$var.names <- var.names
vi.ncor.rt <- vi.ncor.rt[order(-vi.ncor.rt$vi.ncor.rt),]    # Order the data frame such that the covariate with the highest variable importance is in the last row
vi.ncor.rt$var.names <- factor(vi.ncor.rt$var.names, levels=vi.ncor.rt$var.names)
vi.ncor.rt$vi.ncor.rt.100 <- vi.ncor.rt$vi.ncor.rt*100  # Multiply the variable importance measure by 100 for readability

# Number of correct answers
vi.ncorrect <- variable_importance(models.cf[["ncorrect"]])
vi.ncorrect <- data.frame(vi.ncorrect)
vi.ncorrect$var.names <- var.names
vi.ncorrect <- vi.ncorrect[order(-vi.ncorrect$vi.ncorrect),]  # Order the data frame such that the covariate with the highest variable importance is in the last row
vi.ncorrect$var.names <- factor(vi.ncorrect$var.names, levels=vi.ncorrect$var.names)
vi.ncorrect$vi.ncorrect.100 <- vi.ncorrect$vi.ncorrect*100  # Multiply the variable importance measure by 100 for readability

# Total response time
vi.rt <- variable_importance(models.cf[["rt"]])
vi.rt <- data.frame(vi.rt)
vi.rt$var.names <- var.names
vi.rt <- vi.rt[order(-vi.rt$vi.rt),]  # Order the data frame such that the covariate with the highest variable importance is in the last row
vi.rt$var.names <- factor(vi.rt$var.names, levels=vi.rt$var.names)
vi.rt$vi.rt.100 <- vi.rt$vi.rt*100  # Multiply the variable importance measure by 100 for readability


# Create the actual plots and save them

# Number of correct answers per second
vi.plot.ncor.rt <- ggplot(data=vi.ncor.rt, aes(x=var.names, y=vi.ncor.rt.100)) + geom_col(width=0.85, fill="darkblue", colour="darkblue") +
  theme_bw() +
  scale_y_continuous(expand=c(0.007,0), breaks=c(seq(0,30,5)), limits=c(0, 31)) +
  #theme(text=element_text(family="Times New Roman"),
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=23),
        axis.text.x = element_text(size=20, margin = margin(t = 7, r = 0, b = 7, l = 0)),
        axis.title=element_text(size=25),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.background = element_rect(colour = "black", size=0.7),
        axis.ticks.x=element_line(size=0.7)) +
  labs(y="Variable Importance") + coord_flip()
# Appendix Figure 1
ggsave("KP_figures/Figure_A1_vi_plot_ncor_rt.png", vi.plot.ncor.rt, height=300, width=200, units="mm", device="png")

# Number of correct answers
vi.plot.ncorrect <- ggplot(data=vi.ncorrect, aes(x=var.names, y=vi.ncorrect.100)) + geom_col(width=0.85, fill="darkblue", colour="darkblue") +
  theme_bw() +
  scale_y_continuous(expand=c(0.007,0), breaks=c(seq(0,30,5)), limits=c(0, 31)) +
  #theme(text=element_text(family="Times New Roman"),
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=23),
        axis.text.x = element_text(size=20, margin = margin(t = 7, r = 0, b = 7, l = 0)),
        axis.title=element_text(size=25),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.background = element_rect(colour = "black", size=0.7),
        axis.ticks.x=element_line(size=0.7)) +
  labs(y="Variable Importance") + coord_flip()
# Appendix Figure 1
ggsave("KP_figures/Figure_A1_vi_plot_ncorrect.png", vi.plot.ncorrect, height=300, width=200, units="mm", device="png")

# Total response time
vi.plot.rt <- ggplot(data=vi.rt, aes(x=var.names, y=vi.rt.100)) + geom_col(width=0.85, fill="darkblue", colour="darkblue") +
  theme_bw() +
  scale_y_continuous(expand=c(0.007,0), breaks=c(seq(0,30,5)), limits=c(0, 31)) +
  #theme(text=element_text(family="Times New Roman"),
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=23),
        axis.text.x = element_text(size=20, margin = margin(t = 7, r = 0, b = 7, l = 0)),
        axis.title=element_text(size=25),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.background = element_rect(colour = "black", size=0.7),
        axis.ticks.x=element_line(size=0.7)) +
  labs(y="Variable Importance") + coord_flip()
# Appendix Figure 1
ggsave("KP_figures/Figure_A1_vi_plot_rt.png", vi.plot.rt, height=300, width=200, units="mm", device="png")

###
# 7b) Create the heatmaps
###

# Assert that the maximum current income amount corresponds to the median current income in our full sample
stopifnot(quantile(full.kp.mod$total_payamount, probs=0.5)==max(X.full.pred$total_payamount))

# Correct answers per second
# Set the effect magnitude intervals which are the categories for the estimated effects in the heatmaps
# (rename the lowest and highest category such that not too many digits are displayed)
val.cut.ncor.rt <- c(min(X.full.pred$pred.ncor.rt), -0.03, -0.02, -0.01, 0, 0.01, max(X.full.pred$pred.ncor.rt))
X.full.pred$cut.ncor.rt <- cut(X.full.pred$pred.ncor.rt, val.cut.ncor.rt, include.lowest = T, right=F)
levels(X.full.pred$cut.ncor.rt)[1] <- "[-0.045,-0.03)"
levels(X.full.pred$cut.ncor.rt)[6] <- "[0.01,0.020]"

p.heat.ncor.rt <- ggplot(data=X.full.pred, aes(x=total_payamount, y=age)) +
  geom_raster(aes(fill=cut.ncor.rt), interpolate = T) +
  scale_fill_manual(name="",values = rev(diverge_hsv(length(val.cut.ncor.rt)-1+2))[1:6], drop=F,
                    guide = guide_legend(direction = "horizontal", nrow=1, byrow=T, label.position = "bottom",
                                         label.theme = element_text(angle=0, size=15))) +
  theme(legend.position = "bottom", legend.key.width = unit(3.7, "cm"), legend.justification = 'center', legend.text.align = 0.5,
        axis.text.x = element_text(angle=45, vjust=0.75, hjust=0.9, size=29),
        axis.text.y = element_text(size=29),
        axis.title = element_text(size=29)) +
  scale_x_continuous(breaks=c(1,seq(250,1500,250)), expand=c(0,0)) +
  scale_y_continuous(breaks=c(seq(20,90,5)), expand=c(0,0)) +
  labs(x="Current income", y="Age")
# Figure 1
ggsave("KP_figures/Figure_1_p_heat_ncor_rt.png", p.heat.ncor.rt, height=300, width=300, units="mm", device="png")

# Number of correct answers
# Set the effect magnitude intervals which are the categories for the estimated effects in the heatmaps
# (rename the lowest and highest category such that not too many digits are displayed)
val.cut.ncorrect <- c(min(X.full.pred$pred.ncorrect), -1, -0.5, 0, 0.5, 1, max(X.full.pred$pred.ncorrect))
X.full.pred$cut.ncorrect <- cut(X.full.pred$pred.ncorrect, val.cut.ncorrect, include.lowest = T, right=F)
p.heat.ncorrect <- ggplot(data=X.full.pred, aes(x=total_payamount, y=age)) +
  geom_raster(aes(fill=cut.ncorrect), interpolate = T) +
  scale_fill_manual(name="",values = rev(diverge_hsv(length(val.cut.ncorrect)-1)), drop=F,
                    guide = guide_legend(direction = "horizontal", nrow=1, byrow=T, label.position = "bottom",
                                         label.theme = element_text(angle=0, size=15))) +
  theme(legend.position = "bottom", legend.key.width = unit(3.7, "cm"), legend.justification = 'center', legend.text.align = 0.5,
        axis.text.x = element_text(angle=45, vjust=0.75, hjust=0.9, size=29),
        axis.text.y = element_text(size=29),
        axis.title = element_text(size=29)) +
  scale_x_continuous(breaks=c(1,seq(250,1500,250)), expand=c(0,0)) +
  scale_y_continuous(breaks=c(seq(20,90,5)), expand=c(0,0)) +
  labs(x="Current income", y="Age")
# Figure 1
ggsave(paste0("KP_figures/Figure_1_p_heat_ncorrect.png"), p.heat.ncorrect, height=300, width=300, units="mm", device="png")

# Total response time
# Set the effect magnitude intervals which are the categories for the estimated effects in the heatmaps
# (rename the lowest and highest category such that not too many digits are displayed)
val.cut.rt <- c(min(X.full.pred$pred.rt), -2, -1, 0, 1, 2, 3, max(X.full.pred$pred.rt))
X.full.pred$cut.rt <- cut(X.full.pred$pred.rt, val.cut.rt, include.lowest = T, right=F)

p.heat.rt <- ggplot(data=X.full.pred, aes(x=total_payamount, y=age)) +
  geom_raster(aes(fill=cut.rt), interpolate = T) +
  scale_fill_manual(name="",values = (diverge_hsv(length(val.cut.rt)-1+1)[2:8]), drop=F,
                    guide = guide_legend(direction = "horizontal", nrow=1, byrow=T, label.position = "bottom",
                                         label.theme = element_text(angle=0, size=15))) +
  theme(legend.position = "bottom", legend.key.width = unit(3.7, "cm"), legend.justification = 'center', legend.text.align = 0.5,
        axis.text.x = element_text(angle=45, vjust=0.75, hjust=0.9, size=29),
        axis.text.y = element_text(size=29),
        axis.title = element_text(size=29)) +
  scale_x_continuous(breaks=c(1,seq(250,1500,250)), expand=c(0,0)) +
  scale_y_continuous(breaks=c(seq(20,90,5)), expand=c(0,0)) +
  labs(x="Current income", y="Age")
# Figure 1
ggsave("KP_figures/Figure_1_p_heat_rt.png", p.heat.rt, height=300, width=300, units="mm", device="png")

###
# 7c) Create the plots for the typical individuals
###

# Create the y labels that we use for the typical individual plots
change.var <- c("Typical individual: baseline",
                "Live from check to check=1",
                "Live from check to check=0",
                "Caloric crunch=1",
                "Caloric crunch=0",
                "Liquidity constrained=1",
                "Liquidity constrained=0",
                "Hardship=1",
                "Hardship=0",
                "Male=1",
                "Male=0",
                "Married=1",
                "Divorced=1",
                "Widowed=1",
                "Never married=1",
                "White=1",
                "Black=1",
                "Hispanic=1",
                "Other race=1",
                "Working=1",
                "Unemployed=1",
                "Disabled=1",
                "Retired=1",
                "Other empstatus=1",
                "College graduate=1",
                "Some college=1",
                "High school=1",
                "Less than high school=1",
                "An. hh inc less than $5k=1",
                "An. hh inc btw $5k & $10k=1",
                "An. hh inc btw $10k & $15k=1",
                "An. hh inc btw $15k & $20k=1",
                "An. hh inc btw $20k & $25k=1",
                "An. hh inc btw $25k & $30k=1",
                "An. hh inc btw $30k & $35k=1",
                "An. hh inc btw $35k & $40k=1",
                "Hhld head=1",
                "Hhld head=0",
                "Children=1",
                "Children=0",
                "Metropolitan area=1",
                "Metropolitan area=0",
                "Household size=1",
                "Household size=2",
                "Household size=3",
                "Household size=4",
                "Household size=5",
                "Share payday pay amount=0.25",
                "Share payday pay amount=0.5",
                "Share payday pay amount=0.75",
                "Share payday pay amount=1"
  ) 

###
# Typical individual with age=20 and current income=$450
###

# Extract the data which we want to plot from X.typical.list.pred
typical.20450 <- X.typical.list.pred[["tpa450age20"]]

# Add the adapted change.var label vector to the data frame
typical.20450$change.var.orig <- typical.20450$change.var
typical.20450$change.var <- factor(change.var, levels=rev(change.var))

# Create the actual plot and save it
p20450.ncor.rt <- ggplot(data=typical.20450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(seq(-0.1, -0.01,0.01), 0.01, 0.02), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.ncor.rt, ymin = lb.ncor.rt, ymax = ub.ncor.rt),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-0.1, 0.02,0.02)), limits=c(-0.10, 0.02)) +
  coord_flip() + theme_bw() + ggtitle("Panel A. Correct \n Answers per Second") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        axis.text.y = element_text(size=18),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

p20450.ncorrect.nl <- ggplot(data=typical.20450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(seq(-6, -1, 1), 1 , 2), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.ncorrect, ymin = lb.ncorrect, ymax = ub.ncorrect),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-6, 2, 1)), limits=c(-6.2, 2.5)) +
  coord_flip() + theme_bw() + ggtitle("Panel B. Number \n of Correct Answers") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        # The next command switches off the y text, which is switched on for the first graph
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

p20450.rt.nl <- ggplot(data=typical.20450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(-2, -1, seq(1,12)), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.rt, ymin = lb.rt, ymax = ub.rt),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-2, 12, 2)), limits=c(-3.6, 12.6)) +
  coord_flip() + theme_bw() + ggtitle("Panel C. Total \n Response Time") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        # The next command switches off the y text, which is switched on for the first graph
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

# The following three commands return the warning message: "Removed 15 rows containing missing values (geom_pointrange)."
# The warning messages result from the 15 "blank" rows in each data frame, which give information about the covariate values
# for the typical individual in its baseline specification -> warning messages can be ignored
gt1 <- ggplotGrob(p20450.ncor.rt)
gt2 <- ggplotGrob(p20450.ncorrect.nl)
gt3 <- ggplotGrob(p20450.rt.nl)
gt1$widths[3] = gt2$widths[3]
gt3$widths[3] = gt2$widths[3]

# Create the empty gtable
gt = gtable(widths = unit(c(1.9, 1, 1), "null"), heights = unit(1, "null"))
# Put the plots into the empty gtable
gt <- gtable_add_grob(gt, gt1, 1, 1)
gt <- gtable_add_grob(gt, gt2, 1, 2)
gt <- gtable_add_grob(gt, gt3, 1, 3)
# Arrange the plots and save them
# Figure 2
png("KP_figures/Figure_2_typical_20450.png", height=350, width=400, units="mm", res=150)
grid.newpage()
grid.draw(gt)
dev.off()

###
# Typical individual with age=75 and current income=$450
###

# Extract the data which we want to plot from X.typical.list.pred
typical.75450 <- X.typical.list.pred[["tpa450age75"]]

# Add the adapted change.var label vector to the data frame
typical.75450$change.var.orig <- typical.75450$change.var
typical.75450$change.var <- factor(change.var, levels=rev(change.var))

p75450.ncor.rt <- ggplot(data=typical.75450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(seq(-0.1, -0.01,0.01), 0.01, 0.02), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.ncor.rt, ymin = lb.ncor.rt, ymax = ub.ncor.rt),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-0.1, 0.02,0.02)), limits=c(-0.10, 0.02)) +
  coord_flip() + theme_bw() + ggtitle("Panel A. Correct \n Answers per Second") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        axis.text.y = element_text(size=18),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

p75450.ncorrect.nl <- ggplot(data=typical.75450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(seq(-6, -1, 1), 1 , 2), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.ncorrect, ymin = lb.ncorrect, ymax = ub.ncorrect),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-6, 2, 1)), limits=c(-6.2, 2.5)) +
  coord_flip() + theme_bw() + ggtitle("Panel B. Number \n of Correct Answers") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        # The next command switches off the y text, which is switched on for the first graph
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

p75450.rt.nl <- ggplot(data=typical.75450, aes(colour = change.var)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = "dashed") +
  geom_hline(yintercept=c(-2, -1, seq(1,12)), colour="grey92") +
  geom_pointrange(aes(x = change.var, y = pred.rt, ymin = lb.rt, ymax = ub.rt),
                  lwd = 1, position = position_dodge(width = 1), shape = 20, colour="black") +
  scale_y_continuous(name=element_blank(), breaks=c(seq(-2, 12, 2)), limits=c(-3.6, 12.6)) +
  coord_flip() + theme_bw() + ggtitle("Panel C. Total \n Response Time") +
  theme(legend.position="none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle=element_text(size=25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=0, size=18, vjust=0),
        # The next command switches off the y text, which is switched on for the first graph
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title=element_text(size=29),
        plot.caption = element_text(size = 25, hjust = 0))

# The following three commands return the warning message: "Removed 15 rows containing missing values (geom_pointrange)."
# The warning messages result from the 15 "blank" rows in each data frame, which give information about the covariate values
# for the typical individual in its baseline specification -> warning messages can be ignored
gt1 <- ggplotGrob(p75450.ncor.rt)
gt2 <- ggplotGrob(p75450.ncorrect.nl)
gt3 <- ggplotGrob(p75450.rt.nl)
gt1$widths[3] = gt2$widths[3]
gt3$widths[3] = gt2$widths[3]

# Create the empty gtable
gt = gtable(widths = unit(c(1.9, 1, 1), "null"), heights = unit(1, "null"))
# Put the plots into the empty gtable
gt <- gtable_add_grob(gt, gt1, 1, 1)
gt <- gtable_add_grob(gt, gt2, 1, 2)
gt <- gtable_add_grob(gt, gt3, 1, 3)
# Arrange the plots and save them
# Figure 3
png("KP_figures/Figure_3_typical_75450.png", height=350, width=400, units="mm", res=150)
grid.newpage()
grid.draw(gt)
dev.off()

######################################
#### Script end
######################################
