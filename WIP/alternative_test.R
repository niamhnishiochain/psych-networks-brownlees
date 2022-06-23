library("devtools")
install_github("sachaepskamp/psychonetrics")
library("psychonetrics")
library("dplyr")

setwd('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/')
Data <- read.csv("Supplementary2_data.csv")

# Variables to use:
Vars <- c("relaxed", "sad", "nervous",
          "concentration", "tired",
          "rumination", "bodily.discomfort")

# Encode time variable in a way R understands:
Data$time <- as.POSIXct(Data$time, tz = "Europe/Amsterdam", format= '%m/%d/%Y %H:%M' )

# Extract days:
Data$Day <- as.Date(Data$time, tz = "Europe/Amsterdam", format= '%m/%d/%Y %H:%M' )

# Model, using FIML for missing data:
mod <- gvar(Data, vars = Vars,  beta = "full", 
            omega_zeta = "full", estimator = "FIML")

mod <- mod %>% runmodel %>% prune %>% stepup(criterion = "bic")