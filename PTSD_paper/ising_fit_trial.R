library(IsingFit)
library(foreign)
library(bootnet)
library(qgraph)
library(mgm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(fastDummies)


data_full <- read.spss('C:/Users/35387/Dropbox/Documents/Brownlees/final project/git/psych-networks-brownlees/PTSD_paper/PTSD_data_V1.sav', to.data.frame=TRUE)

#subset to only binary variables

#original approach with the four binary variables 
#binary_vars <- c('PPGENDER', 'PM_PCL5_38', 'PTSD_2cat', 'Any_SI')
#binary_df <- data_full[binary_vars]

#want more variables so
#checked manually (lapply(data_full, table))

#plot and see if any we can force to be binary
#data_full %>% select_if(is.numeric) %>%  gather() %>%  
#  ggplot(aes(value)) + facet_wrap(~ key, scales = "free") +
#  geom_histogram()

#use the fast dummies package to make dummies of the ones we see are suitable
not_dummy <- c('CaseID', 'MCS', 'PCS', 'QualityofLife_SUM',
               'weight', 'Subthreshold_or_Higher_PM_PCL', 
               'SUM_PM_PCL5')
dummy_data <- select(data_full, -c(not_dummy)) #31 vars
dummy_data$PPAGE <- as.numeric(dummy_data$PPAGE)
dummy_data$PPGENDER <- as.numeric(dummy_data$PPGENDER)
dummy_data <-  fastDummies::dummy_cols(dummy_data, select_columns = names(dummy_data),
                                       remove_selected_columns = TRUE,
                                       remove_first_dummy = TRUE)

X <- as.matrix(sapply(dummy_data, as.numeric))  

IsingFit(X, AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE,
         lowerbound.lambda = NA, label.cex = 8)