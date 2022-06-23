library(foreign)
library("devtools")
devtools::install_github("sachaepskamp/bootnet")
#install.packages('bootnet')
library(bootnet)
devtools::install_github("sachaepskamp/qgraph")
#install.packages('qgraph')
library(qgraph)
devtools::install_github("jmbh/mgm")
#install.packages('mgm')
library(mgm)
install.packages('NetworkComparisonTest')
library('NetworkComparisonTest')

setwd('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/PTSD_paper')
data_full<-read.spss("PTSD_data_V1.sav", to.data.frame=TRUE)

#-------------------------------------------------------------#
###########Clean Data Frame as per Armour et al
#-------------------------------------------------------------#

# Create data frame of variables that we need for analysis: 
# age, gender, all 20 symptoms (_MONTH extension), Sum_GAD2, Sum_PHQ2, Passive_SI_FINAL, Active_SI_FINAL, PCS, MCS, QualityofLife_SUM)
data<-as.data.frame(data_full[,c(2:3,11:37)]) # data for analysis
data <- subset(data, select=c(Q28_01_MONTH:QualityofLife_SUM, PPAGE, PPGENDER)) #put covariates last
data[, 29]<-ifelse(data$PPGENDER=='Male',1,2) # 1 is male; 2 is female
data<-data[,-c(23,27)] # Remove SI_passive & QoL, we decided later we don't want them
data_ptsd<-as.data.frame(data[,c(1:20)]) # PTSD items only

### A. Make correlation matrix
data_ptsd.cor<-cor_auto(data_ptsd) #only ptsd symptoms
data.cor<-cor_auto(data) #data including covariates


### B. Create names object & group object for graph
names<-c("B1","B2", "B3", "B4", "B5", "C1", "C2", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "E1", "E2", "E3", "E4", "E5", "E6")
names1<-c("B1","B2", "B3", "B4", "B5", "C1", "C2", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "E1", "E2", "E3", "E4", "E5", "E6", "Anx", "Dep", "SI", "PFunc", "Mfunc", "Age", "Sex")
longnames <- c("Intrusive thoughts", "Nightmares", "Flashbacks", "Emotional cue reactivity", "Physiological cue reactivity", "Avoidance of thoughts", "Avoidance of reminders", "Trauma-related amnesia", "Negative beliefs", "Blame of self or others", "Negative trauma-related emotions", "Loss of interest", "Detachment", "Restricted affect", "Irritability/anger", "Self-destructive/reckless behavior", "Hypervigilance", "Exaggerated startle response", "Difficulty concentrating", "Sleep disturbance")
longnames1 <- c("Intrusive thoughts", "Nightmares", "Flashbacks", "Emotional cue reactivity", "Physiological cue reactivity", "Avoidance of thoughts", "Avoidance of reminders", "Trauma-related amnesia", "Negative beliefs", "Blame of self or others", "Negative trauma-related emotions", "Loss of interest", "Detachment", "Restricted affect", "Irritability/anger", "Self-destructive/reckless behavior", "Hypervigilance", "Exaggerated startle response", "Difficulty concentrating", "Sleep disturbance", "Anxiety", "Depression", "Suicidal ideation", "Physical functioning", "Mental functioning", "Age", "Sex")

gr1 <- list(c(21:27), c(1:20))                         #first covariates, second group ptsd symptoms
gr2 <- list(c(21:27), c(1:5), c(6:7),c(8:14),c(15:20)) #first covariates, then PTSD symptoms categories B C D E
gr3 <- list("Intrusions"=c(1:5), "Avoidance"=c(6:7), "Cognition & mood alterations"=c(8:14),"Arousal & reactivity alterations"=c(15:20)) #PTSD symptoms categories B C D E
gr4 <- list("Intrusions"=c(1:5), "Avoidance"=c(6:7), "Cognition & mood alterations"=c(8:14),"Arousal & reactivity alterations"=c(15:20), "Clinical covariates"=c(21:27))           


#-------------------------------------------------------------#
###########Relevant elements for NCT from Armour et al
#-------------------------------------------------------------#

#note that data is 27 covariates, data_ptsd is 20

## EBICglasso outputs
graph_ptsd.m <-EBICglasso(data_ptsd.cor, n = nrow(data_ptsd)) 
graph.m <-EBICglasso(data.cor, n = nrow(data))



## Qgraph objects
graph_ptsd.g<-qgraph(data_ptsd.cor, labels=names, graph="glasso", layout="spring", 
                     vsize=6, cut=0, maximum=.45, sampleSize = nrow(data_ptsd),
                     border.width=1.5, border.color="black", minimum=.03, 
                     groups=gr3, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94'),
                     nodeNames = longnames,legend.cex=.4)

graph.g<-qgraph(data.cor, graph="glasso", layout="spring", sampleSize = nrow(data),
                vsize=6, cut=0, maximum=.45, minimum=.03, border.width=1.5, border.color="black", 
                groups=gr4, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94', '#bbbbbb'),
                labels=names1, nodeNames = longnames1,legend.cex=.4)

## bootnet estimate network object
network_ptsd <-estimateNetwork(data_ptsd, default = "EBICglasso", corMethod = 'cor_auto') 
network <-estimateNetwork(data, default = "EBICglasso", corMethod = 'cor_auto') #note that the last option was added in order to compute the correlations between ordinal and non-ordinal data (polychoric and polyserial)

#now i need to adjust the dimensions of the network in order to have the same elements for comparison

network$graph = network$graph[1:20, 1:20]
network$results$results$w = network$results$results$w[1:20, 1:20, ]
network$results$results$wi = network$results$results$wi[1:20, 1:20, ]
network$results$optnet =network$results$optnet[1:20, 1:20]
network$results$optwi = network$results$optwi[1:20, 1:20]
network$labels = network$labels[1:20]
network$data$Sum_GAD2 = NULL
network$data$Sum_PHQ2 = NULL
network$data$Active_SI_FINAL = NULL
network$data$PCS = NULL
network$data$MCS = NULL
network$data$PPAGE = NULL
network$data$PPGENDER = NULL


# ---------------------------------------------------------------------------------------
########### NCT
# ---------------------------------------------------------------------------------------

NCT_test = NCT(network_ptsd, network)

NCT_and = NCT(network, #takes either data or estimateNetwork object as input. We need the latter
              network_ptsd, #takes either data or estimateNetwork object as input. We need the latter
              gamma = 0.5, #same as that used by the paper (see that they use the default EBIC glasso. See bootnet documentation for Estimate Network - tuning - EBICglasso)
              it = 100, 
              binary.data = FALSE,
              paired=FALSE, #samples are independent
              weighted=TRUE,  #compare network weights
              AND = TRUE, #authors don't say which they use so we run for both
              abs = FALSE, #not sure about this one
              test.edges=TRUE, #test individual edges
              #edges,  #optional list of edges to test
              progressbar=TRUE,
              make.positive.definite=TRUE,
              p.adjust.methods="none", #not controlling for testing of multiple edges
              test.centrality = TRUE,
              centrality = c('betweenness', 'closeness', 'strength'),
              #nodes, #optional list of nodes to select for testing centrality. 
              verbose = TRUE)

NCT_or = NCT(network, #takes either data or estimateNetwork object as input. We need the latter
              network_ptsd, #takes either data or estimateNetwork object as input. We need the latter
              gamma = 0.5, #same as that used by the paper (see that they use the default EBIC glasso. See bootnet documentation for Estimate Network - tuning - EBICglasso)
              it = 100, 
              binary.data = FALSE,
              paired=FALSE, #samples are independent
              weighted=TRUE,  #compare network weights
              AND = TRUE, #authors don't say which they use so we run for both
              abs = FALSE, #not sure about this one
              test.edges=TRUE, #test individual edges
              #edges,  #optional list of edges to test
              progressbar=TRUE,
              make.positive.definite=TRUE,
              p.adjust.methods="none", #not controlling for testing of multiple edges
              test.centrality = TRUE,
              centrality = c('betweenness', 'closeness', 'strength'),
              #nodes, #optional list of nodes to select for testing centrality. 
              verbose = TRUE)









