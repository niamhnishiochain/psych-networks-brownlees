################################################################
#                                                              #
#                 Armour, Fried et al. 2017                    #
#                                                              #
#      A Network Analysis of DSM-5 posttraumatic stress        #
#  disorder symptoms and correlates in U.S. military veterans  #
#                                                              #
#             Journal of Anxiety Disorders                     #                   
#                                                              #
################################################################


# ---------------------------------------------------------------------------------------
# ---------- 1. Load libraries ----------------------------------------------------------
# ---------------------------------------------------------------------------------------

library(foreign)
# library("devtools")
# devtools::install_github("sachaepskamp/bootnet")
library(bootnet)
# devtools::install_github("sachaepskamp/qgraph")
library(qgraph)
# devtools::install_github("jmbh/mgm")
library(mgm)
install.packages('NetworkComparisonTest')
library('NetworkComparisonTest')


# ---------------------------------------------------------------------------------------
# ---------- 2. Data preparation --------------------------------------------------------
# ---------------------------------------------------------------------------------------

setwd('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/PTSD_paper')
data_full<-read.spss("PTSD_data_V1.sav", to.data.frame=TRUE)

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

# ---------------------------------------------------------------------------------------
# ---------- 3. Demographics ------------------------------------------------------------
# ---------------------------------------------------------------------------------------

mean(data$PPAGE)                      #54.97285
sd(data$PPAGE)                        #14.75371
range(data$PPAGE)                     #21-89

mean(data$PPGENDER)                   #1.171946
mean(data_full$SUM_LIFETIME_TRAUMAS)  # 5.954751
sd(data_full$SUM_LIFETIME_TRAUMAS)    #3.207624
range(data_full$SUM_LIFETIME_TRAUMAS) #1-15


# ---------------------------------------------------------------------------------------
# ---------- 4. Network 1: 20 PTSD symptoms ---------------------------------------------
# ---------------------------------------------------------------------------------------

graph_ptsd.m <-EBICglasso(data_ptsd.cor, n = nrow(data_ptsd)) 

### A. Figure 1
pdf("Fig1.pdf", width=9, height=7)
graph_ptsd.g<-qgraph(data_ptsd.cor, labels=names, graph="glasso", layout="spring", 
                     vsize=6, cut=0, maximum=.45, sampleSize = nrow(data_ptsd),
                     border.width=1.5, border.color="black", minimum=.03, 
                     groups=gr3, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94'),
                     nodeNames = longnames,legend.cex=.4)
dev.off()

### B. Figure 2 
pdf("Fig2.pdf")
centralityPlot(graph_ptsd.g)
dev.off()

### C. intercorrelation of different centrality measures
cor(c.ptsd$InDegree, c.ptsd$Closeness, method="spearman")     #0.69 spearman=0.65
cor(c.ptsd$InDegree, c.ptsd$Betweenness, method="spearman")   #0.80 spearman=0.72
cor(c.ptsd$Betweenness, c.ptsd$Closeness, method="spearman")  #0.81 spearman=0.82


### D. Robustness analysis
colnames(data_ptsd)<-c(1:20)
colnames(data)<-c(1:27)

network1<-estimateNetwork(data_ptsd, default = "EBICglasso") 
network2<-estimateNetwork(data, default = "EBICglasso", corMethod = 'cor_auto') #note that the last option was added in order to compute the correlations between ordinal and non-ordinal data (polychoric and polyserial)

# boot1 <- bootnet(network1, nBoots = 1000,  nCores = 4) # edge weights bootstrap network 1 
# boot2 <- bootnet(network2, nBoots = 1000,  nCores = 4) # edge weights bootstrap network 2
# boot3 <- bootnet(network1, nBoots = 1000, type = "case",  nCores = 4) # subsetting bootstrap network 1
# save(boot1, file = "boot1.Rdata")
# save(boot2, file = "boot2.Rdata")
# save(boot3, file = "boot3.Rdata")
load(file = "boot1.Rdata")
load(file = "boot2.Rdata")
load(file = "boot3.Rdata")

pdf("boot1.pdf", useDingbats=FALSE) # edge weights network 1
plot(boot1, labels = F, order = "sample" ) 
dev.off()

pdf("boot2.pdf", useDingbats=FALSE) # edge weights network 2
plot(boot2, labels = F, order = "sample" ) 
dev.off()

pdf("boot3.pdf", useDingbats=FALSE) # centrality network 1
plot(boot3)
dev.off()

corStability(boot3) # "CS-coefficient should not be below 0.25, and preferably above 0.5" (bootnet paper)
# betweenness   closeness    strength 
# 0.04977376    0.12669683   0.36199095  

# Gray boxes indicate nodes or edges that do not differ significantly from one-another 
# and black boxes represent nodes or edges that do differ significantly from one-another
boot4 <- plot(boot1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
boot5 <- plot(boot1, "strength", order="sample", labels=F)
boot6 <- plot(boot2, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")

pdf("boot4.pdf")
plot(boot4, useDingbats=FALSE) # differences among edges in 20-item network
dev.off()

pdf("boot5.pdf")
plot(boot5, useDingbats=FALSE) # differences among centrality in 20-item network
dev.off()

pdf("boot6.pdf")
plot(boot6, useDingbats=FALSE) # differences among centrality in 27-item network
dev.off()


# ---------------------------------------------------------------------------------------
# ---------- 5. Network 2: 20 PTSD symptoms + 7 covariates ------------------------------
# ---------------------------------------------------------------------------------------

graph.m <-EBICglasso(data.cor, n = nrow(data))



pdf("Fig3.pdf", width=10, height=8)
graph_ptsd.g<-qgraph(data.cor, graph="glasso", layout="spring", sampleSize = nrow(data),
                     vsize=6, cut=0, maximum=.45, minimum=.03, border.width=1.5, border.color="black", 
                     groups=gr4, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94', '#bbbbbb'),
                     labels=names1, nodeNames = longnames1,legend.cex=.4)
dev.off()

# to address reviewer comment::
sum1<-sum(abs(graph.m[upper.tri(graph.m,diag=FALSE)])); sum1     # 11.88; sum/2 of full network matrix

# ---------------------------------------------------------------------------------------
# ---------- 6. Examine changes upon adding covariates ----------------------------------
# ---------------------------------------------------------------------------------------

graph.m2<-graph.m[-c(21:27),-c(21:27)] # delete covariates from full network with covariates

delta<-graph_ptsd.m-graph.m2 #delta network
max(delta)
mean(abs(delta[upper.tri(delta,diag=FALSE)]))
colnames(delta)<-rownames(delta)<-c(1:20)

L<-averageLayout(graph_ptsd.m) # use same layout as Figure 1
pdf("Fig4.pdf", width=9, height=7)
qgraph(delta, labels=names, layout=L, vsize=6, cut=0, maximum=.45, 
       border.width=1.5, border.color="black", minimum=.03, 
       groups=gr3, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94'),
       nodeNames = longnames,legend.cex=.4)
dev.off()

cor(as.vector(graph.m2),as.vector(graph_ptsd.m), method="spearman") #0.97

### A. mean edge weight among symptoms=0.04; adding covariates reduces sum of edges among symptoms by 11.5%
mean(abs(delta[upper.tri(delta,diag=FALSE)]))
sum1<-sum(abs(graph_ptsd.m[upper.tri(graph_ptsd.m,diag=FALSE)])); sum1     # 8.68; sum/2 of ptsd network matrix
mean1<-mean(abs(graph_ptsd.m[upper.tri(graph_ptsd.m,diag=FALSE)])); mean1  # 0.05; mean edge strength
sum2<-sum(abs(graph.m2[upper.tri(graph.m2,diag=FALSE)])); sum2             # 7.78; sum/2 of ptsd + covariates matrix, with covariate cells deleted (ptsd network accounting for covariates but without them)
mean2<-mean(abs(graph.m2[upper.tri(graph.m2,diag=FALSE)])); mean2          # 0.04; mean edge strength
1-sum2/sum1 # change in connectivity of PTSD network once symptoms are added 11.5%

### B. mean edge weight among covariates=0.09
graph.m3<-graph.m[-c(1:20),-c(1:20)]   # just covariate connections of full network
sum3<-sum(abs(graph.m3[upper.tri(graph.m3,diag=FALSE)])); sum3            # 1.91
mean3<-mean(abs(graph.m3[upper.tri(graph.m3,diag=FALSE)])); mean3         # 0.09

### C. mean edge weight symptoms-covaraites=0.02
graph.m4<-graph.m[-c(1:20),-c(21:27)]  # just connection betw covariates and symptoms
sum4<-sum(abs(graph.m4[upper.tri(graph.m4,diag=FALSE)])); sum4            # 2.13
mean4<-mean(abs(graph.m4[upper.tri(graph.m3,diag=FALSE)])); mean4         # 0.02

# ---------------------------------------------------------------------------------------
# ---------- 7. Network comparison test  ----------------------------------
# ---------------------------------------------------------------------------------------


NCT_and = NCT(data, 
                         data_ptsd, 
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
              centrality = c('betweeness', 'closeness', 'strength'),
              #nodes, #optional list of nodes to select for testing centrality. 
                         verbose = TRUE)

NCT_or = NCT(network1, 
              network2, 
              gamma = 0.5, #same as that used by the paper (see that they use the default EBIC glasso. See bootnet documentation for Estimate Network - tuning - EBICglasso)
              it = 100, 
              binary.data = FALSE,
              paired=FALSE, #samples are independent
              weighted=TRUE,  #compare network weights
              AND = FALSE,  #authors don't say which they use so we run for both
              test.edges=TRUE, #test individual edges
              #edges,  #optional list of edges to test
              progressbar=TRUE,
              make.positive.definite=TRUE,
              p.adjust.methods="none", #not controlling for testing of multiple edges
              verbose = TRUE)





