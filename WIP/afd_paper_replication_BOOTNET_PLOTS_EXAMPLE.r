#[0] Upload packages
library(dplyr)
library(qgraph)
library(bootnet)
library(psych)
library(polycor)
library(networktools)
library(ggpubr)
library(foreign)
library(mgm) #predictability

#/////////////////////////////////////////////////////////////////////////////
#///////////////////////////[1] Prepare the data//////////////////////////////
#/////////////////////////////////////////////////////////////////////////////.
# Set working directory
setwd("")

#Upload data
data <- read.csv("data Avoidant.csv") 
table(data$Sample) # 718 non-clinical and 354 clinical

AVdata<-subset(data, Sample==1,select=c(AV1:AV7))
AVdataCli<-subset(data, Sample==2,select=c(AV1:AV7))

colnames(AVdata) <- colnames(AVdataCli) <- c("Avo", "Like", "Res", "FCR", "Inh", "Inf", "Risk")

#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////[2] ESTIMATING NETWORKS/////////////////////////////
#/////////////////////////////////////////////////////////////////////////////.
#Identify redundant nodes
goldbricker(AVdata)
goldbricker(AVdataCli)

#Networks
# A_netAV.gen <- estimateNetwork(AVdata, default = "IsingFit") 
# A_netAV.cli <- estimateNetwork(AVdataCli, default = "IsingFit") 

A_netAV.gen <- estimateNetwork(AVdata, default = "IsingSampler") 
A_netAV.cli <- estimateNetwork(AVdataCli, default = "IsingSampler")


#/////////////////////////////////////////////////////////////////////////////
#//////////////////////////[3] PLOTTING NETWORKS//////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
###############Non-clinical
set.seed(100)
fitGen <- mgm(data = na.omit(AVdata),
              type = rep("c",7),
              level = rep(2,7),
              ruleReg = 'AND',
              k = 2, 
              binarySign = F, 
              overparameterize = FALSE)
predGen <- predict(fitGen, na.omit(AVdata))
predGen$errors


#Graph
B_plotAV.gen <- qgraph(A_netAV.gen$graph,labels = colnames(AVdata),
                       nodeNames = c("Avoidance","Certainty of being liked","Restaint in relationships","Fear of criticism and rejection","Inhibition","Inferiority","Reluctance to risk"),
                       title ="Non-clinical",title.cex = 2,
                       # filetype = "png",filename = "B_plot_AV.un",
                       aspect = F,
                       #edge.label.position =0.5,
                       edge.labels = F,edge.label.bg = F,edge.label.color = "black",edge.label.cex = 0.7,
                       legend = FALSE,legend.mode = "names",legend.cex = 0.8,
                       negCol = "blue",vsize = 8,layout = "spring",
                       layoutScale = c(0.9,0.9),GLratio = 1.3,
                       pie=predGen$errors[,3])
x<-averageLayout(B_plotAV.gen)



###############Clinical
set.seed(100)
fitCli <- mgm(data = na.omit(AVdataCli), 
              type = rep("c",7),
              level = rep(2,7),
              ruleReg = 'AND',
              k = 2, 
              binarySign = F, 
              overparameterize = FALSE)
predCli <- predict(fitCli, na.omit(AVdataCli))
predCli$errors


B_plotAV.cli <- qgraph(A_netAV.cli$graph,
                       labels = colnames(AVdata),
                       nodeNames = c("Avoidance","Certainty of being liked","Restaint in relationships","Fear of criticism and rejection","Inhibition","Inferiority","Reluctance to risk"),
                       title ="Clinical",title.cex = 2,
                       # filetype = "png",filename = "B_plot_AV.25",
                       aspect = F,
                       #edge.label.position =0.5,
                       edge.labels = F,edge.label.bg = F,edge.label.color = "black",edge.label.cex = 0.7,
                       legend = F,legend.mode = "names",legend.cex = 0.8,
                       negCol = "blue",vsize = 8,layout = x,
                       layoutScale = c(0.9,0.9),GLratio = 1.3,
                       pie=predCli$errors[,3])


tiff(width=15, height=7, filename="network.tiff",res=300, units="in")
par(mfrow=c(1,2))
plot(B_plotAV.gen)
plot(B_plotAV.cli)
dev.off()
getwd()

##### Gray scale
tiff(width=15, height=7, filename="network gray.tiff",res=300, units="in")
par(mfrow=c(1,2))
qgraph(B_plotAV.gen,theme="gray",edge.color="#333333")
qgraph(B_plotAV.cli,theme="gray",edge.color="#333333")
dev.off()



#/////////////////////////////////////////////////////////////////////////////
#////////////////////////////[4] NETWORK INDICES//////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#[4.1] Smallworldness
set.seed(100)
smallworldness(B_plotAV.gen)
smallworldness(B_plotAV.cli)

#[4.2] Clustering coefficients
clustPlot<-clusteringPlot(list("Non-clinical" = B_plotAV.gen,
                               "Clinical" = B_plotAV.cli), include="Zhang", signed=T)


######## CENTRALITY PLOTS ########
C_cent_AV <- centralityPlot(list("Non-clinical" = B_plotAV.gen,
                                         "Clinical" = B_plotAV.cli),
                            include = c("Strength","Closeness", "Betweenness"),
                            scale="z-scores")


###########SAVE PLOTS###########
p1<-C_cent_AV+theme(legend.position = "none")
p2<-clustPlot+theme(legend.position = "right",legend.title = element_blank())

p3<-ggarrange(p1,p2,ncol = 2, nrow = 1,widths=c(1,0.815))
p3

ggsave("cent_clust.tiff", plot = p3,
       width = 17, height=10, dpi=400, units="cm")

#Gray scale
p3_<-ggarrange(p1+scale_colour_grey(),
              p2+scale_colour_grey(),ncol = 2, nrow = 1,widths=c(1,0.815))

ggsave("cent_clust gray.tiff", plot = p3_,
       width = 17, height=10, dpi=400, units="cm")




#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#///////////////////////////////////////BOOTSTRAPPING///////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
set.seed(100)
#@@@@@@@@@@@@@@@@@(A) EDGE WEIGHT ACCURACY@@@@@@@@@@@@@@@@@.
D_boot1_AV.gen <- bootnet(A_netAV.gen, 2000, default = "IsingSampler", type = "nonparametric", nCores = 6)
D_boot1_AV.cli <- bootnet(A_netAV.cli, 2000, default = "IsingSampler", type = "nonparametric", nCores = 6)
plot(D_boot1_AV.gen, labels = TRUE, order = "sample")
plot(D_boot1_AV.cli, labels = TRUE, order = "sample")


p1<-plot(D_boot1_AV.gen, labels = TRUE, order = "sample")
p2<-plot(D_boot1_AV.cli, labels = TRUE, order = "sample")

p3<-ggarrange(p1,p2,
          labels = c("Non-clinical", "Clinical"),
          ncol = 2, nrow = 1,vjust=1.2)
p3
ggsave("D_boot1_nonparam.tiff", plot = p3,
       width = 17, height=14, dpi=400, units="cm")



#@@@@@@@@@@@@@@@@@(B) CENTRALITY STABILITY@@@@@@@@@@@@@@@@@.
set.seed(100)
G_boot2_AV.gen <- bootnet(A_netAV.gen, nboots=2000, statistics = c("strength", "closeness", "betweenness","edge"),  
                         default = "IsingSampler", type = "case", nCores = 6)
G_boot2_AV.cli <- bootnet(A_netAV.cli, nboots=2000, statistics = c("strength", "closeness", "betweenness","edge"),  
                         default = "IsingSampler", type = "case", nCores = 6)

corStability(G_boot2_AV.gen,statistics = c("strength", "closeness", "betweenness"))
corStability(G_boot2_AV.cli,statistics = c("strength", "closeness", "betweenness"))



p1<-plot(G_boot2_AV.gen,statistics = c("strength", "closeness", "betweenness"))+
  labs(x = "", y = "")+ ggtitle("Non-clinical")+theme(legend.position = "top")

p2<-plot(G_boot2_AV.cli,statistics = c("strength", "closeness", "betweenness"))+
  labs(y = "")+ ggtitle("Clinical")+theme(legend.position = "none")

p3<-ggarrange(p1,p2,
              ncol = 1, nrow = 2,
              heights = c(1,0.836))

p3<-annotate_figure(p3,left = text_grob("Average correlation with original sample", rot = 90,vjust=2))


ggsave("D_boot2_case.tiff", plot = p3,
       width = 17, height=18, dpi=400, units="cm")

#@@@@@@@@@@@@@@@@@(C) EDGE STABILITY@@@@@@@@@@@@@@@@@.
p1_edge<-plot(G_boot2_AV.gen,statistics = c("edge"))+
  labs(x = "", y = "")+ 
  ggtitle("Non-clinical")+
  theme(legend.position = "none")

p2_edge<-plot(G_boot2_AV.cli,statistics = c("edge"))+
  labs(y = "")+ 
  ggtitle("Clinical")+
  theme(legend.position = "none")


p3_edge<-ggarrange(p1_edge,p2_edge, ncol = 1, nrow = 2, heights = c(1,1))

p3_edge<-annotate_figure(p3_edge,left = text_grob("Average correlation with original sample", rot = 90,vjust=2))

tiff("D_boot2_case_edge.tiff",res=300, width = 480*3.5, height = 480*2.7)
p3_edge
dev.off()


#@@@@@@@@@@@@@@@@@(D) Testing for significant differences@@@@@@@@@@@@@@@@@
H_edgediff_AV.gen<-plot(D_boot1_AV.gen, "edge", plot = "difference",onlyNonZero = TRUE, order = "sample")	
H_edgediff_AV.cli<-plot(D_boot1_AV.cli, "edge", plot = "difference",onlyNonZero = TRUE, order = "sample")

p1<-H_edgediff_AV.gen+ ggtitle("Non-clinical")
p2<-H_edgediff_AV.cli+ ggtitle("Clinical")

p3<-ggarrange(p1,p2,
              ncol = 2, nrow = 1)

ggsave("H_edgediff.tiff", plot = p3,
       width = 24, height=12, dpi=400, units="cm")

#between nodes
#if a node strength is significantly different from other nodes strength
I_nodediff_AV.gen<-plot(D_boot1_AV.gen, "strength")
I_nodediff_AV.cli<-plot(D_boot1_AV.cli, "strength")

p1<-I_nodediff_AV.gen+ ggtitle("Non-clinical")
p2<-I_nodediff_AV.cli+ ggtitle("Clinical")

p3<-ggarrange(p1,p2,
              ncol = 2, nrow = 1)

ggsave("H_nodediff.tiff", plot = p3,
       width = 24, height=12, dpi=300, units="cm")



#################################################################
#Network connectivity comparison
#################################################################
x1 <- centrality(A_netAV.gen)$ShortestPathLengths
x2 <- centrality(A_netAV.cli)$ShortestPathLengths
y1 <- x1[upper.tri(x1)]
y2 <- x2[upper.tri(x2)]
mean(y1)
mean(y2)

library(car)
df<-data.frame(con=c(y1,y2),group=c(rep(1,length(y1)),rep(2,length(y1))))
leveneTest(data=df,con~factor(group))
t.test(y1,y2,alternative = "two.sided", var.equal = T)


##########Network comparison test
library(NetworkComparisonTest)
NetCompar<-NCT(A_netAV.gen, A_netAV.cli, it=1000,binary.data= TRUE, 
                 test.edges	 = T,edges="all",
                 test.centrality=T,centrality=c("strength","closeness","betweenness"))
NetCompar$diffcen.real
NetCompar$diffcen.pval







