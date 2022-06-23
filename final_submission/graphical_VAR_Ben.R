library('graphicalVAR')
library('dplyr')
library('igraph')
library(qgraph)

# ---------------------------------------------------------------------------------------
# ----------  Data preparation --------------------------------------------------------
# ---------------------------------------------------------------------------------------

#read in data
#NIAMH's version
df <- read.csv("C:/Users/35387/Dropbox/Documents/Brownlees/final project/git/psych-networks-brownlees/ESM_first.csv")

#BEN's version
#df <- read.csv('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/ESM_first.csv')


##make two sets of variables
cont_vars = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy')
bin_vars = c('Company', 'Safe', 'Cancel', 'Sleep', 'Pleas', 'Unpl')

#create subsets of dfs
df_all = df[c(cont_vars, bin_vars)] 
df_cont = df[cont_vars]
df_bin = df[bin_vars]

# ---------------------------------------------------------------------------------------
# ----------  Run var models --------------------------------------------------------
# ---------------------------------------------------------------------------------------

#gvar for all data
gvar_all <- graphicalVAR(df_all, 
             nLambda = 50, 
             verbose = TRUE, 
             gamma = 0.5, 
             scale = TRUE, 
             #lambda_beta,
             #lambda_kappa,
             #regularize_mat_beta,
             #regularize_mat_kappa,
             maxit.in = 100, 
             maxit.out = 100, 
             deleteMissings = FALSE, 
             penalize.diagonal = TRUE, 
             lambda_min_kappa = 0.05,
             lambda_min_beta = 0.05, 
             #mimic,
             #vars,
             #beepvar =  "beep_var", 
             #dayvar = 'day_var', 
             #idvar = "id_var", 
             lags = 1, 
             centerWithin = FALSE,
             likelihood = "penalized"
             )

#gvar for continuous data
gvar_cont <- graphicalVAR(df_cont, 
                         nLambda = 50, 
                         verbose = TRUE, 
                         gamma = 0.5, 
                         scale = TRUE, 
                         #lambda_beta,
                         #lambda_kappa,
                         #regularize_mat_beta,
                         #regularize_mat_kappa,
                         maxit.in = 100, 
                         maxit.out = 100, 
                         deleteMissings = FALSE, 
                         penalize.diagonal = TRUE, 
                         lambda_min_kappa = 0.05,
                         lambda_min_beta = 0.05, 
                         #mimic,
                         #vars,
                         #beepvar =  "beep_var", 
                         #dayvar = 'day_var', 
                         #idvar = "id_var", 
                         lags = 1, 
                         centerWithin = FALSE,
                         likelihood = "penalized"
)



# ---------------------------------------------------------------------------------------
# ----------  Plots --------------------------------------------------------
# ---------------------------------------------------------------------------------------

#generate plot
#for plot
names_all = c(names(df_all))
names_cont = c(names(df_cont))


plot(gvar_all, include = c('PCC', 'PDC'), titles = TRUE, sameLayout = TRUE
     #, graph = 'glasso'
     , layout="spring", 
     vsize=6, cut=0, maximum=.45, sampleSize = nrow(df), border.width=1.5, border.color="black", minimum=.03
     #, groups=gr3, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94')
     ,nodeNames = names_all ,legend.cex=.4
     )

plot(gvar_cont, include = c('PDC'), titles = TRUE, sameLayout = TRUE
     #, graph = 'glasso'
     , layout="spring", 
     vsize=6, cut=0, maximum=.45, sampleSize = nrow(df), border.width=1.5, border.color="black", minimum=.03
     #, groups=gr3, color=c('#a8e6cf', '#dcedc1', '#ffd3b6', '#ff8b94')
     ,nodeNames = names_cont ,legend.cex=.4, filetype = 'png', filename = 'correct_plot_pdc'
)


# ---------------------------------------------------------------------------------------
# ----------  Community detection --------------------------------------------------------
# ---------------------------------------------------------------------------------------

#we cannot have negative correlations for community detection, so take the absolute value of the adj matrix
cont_adj_matrix = abs(gvar_cont$PCC)

#use igraph to create graph from adjacency matrix. 
#need this because only igraph does community detection in R
i_cont = graph_from_adjacency_matrix(
  cont_adj_matrix,
  mode = 'undirected',
  weighted = TRUE,
  diag = TRUE,
  #add.colnames = names_cont,
  add.rownames = NA
)

#cluster louvain algorithm 
cont_louvain = cluster_louvain(i_cont)

#plot again with communities

##note that need to do mapping of cluster output to nodes to plot communities
#manually
groups_man <- list("C1"=c(1,2, 3, 7, 8), "C2"=c(4, 5, 6, 9, 12, 13, 15), "C3"=c(10, 11, 14))

plot(gvar_cont, include = c('PCC'), 
     #, graph = 'glasso'
      layout="spring", 
     vsize=6, cut=0, maximum=.45, sampleSize = nrow(df), border.width=1.5, border.color="black", minimum=.03
     , groups=groups_man, color=c('gold2', 'skyblue', 'mediumpurple1')
     ,nodeNames = names_cont ,legend.cex=.4,
     filename = 'legend_comms', filetype = 'png', 
     layoutOffset = c(50, 0), title = 'PCC with Louvain communities')





#########
# Making centrality plots
#########

df_cont.cor<-cor_auto(df_cont)
df_all.cor <- cor_auto(df_all)
### A. Figure 1
graph_cont.g<-qgraph(df_cont.cor, graph="glasso", sampleSize = nrow(df_cont))
pdf("Centrality_Cont_ESM.pdf")
centralityPlot(graph_cont.g, include = c('Strength', 'Closeness', 'Betweenness'))
dev.off()
clusteringPlot(graph_cont.g) 


graph_all.g<-qgraph(df_all.cor) graph="glasso", sampleSize = nrow(df_all))
pdf("Centrality_All_ESM.pdf")
centralityPlot(graph_all.g, include = c('Strength', 'Closeness', 'Betweenness'))
dev.off()
