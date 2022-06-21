library('graphicalVAR')
library('dplyr')

df <- read.csv('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/ESM_data_github_issue.csv')

vars = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy', 
         'day_var', 'id_var', 'beep_var')

df = df[vars] #in order to drop 'SentAt'


x <- graphicalVAR(df, 
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
             vars,
             beepvar =  "beep_var", 
             dayvar = 'day_var', 
             idvar = "id_var", 
             lags = 1, 
             centerWithin = FALSE,
             likelihood = c("unpenalized", "penalized"))
