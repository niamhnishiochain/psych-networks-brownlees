library('graphicalVAR')

df <- read.csv('C:/Users/35387/Dropbox/Documents/Brownlees/final project/git/psych-networks-brownlees/ESM_using.csv', stringsAsFactors = TRUE)

df$day_var = as.factor(df$day_var)
df$beep_var = as.factor(df$beep_var)
df$id_var = as.factor(df$id_var)

vars = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy')

keep = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying', 'Out',
         'Eaten', 'Avoid', 'Useful', 'Enjoy',
         'day_var', 'beep_var', 'id_var')


fresh <- df[keep]
fresh_sub <- fresh[1:25,]

dayvar = df[1:25,]$day_var
beepvar = df[1:25,]$beep_var
idvar = df[1:25,]$id_var

lambda_beta = c(0.01)
lambda_kappa = c(0.01)


graphicalVAR(fresh_sub, nLambda = 50, verbose = TRUE, 
             gamma = 0.5, 
             scale = TRUE, 
             lambda_beta,
             lambda_kappa,
             #regularize_mat_beta, 
             #regularize_mat_kappa,
             maxit.in = 100, 
             maxit.out = 100, 
             deleteMissings = TRUE, 
             penalize.diagonal = TRUE,
             lambda_min_kappa = 0.05,
             lambda_min_beta = 0.05, 
             mimic = c("current","0.1.2","0.1.4","0.1.5","0.2"),
             vars,
             beepvar, 
             dayvar, 
             idvar, 
             lags = 1, 
             centerWithin = TRUE,
             likelihood = c("unpenalized", "penalized"))
