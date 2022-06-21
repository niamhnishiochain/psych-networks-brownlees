library('graphicalVAR')

df <- read.csv('C:/Users/35387/Dropbox/Documents/Brownlees/final project/git/psych-networks-brownlees/ESM_using.csv')

vars = c('Control') # 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying', 'Qualitysleep',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy', 
         'Pleasant')

keep = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying', 'Out',
         'Eaten', 'Avoid', 'Useful', 'Enjoy')
         #'day_var', 'beep_var', 'id_var')


fresh <- df[keep]
fresh_sub <- fresh[1:25,]

dayvar = paste(df[1:25,]$day_var)
beepvar = paste(df[1:25,]$beep_var)
idvar = paste(df[1:25,]$id_var)



graphicalVAR(fresh_sub, nLambda = 50, verbose = TRUE, 
             gamma = 0.5, 
             scale = TRUE, 
             maxit.in = 100, 
             maxit.out = 100, 
             deleteMissings = TRUE, 
             penalize.diagonal = TRUE, 
             lambda_min_kappa = 0.05,
             lambda_min_beta = 0.05, 
             vars,
             beepvar, 
             dayvar, 
             idvar, 
             lags = 1, 
             centerWithin = TRUE,
             likelihood = c("unpenalized", "penalized"))


#trying stuff
no_time <- subset(df, select=-c(SentAt,StartedAt))

trial1 = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying', 'Qualitysleep',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy', 
         'Pleasant')
bigger <- df[trial1]


trial2 <- c("Control", "Sad")
small <- df[trial2]


graphicalVAR(df, nLambda = 50, verbose = TRUE, 
             gamma = 0.5, 
             scale = TRUE, 
             maxit.in = 100, 
             maxit.out = 100, 
             deleteMissings = TRUE, 
             penalize.diagonal = TRUE, 
             lambda_min_kappa = 0.05,
             lambda_min_beta = lambda_min_kappa, 
             vars,
             beepvar, 
             dayvar, 
             idvar, 
             lags = 1, 
             centerWithin = TRUE,
             likelihood = c("unpenalized", "penalized"))

graphicalVAR(bigger, gamma = 0.5, nLambda = 5)


