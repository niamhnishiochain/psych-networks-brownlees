library("devtools")
install_github("sachaepskamp/psychonetrics")
library("psychonetrics")
library("dplyr")

df <- read.csv('/Users/benseimon/Documents/Barca GSE/Studies/Term 3/Advanced Techniques/term_paper/psych-networks-brownlees/ESM_using.csv')
keep = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy',
         'day_var', 'beep_var', 'id_var')


vars = c('Control', 'Sad', 'Intrusion', 'Encourage', 
         'FunThings', 'IgotThis', 'Alone', 'Thoughts', 
         'Appointments', 'Laying',
         'Out', 'Eaten', 'Avoid', 'Useful', 'Enjoy'
)

df = df[keep]


test = var1(data, contemporaneous = c("ggm"), beta = "full", omega_zeta = "full", delta_zeta
     = "full", beepvar = 'beep_var', dayvar = 'day_var', idvar = 'id_var',
     vars = vars, estimator =
       "FIML", storedata = TRUE, verbose = TRUE)


