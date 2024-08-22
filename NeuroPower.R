# RD posthoc neuro power
# Dr. Nicholas Judd
# 21-08-24 njudd.com



### notes

# example script: https://github.com/rdpackages/rdpower/blob/master/R/rdpower_illustration.R
# paper: https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf#page=18.40

# https://khakieconomics.github.io/2017/11/26/Bayesian_iv.html
# Bayesian IV-power




if (!require(pacman)){install.packages('pacman')}
pacman::p_load(tidyverse, lubridate, stringr, lm.beta,
               RDHonest, rdrobust, rdpower, data.table)

telomere_set <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/telomere_set.csv")
# telomere_set <- telomere_set[complete.cases(telomere_set[, .(ltl, running_var, EduAge16)]),] # no missingness allows in Y, running, or IV
# ^there already is none!

neuro_set <- fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_fullset.csv")

sa <- neuro_set[, c("SA", "running_var", "EduAge16", "EduAge"), with=FALSE]
ct = neuro_set[, c("CT", "running_var", "EduAge16", "EduAge"), with=FALSE]
WMh = neuro_set[, c("WM_hyper", "running_var", "EduAge16", "EduAge"), with=FALSE]
CSFn = neuro_set[, c("CSF_norm", "running_var", "EduAge16", "EduAge"), with=FALSE]
TBVn = neuro_set[, c("TBV_norm", "running_var", "EduAge16", "EduAge"), with=FALSE]
wFAs = neuro_set[, c("wFA", "running_var", "EduAge16", "EduAge"), with=FALSE] # wFAs; s for subset since you already have a wFA DT

# not allowing any missing in Y, running or instrument
sa <- sa[complete.cases(sa[, .(SA, running_var, EduAge16)]),] # no missingness allows in Y, running, or IV
# removed the fencing arguemnt for outliers
ct <- ct[complete.cases(ct[, .(CT, running_var, EduAge16)]),]
WMh <- WMh[complete.cases(WMh[, .(WM_hyper, running_var, EduAge16)]),]
CSFn <- CSFn[complete.cases(CSFn[, .(CSF_norm, running_var, EduAge16)]),]
TBVn <- TBVn[complete.cases(TBVn[, .(TBV_norm, running_var, EduAge16)]),]
wFAs <- wFAs[complete.cases(wFAs[, .(wFA, running_var, EduAge16)]),]


saEdu_mod <- lm(SA ~ EduAge, data = sa)
saEdu_mod_noCollege <- lm(SA ~ EduAge, data = sa[EduAge<21])

summary(lm.beta(saEdu_mod)); summary(lm.beta(saEdu_mod_noCollege))

rdpower(data=cbind(sa$SA, sa$running_var), fuzzy = sa$EduAge16)
rdpower(data=cbind(sa$SA[sa$EduAge <21], sa$running_var[sa$EduAge <21]), fuzzy = sa$EduAge16[sa$EduAge <21])


# the 2nd-stage is what is fucking you up... too many people before did the exp
rdpower(data=cbind(sa$SA, sa$running_var)) # loads of power when not fuzzy
table(sa$EduAge16, sa$running_var>0)

# it is only 2216 people that didn't go on to education at 16

# artificially lowering the amount before that continued
fake_2ndStage <- sa$EduAge16
fake_2ndStage[sa$running_var <0] <- sample(0:1, length(fake_2ndStage[sa$running_var <0]), replace = T, prob = c(.5,.5))
rdpower(data=cbind(sa$SA, sa$running_var), fuzzy = fake_2ndStage) # power =.998

# illustrating this effect; with standardizing the outcome (Y)

# function to get power, takes the percentage of people before the NatExp not going & a std effect
secound_stage_STDpower <- function(y = NULL, running = NULL, stage_two = NULL, pi_notgoing = .5, std_eff = .5){
  
  # if none of y, running or stage two are null & they have the same dimensions...
  if (is.null(y) + is.null(running) + is.null(stage_two) == 0 && dim(table(lengths(list(y, running, stage_two)))) == 1){
    pi_going <- 1-pi_notgoing # % going before the natexp
    
    # making a fake 2nd stage var
    stage_two[running <0] <- # replace those values before the Nat Exp on stage_two
      sample(0:1, # either false or true
             length(stage_two[running <0]), replace = T, 
             prob = c(pi_notgoing, pi_going)) # with the probabilities we want
    
    y_STD <- as.numeric(scale(y)) # scaling the outcome, Y variable
    
    pwr <-rdpower(data=cbind(y_STD, running), fuzzy = stage_two, tau = std_eff)$power.rbc
    return(pwr)
  } else
    print("Error: either missing an y (outcome), a x (running) or secound stage (stage_two) argument OR they are incompatible dimensions try length()")
}

# I want to get all values of power for all combinations from 0-1 in steps of .05 of pi_notgoing and std_eff
stage2 <- rep(seq(0,1, by = .05), each = 21)
std_eff<- rep(seq(0,1, by = .05), times = 21)

# doing the rdrohust model for Surface Area
rdr_SA <- rdrobust(sa$SA, sa$running_var, fuzzy = sa$EduAge16)
pwr_sa <- data.frame(stage2 = stage2,
                        std_eff = std_eff,
                        pwr = map2_dbl(stage2, std_eff, 
                                              ~secound_stage_STDpower(pi_notgoing = .x, std_eff = .y,
                                                                      y = sa$SA, 
                                                                      running = sa$running_var, 
                                                                      stage_two = sa$EduAge16)))
# fwrite(pwr_sa, "~/projects/EduTelomere/temp_data/pwr_sa.csv")

# there is a significant slow down in the function after adding checks and more args...

# doing the rdrohust model for ltl
rdr_ltl <- rdrobust(telomere_set$ltl, telomere_set$running_var, fuzzy = telomere_set$EduAge16)
pwr_ltl <- data.frame(stage2 = stage2,
                     std_eff = std_eff,
                     pwr = map2_dbl(stage2, std_eff, 
                                    ~secound_stage_STDpower(pi_notgoing = .x, std_eff = .y,
                                                            y = telomere_set$ltl, 
                                                            running = telomere_set$running_var, 
                                                            stage_two = telomere_set$EduAge16)))

# fwrite(pwr_ltl, "~/projects/EduTelomere/temp_data/pwr_ltl.csv")



# graph of neuro power, was gonna do power curves but than Copilot recommended this & I like it!
pwr_plt_neuro <-ggplot(pwrSTD_df, aes(x = stage2, y = std_eff, fill = power_.5sd))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  geom_vline(xintercept = .1, color = "red")+
  geom_vline(xintercept = .16, color = "darkred")+
  labs(title = "RD power graph for SA",
       subtitle = "n = 34010; effobs ~ 10042",
       x = "Proportion of people not going to school at 16",
       y = "Standardized effect size (s.d.)",
       fill = "Power (RBC)",
       caption = "red line: 2nd stage for neuro \n white line: 2nd stage for telomere")

# do the same graph but with ltl as the outcome
pwr_plt_ltl <- ggplot(pwrSTD_df_TELOMERE, aes(x = stage2, y = std_eff, fill = power_.5sd))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  geom_vline(xintercept = .1, color = "red")+
  geom_vline(xintercept = .16, color = "darkred")+
  labs(title = "RD power graph for ltl",
       subtitle = "n = 256243; effobs ~ 61828",
       x = "Proportion of people not going to school at 16",
       y = "Standardized effect size (s.d.)",
       fill = "Power (RBC)",
       caption = "red line: 2nd stage for neuro \n white line: 2nd stage for telomere")

# comparison from neuro to ltl

# mean power across all effect sizes at a .25 stage2
mean(pwrSTD_sa[pwrSTD_sa$stage2==.25,]$pwr)
mean(pwrSTD_df_TELOMERE[pwrSTD_df_TELOMERE$stage2==.25,]$power_.5sd)

# it is the same... possibly because a large n doesn't add more discrete data?


# next I want to explore how discreteness makes sample size stop mattering?



# rdpower without data
rdr_ltl



# option for rdrobust(): if TRUE, x and y are standardized before computing the bandwidths. Default is stdvars=TRUE.


#
rdpower(data=cbind(sa$SA, sa$running_var), fuzzy = sa$EduAge16)













aux <- rdpower(tau=5,
               nsamples=c(rdr_ltl$N[1],rdr_ltl$N_h[1],rdr_ltl$N[2],rdr_ltl$N_h[2]),
               bias=c(rdr_ltl$bias[1], rdr_ltl$bias[2]),
               variance=c(rdr_ltl$V_rb_l[1,1],rdr_ltl$V_rb_r[1,1]),
               sampsi=c(rdr_ltl$N_h[1],rdr_ltl$N_h[2]),
               samph=c(rdr_ltl$bws[1,1],rdr_ltl$bws[1,1]))

# comparing exp-post power across specifications

aux1 <- rdpower(data=Z,tau=5,p=1,h=20,plot=TRUE)
aux2 <- rdpower(data=Z,tau=5,p=2,h=20,plot=TRUE)
aux3 <- rdpower(data=Z,tau=5,p=1,plot=TRUE)
aux4 <- rdpower(data=Z,tau=5,p=2,plot=TRUE)






# default is half the SD (of the untreated group)
rdpower(data=cbind(sa$SA, sa$running_var), fuzzy = sa$EduAge16)
rdpower(data=cbind(ct$CT, ct$running_var), fuzzy = ct$EduAge16)$power.rbc
rdpower(data=cbind(WMh$WM_hyper, WMh$running_var), fuzzy = WMh$EduAge16)$power.rbc
rdpower(data=cbind(CSFn$CSF_norm, CSFn$running_var), fuzzy = CSFn$EduAge16)$power.rbc
rdpower(data=cbind(TBVn$TBV_norm, TBVn$running_var), fuzzy = TBVn$EduAge16)$power.rbc
rdpower(data=cbind(wFAs$wFA, wFAs$running_var), fuzzy = wFAs$EduAge16)$power.rbc




X <- array(rnorm(2000),dim=c(1000,2))
R <- X[,1] + X[,2] + rnorm(1000)

Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(1000)

cor(R,Y)

tmp <- rdpower(data=cbind(Y,R))

length(unique(round(R, 2)))# 497 unique values
tmp <- rdpower(data=cbind(Y,round(R, 2)))
# power stays similar...

length(unique(round(R, 1)))# 91 unique values
tmp <- rdpower(data=cbind(Y,round(R, 1)))

length(unique(floor(R*10)))# 93 unique values

tmp <- rdpower(data=cbind(Y,floor(R*10)))



