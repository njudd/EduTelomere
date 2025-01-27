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

pwr_sa <- fread("~/projects/EduTelomere/temp_data/pwr_sa.csv")
pwr_ltl <- fread("~/projects/EduTelomere/temp_data/pwr_ltl.csv")


# graph of neuro power, was gonna do power curves but than Copilot recommended this & I like it!

# subset if you wish pwr_ltl[std_eff < .5 & stage2 < .5]
pwr_plt_neuro <-ggplot(pwr_sa, aes(x = stage2, y = std_eff, fill = pwr))+
  geom_vline(xintercept = .1, color = "red")+
  geom_vline(xintercept = .16, color = "black")+
  geom_tile()+
  geom_text(aes(label = round(pwr,2)), color = "white", size = 2.5) +
  scale_fill_viridis_c()+
  theme_minimal()+
  labs(title = "RD power graph for SA",
       subtitle = "n = 34010; effobs ~ 10042",
       x = "Proportion of people not going to school at 16",
       y = "Standardized effect size (s.d.)",
       fill = "Power (RBC)",
       caption = "red line: 2nd stage estimate (10%) \n black line: Clark & Royer 2nd stage estimate")

# do the same graph but with ltl as the outcome
pwr_plt_ltl <- ggplot(pwr_ltl, aes(x = stage2, y = std_eff, fill = pwr))+
  geom_vline(xintercept = .16, color = "red")+
  geom_vline(xintercept = .16, color = "black")+
  geom_tile()+
  geom_text(aes(label = round(pwr,2)), color = "white", size = 2.5) +
  scale_fill_viridis_c()+
  theme_minimal()+
  labs(title = "RD power graph for ltl",
       subtitle = "n = 256243; effobs ~ 61828",
       x = "Proportion of people not going to school at 16",
       y = "Standardized effect size (s.d.)",
       fill = "Power (RBC)",
       caption = "red line: 2nd stage estimate (16%) \n black line: Clark & Royer 2nd stage estimate")

# comparison from neuro to ltl

ggsave("~/projects/EduTelomere/temp_plts/pwr_plt_neuro.png", pwr_plt_neuro, 
       bg = "white", width = 9, height =7.3)
ggsave("~/projects/EduTelomere/temp_plts/pwr_plt_ltl.png", pwr_plt_ltl, 
       bg = "white", width = 9, height =7.3)


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





################## -------- ##################
#################-- simulating data --#################
################## -------- ##################

# I want to simulate 7000, 1000 & 240 discrete values

# for telomere the dataset is n=255243
# for neuro n=34010

# example:
# library(arules)
# discretize(iris$Sepal.Length, breaks =3)

# function to donut hole (get rid of the entry over the cutoff)
first_positive <- function(vec){
  for (i in 1:length(vec)){
    if (attr(vec,"discretized:breaks")[i] > 0){return(i)}}}

# Toy dataset
X <- array(rnorm(34000),dim=c(34000,2))
R <- X[,1] + X[,2] + rnorm(34000)

# simulating the relationship contineous n = 34,000
Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(34000)
summary(rdrobust::rdrobust(Y,R))

# now we will discrete it to our largest unit days over 20 years (~7000)
R_7000 <- discretize(R, breaks = 7000, labels = FALSE)
c_7000 <- first_positive(R_7000)
R_7000[R_7000 == c_7000] <- rep(NA, length(R_7000[R_7000 == c_7000]))
rdpower(data=cbind(Y,R_7000), cutoff = c_7000)$power.rbc

# now we will discrete it to weeks over 20 years (~1000)
R_1000 <- discretize(R, breaks = 1000, labels = FALSE)
c_1000 <- first_positive(R_1000)
R_1000[R_1000 == c_1000] <- rep(NA, length(R_1000[R_1000 == c_1000]))
rdpower(data=cbind(Y_1000,R_1000), cutoff = c_1000)$power.rbc

# now we will discrete it to months over 20 years (~240)
R_240 <- discretize(R, breaks = 240, labels = FALSE)
c_240 <- first_positive(R_240)
R_240[R_240 == c_240] <- rep(NA, length(R_240[R_240 == c_240]))
rdpower(data=cbind(Y_240,R_240), cutoff = c_240)$power.rbc

# this increased power which was kinda weird...



#### trying my own simulation...
# if we do a year with 336 days it is easier as each month has 28 days and 4 weeks...
# ^this might artificially inflate the power though; so we will also make it larger

effect_per_day <- 0.001767/30

running_days <- -3360:3359
running_months <- rep(-120:119, each = 28) # there are 28 days in a month
running_years <- rep(-12:19, each = 280) # there are 28 days in a month


# in our data set
# a mean of 1068 objs per month 
# means around 35 per day...

# n_per_day <- rbinom(length(y), 70, p=0.5)

# an equal per day strategy for now...
n_per_day <- 35

# length(running_days)*35 # n ~ 235k

df <- data.frame(running_days = rep(running_days, times = n_per_day), 
                 running_months = rep(running_months, times = n_per_day),
                 running_years = rep(running_years, times = n_per_day))

# sum(df$running_months == -120); sum(df$running_months == 199)
# looks good...

Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(34000)

df$y <- rnorm(length(df$running_days), 0, 1) + (effect_per_day/10)*df$running_days


rdpower(data=cbind(df$y,df$running_days), cutoff = 0)$power.rbc
summary(rdrobust::rdrobust(df$y, df$running_days))


rdpower(data=cbind(df$y,df$running_months), cutoff = 0)$power.rbc
summary(rdrobust::rdrobust(df$y, df$running_months))


RDHonest::RDHonest(y ~ running_months, data = df)
RDHonest::RDHonest(y ~ running_days, data = df)


















# Power against tau = 1
pwr_1000 <- rdpower(data=cbind(Y,R),tau=1)$power.rbc
# Power against tau = 1 including covariates
pwr_1000_cov <- rdpower(data=cbind(Y,R),tau=1,covs=X)$power.rbc


Y_30 <- discretize(Y, breaks = 30)

# finding zero

Y_30[Y_30 == levels(discretize(Y, breaks = 30))[12]] <-
  rep(NA, length(Y_30[Y_30 == levels(discretize(Y, breaks = 30))[12]]))






R_30 <- discretize(R, breaks = 30, labels = FALSE)
R_30[R_30 == first_positive(attr(R_30,"discretized:breaks"))] <- rep(NA, length(R_30[R_30 == first_positive(attr(R_30,"discretized:breaks"))]))
pwr_30 <- rdpower(data=cbind(Y,R_30), cutoff = first_positive(attr(R_30,"discretized:breaks")))$power.rbc


R_300 <- discretize(R, breaks = 300, labels = FALSE)
R_300[R_300 == first_positive(attr(R_300,"discretized:breaks"))] <- rep(NA, length(R_300[R_300 == first_positive(attr(R_300,"discretized:breaks"))]))
pwr_300 <- rdpower(data=cbind(Y,R_300), cutoff = first_positive(attr(R_300,"discretized:breaks")))$power.rbc



pwr_1000; pwr_30; pwr_300


### how you need to loop this a bit & mess around for your context


# length(unique(telomere_set$running_var))*4 # 960 if you have weeks
# length(unique(telomere_set$running_var))*28 # 6720 if you have days


# 1) first show how more data != more power



# is there a difference between simulating discrete and making it discrete?








# 2) then show how more unique data point can really help power






















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



