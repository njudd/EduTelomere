# doing some power calculations
# Nicholas Judd 19/08/24
# njudd.com

if (!require(pacman)){install.packages('pacman')}
pacman::p_load(tidyverse, lubridate, stringr, 
               RDHonest, rdrobust, rdpower)


telomere_set <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/telomere_set.csv")

# linear modele of age & education
ltl_age <- lm(ltl ~ running_var, data = telomere_set)
ltl_edu <- lm(ltl ~ EduAge, data = telomere_set)

# standardized results
ltl_age.s <- lm(ltl ~ running_var.s, data = telomere_set)
ltl_edu.s <- lm(ltl ~ EduAge.s, data = telomere_set)

# could you see the effect of a day with the visit day error?


telomere_set$visit_day_correct.s <- as.numeric(scale(telomere_set$visit_day_correct))

summary(lm(ltl ~ running_var.s + visit_day_correct.s, data = telomere_set))

summary(lm(ltl ~ running_var + visit_day_correct, data = telomere_set))


# there is almost no correlation between visit date & running var
cor(telomere_set$visit_day_correct, telomere_set$running_var)

summary(lm(ltl ~ visit_day_correct, data = telomere_set))
summary(lm(ltl ~ visit_day_correct.s, data = telomere_set))


# visit days might just be a measure of how long they had the study for 
# not the error between the two... but you could make it into months 
# and than see how it helps if it is days...?

#######

m1_pre <- RDHonest(ltl | EduAge16  ~ running_var, data = telomere_set)
m1 <- RDHonest(ltl | EduAge16 ~ running_var, data = telomere_set,
               T0 = m1_pre$coefficients$estimate)
m1


####### try an RD with the college people excluded?
# m1_pre <- RDHonest(ltl | EduAge16  ~ running_var, data = telomere_set[EduAge %in% 14:20])
# m1 <- RDHonest(ltl | EduAge16 ~ running_var, data = telomere_set[EduAge %in% 14:20],
#                T0 = m1_pre$coefficients$estimate)
# m1

# I don't understand why that would work...

telomere_set <- as.data.frame(telomere_set)
Z <- telomere_set[c('ltl','running_var')]
# data table is wierd where you have to order it the other way around

rdpower(Z, fuzzy = telomere_set$EduAge16, tau = .4)
# the fuzzy is really messing with power

rdsampsi(data=Z)
rdsampsi(data=Z, tau = .2)

rdpower(Z, tau = .2)


ggplot(telomere_set, aes(x = running_var, y = ltl)) + 
  geom_point(alpha = .05) + geom_smooth(method = "lm")

ggplot(telomere_set, aes(x = visit_day_correct, y = ltl)) + 
  geom_point(alpha = .05) + geom_smooth(method = "lm")



# trying for UKB neuro

neuro_set <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_fullset.csv")
neuro_set <- as.data.frame(neuro_set)

neuro_set$CT.s <- as.numeric(scale(neuro_set$CT))
neuro_set$SA.s <- as.numeric(scale(neuro_set$SA))

# CT
rdpower(neuro_set[c('CT.s','running_var')], fuzzy = neuro_set$EduAge16, tau = .5)

rdpower(neuro_set[c('CT','running_var')], fuzzy = neuro_set$EduAge16, tau = .5)


# SA
rdpower(neuro_set[c('SA.s','running_var')], fuzzy = neuro_set$EduAge16, tau = 1)

var(neuro_set$CT.s, na.rm = T)
var(neuro_set$SA.s, na.rm = T)


# you need to understand why it changes 
# so dramatically with different neuro vars





# Z2 <- telomere_set[c('ltl','visit_day_correct')]
# rdpower(Z2, fuzzy = telomere_set$EduAge16, tau = .2, cutoff = 1041)


########## simulating data
# 0.001767 per month



# if we do a year with 336 days it is easier as each month has 28 days and 4 weeks...
# ^this might artificially inflate the power though; so we will also make it larger

effect_per_day <- 0.001767/30

running_var <- -3360:3359
running_months <- rep(-120:119, each = 28) # there are 28 days in a month

Z %>% group_by(running_var) %>% 
  summarize(n = n()) %>% 
  summarise(mean = mean(n), min = min(n), max = max(n))

# a mean of 1068 per month 
# means around 35 per day...

# sample(1:70, 30)
# this is an equal sampling strategy which would be a bit wierd...

# making a discrete sampling stragety
# this is something you should try out with more variance

n_per_day <- rbinom(length(y), 70, p=0.5)

# now we will get to a similar number of observations per day


df <- data.frame(running_var = rep(running_var, times = n_per_day), 
                 running_months = rep(running_months, times = n_per_day))

dim(df); dim(telomere_set)


# making it fuzzy, by having 10% before not going to school until 16
df$fuzzy <- sample(c(0,1), length(df$running_var), replace = T, prob = c(.1,.9))
# everything about the cutoff needs to be equal to 1
df$fuzzy <- ifelse(df$running_var < 0, 1, df$fuzzy)


# now do the effect
df$y <- rnorm(length(df$running_var)) + effect_per_day*df$running_var


#### trying the outcome again
lm(y ~ running_var, data = df)
0.00005945*336 # gets close to recovering the parameter...

lm(y ~ running_months, data = df)
0.0016433*12 # also gets close

# lets do a power test with this beeeach

sim_Z <- df[c('y','running_var')]

# data table is wierd where you have to order it the other way around
rdpower(sim_Z, fuzzy = df$fuzzy, tau = .4)

sim_Z_months <- df[c('y','running_months')]
rdpower(sim_Z_months, fuzzy = df$fuzzy, tau = .4)



# now we will group the running var...

-3360:3359

# there are 336 days in a year and each month has 28 days, with 4 weeks
3360/28 # 120 months = 10 years




df$running_var_group <- cut(df$running_var, breaks = seq(-3360, 3360, by = 28))




######### playspace with manual power

rdr <- rdrobust(telomere_set$ltl, telomere_set$running_var, fuzzy = telomere_set$EduAge16)

samph <- rdr$bws[1]

sampsi.l <- rdr$N_h[1]
sampsi.r <- rdr$N_h[2]

bias.l <- rdr$bias[1]/(rdr$bws[1]^2)
bias.r <- rdr$bias[2]/(rdr$bws[2]^2)

VL <- rdr$V_rb_l
VR <- rdr$V_rb_r

N <- sum(rdr$N)

Vl.rb <- N*rdr$bws[1]*VL[1,1]
Vr.rb <- N*rdr$bws[1]*VR[1,1]

aux <- rdpower(data=Z,tau=.2,bias=c(bias.l,bias.r),variance=c(Vl.rb,Vr.rb),samph=samph,sampsi=c(sampsi.l,sampsi.r))


# increasing the variance by 20%
aux <- rdpower(data=Z,tau=.2,variance=c(Vl.rb*1.2,Vr.rb*1.2))

# you want less variance obviously...

rdpower(data=Z,tau=.2,variance=c(Vl.rb*.5,Vr.rb*.5))


# https://stats.stackexchange.com/questions/93303/variance-covariance-matrix-interpretation




# you could make UKB even more discrete to see how that effects power negatively
# for instance by pushing months together (might want to do this in a way were it is equal)
# e.g., Sep gets Aug & Oct, Dec gets Nov & Jan, etc.





###### standard RCT power with only 10% treated impacted by .2 SD


rct_power <- function(n = 100, eff = .2, effectly_treated = .1){
  
  
  y_if_control <- rnorm(n)
  
  ten_index <- sample(c(0,1), n, replace = T, prob = c(1-effectly_treated, effectly_treated))
  
  y_if_treated <- y_if_control
  y_if_treated[as.logical(ten_index)] <- y_if_treated[as.logical(ten_index)] + eff
  
  
  z <- sample(rep(c(0,1), n/2))
  y <- ifelse(z==1, y_if_treated, y_if_control)
  fake <- data.frame(y, z)
  
  summary(lm(y ~ z, data = fake))$coefficients[2,4]
  

}

sum(map_dbl(1:1000, ~rct_power(n = 1000, eff = .4, effectly_treated = .1)) < .05)/1000




n = 50000
effectly_treated = .1
eff = .2
y_if_control <- rnorm(n)

ten_index <- sample(c(0,1), n, replace = T, prob = c(1-effectly_treated, effectly_treated))

y_if_treated <- y_if_control
y_if_treated[as.logical(ten_index)] <- y_if_treated[as.logical(ten_index)] + eff


z <- sample(rep(c(0,1), n/2))
y <- ifelse(z==1, y_if_treated, y_if_control)
fake <- data.frame(y, z)


ggplot(fake, aes(group = z, x = y)) + geom_density()




