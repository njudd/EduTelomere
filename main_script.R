#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-06-04 
# https://njudd.com 

# Looking at if ROSLA causes a difference in telomere length

# first seeing if there is an assosiational effect
# than I want to be powered to a 1/3rd of that effect size

### Table of Contents
# 1.0 Packages & Settings
# 1.1 Cleaning & tidying
# 1.2 Loading cleaned/saved data 
# 1.3 Basic descriptives and lm's
# 1.4 LOCAL RANDOMIZATION FREQUENTIST
# 1.5 LOCAL RANDOMIZATION BAYESIAN
# 1.6 Cont RD
#### ### ----------------- ### ####


# UKB: I just read straight off a copy of /phenotypes/current in the Donders server. 
# These are already funpacked, 
# I doublechecked using the command line fmrib-unpack

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.1 Packages & Settings ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (!require(pacman)){install.packages('pacman')}
pacman::p_load(tidyverse, lubridate, stringr, fastDummies, RDHonest, # using the new RDHonest syntax packageVersion("RDHonest")
               kableExtra, rddensity, patchwork, ggrain, report, ggpointdensity, viridis, ggnewscale,
               rstanarm, bayestestR, insight,bayesplot, rdlocrand)
options(scipen=999) # no scientific notation


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.1 Cleaning & tidying ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# attach to data with command K & the donders server
witte_vars <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_02_06.csv") # he (Ward.deWitte radboudumc.nl) said he doesn't have 196 & 24419
witte_UKbirth <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_02_23.csv")
witte_telomere <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_06_07.csv")

# rm(list = ls()[ls() != "witte_vars"])
fullset <- data.table::copy(witte_vars)

# renaming the UKB id's: https://biobank.ndph.ox.ac.uk/ukb/search.cgi
# UKBID-instance-array (instance 0 is the first visit)
namevar = c(year = "34-0.0", month = "52-0.0", sex = "31-0.0", visit_date = "53-0.0", site = "54-0.0",
            EduAge_1 = "845-0.0", EduAge_2 = "845-1.0", EduAge_3image = "845-2.0")
cols_remove = c("54-1.0", "54-3.0", "53-2.0", "53-1.0", "53-3.0", # I don't care about the other site visits
                "25921-3.0", "25922-3.0", "25928-3.0", "25925-3.0", "26500-3.0", #followup confounds
                "190-0.0", # reasons for lost followup
                "25921-2.0", "25922-2.0", "25928-2.0","25925-2.0", "26500-2.0", "54-2.0", #removed NEURO cols Row1
                "25781-2.0", "25003-2.0", "25009-2.0", "26721-2.0", "26755-2.0", "26822-2.0", "26856-2.0", #removed NEURO cols Row2
                "135-0.0", "135-1.0", "135-2.0", "135-3.0", # mistakenly given to me; instead on an imaging category
                "25781-3.0","25003-3.0", "25009-3.0", "26721-3.0", "26755-3.0", "26822-3.0", "26856-3.0") #followup Y's)

fullset[, (cols_remove):=NULL]
fullset <- rename(fullset, all_of(namevar))

# trying new thing 
# now lets only include 10 years Sept 1947 until Aug 1967
fullset <- fullset %>% 
  mutate(running_var = ym(str_c(year, "-", month))) %>% 
  filter(between(running_var, ym('1947-09'), ym('1967-08'))) # as.date('1947-09-01')

# missing GEOGRAPHICAL area
witte_UKbirth <- witte_UKbirth %>% 
  mutate(UKbirth = coalesce(`1647-2.0`, `1647-1.0`, `1647-0.0`)) %>% 
  mutate(UKbirth = UKbirth %in% c(1,2,3)) %>% #England, Wales, Scotland
  select(eid, UKbirth)
fullset <- witte_UKbirth[fullset, on = "eid"] #joining it to fullset
# must be born in England, Wales or Scotland
fullset <- fullset[UKbirth==TRUE]

# making a running var

fullset <- fullset %>% 
  mutate(birth_quarters = interval(ym("1957-9"), running_var) %/% months(3), # Sept 1957 is the DOB of the effected cohort
         running_var1959 = interval(ym("1959-9"),running_var) %/% months(1), # 2 years up
         running_var1955 = interval(ym("1955-9"),running_var) %/% months(1), # 2 years down
         running_var = interval(ym("1957-9"), running_var) %/% months(1), # making the running var actually a running var centered on Sept 1957
         visit_day_correct = interval(ymd(min(fullset$visit_date, na.rm = T)), ymd(fullset$visit_date)) %/% days(1), # the amount of days away from the min subject
         visit_day_correct2 = visit_day_correct^2) # quad effect of # of scan days

# recoding values to missing
## EduAge
# -2	Never went to school
# -1	Do not know
# -3	Prefer not to answer

# recode was complicated.. so I used data.table instead
fullset[, c("EduAge_1", "EduAge_2", "EduAge_3image") := # assigning these columns the same name
          lapply(.SD, function(x) replace(x, which(x %in% c(-3,-2,-1)), NA)), # replace the values -3,-2,-1 with NA
        .SDcols = c("EduAge_1", "EduAge_2", "EduAge_3image")] # for these columns

fullset <- fullset %>% 
  mutate(EduAge = coalesce(EduAge_3image, EduAge_2, EduAge_1)) %>% # find the first non-missing column in imaging followup, followup1 or T1
  select(-c(EduAge_3image, EduAge_2, EduAge_1)) # remove these cols

# Qualifications (6138-instance-array) are coded in the wrong direction
# https://github.com/margotvandeweijer/EA_causality/blob/main/1.%20Clean%20Data/2_clean_education.R
# instance is the followup 3 is equal to neuroimaging wave 1; array allows them to answer multiple
# 1	College or University degree
# 2	A levels/AS levels or equivalent
# 3	O levels/GCSEs or equivalent
# 4	CSEs or equivalent
# 5	NVQ or HND or HNC or equivalent
# 6	Other professional qualifications eg: nursing, teaching
# -7	None of the above
# -3	Prefer not to answer

# I don't care about their Quals, YET those that went to college weren't asked how many years they went to Education for!!!
fullset <- fullset %>% 
  mutate(EduAge = ifelse(if_any(matches("6138"), ~.x %in% c(1)), # if anyone ever reports college (value =1 ); one is equal to college
                                    21, EduAge)) %>% # than make it 21 years in EduAge (the age left education)
  select(-matches("6138")) # remove all 6138 cols

# last thing is making a dummy if someone has 16 years of Education or more because it is fuzzy
# local linear RD
# stage 1: EduAge16 ~ running_var + ROSLA_treat + running_var:ROSLA_treat
# stage2: Y ~ running_var + EduAge16_fitted + running_var:EduAge16_fitted
fullset$EduAge16 <- fullset$EduAge >= 16
# table(fullset$EduAge)

# making sites dummy coded
# fullset <- fastDummies::dummy_cols(fullset, select_columns = "imaging_center", remove_first_dummy = TRUE)

# a summer dummy, since these people could technically end at 15 (but had the same amount of school)
fullset$summer <- as.numeric(fullset$month %in% c(7,8)) # (July + Aug)

#adding telomeres
witte_telomere <- witte_telomere[, .(eid, `22191-0.0`, `22192-0.0`)]
tel_names <- c("adj_TS_ratio", "Z-adjusted_TS_log")
colnames(witte_telomere)[2:3] <- tel_names

fullset <- witte_telomere[fullset, on = "eid"] #joining it to fullset

# doing what is recommended in the article
fullset$ltl <- as.numeric(scale(log(fullset$adj_TS_ratio)))
fullset[, (tel_names):=NULL]

# data.table::fwrite(fullset, "/Volumes/home/lifespan/nicjud/UKB/proc/telomere_fullset.csv")
# fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/telomere_fullset.csv")

vec_to_fence <- function(vec){
  stats <- boxplot.stats(vec)$stats
  vec[vec < stats[1]] <- stats[1]
  vec[vec > stats[5]] <- stats[5]
  return(vec)}

telomere_set <- data.table::copy(fullset)

# library(ggrain); ggplot(fullset, aes(1, telomerelogSTD)) + geom_rain()
# some serious outliers!
# you can do one with & one without any outlier treatment...
telomere_set <- telomere_set[complete.cases(telomere_set[, .(ltl, running_var, EduAge16)]),][ # no missingness allows in Y, running, or IV
  ,c("ltl_NoOuts") := lapply(.SD, vec_to_fence), .SDcols=c('ltl')] # fencing covs with error & Y

telomere_set$running_var.s <- as.numeric(scale(telomere_set$running_var))
telomere_set$EduAge.s <- as.numeric(scale(telomere_set$EduAge))

# doing a 3SD outlier test
telomere_set$ltl_3SD <- telomere_set$ltl
telomere_set$ltl_3SD[telomere_set$ltl_3SD>3 | telomere_set$ltl_3SD<(-3)] <- rep(NA, length(telomere_set$ltl_3SD[telomere_set$ltl_3SD>3 | telomere_set$ltl_3SD<(-3)]))
sum(is.na(telomere_set$ltl_3SD))/length(telomere_set$ltl_3SD)

# fixing sex

telomere_set$sex <- as.character(telomere_set$sex)
telomere_set$sex[telomere_set$sex ==0] <- rep("Female", length(telomere_set$sex[telomere_set$sex ==0]))
telomere_set$sex[telomere_set$sex ==1] <- rep("Male", length(telomere_set$sex[telomere_set$sex ==1]))

# age at visit in months
telomere_set$age <- interval(ym(str_c(telomere_set$year,"-", telomere_set$month)), ymd(telomere_set$visit_date)) %/% months(1)


# data.table::fwrite(telomere_set, "/Volumes/home/lifespan/nicjud/UKB/proc/telomere_set.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.2 Loading cleaned/saved data   ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# telomere_set <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/telomere_set.csv")
# pacman::p_load(tidyverse, lubridate, stringr, fastDummies, RDHonest, # using the new RDHonest syntax packageVersion("RDHonest")
#                kableExtra, rddensity, patchwork, ggrain, report, ggpointdensity, viridis, ggnewscale,
#                rstanarm, bayestestR, insight,bayesplot, rdlocrand)
# options(scipen=999) # no scientific notation

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.3 Basic descriptives and lm's   ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# make sure this is done with the same data set called "telomere_set"
# have a standard descriptive table of this dataset
# first just uncorrected
# than do outlier treatment
# than check visit date for percision

# cor(telomere_set$running_var, telomere_set$age)
# there is a .99 correlation between age in months & our running var

# for the structured abstract & table
psych::describe(telomere_set$age/12)

# basic SI plot of telomere lenght
# basic_SIplt <- ggplot(telomere_set[ltl > -5 & ltl < 5], aes(running_var, ltl)) + 
#   geom_point(alpha= .08, color = "darkblue") +
#   geom_smooth(aes(color = sex), method = "lm", se = FALSE) +
#   theme_minimal(base_size = 25) +
#   labs(y = "Telomere Length (ltl)", x = "Running Variable (age)", color = "") +
#   scale_color_brewer(palette = 'Dark2')
# 
# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_rawdatafig_MvF.png",
#        basic_SIplt, width = 10, height = 6, bg = "white")

# axises are limited to 5SD's this excluded the illustration of XXX subjects (XX%).
dim(telomere_set)[1] - dim(telomere_set[ltl > -5 & ltl < 5])[1]


# trying geom_point density again... 
basic_SIplt_dens <- ggplot(telomere_set[ltl > -5 & ltl < 5], aes(age/12, ltl)) + 
  geom_pointdensity(alpha= .2, adjust = .1) +
  scale_color_viridis("Density", option = "C") +
  new_scale_color() +
  geom_smooth(method = "lm", color = "black", size = 1.5) +
  # scale_color_manual(values = c("#28B463", "#80DEEA")) +
  theme_minimal(base_size = 30) +
  labs(y = "Telomere Length\n(Std.)", x = "Age", color = "Sex", fill = "Density") + 
  theme(legend.position="none") +
  scale_x_continuous(breaks=c(42, 47, 52, 57, 62),
                     labels=c("Age\n42", "Age\n47", "Age\n52", "Age\n57", "Age\n62")) +
  theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())

# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_rawdatafig_MvF.png",
#        basic_SIplt_dens, width = 10, height = 6, bg = "white")

# making a descp table of relevant vars
options(knitr.kable.NA = '')
# report_table(telomere_set[, .(ltl, ltl_3SD, running_var, running_var.s, visit_day_correct, visit_day_correct2, EduAge, EduAge.s, EduAge16, sex)]) %>% 
#   mutate_if(is.numeric, round, 2) %>% 
#   kbl(caption = "Descriptives Leukocyte Telomere Length (ltl)") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/descp_table.html")

ltl_age <- lm(ltl ~ running_var, data = telomere_set)
ltl_edu <- lm(ltl ~ EduAge, data = telomere_set)
ltl_age$coefficients[2]*12; ltl_edu$coefficients[2]

confint(ltl_age); confint(ltl_edu)

# std results
ltl_age.s <- lm(ltl ~ running_var.s, data = telomere_set)
ltl_edu.s <- lm(ltl ~ EduAge.s, data = telomere_set)
ltl_age.s$coefficients[2]; ltl_edu.s$coefficients[2]


# checking running var vs age
summary(lm(ltl ~ running_var, data = telomere_set))
summary(lm(ltl ~ age, data = telomere_set))


telomere_set$age.s <- as.numeric(scale(telomere_set$age))
summary(lm(ltl ~ running_var.s, data = telomere_set))
summary(lm(ltl ~ age.s, data = telomere_set))


# but there is heterogeneity in recruitment
hist(telomere_set$visit_day_correct/365)
# with very little for to the running var
# cor(telomere_set$visit_day_correct, telomere_set$running_var)

ltl_age_cor <- lm(ltl ~ running_var + visit_day_correct + visit_day_correct2, data = telomere_set)
ltl_edu_cor <- lm(ltl ~ EduAge + visit_day_correct + visit_day_correct2, data = telomere_set)
ltl_age_cor$coefficients[2]*12; ltl_edu_cor$coefficients[2]



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.4 LOCAL RANDOMIZATION FREQUENTIST #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


telomere_set$EduAge16 <- as.numeric(telomere_set$EduAge16)

set.seed(42)


# reporting ITT & TSLS; power for both at 80% and on an effect of .2SD



lr_w1_ITT <- rdrandinf(Y = telomere_set$ltl, R = telomere_set$running_var, 
                   fuzzy = telomere_set$EduAge16,
                   wmasspoints = TRUE, wl=c(-1), wr=c(0), d = .13) #need to manually figure out 80%
lr_w1_TSLS <- rdrandinf(Y = telomere_set$ltl, R = telomere_set$running_var, 
                   fuzzy = c(telomere_set$EduAge16, "tsls"),
                   wmasspoints = TRUE, wl=c(-1), wr=c(0), d = .57)

lr_w5_ITT <- rdrandinf(Y = telomere_set$ltl, R = telomere_set$running_var, 
                   fuzzy = c(telomere_set$EduAge16), #statistic = "tsls",
                   wmasspoints = TRUE, wl=-5, wr=4, d = .06)
lr_w5_TSLS <- rdrandinf(Y = telomere_set$ltl, R = telomere_set$running_var, 
                   fuzzy = c(telomere_set$EduAge16, "tsls"), #statistic = "tsls",
                   wmasspoints = TRUE, wl=-5, wr=4, d = .35)

# as the reviewer said, subseting increases power... as then there is less college people water down my eff
# telomere_setNOcollege <- telomere_set[telomere_set$EduAge<18,]
# rdrandinf(Y = telomere_setNOcollege$ltl, R = telomere_setNOcollege$running_var, 
#           fuzzy = c(telomere_setNOcollege$EduAge16, "tsls"), #statistic = "tsls",
#           wmasspoints = TRUE, wl=-5, wr=4, d = .2)

# it is too simple see https://pmc.ncbi.nlm.nih.gov/articles/PMC5837396/ (altho 2 binary things)
# power.t.test(n = 400, delta = 0.2, sd = sd(telomere_set$ltl),
#              type = "two.sample")
# https://venexia.shinyapps.io/PharmIV/



Rain_1m <- ggplot(telomere_set[running_var %in% c(-1,0)], aes(as.character(running_var),ltl_3SD, fill = as.character(running_var))) +
  geom_rain(rain.side = "f1x1", alpha = .8) +
  scale_fill_brewer(palette = 'Dark2') +
  ggsignif::geom_signif(
    comparisons = list(c("-1", "0")),
    map_signif_level = TRUE) +
  theme_minimal(base_size = 35) +
  theme(legend.position = "none") +
  labs(y = "Telomere Length", x = "") + 
  scale_x_discrete(labels = c("Aug 1957 (- ROSLA)", "Sept 1957 (+ ROSLA)"))

# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_Rain_1m.png",
#        Rain_1m, width = 14, height = 10, bg = "white")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.5 LOCAL RANDOMIZATION BAYESIAN #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# subseting the data
m1 <- telomere_set[running_var %in% c(-1,0)][
  , ROSLA := running_var >= 0
]
# table(m1$ROSLA)


BayesMod_m1 <- stan_glm("ltl ~ ROSLA", data = m1, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
round(mean(get_parameters(BayesMod_m1)[,2]), 4)
hdi(BayesMod_m1)
bayesfactor_parameters(BayesMod_m1)
BF_BayesMod_m1<- bayesfactor_parameters(BayesMod_m1)


quantile(get_parameters(BayesMod_m1)[,2], probs=c(.025,.975))
hdi(BayesMod_m1)



BayesMod_m1_p.5 <- stan_glm("ltl ~ ROSLA", data = m1, iter = 40000, refresh=0, prior = normal(location = 0, scale = .5, autoscale = TRUE))
BayesMod_m1_p1.5 <- stan_glm("ltl ~ ROSLA", data = m1, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1.5, autoscale = TRUE))

ROSLA_posterior_BayesMod_m1_p1 <- get_parameters(BayesMod_m1)[,2]
ROSLA_posterior_BayesMod_m1_p.5 <- get_parameters(BayesMod_m1_p.5)[,2]
ROSLA_posterior_BayesMod_m1_p1.5 <- get_parameters(BayesMod_m1_p1.5)[,2]






round(sd(get_parameters(BayesMod_m1)[,2]), 4)


summary(BayesMod_m1_p.5); summary(BayesMod_m1); summary(BayesMod_m1_p1.5)


rope(BayesMod_m1, range = c(-0.10, 0.10))
rope(BayesMod_m1, range = c(-0.05, 0.05))



1/0.023

m5 <- telomere_set[running_var %in% c(-5,-4, -3, -2, -1, 0, 1, 2, 3, 4)][
  , ROSLA := running_var >= 0
]

# + visit_day_correct + visit_day_correct2

table(m5$ROSLA)
BayesMod_m5 <- stan_glm("ltl ~ ROSLA", data = m5, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
round(mean(get_parameters(BayesMod_m5)[,2]), 4)
hdi(BayesMod_m5)
BF_BayesMod_m5<- bayesfactor_parameters(BayesMod_m5)
1/0.012

BayesMod_m5_p.5 <- stan_glm("ltl ~ ROSLA", data = m5, iter = 40000, refresh=0, prior = normal(location = 0, scale = .5, autoscale = TRUE))
BayesMod_m5_p1.5 <- stan_glm("ltl ~ ROSLA", data = m5, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1.5, autoscale = TRUE))



ROSLA_posterior_BayesMod_m5_p1 <- get_parameters(BayesMod_m5)[,2]
ROSLA_posterior_BayesMod_m5_p.5 <- get_parameters(BayesMod_m5_p.5)[,2]
ROSLA_posterior_BayesMod_m5_p1.5 <- get_parameters(BayesMod_m5_p1.5)[,2]



names <- c(rep('m1_p.5', 80000), rep('m1_p1', 80000), rep('m1_p1.5', 80000), 
  rep('m5_p.5', 80000),rep('m5_p1', 80000), rep('m5_p1.5', 80000))

values <- c(ROSLA_posterior_BayesMod_m1_p.5, ROSLA_posterior_BayesMod_m1_p1, ROSLA_posterior_BayesMod_m1_p1.5,
            ROSLA_posterior_BayesMod_m5_p.5, ROSLA_posterior_BayesMod_m5_p1, ROSLA_posterior_BayesMod_m5_p1.5)

data.frame(names = names, values= values) %>% 
  group_by(names) %>% 
  summarise(mean = mean(values), sd = sd(values), lower = quantile(values, .025), upper = quantile(values, .975))


BF_BayesMod_m5_p.5 <- bayesfactor_parameters(BayesMod_m5_p.5)
BF_BayesMod_m5_p1.5 <- bayesfactor_parameters(BayesMod_m5_p1.5)

BF_BayesMod_m1_p.5 <- bayesfactor_parameters(BayesMod_m1_p.5)
BF_BayesMod_m1_p1.5 <- bayesfactor_parameters(BayesMod_m1_p1.5)

rope(BayesMod_m1_p.5, range = c(-0.1, 0.1))
rope(BayesMod_m1, range = c(-0.1, 0.1))
rope(BayesMod_m1_p1.5, range = c(-0.1, 0.1))
rope(BayesMod_m5_p.5, range = c(-0.1, 0.1))
rope(BayesMod_m5, range = c(-0.1, 0.1))
rope(BayesMod_m5_p1.5, range = c(-0.1, 0.1))

rope(BayesMod_m1_p.5, range = c(-0.05, 0.05))
rope(BayesMod_m1, range = c(-0.05, 0.05))
rope(BayesMod_m1_p1.5, range = c(-0.05, 0.05))
rope(BayesMod_m5_p.5, range = c(-0.05, 0.05))
rope(BayesMod_m5, range = c(-0.05, 0.05))
rope(BayesMod_m5_p1.5, range = c(-0.05, 0.05))



# summary(lm(telomerelogSTD_out  ~ EduAge, data = telomere_set))
rope(BayesMod_m5, range = c(-0.02, 0.02))

######## remeber this is cohen D to SD!
# summary(lm(telomerelogSTD_out  ~ EduAge.s, data = telomere_set))
rope(BayesMod_m5, range = c(-0.05, 0.05))
# so is .05 half...?

# education increases & age decreases...
# so both give more length... eventually age will just show up like it does here!
# m12 <- telomere_set[running_var %in% c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3, -2, -1, 0, 1, 2, 3, 4,5,6,7,8,9,10,11)][
#   , ROSLA := running_var >= 0
# ]
# table(m12$ROSLA)
# BayesMod_m12 <- stan_glm("telomerelogSTD_out ~ ROSLA", data = m12, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
# hdi(BayesMod_m12)


# Bayesian model diags (see brain script for parallel action)

# getting trace plots
color_scheme_set("viridisA")
m1_trace <- mcmc_trace(BayesMod_m1, pars = "ROSLATRUE") + 
  labs(title = "Diagnostic Trace plot BayesLocal 1 month") + theme_minimal(base_size = 20)
m5_trace <- mcmc_trace(BayesMod_m5, pars = "ROSLATRUE") + 
  labs(title = "Diagnostic Trace plot BayesLocal 5 month") + theme_minimal(base_size = 20)

tracePlt <- m1_trace / m5_trace + 
  plot_annotation(tag_levels = 'a')
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_diag_trace.png",
#        tracePlt, width = 14, height = 10)

# plotting and saving posterior draws
color_scheme_set("pink"); 
m1_ppc <- ppc_dens_overlay(y = m1$ltl, yrep = posterior_predict(BayesMod_m1, draws = 300)) + 
  labs(title = "Posterior draws BayesLocal (1 mnth)") + 
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom") +
  xlim(-4,4)
  
m5_ppc <- ppc_dens_overlay(y = m5$ltl, yrep = posterior_predict(BayesMod_m5, draws = 300)) + 
  labs(title = "Posterior draws BayesLocal (5 mnth)") + 
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom") +
  xlim(-4,4)

tracePlt <- m1_ppc + m5_ppc + 
  plot_annotation(tag_levels = 'a')
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_diag_pcc.png",
#        tracePlt, width = 10, height = 5)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.6 Cont RD ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# double check syntax in the early manual
m1_pre <- RDHonest(ltl | EduAge16  ~ running_var, data = telomere_set)
m1 <- RDHonest(ltl | EduAge16 ~ running_var, data = telomere_set,
               T0 = m1_pre$coefficients$estimate)
m1
# outlier treated (very little of the data)
m2_pre <- RDHonest(ltl_3SD | EduAge16  ~ running_var, data = telomere_set)
m2 <- RDHonest(ltl_3SD | EduAge16  ~ running_var, data = telomere_set, 
               T0 = m2_pre$coefficients$estimate)

# check there is no association with visit date
m0_1_pre <- RDHonest(visit_day_correct | EduAge16 ~ running_var, data = telomere_set)
m0_1 <- RDHonest(visit_day_correct | EduAge16 ~ running_var, data = telomere_set, T0 = m0_1_pre$coefficients$estimate)
m0_2_pre <- RDHonest(visit_day_correct2 | EduAge16 ~ running_var, data = telomere_set)
m0_2 <- RDHonest(visit_day_correct2 | EduAge16 ~ running_var, data = telomere_set, T0 = m0_2_pre$coefficients$estimate)
m0_1;m0_2

# with covariates
m3_pre <- RDHonest(ltl | EduAge16  ~ running_var | visit_day_correct + visit_day_correct2, 
                   data = telomere_set)
m3 <- RDHonest(ltl | EduAge16  ~ running_var | visit_day_correct + visit_day_correct2, 
               data = telomere_set,
               T0 = m3_pre$coefficients$estimate)

# all of them clustered
m1_clust <- RDHonest(ltl | EduAge16 ~ running_var, data = telomere_set,
                     T0 = m1_pre$coefficients$estimate, clusterid = site, se.method="EHW")
m2_clust <- RDHonest(ltl_3SD | EduAge16  ~ running_var, data = telomere_set, 
                     T0 = m2_pre$coefficients$estimate, clusterid = site, se.method="EHW")
m3_clust <- RDHonest(ltl | EduAge16  ~ running_var | visit_day_correct + visit_day_correct2, 
                     data = telomere_set,
                     T0 = m3_pre$coefficients$estimate, clusterid = site, se.method="EHW")

# Cont. RD plotting 
telomere_plt <- telomere_set %>% 
  group_by(running_var) %>% 
  summarise(sa = mean(ltl, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, sa)) +
      geom_vline(xintercept = 0,linetype="dashed", size = 1) +
      geom_point(color = "darkblue", alpha = .8, size = 3) +
      # geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "darkblue", alpha = .8, size = 3) +
      # geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "blue", alpha = .3, size = 3) +
      # geom_smooth(data=subset(., running_var > -m1$coefficients$bandwidth & running_var  < m1$coefficients$bandwidth), method='glm',formula=y~poly(x,1),se=F, color = "black",  size = 1.5) +
      # geom_smooth(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < 0), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      # geom_smooth(data=subset(., running_var > 0 & running_var  < m1$coefficients$bandwidth), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      labs(y = "Telomere Length\n(Std. monthly)", x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # scale_y_continuous(breaks=c(164000, 168000, 172000, 176000),
      #                    labels=c(expression("164000" ^mm2 ~ ""), "Sept.\n1952", "Sept.\n1957", "Sept.\n1962"),
      #                    limits = c(164000, 176000))+ 
      ylim(c(-.3, .3)) +
      theme_minimal(base_size = 30) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/Plt1.png",
#        telomere_plt_linear, width = 14, height = 8, bg = "white")

pnas_plt <- basic_SIplt_dens / telomere_plt + 
  plot_annotation(tag_levels = 'a')
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/PNAS_Fig1.png",
#        pnas_plt, width = 14, height = 16)



stage1Plt_inter <- telomere_set %>% 
  group_by(running_var) %>% 
  summarise(piEdu16 = mean(EduAge16, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point(color = "#34495E", alpha = .3, size = 3) +
      geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "black", alpha = .8, size = 3) +
      geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "#34495E", alpha = .3, size = 3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < 0), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      geom_smooth(data=subset(., running_var > 0 & running_var  < m1$coefficients$bandwidth), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      labs(y = bquote('Completed 16 yrs\nof Education (%)'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.6, 1)) +
      theme_minimal(base_size = 30) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

stage2Plt_inter <- telomere_set %>% 
  group_by(running_var) %>% 
  summarise(mltl = mean(ltl, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, mltl)) +
      geom_point(color = "#34495E", alpha = .3, size = 3) +
      geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "black", alpha = .8, size = 3) +
      geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "#34495E", alpha = .3, size = 3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < 0), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      geom_smooth(data=subset(., running_var > 0 & running_var  < m1$coefficients$bandwidth), method='glm',formula=y~poly(x,1),se=F, color = "blue",  size = 1.5) +
      labs(y = "Leukocyte Telomere\nLength (LTL)", x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(-.3, .3)) +
      theme_minimal(base_size = 30) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

StagedPlt <- stage1Plt_inter / stage2Plt_inter + 
  plot_annotation(tag_levels = 'a')

# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_plt1.png",
#        StagedPlt, width = 14, height = 16)

