#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-06-04 
# https://njudd.com 

# Looking at if ROSLA causes a difference in telomere length

# first seeing if there is an assosiational effect
# than I want to be powered to a 1/3rd of that effect size

### Table of Contents
# 1.1 Loading & cleaning
# 1.2 making a running var
# 1.3 plotting
# 1.4 fuzzy RD
# 1.5 assumption tests & misc
#### ### ----------------- ### ####


# UKB: I just read straight off a copy of /phenotypes/current in the Donders server. 
# These are already funpacked, 
# I doublechecked using the command line fmrib-unpack


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.1 Loading & cleaning ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (!require(pacman)){
  install.packages('pacman')
}


pacman::p_load(tidyverse, lubridate, stringr, fastDummies, mice, ggseg, RDHonest, # see notes must be v 0.4.1
               kableExtra, rddensity, patchwork, ggrain, report)

# using an old version of RDHonest. I shouldn't matter in this context, yet is the one I have installed
# from the ROSLA brain paper, which needed covariates
# the RDHonest package added covariates, we will use an old verion from commit aa616f4 where p-vals were added
# this is because we preregisterd, otherwise the new version with automatic covariate correction is prefered
## new RD package with covs SA | EduAge16 ~ running_var | covs  
# packageVersion("RDHonest") # 0.4.1 -----> 0.9
# 
# RDHonest(SA | EduAge16 ~ running_var, data = sa)
# devtools::install_github("kolesarm/RDHonest@aa616f4", auth_token = auth)

# attach to data with command K
witte_vars <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_02_06.csv") # he (Ward.deWitte radboudumc.nl) said he doesn't have 196 & 24419
witte_UKbirth <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_02_23.csv")
witte_telomere <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/N.Judd_2024_06_07.csv")


# rm(list = ls()[ls() != "witte_vars"])
fullset <- data.table::copy(witte_vars)

# renaming the UKB id's: https://biobank.ndph.ox.ac.uk/ukb/search.cgi
# UKBID-instance-array
# instance 2 is the imaging visit!

namevar = c(year = "34-0.0", month = "52-0.0", sex = "31-0.0", visit_date = "53-0.0",
            EduAge_1 = "845-0.0", EduAge_2 = "845-1.0", EduAge_3image = "845-2.0")
cols_remove = c("54-0.0", "54-1.0", "54-3.0", "53-2.0", "53-1.0", "53-3.0", # I don't care about the other site visits
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

# adding telomere IDs

#
# insert telomeres
#



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.2 making a running var ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

fullset <- fullset %>% 
  # mutate(running_var = ym(str_c(year, "-", month))) %>% # right now running_var is DOB in Year-Month (done above)
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

# data.table::fwrite(fullset, "/Volumes/home/lifespan/nicjud/UKB/proc/xx")
# fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/xx")


#adding telomeres
witte_telomere <- witte_telomere[, .(eid, `22191-0.0`, `22192-0.0`)]
tel_names <- c("adj_TS_ratio", "Z-adjusted_TS_log")
colnames(witte_telomere)[2:3] <- tel_names

fullset <- witte_telomere[fullset, on = "eid"] #joining it to fullset

# doing what is recommended in the article
fullset$ltl <- as.numeric(scale(log(fullset$adj_TS_ratio)))
fullset[, (tel_names):=NULL]

# data.table::fwrite(fullset, "~/Desktop/fullset_hold.csv")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.4 Analysis fuzzy local-linear RD   ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# make sure this is done with the same data set called "telomere_set"
# have a standard descriptive table of this dataset
# first just uncorrected
# than do outlier treatment
# than check visit date for percision

# fullset <- data.table::fread("~/Desktop/fullset_hold.csv")

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


# making a descp table of relevant vars
options(knitr.kable.NA = '')

report_table(telomere_set[, .(ltl, ltl_NoOuts, running_var, running_var.s, EduAge, EduAge.s, EduAge16)]) %>% 
  kbl(caption = "Descriptives Leukocyte Telomere Length (ltl)") %>%
  kable_styling("hover", full_width = F) %>%
  save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/descp_table.html")

ltl_age <- lm(ltl ~ running_var, data = telomere_set)
ltl_edu <- lm(ltl ~ EduAge, data = telomere_set)
ltl_age$coefficients[2]*12; ltl_edu$coefficients[2]

# std results
ltl_age.s <- lm(ltl ~ running_var.s, data = telomere_set)
ltl_edu.s <- lm(ltl ~ EduAge.s, data = telomere_set)
ltl_age.s$coefficients[2]; ltl_edu.s$coefficients[2]

# but there is heterogeneity in recruitment
hist(telomere_set$visit_day_correct/365)
# with very little for to the running var
# cor(telomere_set$visit_day_correct, telomere_set$running_var)

ltl_age_cor <- lm(ltl ~ running_var + visit_day_correct + visit_day_correct2, data = telomere_set)
ltl_edu_cor <- lm(ltl ~ EduAge + visit_day_correct + visit_day_correct2, data = telomere_set)
ltl_age_cor$coefficients[2]*12; ltl_edu_cor$coefficients[2]

# double check syntax in the early manual***
m1_pre <- RDHonest(ltl ~ EduAge16 |running_var, data = telomere_set)
m1 <- RDHonest(ltl ~ EduAge16 |running_var, data = telomere_set,
         T0 = m1_pre$coefficients$estimate)


# getting the time between ROSLA & scanning :)
visit <- fullset[!is.na(visit_date),]
visit$DOB <- ym(str_c(visit$year,"-", visit$month))
#age of vist AoV
visit$AoV <- interval(visit$DOB, visit$visit_date) %/% months(1)/12



r <- RDHonest(AoV ~ EduAge16 |running_var, data = visit)

RDHonest(AoV ~ EduAge16 |running_var, data = visit,
         T0 = r$coefficients$estimate)




summary(lm(telomerelogSTD ~ running_var.s + AoV, data = visit))
summary(lm(telomerelogSTD ~ EduAge.s + AoV, data = visit))



# SI_descpAGEplt <- ggplot(scanage, aes(1, AOS)) + 
#   geom_rain(point.args = list(alpha = .08), fill = "#EC7063",
#                     point.args.pos = list(position = position_jitter(width = 0.06, height = 0, seed = 42))) +
#   theme_minimal(base_size = 25) +
#   labs(x = "", y = "Age at neuroimaging") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) +
#   scale_y_continuous(breaks = c(50,55,60,65,70,75,80))
# 
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_FigScanAGE.png", 
#        SI_descpAGEplt, bg = "white")

# looks normal take the mean
mean(scanage$AOS) # 61.89 ~ 62

62-16 # 16 because it is the years AFTER the intervention

# 46 years

# issue of visit days vs Running var & DOB

# cor(fullset$visit_day_correct, fullset$running_var, use = "pairwise.complete.obs")
# cor(fullset$visit_day_correct, fullset$year, use = "pairwise.complete.obs")
# 
# ggplot(fullset, aes(visit_day_correct, running_var)) +
#   geom_point(alpha = .1) +
#   geom_smooth()
# 
# 
# visit_day_correct ~ EduAge16 | running_var
# 
# cor(sa$visit_day_correct, sa$running_var, use = "pairwise.complete.obs")
# 
# sa$ROSLA <- sa$running_var >0 #this is the NatEx grouped
# 
# summary(lm(ROSLA ~ visit_day_correct, data =sa))
# cor(sa$visit_day_correct, sa$running_var, use = "pairwise.complete.obs")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.3 plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#moving plotting below so I can put the MSE derived bounds

# RDHonest::RDScatter(SA~birth_quarters, data = fullimage, avg = Inf, propdotsize = T, vert = T)
# rdrobust::rdplot(fullimage$SA, fullimage$running_var, p = 3)


Edu16_plt <- fullset[!is.na(fullset$visit_date),] %>% #making sure it is imaging subjects
  group_by(running_var) %>% 
  summarise(piEdu16 = sum(EduAge16, na.rm = T)/n()) %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point(color = "blue", alpha = .3) +
      geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Percent staying until 16'), x = "Date of Birth in Months", title = "First Stage ROSLA") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.6, 1)) +
      theme_minimal(base_size = 22) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())
  }


# telomere plt
telomere_plt <- telomere_set %>% 
  group_by(running_var) %>% 
  summarise(sa = mean(telomerelogSTD, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, sa)) +
      geom_point(color = "blue", alpha = .3) +
      geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = "Leukocyte Telomere Length", x = "Date of Birth in Months",  title = "Secound Stage ROSLA") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # scale_y_continuous(breaks=c(164000, 168000, 172000, 176000),
      #                    labels=c(expression("164000" ^mm2 ~ ""), "Sept.\n1952", "Sept.\n1957", "Sept.\n1962"),
      #                    limits = c(164000, 176000))+ 
      ylim(c(-.3, .3)) +
      theme_minimal(base_size = 22) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

SI_plt1 <- Edu16_plt / telomere_plt + 
  plot_annotation(tag_levels = 'a')

ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/SI_plt1.png", 
       SI_plt1, width = 12, height = 10)


telomere_plt_linear <- telomere_set %>% 
  group_by(running_var) %>% 
  summarise(sa = mean(telomerelogSTD, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, sa)) +
      geom_point(color = "blue", alpha = .3) +
      geom_point(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < m1$coefficients$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m1$coefficients$bandwidth | running_var  > m1$coefficients$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var > -m1$coefficients$bandwidth & running_var < 0), method='glm',formula=y~poly(x,1),se=F, color = "red") +
      geom_smooth(data=subset(., running_var > 0 & running_var  < m1$coefficients$bandwidth), method='glm',formula=y~poly(x,1),se=F, color = "red") +
      labs(y = "Leukocyte Telomere Length", x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # scale_y_continuous(breaks=c(164000, 168000, 172000, 176000),
      #                    labels=c(expression("164000" ^mm2 ~ ""), "Sept.\n1952", "Sept.\n1957", "Sept.\n1962"),
      #                    limits = c(164000, 176000))+ 
      ylim(c(-.3, .3)) +
      theme_minimal(base_size = 22) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}


ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/10.02.02 ROSLA Telomere/figs/Plt1.png", 
       telomere_plt_linear, width = 12, height = 8, bg = "white")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.X LOCAL RANDOMIZATION #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

library(rstanarm); library(bayestestR); library(insight)


m1 <- telomere_set[running_var %in% c(-1,0)][
  , ROSLA := running_var >= 0
]
table(m1$ROSLA)
BayesMod_m1 <- stan_glm("telomerelogSTD_out ~ ROSLA", data = m1, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
round(mean(get_parameters(BayesMod_m1)[,2]), 4)
hdi(BayesMod_m1)
bayesfactor_parameters(BayesMod_m1)

m5 <- telomere_set[running_var %in% c(-5,-4, -3, -2, -1, 0, 1, 2, 3, 4)][
  , ROSLA := running_var >= 0
]

table(m5$ROSLA)
BayesMod_m5 <- stan_glm("telomerelogSTD_out ~ ROSLA + visit_day_correct + visit_day_correct2", data = m5, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
round(mean(get_parameters(BayesMod_m5)[,2]), 4)
hdi(BayesMod_m5)
bayesfactor_parameters(BayesMod_m5)

# summary(lm(telomerelogSTD_out  ~ EduAge, data = telomere_set))
rope(BayesMod_m5, range = c(-0.02, 0.02))

######## remeber this is cohen D to SD!
# summary(lm(telomerelogSTD_out  ~ EduAge.s, data = telomere_set))
rope(BayesMod_m5, range = c(-0.05, 0.05))
# so is .05 half...?



# education increases & age decreases...
# so both give more length... eventually age will just show up like it does here!
m12 <- telomere_set[running_var %in% c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3, -2, -1, 0, 1, 2, 3, 4,5,6,7,8,9,10,11)][
  , ROSLA := running_var >= 0
]

table(m12$ROSLA)
BayesMod_m12 <- stan_glm("telomerelogSTD_out ~ ROSLA", data = m12, iter = 40000, refresh=0, prior = normal(location = 0, scale = 1, autoscale = TRUE))
hdi(BayesMod_m12)


# you could totally, plot each year against each other in stem & leaf...


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.5 Assumptions - Checking the RD design #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# for these tests, I obviously want to use the sample I analyzed (the neuroimaging followup; instance 2)
# I select these subjects by filtering by visit date, which is only instance 2
# fullset[!is.na(fullset$visit_date),]

# Density of the running variable, silly test for ROSLA
# rddensity: all statistics & default options
# it is assuming something not great about zero; redoing the running var to fix this
scanage$running_var_SHIFT <- scanage$running_var # shifting positive values by 1; so 0 -> 1, etc etc...
scanage$running_var_SHIFT[scanage$running_var_SHIFT >=0] <- scanage$running_var_SHIFT[scanage$running_var_SHIFT >=0] +1

# scanage is a better representation of the included subjects!
dens_test <- rddensity(X = scanage$running_var_SHIFT) 

# testing the covariates
simple_fuzzy <- function(f, dt){
  p0 <- RDHonest(f, data = dt)
  p1 <- RDHonest(f, data = dt, T0 = p0$coefficients$estimate) #as seen in the vinjette; giving a starting val
  out <- p1$coefficients[c("term", "bandwidth", "eff.obs", "estimate", "conf.low", "conf.high", "p.value")]
  out[2:6] <- round(out[2:6], 2) #rounding; not the p because I will FDR
  out$term = f # adding the forumula to the output
  return(out)
}

# dummy coding cols means you miss one center on the test
fullset$imaging_center_11025 <- fullset$imaging_center == 11025
covsT <- c(covs, "imaging_center_11025")

covsT %>%
  str_c(., " ~ EduAge16 | running_var") %>% #making a forumula from the covs
  map_dfr(~simple_fuzzy(., dt = fullset[!is.na(fullset$visit_date),])) %>% #RDHonest automatically listwise deletes
  mutate(pFDR = round(p.adjust(p.value, method = 'fdr'), 3),
         p.value = round(p.value, 3)) %>%
  kbl(caption = "Placebo Outcomes from the covs of no interest") %>%
  kable_styling("hover", full_width = F) %>%
  save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/SI_Table1_placebo_outcomes2.html")



# looking closer into summer, need to be donut holed removing running var 7,8,9,10 of 1957
# summer is deterministically related to ROSLA for a 2 mnth window around the cutoff
summer0 <- RDHonest(summer ~ EduAge16 | running_var, data = fullset[!is.na(fullset$visit_date),])
summer1 <- RDHonest(summer ~ EduAge16 | running_var, data = fullset[!is.na(fullset$visit_date),], T0 = summer0$coefficients$estimate)

# I think it wants to running var remapped
twoDonuts <- fullset[!fullset$running_var %in% c(-2,-1, 0, 1),]
twoDonuts$running_var[twoDonuts$running_var >=0] <- twoDonuts$running_var[twoDonuts$running_var >=0] -2
twoDonuts$running_var[twoDonuts$running_var <0] <- twoDonuts$running_var[twoDonuts$running_var <0] +2
table(twoDonuts$running_var)

summerD0 <- RDHonest(summer ~ EduAge16 | running_var, data = twoDonuts[!is.na(twoDonuts$visit_date),])
summerD1 <- RDHonest(summer ~ EduAge16 | running_var, data = twoDonuts[!is.na(twoDonuts$visit_date),], T0 = summerD0$coefficients$estimate)


# moving the running var is used commonly; yet we have no sig result to test...


