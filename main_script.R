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

# grab functions
source("~/projects/roslaUKB/main_analysis/0_functions.R") 
# function: vec_to_fence() takes a vector and puts outliers to the fence (i.e., boxplot.stats)
# function: RDjacked() for cov corrected fuzzy RD

pacman::p_load(tidyverse, lubridate, stringr, fastDummies, mice, ggseg, kableExtra, rddensity, patchwork, ggrain)

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
witte_UKbirth <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/XXX")


# rm(list = ls()[ls() != "witte_vars"])
fullset <- data.table::copy(witte_vars)

# renaming the UKB id's: https://biobank.ndph.ox.ac.uk/ukb/search.cgi
# UKBID-instance-array
# instance 2 is the imaging visit!

namevar = c(year = "34-0.0", month = "52-0.0", sex = "31-0.0", visit_date = "53-2.0",
            EduAge_1 = "845-0.0", EduAge_2 = "845-1.0", EduAge_3image = "845-2.0")
cols_remove = c("54-0.0", "54-1.0", "54-3.0", "53-0.0", "53-1.0", "53-3.0", # I don't care about site visits that aren't imaging
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
fullset <- witte_UK_headmotion[fullset, on = "eid"] #joining it to fullset
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
table(fullset$EduAge)

# making sites dummy coded
fullset <- fastDummies::dummy_cols(fullset, select_columns = "imaging_center", remove_first_dummy = TRUE)

# a summer dummy, since these people could technically end at 15 (but had the same amount of school)
fullset$summer <- as.numeric(fullset$month %in% c(7,8)) # (July + Aug)

# making global vars
fullset$SA <- fullset$lSA + fullset$rSA
fullset$CT <- fullset$lCT + fullset$rCT

# adding mean global weighted FA (comes from 2_ContinuityROI.R)
fullset <- wFA[fullset, on ="eid"]

### saving the result
# data.table::fwrite(fullset, "/Volumes/home/lifespan/nicjud/UKB/proc/20240223_fullset.csv")
# fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240223_fullset.csv")

fullset$CT <- fullset$CT/2 # you summed the hemisphere's for CT, when you want average CT

# getting the time between ROSLA & scanning :) for neuroimaging peeepz
scanage <- fullset[!is.na(visit_date),]

# lets make it a bit more representative of the sample by making sure they have at least 1 global neuro cov

scanage$nodata <- is.na(scanage$TBV_norm) + is.na(scanage$CT) + is.na(scanage$SA) + is.na(scanage$wFA) + is.na(scanage$CSF_norm) + is.na(scanage$WM_hyper)
scanage <- scanage[!nodata == 6] # CANNOT have all 6 missing

# visit_data_correct is just the error (first sub minus everyone)

scanage$DOB <- ym(str_c(scanage$year,"-", scanage$month))

# ggplot(scanage, aes(1, DOB)) + ggrain::geom_rain(point.args = list(alpha = .01))

#age of scan
scanage$AOS <- interval(scanage$DOB, scanage$visit_date) %/% months(1)/12
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
#### 1.4 fuzzy RD for global vars   ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# not allowing missing in SA (Y), running_var or EduAge16

covs <- c("sex", "t2_FLAIR", "visit_day_correct", "visit_day_correct2", "headmotion", "summer",
          "imaging_center_11026", "imaging_center_11027", "imaging_center_11028", 
           "dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1")
covs_fence <- c("dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1")
cols = c("EduAge16", "running_var", covs)

# ct <- ct[,(col_fence.s) := lapply(.SD, vec_to_fence), .SDcols=col_fence.s]


dim(fullset)
# effect of covs
# summary(lm(paste(c("SA ~ acq_time", covs), collapse=" + "), data = fullset))
# summary(lm(paste(c("CT ~ acq_time", covs), collapse=" + "), data = fullset))
# summary(lm(paste(c("WM_hyper ~ acq_time", covs), collapse=" + "), data = fullset))
# summary(lm(paste(c("CSF_norm ~ acq_time", covs), collapse=" + "), data = fullset))
# summary(lm(paste(c("TBV_norm ~ acq_time", covs), collapse=" + "), data = fullset))

sa <- fullset[, c("SA", cols), with=FALSE]
ct = fullset[, c("CT", cols), with=FALSE]
WMh = fullset[, c("WM_hyper", cols), with=FALSE]
CSFn = fullset[, c("CSF_norm", cols), with=FALSE]
TBVn = fullset[, c("TBV_norm", cols), with=FALSE]
wFAs = fullset[, c("wFA", cols), with=FALSE] # wFAs; s for subset since you already have a wFA DT

# not allowing any missing in Y, running or instrument
sa <- sa[complete.cases(sa[, .(SA, running_var, EduAge16)]),][ # no missingness allows in Y, running, or IV
  ,c("SA", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("SA", covs_fence)] # fencing covs with error & Y
ct <- ct[complete.cases(ct[, .(CT, running_var, EduAge16)]),][ 
  ,c("CT", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("CT", covs_fence)]
WMh <- WMh[complete.cases(WMh[, .(WM_hyper, running_var, EduAge16)]),][ 
  ,c("WM_hyper", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("WM_hyper", covs_fence)]
CSFn <- CSFn[complete.cases(CSFn[, .(CSF_norm, running_var, EduAge16)]),][ 
  ,c("CSF_norm", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("CSF_norm", covs_fence)]
TBVn <- TBVn[complete.cases(TBVn[, .(TBV_norm, running_var, EduAge16)]),][ 
  ,c("TBV_norm", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("TBV_norm", covs_fence)]
wFAs <- wFAs[complete.cases(wFAs[, .(wFA, running_var, EduAge16)]),][ 
  ,c("wFA", covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c("wFA", covs_fence)]

#you should also state the amount of missing.

# sa_imp <- MICE_imp(sa, n = 10)
# ct_imp <- MICE_imp(ct, n = 10)
# CSFn_imp <- MICE_imp(CSFn, n = 10)
# TBVn_imp <- MICE_imp(TBVn, n = 10)
# wFAs_imp <- MICE_imp(wFAs, n = 10)
# 
# sa_imp_r <- with(sa_imp, RDjacked("SA", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs))
# ct_imp_r <- with(ct_imp, RDjacked("CT", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs))
# CSFn_imp_r <- with(CSFn_imp, RDjacked("CSF_norm", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs))
# TBVn_imp_r <- with(TBVn_imp, RDjacked("TBV_norm", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs))
# wFAs_imp_r <- with(wFAs_imp, RDjacked("wFA", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs))
# 
# # covs[covs != "t2_FLAIR"]
# WMh_imp <- MICE_imp(WMh[,!"t2_FLAIR"], n = 10)
# WMh_imp_r <- with(WMh_imp, RDjacked("WM_hyper", "running_var", fuzzy = 'EduAge16', df = data.frame(mget(ls())), covs = covs[covs != "t2_FLAIR"]))

# checking missingness in covs
# are they complete?
1-sum(complete.cases(sa))/dim(sa)[1]
1-sum(complete.cases(ct))/dim(ct)[1]
1-sum(complete.cases(WMh))/dim(WMh)[1]
1-sum(complete.cases(CSFn))/dim(CSFn)[1]
1-sum(complete.cases(TBVn))/dim(TBVn)[1]
1-sum(complete.cases(wFAs))/dim(wFAs)[1]
# 0.04126855

# missing dMRI covs
sa = sa[complete.cases(sa),] #34010 to 32102
ct = ct[complete.cases(ct),] #same as SA
# Fuck these 3 have missingness in other Covs
WMh = WMh[complete.cases(WMh),] #32676 to 31698
CSFn = CSFn[complete.cases(CSFn),] #33634 to 32102
TBVn = TBVn[complete.cases(TBVn),] #33634 to 32102
wFAs = wFAs[complete.cases(wFAs),]

m1_sa_fuzzy = RDjacked("SA", "running_var", fuzzy = 'EduAge16', df = sa, covs = covs) #352 t2_FLAIR
m2_ct_fuzzy = RDjacked("CT", "running_var", fuzzy = 'EduAge16', df = ct, covs = covs) #same as above
m3_WMh_fuzzy = RDjacked("WM_hyper", "running_var", fuzzy = 'EduAge16', df = WMh, covs = covs[covs != "t2_FLAIR"]) # only 2 occassions
m4_CSFn_fuzzy = RDjacked("CSF_norm", "running_var", fuzzy = 'EduAge16', df = CSFn, covs = covs) #352 t2_FLAIR
m5_TBVn_fuzzy = RDjacked("TBV_norm", "running_var", fuzzy = 'EduAge16', df = TBVn, covs = covs) #352 t2_FLAIR
m6_wFA_fuzzy = RDjacked("wFA", "running_var", fuzzy = 'EduAge16', df = wFAs, covs = covs) # only 37 occassions

global_results <- rbind(m1_sa_fuzzy, m2_ct_fuzzy, m3_WMh_fuzzy, m4_CSFn_fuzzy, m5_TBVn_fuzzy, m6_wFA_fuzzy)

# global_results %>% #saving the results
#   kbl(caption = "Global neuroimaging results") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/global_results.html")

# making a table without covs

# m1_sa_fuzzy_uY = RDjacked("SA", "running_var", fuzzy = 'EduAge16', df = sa) #352 t2_FLAIR
# m2_ct_fuzzy_uY = RDjacked("CT", "running_var", fuzzy = 'EduAge16', df = ct) #same as above
# m3_WMh_fuzzy_uY = RDjacked("WM_hyper", "running_var", fuzzy = 'EduAge16', df = WMh) # only 2 occassions
# m4_CSFn_fuzzy_uY = RDjacked("CSF_norm", "running_var", fuzzy = 'EduAge16', df = CSFn) #352 t2_FLAIR
# m5_TBVn_fuzzy_uY = RDjacked("TBV_norm", "running_var", fuzzy = 'EduAge16', df = TBVn) #352 t2_FLAIR
# m6_wFA_fuzzy_uY = RDjacked("wFA", "running_var", fuzzy = 'EduAge16', df = wFAs) # only 37 occassions
# 
# global_results_uY <- rbind(m1_sa_fuzzy_uY, m2_ct_fuzzy_uY, m3_WMh_fuzzy_uY, m4_CSFn_fuzzy_uY, m5_TBVn_fuzzy_uY, m6_wFA_fuzzy_uY)
# 
# global_results_uY %>% #saving the results
#   kbl(caption = "Global neuroimaging results uncorrected") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/global_results_uY.html")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 1.3 plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#moving plotting below so I can put the MSE derived bounds

# RDHonest::RDScatter(SA~birth_quarters, data = fullimage, avg = Inf, propdotsize = T, vert = T)
# rdrobust::rdplot(fullimage$SA, fullimage$running_var, p = 3)


Edu16plt <- fullset[!is.na(fullset$visit_date),] %>% #making sure it is imaging subjects
  group_by(running_var) %>% 
  summarise(piEdu16 = sum(EduAge16, na.rm = T)/n()) %>% 
  {ggplot(., aes(running_var, piEdu16)) +
      geom_point(color = "blue", alpha = .3) +
      # geom_point(data=subset(., running_var > -m1_sa_fuzzy$bandwidth & running_var < m1_sa_fuzzy$bandwidth), color = "darkblue") +
      # geom_point(data=subset(., running_var < -m1_sa_fuzzy$bandwidth | running_var  > m1_sa_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,2),se=F, color = "darkgreen") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Percent staying until 16'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(.75, 1)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())
  }
# ggsave("~/Google Drive/My Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_Fig1.png", Edu16plt, width = 10, height = 8, bg = "white")


# surface area
SAplt <- sa %>% 
  group_by(running_var) %>% 
  summarise(sa = mean(SA, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, sa)) +
      geom_point(data=subset(., running_var > -m1_sa_fuzzy$bandwidth & running_var < m1_sa_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m1_sa_fuzzy$bandwidth | running_var  > m1_sa_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = expression("Total Surface Area mm" ^2 ~ ""), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # scale_y_continuous(breaks=c(164000, 168000, 172000, 176000),
      #                    labels=c(expression("164000" ^mm2 ~ ""), "Sept.\n1952", "Sept.\n1957", "Sept.\n1962"),
      #                    limits = c(164000, 176000))+ 
      ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

SAplt0 <- sa %>% 
  group_by(running_var) %>% 
  summarise(sa = mean(SA, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, sa)) +
      geom_point(color = "darkblue") +
      # geom_smooth(method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Total\nSurface Area'), x = "Date of Birth in Months") + # bquote('Total Surface Area'~(mm^3))
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}


# cortical thickness
CTplt <- ct %>% 
  group_by(running_var) %>% 
  summarise(ct = mean(CT, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, ct)) +
      geom_point(data=subset(., running_var > -m2_ct_fuzzy$bandwidth & running_var < m2_ct_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m2_ct_fuzzy$bandwidth | running_var  > m2_ct_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = 'Average Cortical\nThickness mm', x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

CTplt0 <- ct %>% 
  group_by(running_var) %>% 
  summarise(ct = mean(CT, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, ct)) +
      geom_point(color = "darkblue") +
      # geom_smooth(method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = 'Average Cortical\nThickness', x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}
 
# WMh
WMhplt <- WMh %>% 
  group_by(running_var) %>% 
  summarise(WM_hyper = mean(WM_hyper, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, WM_hyper)) +
      geom_point(data=subset(., running_var > -m3_WMh_fuzzy$bandwidth & running_var < m3_WMh_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m3_WMh_fuzzy$bandwidth | running_var  > m3_WMh_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = 'White Matter\nHyperintensities', x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 15) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

WMhplt0 <- WMh %>% 
  group_by(running_var) %>% 
  summarise(WM_hyper = mean(WM_hyper, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, WM_hyper)) +
      geom_point(color = "darkblue") +
      labs(y = 'White Matter\nHyperintensities', x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}


# CSF
CSFnplt <- CSFn %>% 
  group_by(running_var) %>% 
  summarise(CSF_norm = mean(CSF_norm, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, CSF_norm)) +
      geom_point(data=subset(., running_var > -m4_CSFn_fuzzy$bandwidth & running_var < m4_CSFn_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m4_CSFn_fuzzy$bandwidth | running_var  > m4_CSFn_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = "Cerebrospinal\nFluid Volume", x = "Date of Birth in Months") + #this is in ml's
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 15) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

CSFnplt0 <- CSFn %>% 
  group_by(running_var) %>% 
  summarise(CSF_norm = mean(CSF_norm, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, CSF_norm)) +
      geom_point(color = "darkblue") +
      labs(y = "Cerebrospinal\nFluid Volume", x = "Date of Birth in Months") + #this is in ml's
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

# TBV
TBVnplt <- TBVn %>% 
  group_by(running_var) %>% 
  summarise(TBV_norm = mean(TBV_norm, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, TBV_norm)) +
      geom_point(data=subset(., running_var > -m5_TBVn_fuzzy$bandwidth & running_var < m5_TBVn_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m5_TBVn_fuzzy$bandwidth | running_var  > m5_TBVn_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = expression("Total Brain Volume mm" ^3 ~ ""), x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

TBVnplt0 <- TBVn %>% 
  group_by(running_var) %>% 
  summarise(TBV_norm = mean(TBV_norm, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, TBV_norm)) +
      geom_point(color = "darkblue") +
      # geom_smooth(method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = 'Total Brain\nVolume', x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}


# wFA
wFAplt <- wFAs %>% 
  group_by(running_var) %>% 
  summarise(wFA = mean(wFA, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, wFA)) +
      geom_point(data=subset(., running_var > -m6_wFA_fuzzy$bandwidth & running_var < m6_wFA_fuzzy$bandwidth), color = "darkblue") +
      geom_point(data=subset(., running_var < -m6_wFA_fuzzy$bandwidth | running_var  > m6_wFA_fuzzy$bandwidth), color = "blue", alpha = .3) +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_smooth(data=subset(., running_var < 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      geom_smooth(data=subset(., running_var > 0), method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Weighted\nFractional Anisotropy'), x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      scale_y_continuous(breaks=c(.440, .445, .450)) + #.442,.444, .446, .448
      # ylim(c(164000, 176000)) +
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

wFAplt0 <- wFAs %>% 
  group_by(running_var) %>% 
  summarise(wFA = mean(wFA, na.rm = T), n =n()) %>% 
  {ggplot(., aes(running_var, wFA)) +
      geom_point(color = "darkblue") +
      # geom_smooth(method='glm',formula=y~poly(x,3),se=F, color = "darkred") +
      labs(y = bquote('Weighted\nFractional Anisotropy'), x = "Date of Birth in Months") +
      scale_x_continuous(breaks=c(-120,-60,0,60,120),
                         labels=c("Sept.\n1947", "Sept.\n1952", "Sept.\n1957", "Sept.\n1962", "Sept.\n1967")) +
      scale_y_continuous(breaks=c(.440, .445, .450)) + #.442,.444, .446, .448
      theme_minimal(base_size = 20) +
      theme(axis.text.x= element_text(angle=45), axis.title.x = element_blank())}

SAplt + CSFnplt + CTplt + WMhplt + TBVnplt + wFAplt + 
  plot_annotation(tag_levels = 'a')

# ggsave("~/Google Drive/My Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/Fig1.png")

SAplt +  CTplt + TBVnplt + wFAplt + 
  plot_annotation(tag_levels = 'a')
# did base size 20
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/Fig1_4_take2.png", width = 16, height = 9)

# plots without RD effect
SAplt0 + CSFnplt0 + CTplt0 + WMhplt0 + TBVnplt0 + wFAplt0 + 
  plot_annotation(tag_levels = 'a')
# did base size 20
# ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/Fig1_6_0Blank.png", width = 20, height = 8)


CSFnplt / WMhplt + 
  plot_annotation(tag_levels = 'a')

# ggsave("~/Google Drive/My Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/SIfig2.png", width = 10, height = 8)

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


