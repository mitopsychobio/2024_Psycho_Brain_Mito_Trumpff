#Load libraries
library(broom)
library(tidyverse)
library(here)
library(stringr)
library(stringi)
library(DescTools)
library(purrr)
library(broom)
library(modelr) 


####### TMT import and cleaning------

#Import TMT and merge------
key <- read.csv(here::here("data", "Key_TMT400.csv")) %>% #import key id to merge proteomics data
mutate(proteomicsid = (chartr('_', '.', proteomicsid)))
  
tmt <- read.csv(here::here("data", "TMT400logr_v2.csv"))%>%
mutate(Unique.ID = (chartr('|', '.', Unique.ID)))

tmt$proteomicsid 

write.csv(tmt, here::here("data","tmt.list_generow.csv"))
 

#Transpose data
tmt1 <- select(tmt, -c(Gene))
tmt2 <- tmt1 %>%
  unique(.)%>%
  remove_rownames %>%
  column_to_rownames(var = "Unique.ID")%>%
  t(.)%>%
  as.data.frame(.)%>%
  rownames_to_column(var = "proteomicsid")

tmt3 <- tmt2 %>%
  right_join(key,., by ="proteomicsid")%>%
  select(-c(proteomicsid, X))
  
#tmt3$projid  <- stri_sub(tmt3$projid, 2)  
write.csv(tmt3, here::here("data","tmt.list_idrow.csv"))

####### Compute residuals------------
pheno <- read.csv(here::here("data", "datarosmap_all_2020-03-02.csv"))%>% 
  select(projid, pmi, study, age_death, cogs,  msex)%>% 
  unique()
#nest the data
nested <- tmt3 %>% 
  arrange(projid)  %>%
  pivot_longer(-projid, names_to = "gene_id", values_to = "abundance") %>%
  filter(!is.na(abundance)) %>%
  left_join(pheno, by = "projid") %>% 
  group_by(gene_id) %>%
  nest()
# Setup the model 
model <- nested %>%
  mutate(model = map(data, ~ lm(abundance ~ pmi + study + age_death, data=.x))) %>%
  mutate(resids = map2(data, model, add_residuals))

# get residuals
resids <- unnest(model, resids) %>%
  select(projid, gene_id, resid) %>%
  arrange(resid) %>%
  #  select(-projid) %>% ##Anna: I commented this line out
  ungroup() %>%
  unique() %>% ##Anna: You've had some duplicates in, not sure why
  pivot_wider(names_from = "projid", values_from = "resid") ##For a wide format

#export ID row
write.csv(resids, here::here("data", "tmt.residuals_generow.csv"))

#export gene row
tmtid <- resids %>%
  unique(.)%>%
  remove_rownames %>%
  column_to_rownames(var = "gene_id")%>%
  t(.)%>%
  as.data.frame(.)%>%
  rownames_to_column(var = "projid")%>%
  write.csv(here::here("data", "tmt.residuals_idrow.csv"))

####### Mitocarta score--------------
#Import data
tmtg <- read.csv(here::here("data", "tmt.residuals_idrow.csv"))%>%
select(-c(X))
list<- read.csv(here::here("data", "gene_list_2021-03-04.csv"))
list$Gene <- str_replace_all(list$Gene, "-", ".") #change - to . since character changes when transposed

tmt_list <- read.csv(here::here("data", "TMT400logr_v2.csv"))%>%
        select(Unique.ID, Gene)%>%
        merge(list, by=c("Gene"))%>%
        mutate(Unique.ID = (chartr('|', '.', Unique.ID)))%>%
        mutate(Unique.ID = (chartr('-', '.', Unique.ID)))%>%
        #select gene list items 
        select(Unique.ID, Gene, gene_cate, mitocarta,	mt_cate_1,	mt_cate_2,	mt_cate_3, RC, RC_gene_cate, RC1_gene_cate_loc, RC1_gene_cate_modules, RC1_gene_cate_modules_2, core_sub, accessory_sub, assembly_factors, housekeeping)

###Data preparation------------------------------------------------------------------
#create t score
#create function

center_scale <- function(x) {
  scale(x)
}
center_scale_2 <- function(x) {
  (x*10)+100
}
#run the function
data <- tmtg
datac1 <- data[,-1]
datac2 <- center_scale(datac1)
datac2 <- center_scale_2 (datac2)
dataid <- select(data, 1)
datas <- cbind(dataid, datac2)
data <- datas
rm(datac1, datac2, dataid, datas, center_scale, center_scale_2)

#merge and pivot
datal <- data %>%
  pivot_longer(-projid, names_to = "Unique.ID", values_to = "abundance") %>%
  left_join(tmt_list, by = "Unique.ID", all) %>%
  arrange(projid) %>%
  unique()
###Compute the scores------------------------
library(dplyr)
###housekeeping protein score---------
data <- datal
datah <- data %>% 
  filter(housekeeping == 1) %>% 
  group_by(projid) %>%
  mutate(hk_score = mean(abundance, na.rm = TRUE)) 
#select score to be merged at the end  
hk_score <-datah %>% 
  select(projid, hk_score) 
hk_score <- unique(hk_score)
hk_score
#check
print <- datah %>%
  select(projid, Gene, abundance, hk_score ) %>%
  arrange(Gene)
print
rm(print, datah)
###Mitocarta score----------
datah <- data %>% 
  filter(mitocarta == 1) %>% 
  group_by(projid) %>%
  mutate(mitocarta = mean(abundance, na.rm = TRUE)) 
datah
#select score to be merged at the end  
mitocarta_score <-datah %>% 
  select(projid, mitocarta) 
mitocarta_score <- unique(mitocarta_score)
mitocarta_score
#check
print <-select(datah,	 projid, Gene, abundance, mitocarta)%>%
  arrange(Gene)
print
rm(print, datah)
##Average per categories mt_cate_1-------
data2 <- data %>% 
  filter(mitocarta == 1) %>% 
  group_by(mt_cate_1, projid) %>%
  mutate(ab_mt_cate_1 = mean(abundance, na.rm = TRUE))
data2
ungroup(data2, mt_cate_1, projid)
#pivot the data 
mitocarta_cate1 <-  data2  %>% 
  select(projid,  mt_cate_1, ab_mt_cate_1)%>%
  unique()%>%
  pivot_wider(names_from = mt_cate_1, values_from =  ab_mt_cate_1) %>%
  select(-c("NA"))
rm(data2)
#Average per categories mt_cate_2-------
data2 <- data %>% 
  group_by(mt_cate_2, projid) %>%
  mutate(ab_mt_cate_2 = mean(abundance, na.rm = TRUE))
data2
ungroup(data2, mt_cate_2, projid)
#pivot the data 
mitocarta_cate2 <-  data2  %>% 
  select(projid,  mt_cate_2, ab_mt_cate_2)%>%
  unique()%>%
  pivot_wider(names_from = mt_cate_2, values_from =  ab_mt_cate_2) %>%
  select(-c("NA"))
mitocarta_cate2
#Average per categories mt_cate_3-----------
data2 <- data %>% 
  group_by(mt_cate_3, projid) %>%
  mutate(ab_mt_cate_3 = mean(abundance, na.rm = TRUE))
data2
ungroup(data2, mt_cate_3, projid)
#pivot the data 
mitocarta_cate3 <-  data2  %>% 
  select(projid,  mt_cate_3, ab_mt_cate_3)%>%
  unique()%>%
  pivot_wider(names_from = mt_cate_3, values_from =  ab_mt_cate_3) %>%
  select(-c("NA"))
mitocarta_cate3
#Average per categories abundance_gene_cate-------
data2 <- subset(data, assembly_factors != "1"|is.na(assembly_factors))
data3 <- data2 %>% 
  group_by(gene_cate, projid) %>%
  mutate(ab_gene_cate = mean(abundance, na.rm = TRUE))
ungroup(data3, gene_cate, projid)

write.csv(data3, here::here("data", "test.csv"))

#pivot the data 
mhi_cate <-  data3  %>% 
  select(projid, gene_cate, ab_gene_cate)%>%
  unique()%>%
  pivot_wider(names_from = gene_cate, values_from =  ab_gene_cate) %>%
  select(-c("NA"))
mhi_cate

write.csv( mhi_cate, here::here("data", "test.csv"))

### merge scores----------
datat <- merge(hk_score, mitocarta_score, by=c("projid"), all=TRUE)
datat <- merge(datat, mitocarta_cate1, by=c("projid"), all=TRUE)
datat <- merge(datat, mitocarta_cate2, by=c("projid"), all=TRUE)
datat <- merge(datat, mitocarta_cate3, by=c("projid"), all=TRUE)
datat <- merge(datat, mhi_cate, by=c("projid"), all=TRUE)
### compute score adjusted for mitocontent--------
datat <- datat %>%
  mutate(nc1_scorem = (nc1/mitocarta))%>%
  mutate(nc2_scorem = (nc2/mitocarta))%>%
  mutate(nc3_scorem = (nc3/mitocarta))%>%
  mutate(nc4_scorem = (nc4/mitocarta))%>%
  mutate(nc5_scorem = (nc5/mitocarta))%>%
  mutate(mtc1_scorem = (mtc1/mitocarta))%>% 
  mutate(mtc5_scorem = (mtc5/mitocarta))
###export-------
write.csv(datat, here::here("data", "mitocarta_score_tmt.csv"))
######################################
###### Clinical data processing------
#import
d678l <- read.csv(here::here("data", "dataset_678_long_09-16-2018.csv"))
d707l <- read.csv(here::here("data", "dataset_707_long_12-03-2019.csv"))
d157l <- read.csv(here::here("data", "dataset_157_long_06-19-2018.csv"))%>%
  select(-c(scaled_to, neglifeevents, panas, perceivedstress))
d678b <- read.csv(here::here("data", "dataset_678_basic_09-16-2018.csv"))%>%
  select (-c(angerin, angerout, apoe_genotype, educ, msex, spanish, age_death))
d707b <- read.csv(here::here("data", "dataset_707_basic_12-11-2018.csv"))
ctype <- read.csv(here::here("data", "Cell_proportion_2019-11-04.csv"))
mtcn_new <- read.csv(here::here("data", "mtDNAcn_Nov2020.csv"))
#longitudinal
#merge
dbb <-  merge(d157l, d707l, by=c("projid", "fu_year", "study"), all=TRUE)%>%
        merge(d678l, by=c("projid", "fu_year", "study", "age_at_visit"), all=TRUE)%>%
        merge(d707b, by=c("projid", "study", "scaled_to"), all=TRUE)
#cross-sectional
cross <- d707b %>%
        merge( d678b, by=c("projid", "study", "scaled_to"), all=TRUE)%>%
        merge(ctype, by=c("projid"), all=TRUE)%>%
        merge(mtcn_new, by=c("projid"), all=TRUE)
#data cleaning 
###AD status
###recode NIA-reagan diagnosis
cross$ad_reagan <- ifelse(cross$niareagansc==1, 0, ifelse(cross$niareagansc==2, 0, ifelse(cross$niareagansc==3, 1, ifelse(cross$niareagansc==4, 1, NA))))
###recode cognitive status cogs NCI=1, MCI=2, AD=3
cross$cogs <- ifelse(cross$cogdx==1, 1, ifelse(cross$cogdx==2, 2, ifelse(cross$cogdx==3, 2, ifelse(cross$cogdx==4, 3, ifelse(cross$cogdx==5, 3, NA)))))
###year before death
dbb <- mutate(dbb, visit_yrs_death = (age_death - age_at_visit))
dbb <- mutate(dbb, visit_yrs_deathr = round(visit_yrs_death))
#export
write.csv(dbb, here::here("data","rosmap_clinical_longitudinal.csv"))
write.csv(cross, here::here("data","rosmap_clinical_cross.csv"))
###### Psychosocial score------------
# open library
library(dplyr)
library(readr)
library(tidyverse)
library(pastecs)
library(frequency)
library(ggplot2)
library(devtools)
library(easyGgplot2)
library(gridExtra)
library(frequency)
library(Rmisc) 
library(magrittr)
#Import data
#dbb <- read.csv(here::here("data", "rosmap_clinical.csv"))
#store data not to score and merge later
cross <-select(data, projid, study, age_death, msex, tot_adverse_exp,	haanticipatoryworry,  neuroticism_12, angerin, angerout, angertrait, extraversion_6, openness, conscientiousness, agreeableness)%>%
  unique()
#create reverse score for tot_adverse_exp
#tot_adverse_exp_r
#stat.desc(cross$tot_adverse_exp)
cross$tot_adverse_exp_r = 45-cross$tot_adverse_exp

cross$tot_adverse_exp_r
#
cross2 <-select(data, projid, study, msex)
#select longitudinal variables
long <-select(dbb, projid, fu_year, visit_yrs_death,	cesdsum,	late_life_soc_act,	soc_net,	social_isolation,	purpose_total,	satisfaction,	timehoriz,	wellbeing,	wellbeing_hedonic,	wellbeing_eudaimonic,	wellbeing_auton,	wellbeing_envmas,	wellbeing_growth,	wellbeing_posrel,	wellbeing_purpose,	wellbeing_accept,	neglifeevents,	panas,	perceivedstress,	discrim_cnt, negsocexchange, rejection)
#make all the variables numeric
    long1 <- long[,-1]
    long2 <- sapply( long1, as.numeric )
    #put back ID column
    longid <- select(long, 1)
    long <- cbind(longid, long2)
    long
### create reverse variable for positive questionnaire to create a overall negative score
    #satisfaction
    long$satisfaction_r = 7-long$satisfaction
    #late_life_soc_act
    long$late_life_soc_act_r = 6-long$late_life_soc_act
    #wellbeing
    long$wellbeing_r = 8-long$wellbeing
    #timehoriz
    long$timehoriz_r = 8-long$timehoriz
    #purpose_total
    long$purpose_total_r = 6-long$purpose_total
    #cesdsum_r_avg
    long$cesdsum_r = 11-long$cesdsum
    #social_isolation_r_avg	
    long$social_isolation_r = 6-long$social_isolation
    #neglifeevents_r_avg
    long$neglifeevents_r = 18-long$neglifeevents
    #panas_r_avg
    long$panas_r = 6-long$panas
    #perceivedstress_r_avg
    long$perceivedstress_r = 5-long$perceivedstress

#count the number of years of follow up
longid <- select(long, projid, fu_year)
longid %>% group_by(fu_year) %>% tally()
nbyear <- longid %>% group_by(projid) %>% tally()

#Score avg life
# Average all years
long_avg <-long %>% 
  group_by(projid) %>% 
  summarize_all(mean, na.rm = TRUE)
ungroup(long_avg)
#add a suffix
    #remove ID column
    long_avg1 <- long_avg[,-1]
    colnames(long_avg1) <- paste(colnames(long_avg1), "avg", sep = "_")
    #put back ID column
    long_avgid <- select(long_avg, 1)
    long_avg <- cbind(long_avgid, long_avg1)
#merge all new data
#merge with all the data
psy<- merge(cross, long_avg, by=c("projid" ), all =TRUE)

#POSITIVE & NEGATIVE PSYCHOSOCIAL SCORES####################
#######SCALE############
#Remove the variable not to scale
psys <-select(psy, projid, cesdsum_r_avg,	social_isolation_r_avg,	neglifeevents_r_avg,	panas_r_avg,	perceivedstress_r_avg,	tot_adverse_exp_r, late_life_soc_act_r_avg,	purpose_total_r_avg,	satisfaction_r_avg,	timehoriz_r_avg,	wellbeing_r_avg, tot_adverse_exp, neuroticism_12,	haanticipatoryworry,	fu_year_avg,	visit_yrs_death_avg,	cesdsum_avg,	late_life_soc_act_avg,	soc_net_avg,	social_isolation_avg,	purpose_total_avg,	satisfaction_avg,	timehoriz_avg,	wellbeing_avg,	wellbeing_hedonic_avg,	wellbeing_eudaimonic_avg,	wellbeing_auton_avg,	wellbeing_envmas_avg,	wellbeing_growth_avg,	wellbeing_posrel_avg,	wellbeing_purpose_avg,	wellbeing_accept_avg,	neglifeevents_avg,	panas_avg,	perceivedstress_avg,	discrim_cnt_avg)
#to merge later
notscore <- select(psy, projid,	study, age_death,	msex)
#check that all the variables are numerics
dt1 <- psys[,-1]
dt <- sapply( dt1, as.numeric)

#scale the data
center_scale <- function(x) {
  scale(x)
}
dt1 <- center_scale(dt)
dtid <- select(psys, 1)
dts <- cbind(dtid, dt1)
#create a t-score
#tscore
#transfert all but not the first column
dtt <- dts
dtt1 <- dtt [,-1]
center_scale <- function(x) {
  (x*10)+100
}
dtt2 <- center_scale(dtt1)
dttid <- select(dtt, 1)
dtt <- cbind(dttid, dtt2)
#create positive and negative exposure score
##########Avg all years - No missing allowed in the avg
#POSITIVE
dtt$pos_score_avg_nomiss<- rowMeans(subset(dtt, select = c(late_life_soc_act_avg,	purpose_total_avg,	satisfaction_avg,	timehoriz_avg,	wellbeing_avg)), na.rm=F) 
#NEGATIVE
dtt$neg_score_avg_nomiss<- rowMeans(subset(dtt, select = c(cesdsum_avg,	social_isolation_avg,	neglifeevents_avg,	panas_avg,	perceivedstress_avg,	tot_adverse_exp)), na.rm=F) 
########## COMBINED SCORE NEG ##HERE
dtt$combined_neg_score_avg_nomiss<- rowMeans(subset(dtt, select = c(cesdsum_avg,	social_isolation_avg,	neglifeevents_avg,	panas_avg,	perceivedstress_avg,	tot_adverse_exp, late_life_soc_act_r_avg,	purpose_total_r_avg,	satisfaction_r_avg,	timehoriz_r_avg,	wellbeing_r_avg)), na.rm=F) 
########## COMBINED POSITIVE ##HERE
dtt$combined_pos_score_avg_nomiss<- rowMeans(subset(dtt, select = c(cesdsum_r_avg,	social_isolation_r_avg,	neglifeevents_r_avg,	panas_r_avg,	perceivedstress_r_avg,	tot_adverse_exp_r, late_life_soc_act_avg,	purpose_total_avg,	satisfaction_avg,	timehoriz_avg,	wellbeing_avg)), na.rm=F) 
#select the data to merge
psyscore <- select(dtt, projid, pos_score_avg_nomiss,	neg_score_avg_nomiss, combined_neg_score_avg_nomiss, combined_pos_score_avg_nomiss)
#merge
psyf <- merge(psy, psyscore, by=c("projid"), all=TRUE)
#export
write.csv(psyf, here::here("data","rosmap_psy_score.csv"))

############# Descriptives------------------

pheno <- read.csv(here::here("data", "datarosmap_all_2020-03-02.csv"))%>% 
  select(projid, pmi, study, age_death, cogs,  msex, study, amyloid, tangles)%>% 
  unique()

psy <- read.csv(here::here("data", "rosmap_psy_score.csv"))%>%
  select(projid, pos_score_avg_nomiss,  neg_score_avg_nomiss)

mito <- read.csv(here::here("data","mitocarta_score_tmt.csv"))%>% 
  select(projid, mitocarta)

nested <- mito %>%
  left_join(psy, by = "projid")%>%
  left_join(pheno, by = "projid")

#Continuous variables
library (pastecs)
desc <- stat.desc(nested)
write.csv(desc, here::here("results", "descriptive_continuous.csv"))

#Categorical variables

nested%>%
  group_by(cogs) %>%
  summarise(count=n())

nested%>%
  group_by(study) %>%
  summarise(count=n())

nested%>%
  group_by(msex) %>%
  summarise(count=n())

############# Spearman rho------------------

#1. all mitocarta scores--------
#Import psy data
psy <- read.csv(here::here("data", "rosmap_psy_score.csv"))%>%
  select(projid, pos_score_avg_nomiss,  neg_score_avg_nomiss)%>%
  pivot_longer(-projid, names_to = "psy_names", values_to = "psy_score")#pivot in long format

#Import tmt and nest
datat <- read.csv(here::here("data","mitocarta_score_tmt.csv"))
nested <- datat %>% 
  arrange(projid)  %>%
  pivot_longer(-projid, names_to = "mhi_index", values_to = "mhi_score") %>%
  left_join(psy, by = "projid") %>% 
  group_by(mhi_index, psy_names) %>%
  nest()
head(nested)
#model
model <- nested %>% 
  mutate(cor_test3 = map(data, ~SpearmanRho(.x$psy_score, .x$mhi_score, use = c("pairwise.complete.obs"), conf.level = 0.95)), 
         tidied = map(cor_test3, tidy)) %>% 
  unnest(tidied) %>% 
  select(mhi_index, psy_names, names, x) %>% 
  pivot_wider(names_from = "names", values_from = "x")

model <- model%>% 
  select(mhi_index, psy_names, lwr.ci, rho, upr.ci)      
#export the results           
write.csv(model, here::here("results","spearman_mitocarta.csv"))
rm(model)

#2. Individual questionnaires and oxphos scores----------
psy <- read.csv(here::here("data", "rosmap_psy_score.csv"))%>%
select(projid, soc_net_avg, late_life_soc_act_avg,  purpose_total_avg, satisfaction_avg, timehoriz_avg, wellbeing_avg, pos_score_avg_nomiss, tot_adverse_exp, neglifeevents_avg, social_isolation_avg, cesdsum_avg, panas_avg, perceivedstress_avg, neg_score_avg_nomiss)%>%
pivot_longer(-projid, names_to = "psy_names", values_to = "psy_score")

#Import tmt and nest
datat <- read.csv(here::here("data","mitocarta_score_tmt.csv"))%>%
 select(projid, nc1_scorem, 	nc2_scorem, 	nc3_scorem, 	nc4_scorem, 	nc5_scorem, 	mtc1_scorem, 	mtc5_scorem)
  
nested <- datat %>% 
  arrange(projid)  %>%
  pivot_longer(-projid, names_to = "mhi_index", values_to = "mhi_score") %>%
  left_join(psy, by = "projid") %>% 
  group_by(mhi_index, psy_names) %>%
  nest()
head(nested)
#model
model <- nested %>% 
  mutate(cor_test3 = map(data, ~SpearmanRho(.x$psy_score, .x$mhi_score, use = c("pairwise.complete.obs"), conf.level = 0.95)), 
         tidied = map(cor_test3, tidy)) %>% 
  unnest(tidied) %>% 
  select(mhi_index, psy_names, names, x) %>% 
  pivot_wider(names_from = "names", values_from = "x")

model <- model%>% 
  select(mhi_index, psy_names, lwr.ci, rho, upr.ci)      
#export the results           
write.csv(model, here::here("results","spearman_oxphos_ind_questionnaires.csv"))
rm(model)

###### Multivariate regression analysis---------


#import
all <- read.csv(here::here("data", "cleaned_tmt400_residuals_dataset_ALL_2021-10-21.csv"))%>%
  unique()
psy <- read.csv(here::here("data", "psycho_dataset_2023-02-28.csv"))%>%
  select(projid, pos_score_avg_nomiss,	neg_score_avg_nomiss, combined_neg_score_avg_nomiss, combined_pos_score_avg_nomiss)%>%
  unique()
data <- merge(all, psy, by=c("projid"), all=TRUE)


#compute the residuals of nc1_scorem adjusted for  cogs + msex
sel <- select(data, projid, nc1_scorem, cogs, msex)%>%
  unique()%>%
  filter(nc1_scorem > 0)%>%
  filter(cogs > 0)

model <- lm(nc1_scorem ~  cogs + msex, sel)
sel$nc1_scorem_r<-model$resid
sel$nc1_scorem_r

sel <- select(data, projid, combined_neg_score_avg_nomiss, combined_pos_score_avg_nomiss)%>%
  merge(sel, by=c("projid"), all=TRUE)


######multivariate linear regression
#combined pos

modeln <- lm(nc1_scorem  ~ combined_pos_score_avg_nomiss + cogs + msex, data)
summary(modeln)


modeln <- lm(nc1_scorem  ~ combined_neg_score_avg_nomiss + cogs + msex, data)
summary(modeln)


# pos and neg
modeln <- lm(nc1_scorem  ~ neg_score_avg_nomiss + pos_score_avg_nomiss + cogs + msex, data)
summary(modeln)

#pos
modeln <- lm(nc1_scorem  ~ pos_score_avg_nomiss + cogs + msex, data)
summary(modeln)

#neg
modeln <- lm(nc1_scorem  ~ neg_score_avg_nomiss + cogs + msex, data)
summary(modeln)


####. Combined positive score -----

#result used in Fig 2.

sel <- select(data, projid, combined_pos_score_avg_nomiss, nc1_scorem, cogs, msex)%>%
  unique()%>%
  filter(nc1_scorem > 0)%>%
  filter(cogs > 0)%>%
  filter(combined_pos_score_avg_nomiss > 0)

model <- lm(nc1_scorem  ~ combined_pos_score_avg_nomiss + cogs + msex, data)
model
summary(model)

#proportion of variance
model <- lm(nc1_scorem ~  cogs + msex, sel)
sel$nc1_scorem_r<-model$resid
sel$nc1_scorem_r

model <- lm(nc1_scorem_r  ~ combined_pos_score_avg_nomiss , sel)
model
summary(model)

rm(model)


