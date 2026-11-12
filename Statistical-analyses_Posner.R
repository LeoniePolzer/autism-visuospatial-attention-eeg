# =============================================================
# Posner paradigm EEG time-frequency - statistical analyses
# -------------------------------------------------------------
# Author: <Leonie Polzer>
# this script uses R version 4.5.0

## Description:
# data information: Posner task (70% valid cues, 30% invalid cues) in autistic and neurotypical children and youth
# analysis:
# 1) group differences in reaction times and neural (time-frequency) data
# 2) associations between neural data and reaction times
# 3) associations between neural measures
# 4) prediction of clinical symptoms with neural measures
# EEG data have been preprocessed in MATLAB
# =============================================================

# TRANSFORM DATA INTO REQUIRED DATA FORMAT #############################################################################################################

#### load required packages ----------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(flextable)
library(officer)

options(scipen=999)

#### Load data -----------------------------------------------------------------------------------------------------------------------------

## Define base output folders 
base_dir <- "...\\TF_Final_fieldtrip"
trial_dir  <- file.path(base_dir, "trialwise_fieldtrip")
itc_dir    <- file.path(base_dir, "itc_subject_fieldtrip")
ispc_dir   <- file.path(base_dir, "ispc_subject_fieldtrip")

## Read all trialwise CSVs 
trial_files <- list.files(trial_dir, pattern = "_trialwise\\.csv$", full.names = TRUE)

trial_df <- trial_files |>
  map_dfr(read_csv, show_col_types = FALSE, .id = "file_id") |>
  mutate(
    subj_id = as.factor(subj_id),
    which   = factor(which, levels = c("cue", "targ")),
    side    = factor(side, levels = c("L", "R")),
    hemi    = factor(hemi, levels = c("ipsi", "contra")),
    band    = factor(band, levels = c("alpha", "beta", "gamma_bb")),
    time_window = as.factor(time_window),
    ROI     = as.factor(ROI)
  )

cat(sprintf("Loaded %d trialwise rows from %d subjects\n",
            nrow(trial_df), length(unique(trial_df$subj_id))))

## Read all ITC CSVs - split into raw and trial-number-matched
itc_all <- list.files(itc_dir, pattern = "_itc_subject.*\\.csv$", full.names = TRUE)

itc_files_raw      <- itc_all[!grepl("_itc_subject_matchedN\\.csv$", itc_all)]
itc_files_matchedN <- itc_all[ grepl("_itc_subject_matchedN\\.csv$", itc_all)]

itc_df <- itc_files_raw |>
  map_dfr(read_csv, show_col_types = FALSE, .id = "file_id") |>
  mutate(
    subj_id = as.factor(subj_id),
    which   = factor(which, levels = c("cue", "targ")),
    side    = factor(side, levels = c("L", "R")),
    hemi    = factor(hemi, levels = c("ipsi", "contra")),
    band    = factor(band, levels = c("alpha", "beta", "gamma_bb")),
    ROI     = as.factor(ROI),
    validity = as.factor(validity),
    SOA_label = as.factor(SOA_label),
    time_window = as.factor(time_window)
  )

itc_df_matchedN <- itc_files_matchedN |>
  map_dfr(read_csv, show_col_types = FALSE, .id = "file_id") |>
  mutate(
    subj_id = as.factor(subj_id),
    which   = factor(which, levels = c("cue", "targ")),
    side    = factor(side, levels = c("L", "R")),
    hemi    = factor(hemi, levels = c("ipsi", "contra")),
    band    = factor(band, levels = c("alpha", "beta", "gamma_bb")),
    ROI     = as.factor(ROI),
    validity = as.factor(validity),
    SOA_label = as.factor(SOA_label),
    time_window = as.factor(time_window)
  )

cat(sprintf("Loaded ITC raw: %d rows from %d files\n", nrow(itc_df), length(itc_files_raw)))
cat(sprintf("Loaded ITC matchedN: %d rows from %d files\n", nrow(itc_df_matchedN), length(itc_files_matchedN)))


### Read all ISPC CSVs - split into raw and trial-number-matched
ispc_all <- list.files(ispc_dir, pattern = "_ispc_subject.*\\.csv$", full.names = TRUE)

ispc_files_raw      <- ispc_all[!grepl("_ispc_subject_matchedN\\.csv$", ispc_all)]
ispc_files_matchedN <- ispc_all[ grepl("_ispc_subject_matchedN\\.csv$", ispc_all)]

# Read raw ISPC
ispc_df <- ispc_files_raw |>
  map_dfr(read_csv, show_col_types = FALSE, .id = "file_id") |>
  rename(hemi = laterality) |>
  mutate(
    subj_id = as.factor(subj_id),
    which   = factor(which, levels = c("cue", "targ")),
    side    = factor(side, levels = c("L", "R")),
    hemi    = factor(hemi, levels = c("ipsi", "contra")),
    band    = factor(band, levels = c("alpha", "beta", "gamma_bb")),
    roi1    = as.factor(roi1),
    roi2    = as.factor(roi2),
    roi_pair = factor(paste(roi1, roi2, sep = "_")),
    validity = as.factor(validity),
    SOA_label = as.factor(SOA_label),
    time_window = as.factor(time_window)
  )

# Read matched-N ISPC
ispc_df_matchedN <- ispc_files_matchedN |>
  map_dfr(read_csv, show_col_types = FALSE, .id = "file_id") |>
  rename(hemi = laterality) |>
  mutate(
    subj_id = as.factor(subj_id),
    which   = factor(which, levels = c("cue", "targ")),
    side    = factor(side, levels = c("L", "R")),
    hemi    = factor(hemi, levels = c("ipsi", "contra")),
    band    = factor(band, levels = c("alpha", "beta", "gamma_bb")),
    roi1    = as.factor(roi1),
    roi2    = as.factor(roi2),
    roi_pair = factor(paste(roi1, roi2, sep = "_")),
    validity = as.factor(validity),
    SOA_label = as.factor(SOA_label),
    time_window = as.factor(time_window)
  )

cat(sprintf("Loaded ISPC raw: %d rows from %d files\n", nrow(ispc_df), length(ispc_files_raw)))
cat(sprintf("Loaded ISPC matchedN: %d rows from %d files\n", nrow(ispc_df_matchedN), length(ispc_files_matchedN)))

# sanity check
setdiff(ispc_all, c(ispc_files_raw, ispc_files_matchedN))  # should be character(0)


# Sanity summaries
trial_df        |> count(which, band, ROI, hemi)
itc_df          |> count(which, band, ROI, hemi)
itc_df_matchedN  |> count(which, band, ROI, hemi)
ispc_df         |> count(which, band, roi1, roi2, hemi)
ispc_df_matchedN |> count(which, band, roi1, roi2, hemi)

# clean ID names: delete dates, testing times etc.
colnames(itc_df)[2]<-"id"
itc_df$id<-gsub("\\_T2.*","",itc_df$id)
itc_df$id<-toupper(itc_df$id)
itc_df$id<-gsub("\\_POSNER.*","",itc_df$id)
itc_df$id<-gsub("-", "_", itc_df$id)
itc_df$id<-gsub("ASS", "AUT", itc_df$id)
itc_df$id<-gsub("F_AUT", "AUT_F", itc_df$id)
itc_df$id<-gsub("F_AUT", "AUT_F", itc_df$id)
itc_df$id[itc_df$id=="QS_A_001"]<-"SOKO_AUT_001"

# clean ID names: delete dates, testing times etc.
colnames(trial_df)[2]<-"id"
trial_df$id<-gsub("\\_T2.*","",trial_df$id)
trial_df$id<-toupper(trial_df$id)
trial_df$id<-gsub("\\_POSNER.*","",trial_df$id)
trial_df$id<-gsub("-", "_", trial_df$id)
trial_df$id<-gsub("ASS", "AUT", trial_df$id)
trial_df$id<-gsub("F_AUT", "AUT_F", trial_df$id)
trial_df$id<-gsub("F_AUT", "AUT_F", trial_df$id)
trial_df$id[trial_df$id=="QS_A_001"]<-"SOKO_AUT_001"

# clean ID names: delete dates, testing times etc.
colnames(itc_df_matchedN)[2]<-"id"
itc_df_matchedN$id<-gsub("\\_T2.*","",itc_df_matchedN$id)
itc_df_matchedN$id<-toupper(itc_df_matchedN$id)
itc_df_matchedN$id<-gsub("\\_POSNER.*","",itc_df_matchedN$id)
itc_df_matchedN$id<-gsub("-", "_", itc_df_matchedN$id)
itc_df_matchedN$id<-gsub("ASS", "AUT", itc_df_matchedN$id)
itc_df_matchedN$id<-gsub("F_AUT", "AUT_F", itc_df_matchedN$id)
itc_df_matchedN$id<-gsub("F_AUT", "AUT_F", itc_df_matchedN$id)
itc_df_matchedN$id[itc_df_matchedN$id=="QS_A_001"]<-"SOKO_AUT_001"

# clean ID names: delete dates, testing times etc.
colnames(ispc_df)[2]<-"id"
ispc_df$id<-gsub("\\_T2.*","",ispc_df$id)
ispc_df$id<-toupper(ispc_df$id)
ispc_df$id<-gsub("\\_POSNER.*","",ispc_df$id)
ispc_df$id<-gsub("-", "_", ispc_df$id)
ispc_df$id<-gsub("ASS", "AUT", ispc_df$id)
ispc_df$id<-gsub("F_AUT", "AUT_F", ispc_df$id)
ispc_df$id<-gsub("F_AUT", "AUT_F", ispc_df$id)
ispc_df$id[ispc_df$id=="QS_A_001"]<-"SOKO_AUT_001"

# clean ID names: delete dates, testing times etc.
colnames(ispc_df_matchedN)[2]<-"id"
ispc_df_matchedN$id<-gsub("\\_T2.*","",ispc_df_matchedN$id)
ispc_df_matchedN$id<-toupper(ispc_df_matchedN$id)
ispc_df_matchedN$id<-gsub("\\_POSNER.*","",ispc_df_matchedN$id)
ispc_df_matchedN$id<-gsub("-", "_", ispc_df_matchedN$id)
ispc_df_matchedN$id<-gsub("ASS", "AUT", ispc_df_matchedN$id)
ispc_df_matchedN$id<-gsub("F_AUT", "AUT_F", ispc_df_matchedN$id)
ispc_df_matchedN$id<-gsub("F_AUT", "AUT_F", ispc_df_matchedN$id)
ispc_df_matchedN$id[ispc_df_matchedN$id=="QS_A_001"]<-"SOKO_AUT_001"


#### Load demographics -----------------------------------------------------------------------------------------------------------------------------
## Read & preprocess demographics data ##
setwd(".../Re_ Datenbank-Export STIPED_final")
# clean ID names: delete dates, testing times etc.
colnames(itc_df)[2]<-"id"
itc_df$id<-gsub("\\_T2.*","",itc_df$id)
itc_df$id<-toupper(itc_df$id)
itc_df$id<-gsub("\\_POSNER.*","",itc_df$id)
itc_df$id<-gsub("-", "_", itc_df$id)
itc_df$id<-gsub("ASS", "AUT", itc_df$id)
itc_df$id<-gsub("F_AUT", "AUT_F", itc_df$id)
itc_df$id<-gsub("F_AUT", "AUT_F", itc_df$id)
itc_df$id[itc_df$id=="QS_A_001"]<-"SOKO_AUT_001"

# clean ID names: delete dates, testing times etc.
colnames(trial_df)[2]<-"id"
trial_df$id<-gsub("\\_T2.*","",trial_df$id)
trial_df$id<-toupper(trial_df$id)
trial_df$id<-gsub("\\_POSNER.*","",trial_df$id)
trial_df$id<-gsub("-", "_", trial_df$id)
trial_df$id<-gsub("ASS", "AUT", trial_df$id)
trial_df$id<-gsub("F_AUT", "AUT_F", trial_df$id)
trial_df$id<-gsub("F_AUT", "AUT_F", trial_df$id)
trial_df$id[trial_df$id=="QS_A_001"]<-"SOKO_AUT_001"
ADI<-read_xlsx("ADI-R.xlsx") # ADI assessed within STIPED
ADOS2_M2<-read_xlsx("ADOS_2_M2_WPS.xlsx") 
ADOS2_M3<-read_xlsx("ADOS_2_M3_G_WPS.xlsx") 
ADOS2_M4<-read_xlsx("ADOS_2_M4_WPS.xlsx")
ADOS1_M1<-read_xlsx("ADOS_M1_G_WPS.xlsx") 
ADOS1_M2<-read_xlsx("ADOS_M2_WPS.xlsx") 
ADOS1_M3<-read_xlsx("ADOS_M3_G_WPS.xlsx") 
CBCL<-read_xlsx("CBCL_4_18.xlsx")
master_file<-read_xlsx("Checkliste.xlsx")
ADHS<-read_xlsx("FBB-HKS.xlsx")
handedness<-read_xlsx("Haendigkeit.xlsx")
HAWIK<-read_xlsx("HAWIK-III.xlsx")
SRS<-read_xlsx("SRS.xlsx")
Wais<-read_xlsx("Wais_III.xlsx")
Wie<-read_xlsx("WIE.xlsx")

## master_file - rename ids ##
names(master_file)[names(master_file)=="ID_Studie"]<-"id"
master_file$id<-toupper(master_file$id)
master_file$id<-gsub("-", "_", master_file$id)
master_file$id<-gsub("ASS", "AUT", master_file$id)
master_file$id<-gsub("F_AUT", "AUT_F", master_file$id)
master_file<-master_file[,c("ID_Bado", "id", "Geburt_Index", "Geschlecht_Index")]
names(master_file)[names(master_file)=="Geschlecht_Index"]<-"sex"


## ADI ##
ADI<-ADI[,c("ID_Bado", "ADI_alg_SOZ_INT", "ADI_alg_KOM", "ADI_alg_RITUALE", "ADI_alg_ABNORM_ENTW")]

## ADOS ##
# rename severity
# ADOS-1
names(ADOS1_M1)[names(ADOS1_M1)=="Severity_Gesamt"]<-"Severity_Gesamt_ADOS1_M1"
names(ADOS1_M2)[names(ADOS1_M2)=="Severity_Gesamt"]<-"Severity_Gesamt_ADOS1_M2"
names(ADOS1_M3)[names(ADOS1_M3)=="Severity_Gesamt"]<-"Severity_Gesamt_ADOS1_M3"
# ADOS-2
names(ADOS2_M2)[names(ADOS2_M2)=="ADOS_2_M2_comp_Score"]<-"Severity_Gesamt_ADOS2_M2"
names(ADOS2_M3)[names(ADOS2_M3)=="Severity_Gesamt"]<-"Severity_Gesamt_ADOS2_M3"
names(ADOS2_M4)[names(ADOS2_M4)=="Severity_Gesamt"]<-"Severity_Gesamt_ADOS2_M4"
# rename date and age
# ADOS-1
names(ADOS1_M1)[names(ADOS1_M1)=="DatumADOS_M1_G_WPS"]<-"date_ADOS1_M1"
names(ADOS1_M1)[names(ADOS1_M1)=="ADOS_M1_G_WPS_alter"]<-"age_ADOS1_M1"
names(ADOS1_M2)[names(ADOS1_M2)=="DatumADOS_M2"]<-"date_ADOS1_M2"
names(ADOS1_M2)[names(ADOS1_M2)=="ADOS_M2alter"]<-"age_ADOS1_M2"
names(ADOS1_M3)[names(ADOS1_M3)=="DatumADOS_M3"]<-"date_ADOS1_M3"
names(ADOS1_M3)[names(ADOS1_M3)=="ADOS_M3alter"]<-"age_ADOS1_M3"
# ADOS-2
names(ADOS2_M2)[names(ADOS2_M2)=="DatumADOS_2_M2"]<-"date_ADOS2_M2"
names(ADOS2_M2)[names(ADOS2_M2)=="ADOS_2_M2alter"]<-"age_ADOS2_M2"
names(ADOS2_M3)[names(ADOS2_M3)=="DatumADOS_2_M3"]<-"date_ADOS2_M3"
names(ADOS2_M3)[names(ADOS2_M3)=="ADOS_2_M3alter"]<-"age_ADOS2_M3"
names(ADOS2_M4)[names(ADOS2_M4)=="DatumADOS_2_M4"]<-"date_ADOS2_M4"
names(ADOS2_M4)[names(ADOS2_M4)=="ADOS_2_M4alter"]<-"age_ADOS2_M4"
# exclude unnecessary variables
ADOS1_M1<-ADOS1_M1[,c("ID_Bado", "date_ADOS1_M1", "age_ADOS1_M1", "Severity_Gesamt_ADOS1_M1", "Severity_SA", "Severity_RRB")]
ADOS1_M2<-ADOS1_M2[,c("ID_Bado", "date_ADOS1_M2", "age_ADOS1_M2", "Severity_Gesamt_ADOS1_M2", "Severity_SA", "Severity_RRB")]
ADOS1_M3<-ADOS1_M3[,c("ID_Bado", "date_ADOS1_M3", "age_ADOS1_M3", "Severity_Gesamt_ADOS1_M3", "Severity_SA", "Severity_RRB")]
ADOS2_M2<-ADOS2_M2[,c("ID_Bado", "date_ADOS2_M2", "age_ADOS2_M2", "Severity_Gesamt_ADOS2_M2", "Severity_SA", "Severity_RRB")]
ADOS2_M3<-ADOS2_M3[,c("ID_Bado", "date_ADOS2_M3", "age_ADOS2_M3", "Severity_Gesamt_ADOS2_M3", "Severity_SA", "Severity_RRB")]
ADOS2_M4<-ADOS2_M4[,c("ID_Bado", "date_ADOS2_M4", "age_ADOS2_M4", "Severity_Gesamt_ADOS2_M4", "Severity_SA", "Severity_RRB")]

# exclude first ADOS by ID_Bado 11199-2015 (participant had a diagnostic ADOS-1, and a later ADOS-2 -> take only ADOS-2)
ADOS1_M2<-ADOS1_M2[!(ADOS1_M2$ID_Bado=="11199-2015"),]

# merge
df_list <- list(ADOS1_M1, ADOS1_M2, ADOS1_M3, ADOS2_M2, ADOS2_M3, ADOS2_M4)
# merge all data frames in list
ADOS<- df_list %>% reduce(full_join, by=c('ID_Bado',"Severity_SA", "Severity_RRB"))
# dates into long format
ADOS<-as.data.frame(pivot_longer(ADOS, cols=c("date_ADOS1_M1", "date_ADOS1_M2",  "date_ADOS1_M3",  
                                              "date_ADOS2_M2", "date_ADOS2_M3", "date_ADOS2_M4"),
                                 names_to='version_module_date',values_to='date', values_drop_na = T))
# age into long format
ADOS<-as.data.frame(pivot_longer(ADOS, cols=c("age_ADOS1_M1", "age_ADOS1_M2",  "age_ADOS1_M3",  
                                              "age_ADOS2_M2", "age_ADOS2_M3", "age_ADOS2_M4"),
                                 names_to='version_module_age',values_to='age_ADOS', values_drop_na = T))
# severity into long format
ADOS<-as.data.frame(pivot_longer(ADOS, cols=c("Severity_Gesamt_ADOS1_M1", "Severity_Gesamt_ADOS1_M2",  "Severity_Gesamt_ADOS1_M3",  
                                              "Severity_Gesamt_ADOS2_M2", "Severity_Gesamt_ADOS2_M3", "Severity_Gesamt_ADOS2_M4"),
                                 names_to='version_module_severity',values_to='Severity', values_drop_na = T))


## IQ ##
# extract relevant variables
HAWIK<-HAWIK[,c("ID_Bado", "Hawik_III_IQ_P", "HAWIK_III_alter")]
Wais<-Wais[,c("ID_Bado", "Wais_III_IQ_P", "Wais_III_alter")]
Wie<-Wie[,c("ID_Bado", "WIE_III_IQ_P", "WIE_alter")]

# exclude 14983-2019 in WAIS, had HAWIK done
Wais<-Wais[!(Wais$ID_Bado=="14983-2019"),]

# merge HAWIK, Wais and Wie
df_list <- list(HAWIK, Wais, Wie)
# merge all data frames in list
IQ<-df_list %>% reduce(full_join, by='ID_Bado')
IQ<-as.data.frame(pivot_longer(IQ, cols=c("Hawik_III_IQ_P", "Wais_III_IQ_P",  "WIE_III_IQ_P"),
                               names_to='IQ_test',values_to='IQ', values_drop_na = T))
IQ<-as.data.frame(pivot_longer(IQ, cols=c("HAWIK_III_alter", "Wais_III_alter",  "WIE_alter"),
                               names_to='IQ_test_age',values_to='age', values_drop_na = T))
IQ<-IQ[,c("ID_Bado", "IQ_test", "IQ", "age")]

## CBCL ##
CBCL<-CBCL[!(CBCL$ID_Bado=="11014-2015" & (!CBCL$Meßzeit_CBCL=="T2")),]
CBCL<-CBCL[,c("ID_Bado","CBCL_T_INT", "CBCL_T_EXT", "CBCL_T_GES")]

## ADHD ##
ADHS<-ADHS[!(ADHS$ID_Bado=="11014-2015" & (!ADHS$Meßzeit_FBB_HKS=="T2")),]
ADHS<-ADHS[!(ADHS$ID_Bado=="11425-2015" & (!ADHS$Meßzeit_FBB_HKS=="T2")),]
ADHS<-ADHS[,c("ID_Bado", "ADHS_aufm", "ADHS_hyperak", "ADHS_impul","ADHS_Hinweis_Diag_ICD","ADHS_Hinweis_Diag_DSM","ADHs_ges")]

## handedness ##
handedness<-handedness[,c("ID_Bado","EHI_Score_L", "EHI_Score_R","LQ")]

## SRS ##
SRS<-SRS[!(SRS$ID_Bado=="11014-2015" & (!SRS$Meßzeit_SRS=="T2")),] #exlude time points != T2 for 11014-2015
SRS<-SRS[!(SRS$ID_Bado=="13706-2017" & (!SRS$Datum_SRS=="11.07.2019")),] #exlude time points != T2 for 11014-2015
SRS<-SRS[,c("ID_Bado","Gesamtwert_N_k_RW","Gesamt_TW_AB")]

## merge datasets into demographics file ##
df_list <- list(master_file, ADI, ADOS, IQ, SRS, ADHS, CBCL, handedness)
# merge all data frames in list
demogr<-df_list %>% reduce(full_join, by='ID_Bado')
#Note: 13404-2017 no ADOS-severity, but confirmed diagnosis by ADOS (no CSS available for module 3 for 18 year olds)

# check ASD participants with ADOS score below 4 
ADOS$ID_Bado[ADOS$Severity<4]
ADOS_excl<-demogr$id[demogr$Severity<4]

## exclude participants without EEG from demographics df ##
# check if all EEG participants are in demographics file
which(itc_df$id %in% demogr$id == FALSE)

# exclude participants from demographics file for whom there are no EEG data
demogr<-demogr[demogr$id %in% itc_df$id,]
# exclude ASD participants with too low ADOS score
demogr<-demogr[!demogr$id %in% ADOS_excl,]

# define group variable
demogr$group<-"ASD"
demogr$group[grep("HC",demogr$id)]<-"NT"

# exclude NT participants with SRS T-value >=60 and/or CBCL > 60
NT_excl<-demogr$id[demogr$group=="NT" & (demogr$CBCL_T_INT>=60 | demogr$CBCL_T_EXT>=60 | demogr$CBCL_T_INT>=60 | demogr$Gesamt_TW_AB>=60)]
demogr<-demogr[!demogr$id %in% NT_excl,]


#### Matching ---------------------------------------------------------------------------------------------------------------
require(MatchIt)
all.match<-matchit(as.factor(group)~age,
                   data=demogr,
                   method='full', caliper=0.35)

all.match<-match.data(all.match)
table(all.match$group)
# check differences
t.test(all.match$age[all.match$group=="NT"], all.match$age[all.match$group=="ASD"], paired=F, alternative = "two.sided")
# exclude non-matched participants from demographics
demogr<-demogr[demogr$id %in% all.match$id,]
dir <- "...\\input_files"
save(demogr, file = file.path(dir, "demogr"))

## Demographics: test group differences ##
# number participants
table(demogr$group)
# sex
table(demogr$group, demogr$sex)
chisq.test(table(demogr$group, demogr$sex))
# ADHD
table(demogr$ADHS_Hinweis_Diag_DSM, demogr$group)
table(demogr$ADHS_Hinweis_Diag_ICD, demogr$group)
adhd_inattention<-demogr$id[demogr$ADHS_Hinweis_Diag_DSM=="314.00"]
adhd_combined<-demogr$id[demogr$ADHS_Hinweis_Diag_DSM=="314.01"]
adhd_id<-c(adhd_inattention, adhd_combined)
adhd_id<-adhd_id[!is.na(adhd_id)]
save(adhd_id, file = file.path(dir, "adhd_id"))

# get statistics for dimensional variables
by(demogr, demogr$group, psych::describe)
# age
t.test(demogr$age[demogr$group=="ASD"], demogr$age[demogr$group=="NT"], paired=F, alternative = "two.sided")
# IQ
t.test(demogr$IQ[demogr$group=="ASD"], demogr$IQ[demogr$group=="NT"], paired=F, alternative = "two.sided")
# CBCL
t.test(demogr$CBCL_T_INT[demogr$group=="ASD"], demogr$CBCL_T_INT[demogr$group=="NT"], paired=F, alternative = "two.sided")
t.test(demogr$CBCL_T_EXT[demogr$group=="ASD"], demogr$CBCL_T_EXT[demogr$group=="NT"], paired=F, alternative = "two.sided")
t.test(demogr$CBCL_T_GES[demogr$group=="ASD"], demogr$CBCL_T_GES[demogr$group=="NT"], paired=F, alternative = "two.sided")
# handedness
t.test(demogr$LQ[demogr$group=="ASD"], demogr$LQ[demogr$group=="NT"], paired=F, alternative = "two.sided")


#### Clean EEG data frames --------------------------------------------------------------------------------------------------
## Merge demographics to EEG data ##
itc_df<-merge(itc_df, demogr, by="id")
trial_df<-merge(trial_df, demogr, by="id")
ispc_df<-merge(ispc_df, demogr, by="id")
itc_df_matchedN<-merge(itc_df_matchedN, demogr, by="id")
ispc_df_matchedN<-merge(ispc_df_matchedN, demogr, by="id")

## remove gamma activity 
#-> no clear task-evoked response
itc_df<-itc_df[!itc_df$band=="gamma_bb",]
trial_df<-trial_df[!trial_df$band=="gamma_bb",]
ispc_df<-ispc_df[!ispc_df$band=="gamma_bb",]
itc_df_matchedN<-itc_df_matchedN[!itc_df_matchedN$band=="gamma_bb",]
ispc_df_matchedN<-ispc_df_matchedN[!ispc_df_matchedN$band=="gamma_bb",]

## remove parietal ROI 
#-> very high correlation with occipital in each neural measure
itc_df<-itc_df[!itc_df$ROI=="pari",]
trial_df<-trial_df[!trial_df$ROI=="pari",]
ispc_df<-ispc_df[ispc_df$roi_pair=="occi_front",] # keep only occi-front pair
itc_df_matchedN<-itc_df_matchedN[!itc_df_matchedN$ROI=="pari",]
ispc_df_matchedN<-ispc_df_matchedN[ispc_df_matchedN$roi_pair=="occi_front",]

## remove non-cued targets
itc_df<-itc_df[!itc_df$validity=="no_cue",]
trial_df<-trial_df[!trial_df$validity=="no_cue",]
ispc_df<-ispc_df[!ispc_df$validity=="no_cue",]
itc_df_matchedN<-itc_df_matchedN[!itc_df_matchedN$validity=="no_cue",]
ispc_df_matchedN<-ispc_df_matchedN[!ispc_df_matchedN$validity=="no_cue",]

## turn character into factor variables
itc_df[sapply(itc_df, is.character)] <- lapply(itc_df[sapply(itc_df, is.character)], as.factor)
trial_df[sapply(trial_df, is.character)] <- lapply(trial_df[sapply(trial_df, is.character)], as.factor)
ispc_df[sapply(ispc_df, is.character)] <- lapply(ispc_df[sapply(ispc_df, is.character)], as.factor)
itc_df_matchedN[sapply(itc_df_matchedN, is.character)] <- lapply(itc_df_matchedN[sapply(itc_df_matchedN, is.character)], as.factor)
ispc_df_matchedN[sapply(ispc_df_matchedN, is.character)] <- lapply(ispc_df_matchedN[sapply(ispc_df_matchedN, is.character)], as.factor)

## remove cue rows if available
trial_df<-trial_df[trial_df$which=="targ",] 
itc_df<-itc_df[itc_df$which=="targ",]
ispc_df<-ispc_df[ispc_df$which=="targ",] 
itc_df_matchedN_df<-itc_df_matchedN[itc_df_matchedN$which=="targ",] 
ispc_df_matchedN<-ispc_df_matchedN[ispc_df_matchedN$which=="targ",] 

# remove fine-grained time windows to create aggregated df
itc_df_tw<-itc_df[itc_df$time_window=="101" | itc_df$time_window=="102",]
trial_df_tw<-trial_df[trial_df$time_window=="101" | trial_df$time_window=="102",]
ispc_df_tw<-ispc_df[ispc_df$time_window=="101" | ispc_df$time_window=="102",]
itc_df_matchedN_tw<-itc_df_matchedN[itc_df_matchedN$time_window=="101" | itc_df_matchedN$time_window=="102",]
ispc_df_matchedN_tw<-ispc_df_matchedN[ispc_df_matchedN$time_window=="101" | ispc_df_matchedN$time_window=="102",]


#### aggregate trial-wise dfs over unnecessary variables ---------------------------------
trial_full_df <- trial_df_tw |>
  group_by(id, validity, ROI, hemi, band, epoch_idx) |>
  summarise(
    dB_total = mean(dB_total, na.rm = TRUE),
    RT_ms = mean(RT_ms, na.rm=TRUE),
    .groups = "drop"
  )
# merge with demographics
trial_full_df<-merge(trial_full_df, demogr, by="id")

# aggregate trial-wise dfs over unnecessary variables for fine-grained tw dfs
trial_full_df_tw <- trial_df_tw |>
  group_by(id, validity, ROI, hemi, band, epoch_idx, time_window) |>
  summarise(
    dB_total = mean(dB_total, na.rm = TRUE),
    RT_ms = mean(RT_ms, na.rm=TRUE),
    .groups = "drop"
  )
# merge with demographics
trial_full_df_tw<-merge(trial_full_df_tw, demogr, by="id")


## create df with aggregated power over trials -> will be merged with ITPC
# power - target
targ_power <- trial_full_df |>
  group_by(id, validity, ROI, hemi, band) |>
  summarise(
    dB_total_mean = mean(dB_total, na.rm = TRUE),
    RT_mean       = mean(RT_ms, na.rm = TRUE),
    n_trials      = n(),
    .groups = "drop"
  )
# merge with demographics
targ_power<-merge(targ_power, demogr, by="id")

## ITPC
# compute number of trials
N_itc <- itc_df_tw %>%
  distinct(id, validity, side, SOA_label, n_trials_used) %>%   # de-duplicate repeats
  group_by(id, validity) %>%
  summarise(
    n_trials_used_itc = sum(n_trials_used),
    .groups = "drop"
  )
# aggregate
df_itc <- itc_df_tw |>
  group_by(id, validity, ROI, hemi, band) |>
  summarise(
    ITC_mean = mean(ITC, na.rm = TRUE),
    ITPCz_mean = mean(ITPCz, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(N_itc, by = c("id", "validity"))


## ITPC n_matched
# compute number of trials
N_itc_matchedN <- itc_df_matchedN_tw %>%
  distinct(id, validity, side, SOA_label, n_trials_used) %>%   # de-duplicate repeats
  group_by(id, validity) %>%
  summarise(
    n_trials_used_itc = sum(n_trials_used),
    .groups = "drop"
  )

# aggregate
df_itc_matchedN <- itc_df_matchedN_tw |>
  group_by(id, validity, ROI, hemi, band) |>
  summarise(
    ITC_mean = mean(ITPC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(N_itc_matchedN, by = c("id", "validity"))


## ISPC
# compute number of trials
N_ispc <- ispc_df_tw %>%
  distinct(id, validity, side, SOA_label, n_trials) %>%   # de-duplicate repeats
  group_by(id, validity) %>%
  summarise(
    n_trials_used_ispc = sum(n_trials),
    .groups = "drop"
  )
# aggregate
df_ispc <- ispc_df_tw |>
  group_by(id, validity, roi_pair, hemi, band) |>
  summarise(
    ISPC_mean = mean(ISPC, na.rm = TRUE),
    ISPCz_mean = mean(ISPCz, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(N_ispc, by = c("id", "validity"))

## ISPC - n matched
# compute number of trials
N_ispc_matched_N <- ispc_df_matchedN_tw %>%
  distinct(id, validity, side, SOA_label, n_trials_used) %>%   # de-duplicate repeats
  group_by(id, validity) %>%
  summarise(
    n_trials_used_ispc = sum(n_trials_used),
    .groups = "drop"
  )
# aggregate
df_ispc_matchedN <- ispc_df_matchedN_tw |>
  group_by(id, validity, roi_pair, hemi, band) |>
  summarise(
    ISPC_mean = mean(ISPC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(N_ispc_matched_N, by = c("id", "validity"))


# merge aggregated power and itpc together into dfs
targ_df_agg<-merge(targ_power, df_itc, by=c("id", "validity", "ROI", "hemi", "band"), all.x=T, all.y=T)
targ_itc_matchedN<-merge(targ_power, df_itc_matchedN, by=c("id", "validity", "ROI", "hemi", "band"), all.x=T, all.y=T)

# ISPC: merge RT and demogr to df
# extract reaction times
RT<-aggregate(RT_mean ~ validity + id, mean, data=targ_df_agg)
targ_ispc<-merge(df_ispc, RT, by=c("id", "validity"), all.x=T, all.y=T)
targ_ispc_matchedN<-merge(df_ispc_matchedN, RT, by=c("id", "validity"), all.x=T, all.y=T)
# add demographics to ISPC and matched n data frames
targ_ispc<-merge(targ_ispc, demogr, by="id")
targ_ispc_matchedN<-merge(targ_ispc_matchedN, demogr, by="id")

# turn character into factor variables
targ_df_agg[sapply(targ_df_agg, is.character)] <- lapply(targ_df_agg[sapply(targ_df_agg, is.character)], as.factor)

# drop unused levels
targ_df_agg <- droplevels(targ_df_agg)

#### compute lateralization indices based on trial-wise dfs -------------------------------------------------
targ_li_df <- trial_full_df %>%
  group_by(
    id, epoch_idx,
    ROI, band, validity
    #time_window
  ) %>%
  pivot_wider(
    names_from  = hemi,       # "ipsi", "contra" → columns
    values_from = dB_total    # value to spread
  ) %>%
  ungroup() %>%
  mutate(
    LI_ipsi_minus_contra = ipsi - contra,  
    LI_contra_minus_ipsi = contra - ipsi,  
    
  )

# with time window
targ_li_df_tw <- trial_full_df_tw %>%
  group_by(
    id, epoch_idx,
    ROI, band, validity, time_window
    #time_window
  ) %>%
  pivot_wider(
    names_from  = hemi,       # "ipsi", "contra" → columns
    values_from = dB_total    # value to spread
  ) %>%
  ungroup() %>%
  mutate(
    LI_ipsi_minus_contra = ipsi - contra,  
    LI_contra_minus_ipsi = contra - ipsi,  
    
  )

#### Compute Power Connectivity --------------------------------------------------------------------------------
# connectivity defined as long-range co-fluctuation of averaged power within 500 ms target over trials within frequency bands
# -> data doesn't show normal distribution for invalid alpha: use Spearman rank correlations

#### Calculate power connectivity ####
# Subject-level variables (one row per id)
subject_info <- trial_full_df %>%
  select(id, age, group) %>%
  distinct()

# df for correlations
df_corr <- trial_full_df %>%
  select(id, band, hemi, validity, epoch_idx, ROI, dB_total) %>%
  filter(ROI %in% c("front", "occi"))


# reshape to wide format
df_corr_wide <- df_corr %>%
  pivot_wider(
    names_from = ROI,
    values_from = dB_total
  )

# compute Spearman rank correlations 
df_powconnect <- df_corr_wide %>%
  group_by(id, band, hemi, validity) %>%
  group_modify(~ {
    
    dat <- .x %>%
      select(front, occi) %>%
      na.omit()
    
    n_trials <- nrow(dat)
    
    if(n_trials < 5) return(tibble())
    
    rho <- cor(
      dat$front,
      dat$occi,
      method = "spearman"
    )
    
    tibble(
      ROI1 = "front",
      ROI2 = "occi",
      rho  = rho,
      n_trials = n_trials
    )
  }) %>%
  ungroup() %>%
  mutate(
    rho_clip = pmin(pmax(rho, -0.999999), 0.999999), # clipped to prevent z=infinite values (happens if r = 1, the denominator becomes zero)
    pow_cor = atanh(rho_clip) # calculate fisher z to loosen correlations from -1 to 1 range
  )

# join subject-level variables
df_powconnect <- df_powconnect %>%
  left_join(subject_info, by = "id")
df_powconnect$ROI_interaction<-interaction(df_powconnect$ROI1,df_powconnect$ROI2)

# merge df with RT
df_powconnect<-merge(df_powconnect, RT, by=c("id", "validity"), all.x=T, all.y=T)
# merge with demogr
df_powconnect<-merge(df_powconnect, demogr, by=c("id", "group", "age"), all.x=, all.y=T)


#### Create df with condition-level neural levels (wide format) for neural clusters ----------------------------------------
# 0) SETTINGS
# Keep only the bands that are wanted
bands_keep <- c("alpha", "beta")

# Participant-level metadata to keep in the final wide df.
meta_vars <- c("id", "group", "age", "IQ", "sex", 
               "ADHS_aufm", "ADHS_hyperak", "ADHS_impul", "ADHs_ges", 
               "CBCL_T_INT", "CBCL_T_EXT",
               "Severity", "Severity_SA", "Severity_RRB",
               "Gesamtwert_N_k_RW", "Gesamt_TW_AB")

# 1) participant-level meta table
meta_df <- trial_full_df %>%
  select(any_of(meta_vars)) %>%
  distinct()

# 2) total power (trial-wise -> aggregate)
#    Use dB_total
#    Keep ipsi + contra separately
power_agg <- trial_full_df %>%
  filter(
    band %in% bands_keep,
    validity %in% c("valid", "invalid"),
    hemi %in% c("ipsi", "contra")
  ) %>%
  group_by(id, band, ROI, hemi, validity) %>%
  summarise(
    value = mean(dB_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = validity,
    values_from = value
  ) %>%
  mutate(
    diff_valid_minus_invalid = valid - invalid
  ) %>%
  pivot_longer(
    cols = c(valid, invalid, diff_valid_minus_invalid),
    names_to = "contrast",
    values_to = "value"
  ) %>%
  mutate(
    feature = paste("power", band, ROI, hemi, contrast, sep = "_")
  ) %>%
  select(id, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value
  )

# 3) lateralization (trial-wise -> aggregate)
#    Use LI_ipsi_minus_contra
#    No hemi split here
li_agg <- targ_li_df %>%
  filter(
    band %in% bands_keep,
    validity %in% c("valid", "invalid")
  ) %>%
  group_by(id, band, ROI, validity) %>%
  summarise(
    value = mean(LI_ipsi_minus_contra, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = validity,
    values_from = value
  ) %>%
  mutate(
    diff_valid_minus_invalid = valid - invalid
  ) %>%
  pivot_longer(
    cols = c(valid, invalid, diff_valid_minus_invalid),
    names_to = "contrast",
    values_to = "value"
  ) %>%
  mutate(
    feature = paste("LI", band, ROI, contrast, sep = "_")
  ) %>%
  select(id, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value
  )

# 4) ITPC
#    Use ITC_mean
#    Contra only
#    Already condition-level
itpc_agg <- targ_itc_matchedN %>%
  filter(
    band %in% bands_keep,
    validity %in% c("valid", "invalid"),
    hemi == "contra"
  ) %>%
  group_by(id, band, ROI, validity) %>%
  summarise(
    value = mean(ITC_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = validity,
    values_from = value
  ) %>%
  mutate(
    diff_valid_minus_invalid = valid - invalid
  ) %>%
  pivot_longer(
    cols = c(valid, invalid, diff_valid_minus_invalid),
    names_to = "contrast",
    values_to = "value"
  ) %>%
  mutate(
    feature = paste("ITPC", band, ROI, "contra", contrast, sep = "_")
  ) %>%
  select(id, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value
  )

# 5) Power connectivity
#    Use fisher_z -> pow_cor
#    Contra only
#    Already condition-level
powconn_agg <- df_powconnect %>%
  filter(
    band %in% bands_keep,
    validity %in% c("valid", "invalid"),
    hemi == "contra"
  ) %>%
  mutate(
    roi_pair = paste(pmin(ROI1, ROI2), pmax(ROI1, ROI2), sep = "_")
  ) %>%
  group_by(id, band, roi_pair, validity) %>%
  summarise(
    value = mean(pow_cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = validity,
    values_from = value
  ) %>%
  mutate(
    diff_valid_minus_invalid = valid - invalid
  ) %>%
  pivot_longer(
    cols = c(valid, invalid, diff_valid_minus_invalid),
    names_to = "contrast",
    values_to = "value"
  ) %>%
  mutate(
    feature = paste("powconn", band, roi_pair, "contra", contrast, sep = "_")
  ) %>%
  select(id, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value
  )

# 6) Phase-based connectivity (ISPC)
#    Use ISPC_mean
#    Contra only
ispc_agg <- targ_ispc_matchedN %>%
  filter(
    band %in% bands_keep,
    validity %in% c("valid", "invalid"),
    hemi == "contra",
  ) %>%
  group_by(id, band, roi_pair, validity) %>%
  summarise(
    value = mean(ISPC_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = validity,
    values_from = value
  ) %>%
  mutate(
    diff_valid_minus_invalid = valid - invalid
  ) %>%
  pivot_longer(
    cols = c(valid, invalid, diff_valid_minus_invalid),
    names_to = "contrast",
    values_to = "value"
  ) %>%
  mutate(
    feature = paste("ispc", band, roi_pair, "contra", contrast, sep = "_")
  ) %>%
  select(id, feature, value) %>%
  pivot_wider(
    names_from = feature,
    values_from = value
  )

# 7) RT aggregated
rt_agg <- trial_full_df %>%
  filter(validity %in% c("valid", "invalid")) %>%
  group_by(id, validity) %>%
  summarise(
    RT_ms = mean(RT_ms, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = validity,
    values_from = RT_ms
  ) %>%
  mutate(
    RT_diff_valid_minus_invalid = valid - invalid
  ) %>%
  rename(
    RT_valid = valid,
    RT_invalid = invalid
  )

# 8) Merge everything into wide df
shared_neural_df <- meta_df %>%
  left_join(rt_agg,      by = "id") %>%
  left_join(power_agg,   by = "id") %>%
  left_join(li_agg,      by = "id") %>%
  left_join(itpc_agg,    by = "id") %>%
  left_join(powconn_agg, by = "id") %>%
  left_join(ispc_agg,    by = "id")


# 9) create equally weighted valid/invalid means
add_condition_means <- function(df) {
  valid_vars <- grep("_valid$", names(df), value = TRUE)
  
  for (valid_var in valid_vars) {
    invalid_var <- sub("_valid$", "_invalid", valid_var)
    mean_var    <- sub("_valid$", "_mean_valid_invalid", valid_var)
    
    if (invalid_var %in% names(df)) {
      df[[mean_var]] <- (df[[valid_var]] + df[[invalid_var]]) / 2
      df[[mean_var]][!is.finite(df[[mean_var]])] <- NA_real_
    }
  }
  
  df
}

shared_neural_df <- add_condition_means(shared_neural_df)

# Quick check
dim(shared_neural_df)
names(shared_neural_df)[1:30]


## save data frames
out_folder<-".../input_files"
save(trial_full_df, file=file.path(out_folder, "trial_full_df"))
save(targ_li_df, file=file.path(out_folder, "targ_li_df"))
save(targ_df_agg, file=file.path(out_folder, "targ_df_agg"))
save(targ_ispc, file=file.path(out_folder, "targ_ispc"))
save(targ_itc_matchedN, file=file.path(out_folder, "targ_itc_matchedN"))
save(targ_ispc_matchedN, file=file.path(out_folder, "targ_ispc_matchedN"))
save(df_powconnect, file=file.path(out_folder, "df_powconnect"))
save(trial_full_df_tw, file=file.path(out_folder, "trial_full_df_tw"))
save(targ_li_df_tw, file=file.path(out_folder, "targ_li_df_tw"))
save(shared_neural_df, file=file.path(out_folder, "shared_neural_df"))


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# ANALYSES #################################################################################################################################################################

# load data
setwd("...\\input_files")
load("targ_df_agg")
load("trial_full_df")
load("targ_li_df")
load("targ_itc_matchedN")
load("targ_ispc_matchedN")
load("df_powconnect")
load("shared_neural_df")
load("demogr"))

### Trial availability ###
trials_targ<-aggregate(n_trials ~ id + group + validity, mean, na.rm=T, data=targ_df_agg)
# add percentages
max_valid<-356
max_invalid<-156
trials_targ$percentage<-NA
trials_targ$percentage[trials_targ$validity=="valid"]<-trials_targ$n_trials[trials_targ$validity=="valid"]/max_valid*100
trials_targ$percentage[trials_targ$validity=="invalid"]<-trials_targ$n_trials[trials_targ$validity=="invalid"]/max_invalid*100
by(trials_targ, interaction(trials_targ$group, trials_targ$validity), psych::describe)
t.test(trials_targ$n_trials[trials_targ$group=="ASD" & trials_targ$validity=="valid"], trials_targ$n_trials[trials_targ$group=="NT" & trials_targ$validity=="valid"])
t.test(trials_targ$n_trials[trials_targ$group=="ASD" & trials_targ$validity=="invalid"], trials_targ$n_trials[trials_targ$group=="NT" & trials_targ$validity=="invalid"])


# 1 BEHAVIORAL ANALYSES ################################################################################################################################################################

#### Calculate number of excluded trials due to incorrect behavioral response ------------------
beh_dir <- ".../BehaviorClean"

beh_files <- list.files(beh_dir,
                        pattern = "_target_behavior_kept\\.csv$",
                        full.names = TRUE)

beh_df <- map_dfr(beh_files, function(f) {
  read_delim(f, delim = ";", show_col_types = FALSE) %>%
    mutate(
      file_id = f,
      id = basename(f),
      id = gsub("_target_behavior_kept\\.csv$", "", id)
    )
})

beh_df$id <- gsub("\\_T2.*", "", beh_df$id)
beh_df$id <- toupper(beh_df$id)
beh_df$id <- gsub("\\_POSNER.*", "", beh_df$id)
beh_df$id <- gsub("-", "_", beh_df$id)
beh_df$id <- gsub("ASS", "AUT", beh_df$id)
beh_df$id <- gsub("F_AUT", "AUT_F", beh_df$id)
beh_df$id <- gsub("F_AUT", "AUT_F", beh_df$id)
beh_df$id[beh_df$id == "QS_A_001"] <- "SOKO_AUT_001"

beh_df <- beh_df %>%
  select(-any_of("group")) %>%
  left_join(demogr %>% select(id, group), by = "id")

beh_df <- beh_df %>%
  filter(id %in% demogr$id)

beh_df <- beh_df %>%
  mutate(
    category = case_when(
      kept == TRUE ~ "kept",
      str_detect(drop_reason, "wrong_side") ~ "wrong_side",
      str_detect(drop_reason, "miss") ~ "miss",
      str_detect(drop_reason, "rt<") ~ "too_fast",
      str_detect(drop_reason, "rt>") ~ "too_slow",
      TRUE ~ "other"
    )
  )

group_props <- beh_df %>%
  count(group, category) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

group_props_wide <- group_props %>%
  select(group, category, prop) %>%
  pivot_wider(names_from = category, values_from = prop)

group_props_wide

## overview of reaction time characteristics per group
bhv_agg <- trial_full_df %>%
  select(id, group, validity, RT_ms) %>%
  distinct() %>%
  mutate(
    group = factor(group, levels = c("ASD", "NT")),
    validity = factor(validity, levels = c("invalid", "valid")),
    id = factor(id)
  )

by(bhv_agg, interaction(bhv_agg$group, bhv_agg$validity), function(x) {psych::describe(x$RT_ms)})


#### Group differences: Reaction time -----------------------------------------
## compute trial level RT 
bhv_df <- trial_full_df %>%
  select(id, group, validity, epoch_idx, RT_ms, age, IQ, sex) %>%
  distinct() %>%
  mutate(
    group = factor(group, levels = c("ASD", "NT")),
    validity = factor(validity, levels = c("invalid", "valid")),
    sex = factor(sex),
    id = factor(id)
  )
# plot
ggplot(bhv_df, aes(x=RT_ms, fill=group)) + geom_density(alpha=0.2) + facet_wrap(~validity)

## model:
m1_bhv<-lmer(scale(RT_ms) ~ group*validity + scale(age) + scale(IQ) + sex + (1|id), data=bhv_df)
anova(m1_bhv)

## Post-hoc  behavior
m1_bhv_emm <- emmeans(m1_bhv, ~ group * validity)

bhv_contrasts_df <- contrast(
  m1_bhv_emm, method = "consec",
  simple = "each", combine = TRUE, adjust = "none"
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as.data.frame() %>%
  mutate(p_FDR = p.adjust(p.value, method = "BH"))

bhv_contrasts_df

# contrast of contrasts
bhv_emm_all <- emmeans(m1_bhv, ~ group * validity)
interaction_contrast <- contrast(
  bhv_emm_all,
  interaction = "pairwise",
  adjust = "none"
)
summary(interaction_contrast, infer = c(TRUE, TRUE))

# save it as table
fmt_p <- function(x) {
  ifelse(x < .001, "< .001", sub("^0", "", sprintf("%.3f", x)))
}
## Omnibus table
bhv_omnibus <- anova(
  m1_bhv,
  type = 3,
  ddf = "Satterthwaite"
) %>%
  as.data.frame() %>%
  rownames_to_column("effect") %>%
  transmute(
    Effect = recode(
      effect,
      group = "Group",
      validity = "Validity",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      `F value`, NumDF, DenDF
    ),
    p = fmt_p(`Pr(>F)`)
  )
print(as.data.frame(bhv_omnibus), row.names = FALSE, right = TRUE)
## Export to Word
output_file_bhv <- paste0(
  ".../Tables and Figures/",
  "Table_Behavioral_RT_Omnibus_Tests.docx"
)
save_as_docx(
  "Behavioral RT omnibus tests" =
    flextable(bhv_omnibus) %>% autofit(),
  path = output_file_bhv
)


# 2 EEG ANALYSES: Group differences #################################################################################################################################################################

#### 2.0 SIMULATION-BASED POWER ANALYSES -----------------------------------------------------------------------------------------------------------------------------------------

library(simr)
# script runs with:
# lme4 1.1-37
# simr 1.0.8

set.seed(1234)

## SETTINGS
cat("\nPackage versions:\n")
cat("R:        ", R.version.string, "\n")
cat("lme4:     ", as.character(packageVersion("lme4")), "\n")
cat("lmerTest: ", as.character(packageVersion("lmerTest")), "\n")
cat("simr:     ", as.character(packageVersion("simr")), "\n\n")

if(as.character(packageVersion("simr")) != "1.0.8") {
  stop(
    "This script was developed with simr 1.0.8. ",
    "Current version: ", as.character(packageVersion("simr")),
    ". Check package versions before starting the overnight run."
  )
}

# Match the actual analyses:
# Type III ANOVA + Satterthwaite denominator df
simrOptions(
  lmerTestType = 3,
  lmerTestDdf = "Satterthwaite"
)

# Significance level and desired power
alpha <- .05
target_power <- .80

# Outcomes are placed on a standardized scale:
# total random + residual variance = 1.
# ICC=.50 is the central planning assumption.
# ICC=.30 and .70 test robustness to lower/higher within-subject correlation.

total_sd <- 1
icc_values <- c(.30, .50, .70)

# Trial-level power:
# finer/lower range because large number of repeated observations should make much smaller effects detectable.
effect_grid_subject <- seq(0, 1.20, by = .10)
effect_grid_trial   <- seq(0, 0.50, by = .05)

# Fine resolution around the 80%-power threshold
fine_step_subject <- .025
fine_step_trial   <- .01

# Number of simulations
nsim_coarse_subject <- 1000  
nsim_fine_subject   <- 2000

nsim_coarse_trial <- 500    
nsim_fine_trial   <- 1000

# Output directory
out_dir <- file.path(
  getwd(),
  paste0("power_analysis_final_", Sys.Date())
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("Results will be saved to:\n", out_dir, "\n\n")

# Save software environment immediately
capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

package_versions <- data.frame(
  package = c("R", "lme4", "lmerTest", "simr"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("lme4")),
    as.character(packageVersion("lmerTest")),
    as.character(packageVersion("simr"))
  )
)

write.csv(
  package_versions,
  file.path(out_dir, "package_versions.csv"),
  row.names = FALSE
)

#### Prepare ITPC design
# Create itpd_df
itpc_df <- targ_itc_matchedN %>%
  filter(
    ROI != "pari",
    band %in% c("alpha", "beta"),
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(
      ROI,
      levels = c("front", "occi")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()
stopifnot(
  nlevels(itpc_df$group) == 2,
  nlevels(itpc_df$validity) == 2,
  nlevels(itpc_df$ROI) == 2
)


itpc_design <- itpc_df %>%
  filter(
    band == "alpha"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(ROI),
    sex = factor(sex),
    id = factor(id),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ))
  ) %>%
  select(
    id,
    group,
    validity,
    ROI,
    age_z,
    IQ_z,
    sex
  ) %>%
  na.omit() %>%
  droplevels()

stopifnot(
  nlevels(itpc_design$group) == 2,
  nlevels(itpc_design$validity) == 2,
  nlevels(itpc_design$ROI) == 2
)


#### Peprare trial-level power design
power_design <- trial_full_df %>%
  filter(
    band == "alpha",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(ROI),
    sex = factor(sex),
    id = factor(id),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ))
  ) %>%
  select(
    id,
    group,
    validity,
    ROI,
    age_z,
    IQ_z,
    sex
  ) %>%
  na.omit() %>%
  droplevels()

stopifnot(
  nlevels(power_design$group) == 2,
  nlevels(power_design$validity) == 2,
  nlevels(power_design$ROI) == 2
)

#### Prepare connectivity design
# One representative phase-connectivity model:
# alpha, occipito-frontal, contralateral
# Power-connectivity has the same model structure, so run only one representative connectivity power analysis
phasecon_design <- targ_ispc_matchedN %>%
  filter(
    band == "alpha",
    roi_pair == "occi_front",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    sex = factor(sex),
    id = factor(id),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ))
  ) %>%
  select(
    id,
    group,
    validity,
    age_z,
    IQ_z,
    sex
  ) %>%
  na.omit() %>%
  droplevels()

stopifnot(
  nlevels(phasecon_design$group) == 2,
  nlevels(phasecon_design$validity) == 2
)


#### Design checks and summary
summarise_design <- function(data, model_name) {
  
  n_per_id <- data %>%
    count(id, name = "n_rows")
  
  group_counts <- data %>%
    distinct(id, group) %>%
    count(group) %>%
    summarise(
      groups = paste0(group, "=", n, collapse = "; ")
    ) %>%
    pull(groups)
  
  data.frame(
    model = model_name,
    participants = n_distinct(data$id),
    rows = nrow(data),
    group_counts = group_counts,
    median_rows_per_id = median(n_per_id$n_rows),
    min_rows_per_id = min(n_per_id$n_rows),
    max_rows_per_id = max(n_per_id$n_rows)
  )
}

design_summary <- bind_rows(
  summarise_design(
    itpc_design,
    "ITPC"
  ),
  summarise_design(
    power_design,
    "Power"
  ),
  summarise_design(
    phasecon_design,
    "Connectivity"
  )
)

print(design_summary)

write.csv(
  design_summary,
  file.path(out_dir, "design_summary.csv"),
  row.names = FALSE
)

# Additional level checks
cat("\nITPC levels:\n")
print(table(itpc_design$group))
print(table(itpc_design$validity))
print(table(itpc_design$ROI))

cat("\nPower levels:\n")
print(table(power_design$group))
print(table(power_design$validity))
print(table(power_design$ROI))

cat("\nConnectivity levels:\n")
print(table(phasecon_design$group))
print(table(phasecon_design$validity))


#### Helper
## Identify interaction coefficient
find_interaction_coef <- function(model, variables) {
  
  nms <- names(fixef(model))
  parts <- strsplit(nms, ":")
  
  keep <- vapply(
    parts,
    function(x) {
      
      length(x) == length(variables) &&
        all(
          vapply(
            variables,
            function(v) any(startsWith(x, v)),
            logical(1)
          )
        )
    },
    logical(1)
  )
  
  result <- nms[keep]
  
  if(length(result) != 1) {
    stop(
      "Expected exactly one coefficient for ",
      paste(variables, collapse = " x "),
      ", but found: ",
      paste(result, collapse = ", ")
    )
  }
  
  result
}


## Helper: Build model template
build_template <- function(design, formula, icc) {
  
  # Standardized outcome variance:
  # total variance = 1
  rand_var <- icc * total_sd^2
  resid_sd <- sqrt(
    (1 - icc) * total_sd^2
  )
  

  fixed_formula <- delete.response(
    terms(
      lme4::nobars(formula)
    )
  )
  
  X <- model.matrix(
    fixed_formula,
    data = design
  )
  
  template <- makeLmer(
    formula,
    fixef = rep(0, ncol(X)),
    VarCorr = list(
      id = rand_var
    ),
    sigma = resid_sd,
    data = design
  )
  
  template
}


## Helper: Run arbitrary ICC × effect size conditions
# - runs one powerSim for every requested ICC/effect-size combination;
# - saves the result after every completed condition;
# - can resume from an existing checkpoint file;
# - records Monte Carlo confidence intervals and simulation diagnostics.
run_power_conditions <- function(
    design,
    formula,
    coef_variables,
    test_term,
    conditions,
    nsim,
    model_label,
    interaction_label,
    stage,
    checkpoint_file,
    seed_base
) {
  
  # Load existing checkpoint if present
  if(file.exists(checkpoint_file)) {
    
    results <- readRDS(checkpoint_file)
    
    message(
      "Resuming from existing checkpoint: ",
      checkpoint_file
    )
    
  } else {
    
    results <- data.frame()
    
  }
  
  # Ensure deterministic order
  conditions <- conditions %>%
    arrange(ICC, effect_size)
  
  for(j in seq_len(nrow(conditions))) {
    
    this_icc <- conditions$ICC[j]
    beta <- conditions$effect_size[j]
    
    # Skip conditions already completed
    already_done <- FALSE
    
    if(nrow(results) > 0) {
      
      already_done <- any(
        abs(results$ICC - this_icc) < 1e-12 &
          abs(results$effect_size - beta) < 1e-12 &
          results$nsim == nsim
      )
    }
    
    if(already_done) {
      
      message(
        "Skipping completed: ",
        model_label,
        " | ",
        interaction_label,
        " | ICC=",
        this_icc,
        " | effect=",
        beta
      )
      
      next
    }
    
    message(
      "\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | ",
      model_label,
      " | ",
      interaction_label,
      " | ",
      stage,
      " | ICC=",
      this_icc,
      " | effect=",
      beta,
      " | nsim=",
      nsim
    )
    
    template <- build_template(
      design = design,
      formula = formula,
      icc = this_icc
    )
    
    coef_name <- find_interaction_coef(
      template,
      coef_variables
    )
    
    # Isolate the interaction of interest:
    m <- template
    
    fixef(m)[] <- 0
    fixef(m)[coef_name] <- beta
    
    # Stable, reproducible seed for this condition
    seed_this <- seed_base +
      round(this_icc * 10000) +
      round(beta * 100000)
    
    start_time <- Sys.time()
    
    ps <- powerSim(
      m,
      test = fixed(
        test_term,
        "anova"
      ),
      nsim = nsim,
      seed = seed_this,
      alpha = alpha,
      progress = FALSE
    )
    
    end_time <- Sys.time()
    
    ci <- binom.test(
      ps$x,
      ps$n,
      conf.level = .95
    )$conf.int
    
    new_row <- data.frame(
      model = model_label,
      interaction = interaction_label,
      stage = stage,
      ICC = this_icc,
      random_SD = sqrt(this_icc),
      residual_SD = sqrt(1 - this_icc),
      effect_size = beta,
      power = ps$x / ps$n,
      CI_low = ci[1],
      CI_high = ci[2],
      significant = ps$x,
      nsim = ps$n,
      alpha = alpha,
      warnings = NROW(ps$warnings),
      errors = NROW(ps$errors),
      elapsed_minutes = as.numeric(
        difftime(
          end_time,
          start_time,
          units = "mins"
        )
      )
    )
    
    results <- bind_rows(
      results,
      new_row
    )
    
    # Sort before saving
    results <- results %>%
      arrange(
        ICC,
        effect_size
      )
    
    # Checkpoint after every completed condition
    saveRDS(
      results,
      checkpoint_file
    )
    
    write.csv(
      results,
      sub(
        "\\.rds$",
        ".csv",
        checkpoint_file
      ),
      row.names = FALSE
    )
    
    # Clear temporary objects periodically
    rm(template, m, ps)
    
    if(j %% 3 == 0) {
      invisible(gc())
    }
  }
  
  results
}


## Helper: create coarse conditions
make_conditions <- function(
    icc_values,
    effect_sizes
) {
  
  expand.grid(
    ICC = icc_values,
    effect_size = effect_sizes
  ) %>%
    arrange(
      ICC,
      effect_size
    )
}


## Coarse conditions
coarse_subject_conditions <- make_conditions(
  icc_values = icc_values,
  effect_sizes = effect_grid_subject
)

coarse_trial_conditions <- make_conditions(
  icc_values = icc_values,
  effect_sizes = effect_grid_trial
)

##### Coarse power analyses ####
## ITPC two-way
coarse_itpc_2way <- run_power_conditions(
  design = itpc_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    coarse_subject_conditions,
  nsim =
    nsim_coarse_subject,
  model_label =
    "ITPC",
  interaction_label =
    "group x validity",
  stage =
    "coarse",
  checkpoint_file =
    file.path(
      out_dir,
      "coarse_ITPC_2way.rds"
    ),
  seed_base =
    100000
)

## ITPC three-way
coarse_itpc_3way <- run_power_conditions(
  design = itpc_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity", "ROI"),
  test_term =
    "group:validity:ROI",
  conditions =
    coarse_subject_conditions,
  nsim =
    nsim_coarse_subject,
  model_label =
    "ITPC",
  interaction_label =
    "group x validity x ROI",
  stage =
    "coarse",
  checkpoint_file =
    file.path(
      out_dir,
      "coarse_ITPC_3way.rds"
    ),
  seed_base =
    200000
)


## Connectivity
coarse_connectivity <- run_power_conditions(
  design = phasecon_design,
  formula =
    y ~ group * validity +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    coarse_subject_conditions,
  nsim =
    nsim_coarse_subject,
  model_label =
    "Connectivity",
  interaction_label =
    "group x validity",
  stage =
    "coarse",
  checkpoint_file =
    file.path(
      out_dir,
      "coarse_connectivity.rds"
    ),
  seed_base =
    300000
)


## Trial-wise power two way
coarse_power_2way <- run_power_conditions(
  design = power_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    coarse_trial_conditions,
  nsim =
    nsim_coarse_trial,
  model_label =
    "Power",
  interaction_label =
    "group x validity",
  stage =
    "coarse",
  checkpoint_file =
    file.path(
      out_dir,
      "coarse_power_2way.rds"
    ),
  seed_base =
    400000
)


## Trial-wise power-three way
coarse_power_3way <- run_power_conditions(
  design = power_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity", "ROI"),
  test_term =
    "group:validity:ROI",
  conditions =
    coarse_trial_conditions,
  nsim =
    nsim_coarse_trial,
  model_label =
    "Power",
  interaction_label =
    "group x validity x ROI",
  stage =
    "coarse",
  checkpoint_file =
    file.path(
      out_dir,
      "coarse_power_3way.rds"
    ),
  seed_base =
    500000
)


## Combine coarse results
coarse_results <- bind_rows(
  coarse_itpc_2way,
  coarse_itpc_3way,
  coarse_connectivity,
  coarse_power_2way,
  coarse_power_3way
) %>%
  arrange(
    model,
    interaction,
    ICC,
    effect_size
  )

write.csv(
  coarse_results,
  file.path(
    out_dir,
    "ALL_coarse_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  coarse_results,
  file.path(
    out_dir,
    "ALL_coarse_results.rds"
  )
)


#### Fine grid power analyses ----------------------------------------- 
## Helper: Build fine grid around 80% threshold 
# identifies the effect-size region where the simulated
# power curve crosses 80%, then samples that region more densely

make_refinement_conditions <- function(
    coarse_df,
    fine_step,
    max_effect
) {
  
  refined <- list()
  k <- 1
  
  for(this_icc in sort(unique(coarse_df$ICC))) {
    
    d <- coarse_df %>%
      filter(
        ICC == this_icc,
        effect_size > 0
      ) %>%
      arrange(effect_size)
    
    above <- which(
      d$power >= target_power
    )
    
    if(length(above) == 0) {
      # No 80% threshold within the coarse range.
      # Refine the upper end so the result is still informative.
      upper <- max(d$effect_size)
      lower <- max(
        0,
        upper - 2 * median(diff(d$effect_size))
      )
      
    } else {
      first_above <- above[1]
      upper <- d$effect_size[first_above]

      if(first_above == 1) {
        lower <- 0
        
      } else {
        lower <- d$effect_size[first_above - 1]
        
      }
      
      # Add a small buffer around the identified interval
      coarse_step <- median(
        diff(d$effect_size)
      )
      lower <- max(
        0,
        lower - coarse_step
      )
      upper <- min(
        max_effect,
        upper + coarse_step
      )
    }
    effects <- seq(
      lower,
      upper,
      by = fine_step
    )
    effects <- unique(
      round(effects, 4)
    )
    refined[[k]] <- data.frame(
      ICC = this_icc,
      effect_size = effects
    )
    
    k <- k + 1
  }
  
  bind_rows(refined) %>%
    arrange(
      ICC,
      effect_size
    )
}


#### Create automatic fine grids
fine_conditions_itpc_2way <- make_refinement_conditions(
  coarse_itpc_2way,
  fine_step = fine_step_subject,
  max_effect = max(effect_grid_subject)
)

fine_conditions_itpc_3way <- make_refinement_conditions(
  coarse_itpc_3way,
  fine_step = fine_step_subject,
  max_effect = max(effect_grid_subject)
)

fine_conditions_connectivity <- make_refinement_conditions(
  coarse_connectivity,
  fine_step = fine_step_subject,
  max_effect = max(effect_grid_subject)
)

fine_conditions_power_2way <- make_refinement_conditions(
  coarse_power_2way,
  fine_step = fine_step_trial,
  max_effect = max(effect_grid_trial)
)

fine_conditions_power_3way <- make_refinement_conditions(
  coarse_power_3way,
  fine_step = fine_step_trial,
  max_effect = max(effect_grid_trial)
)

# Save automatically selected refinement grids
write.csv(
  fine_conditions_itpc_2way,
  file.path(
    out_dir,
    "fine_grid_ITPC_2way.csv"
  ),
  row.names = FALSE
)

write.csv(
  fine_conditions_itpc_3way,
  file.path(
    out_dir,
    "fine_grid_ITPC_3way.csv"
  ),
  row.names = FALSE
)

write.csv(
  fine_conditions_connectivity,
  file.path(
    out_dir,
    "fine_grid_connectivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  fine_conditions_power_2way,
  file.path(
    out_dir,
    "fine_grid_power_2way.csv"
  ),
  row.names = FALSE
)

write.csv(
  fine_conditions_power_3way,
  file.path(
    out_dir,
    "fine_grid_power_3way.csv"
  ),
  row.names = FALSE
)


## Refined ITPC two-way
fine_itpc_2way <- run_power_conditions(
  design = itpc_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    fine_conditions_itpc_2way,
  nsim =
    nsim_fine_subject,
  model_label =
    "ITPC",
  interaction_label =
    "group x validity",
  stage =
    "refined",
  checkpoint_file =
    file.path(
      out_dir,
      "refined_ITPC_2way.rds"
    ),
  seed_base =
    600000
)


## Refined ITPC three-way
fine_itpc_3way <- run_power_conditions(
  design = itpc_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity", "ROI"),
  test_term =
    "group:validity:ROI",
  conditions =
    fine_conditions_itpc_3way,
  nsim =
    nsim_fine_subject,
  model_label =
    "ITPC",
  interaction_label =
    "group x validity x ROI",
  stage =
    "refined",
  checkpoint_file =
    file.path(
      out_dir,
      "refined_ITPC_3way.rds"
    ),
  seed_base =
    700000
)


## Refined connectivity
fine_connectivity <- run_power_conditions(
  design = phasecon_design,
  formula =
    y ~ group * validity +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    fine_conditions_connectivity,
  nsim =
    nsim_fine_subject,
  model_label =
    "Connectivity",
  interaction_label =
    "group x validity",
  stage =
    "refined",
  checkpoint_file =
    file.path(
      out_dir,
      "refined_connectivity.rds"
    ),
  seed_base =
    800000
)


## Refined trial-wise power two-way
fine_power_2way <- run_power_conditions(
  design = power_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity"),
  test_term =
    "group:validity",
  conditions =
    fine_conditions_power_2way,
  nsim =
    nsim_fine_trial,
  model_label =
    "Power",
  interaction_label =
    "group x validity",
  stage =
    "refined",
  checkpoint_file =
    file.path(
      out_dir,
      "refined_power_2way.rds"
    ),
  seed_base =
    900000
)


## Refined trial-wise power three-way
fine_power_3way <- run_power_conditions(
  design = power_design,
  formula =
    y ~ group * validity * ROI +
    age_z + IQ_z + sex +
    (1 | id),
  coef_variables =
    c("group", "validity", "ROI"),
  test_term =
    "group:validity:ROI",
  conditions =
    fine_conditions_power_3way,
  nsim =
    nsim_fine_trial,
  model_label =
    "Power",
  interaction_label =
    "group x validity x ROI",
  stage =
    "refined",
  checkpoint_file =
    file.path(
      out_dir,
      "refined_power_3way.rds"
    ),
  seed_base =
    1000000
)


## combine results
all_results <- bind_rows(
  coarse_results,
  fine_itpc_2way,
  fine_itpc_3way,
  fine_connectivity,
  fine_power_2way,
  fine_power_3way
) %>%
  arrange(
    model,
    interaction,
    ICC,
    effect_size,
    nsim
  )

write.csv(
  all_results,
  file.path(
    out_dir,
    "FINAL_all_simulations.csv"
  ),
  row.names = FALSE
)

saveRDS(
  all_results,
  file.path(
    out_dir,
    "FINAL_all_simulations.rds"
  )
)


## Keep best estimate for each effect-size point
# Some effect sizes occur in both coarse and refined analyses
# For those points, retain the estimate with the larger number of simulations.
final_results <- all_results %>%
  group_by(
    model,
    interaction,
    ICC,
    effect_size
  ) %>%
  slice_max(
    nsim,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    model,
    interaction,
    ICC,
    effect_size
  )

write.csv(
  final_results,
  file.path(
    out_dir,
    "FINAL_power_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  final_results,
  file.path(
    out_dir,
    "FINAL_power_results.rds"
  )
)


## Estimate 80% power thresholds
threshold_summary <- final_results %>%
  filter(
    effect_size > 0
  ) %>%
  group_by(
    model,
    interaction,
    ICC
  ) %>%
  summarise(
    effect_size_80 = if(
      any(power >= target_power)
    ) {
      min(
        effect_size[
          power >= target_power
        ]
      )
    } else {
      NA_real_
    },
    maximum_tested_effect =
      max(effect_size),
    maximum_observed_power =
      max(power),
    .groups = "drop"
  )

threshold_details <- threshold_summary %>%
  left_join(
    final_results %>%
      select(
        model,
        interaction,
        ICC,
        effect_size,
        power,
        CI_low,
        CI_high,
        nsim
      ),
    by = c(
      "model",
      "interaction",
      "ICC",
      "effect_size_80" = "effect_size"
    )
  ) %>%
  arrange(
    model,
    interaction,
    ICC
  )

print(threshold_details)

write.csv(
  threshold_details,
  file.path(
    out_dir,
    "FINAL_80percent_thresholds.csv"
  ),
  row.names = FALSE
)

#### Primary ICC=.50 summary 
primary_thresholds <- threshold_details %>%
  filter(
    ICC == .50
  )

print(primary_thresholds)

write.csv(
  primary_thresholds,
  file.path(
    out_dir,
    "FINAL_primary_ICC50_thresholds.csv"
  ),
  row.names = FALSE
)


## Benchmark table 
benchmark_effects <- c(
  .10, .20, .30, .40, .50,
  .60, .70, .80, .90, 1.00
)

benchmark_results <- final_results %>%
  filter(
    effect_size %in% benchmark_effects
  ) %>%
  mutate(
    power_percent = round(
      power * 100,
      1
    ),
    CI = paste0(
      round(CI_low * 100, 1),
      "-",
      round(CI_high * 100, 1),
      "%"
    )
  ) %>%
  select(
    model,
    interaction,
    ICC,
    effect_size,
    power,
    power_percent,
    CI,
    nsim
  ) %>%
  arrange(
    model,
    interaction,
    ICC,
    effect_size
  )

write.csv(
  benchmark_results,
  file.path(
    out_dir,
    "FINAL_benchmark_power_table.csv"
  ),
  row.names = FALSE
)

## Sanity check
null_results <- final_results %>%
  filter(
    effect_size == 0
  ) %>%
  select(
    model,
    interaction,
    ICC,
    power,
    CI_low,
    CI_high,
    nsim,
    warnings,
    errors
  ) %>%
  arrange(
    model,
    interaction,
    ICC
  )

print(null_results)

write.csv(
  null_results,
  file.path(
    out_dir,
    "FINAL_null_effect_sanity_check.csv"
  ),
  row.names = FALSE
)


## Check warnings and errors
simulation_problems <- all_results %>%
  filter(
    warnings > 0 |
      errors > 0
  ) %>%
  arrange(
    model,
    interaction,
    ICC,
    effect_size
  )

print(simulation_problems)

write.csv(
  simulation_problems,
  file.path(
    out_dir,
    "FINAL_simulation_warnings_errors.csv"
  ),
  row.names = FALSE
)

## Save workspace objects needed for reporting
save(
  design_summary,
  coarse_results,
  all_results,
  final_results,
  threshold_details,
  primary_thresholds,
  benchmark_results,
  null_results,
  simulation_problems,
  file = file.path(
    out_dir,
    "FINAL_power_analysis_objects.RData"
  )
)

## Console summary
cat("\n\n")
cat("============================================================\n")
cat("FINAL POWER ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Primary ICC = .50 thresholds:\n\n")
print(
  primary_thresholds %>%
    select(
      model,
      interaction,
      effect_size_80,
      power,
      CI_low,
      CI_high,
      nsim
    )
)

cat("\nNull-effect sanity check:\n\n")
print(null_results)

cat("\nSimulation conditions with warnings/errors:\n\n")
print(simulation_problems)

cat("\n============================================================\n")
cat("END\n")
cat("============================================================\n")




#### 2.1 POWER -----------------------------------------------------------------------------------------------------------------------------------------
# Total power, dB-corrected 
# FDR correction across all four models (ipsi/contra, alpha/beta) 
# correction within effect families

## Prepare data
power_df <- trial_full_df %>%
  filter(
    band %in% c("alpha", "beta"),
    hemi %in% c("contra", "ipsi")
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI  = factor(ROI),
    sex  = factor(sex),
    id   = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    hemi = factor(
      hemi,
      levels = c("contra", "ipsi")
    )
  ) %>%
  droplevels()
# check levels
levels(power_df$group)
levels(power_df$validity)
levels(power_df$ROI)

###### Models ####
# Alpha contralateral
m_alpha_contra <- lmer(scale(dB_total) ~ group * validity * ROI +
                         scale(age) + scale(IQ) + sex + (1 | id), data = power_df %>% filter(band == "alpha", hemi == "contra"))
# Alpha ipsilateral
m_alpha_ipsi <- lmer(scale(dB_total) ~ group * validity * ROI +
                       scale(age) + scale(IQ) + sex + (1 | id), data = power_df %>% filter(band == "alpha", hemi == "ipsi"))
# Beta contralateral
m_beta_contra <- lmer(scale(dB_total) ~ group * validity * ROI +
                        scale(age) + scale(IQ) + sex + (1 | id), data = power_df %>% filter(band == "beta", hemi == "contra"))
# Beta ipsilateral
m_beta_ipsi <- lmer(scale(dB_total) ~ group * validity * ROI +
                      scale(age) + scale(IQ) + sex + (1 | id), data = power_df %>% filter(band == "beta", hemi == "ipsi"))

## Define & extract omnibus tests 
effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)
# Function to extract Type III omnibus tests 
extract_power_tests <- function(models) {
  
  imap_dfr(
    models,
    function(model, hemisphere_name) {
      
      anova(
        model,
        type = 3,
        ddf  = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          hemisphere = hemisphere_name,
          effect      = effect,
          df_num      = NumDF,
          df_den      = DenDF,
          F_value     = `F value`,
          p_raw       = `Pr(>F)`
        ) %>%
        filter(
          effect %in% effects_of_interest
        )
    }
  )
}
# Extract alpha omnibus tests
alpha_power_models <- list(
  contra = m_alpha_contra,
  ipsi   = m_alpha_ipsi
)
alpha_power_tests <- extract_power_tests(
  alpha_power_models
)
alpha_power_tests
# Extract beta omnibus tests
beta_power_models <- list(
  contra = m_beta_contra,
  ipsi   = m_beta_ipsi
)

beta_power_tests <- extract_power_tests(
  beta_power_models
)
beta_power_tests

## Run omnibus test
# merge tests
all_power_tests <- bind_rows(
  alpha_power_tests %>%
    mutate(
      band = "alpha",
      .before = hemisphere
    ),
  beta_power_tests %>%
    mutate(
      band = "beta",
      .before = hemisphere
    )
)
# For each omnibus effect, apply BH-FDR across all four models
all_power_tests_FDR <- all_power_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    hemisphere = factor(
      hemisphere,
      levels = c("contra", "ipsi")
    ),
    effect = factor(
      effect,
      levels = effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    effect,
    band,
    hemisphere
  )

all_power_tests_FDR
# Format p-values
fmt_p <- function(x) {
  ifelse(x < .001, "< .001", sub("^0", "", sprintf("%.3f", x)))
}
# Arrange rows by model, then by omnibus effect
power_table <- all_power_tests_FDR %>%
  arrange(
    match(
      paste(band, hemisphere),
      c(
        "alpha ipsi",
        "alpha contra",
        "beta ipsi",
        "beta contra"
      )
    ),
    effect
  ) %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Hemi = recode(
      as.character(hemisphere),
      ipsi   = "Ipsi",
      contra = "Contra"
    ),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value, df_num, df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(as.data.frame(power_table), row.names = FALSE, right = TRUE)

# save as word table
output_file <- paste0(
  ".../Tables and Figures/",
  "Table_Power_Omnibus_Tests_FDR.docx"
)
ft <- flextable(power_table) %>% autofit()
save_as_docx(
  "Power omnibus tests" = ft,
  path = output_file
)


###### Alpha power post-hoc ####
## Ipsilateral alpha power: ROI × validity
# Estimated marginal means
em_alpha_pow_valid_ROI <- emmeans(
  m_alpha_ipsi,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)
# ROI differences within each validity condition
contr_alpha_pow_valid_ROI <- contrast(
  em_alpha_pow_valid_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs + FDR-adjusted p-values
df_alpha_pow_valid_ROI <- summary(
  contr_alpha_pow_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_alpha_pow_valid_ROI

## Contrast of contrasts: ROI × validity
diffs <- contrast(
  em_alpha_pow_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Difference between the two validity contrasts
interaction_contrasts <- contrast(
  diffs,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# CI and raw p-value
df_interaction_contrasts <- summary(
  interaction_contrasts,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the omnibus interaction
interaction_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "alpha",
    hemisphere == "ipsi",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_p_FDR) == 1)
df_interaction_contrasts <- df_interaction_contrasts %>%
  mutate(
    p_FDR = interaction_p_FDR
  )
df_interaction_contrasts

## Ipsilateral alpha power: group × validity
# Estimated marginal means
em_alpha_pow_valid <- emmeans(
  m_alpha_ipsi,
  ~ group | validity,
  lmer.df = "asymptotic"
)
# group differences within validity and validity effects within group
contr_alpha_pow_valid <- contrast(
  em_alpha_pow_valid,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs + FDR-adjusted p-values
df_alpha_pow_valid <- summary(
  contr_alpha_pow_valid,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_alpha_pow_valid

## Contrast of contrasts: group × validity
# Validity contrast separately within ASD and NT
diffs <- contrast(
  em_alpha_pow_valid,
  method = "revpairwise",
  by = "group",
  adjust = "none"
)
# Difference between the two validity contrasts
interaction_contrasts <- contrast(
  diffs,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# CI and raw p-value
df_interaction_contrasts <- summary(
  interaction_contrasts,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the omnibus interaction
interaction_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "alpha",
    hemisphere == "ipsi",
    effect == "group:validity"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_p_FDR) == 1)
df_interaction_contrasts <- df_interaction_contrasts %>%
  mutate(
    p_FDR = interaction_p_FDR
  )
df_interaction_contrasts

## Ipsilateral alpha power: group × ROI
# EMMs for group within ROI
em_alpha_pow_ROI <- emmeans(
  m_alpha_ipsi,
  ~ group | ROI,
  lmer.df = "asymptotic"
)
# group differences within each ROI and ROI differences within each group
contr_alpha_pow_ROI <- contrast(
  em_alpha_pow_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs and BH-FDR-adjusted p-values
df_alpha_pow_ROI <- summary(
  contr_alpha_pow_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_alpha_pow_ROI

## Contrast of contrasts: group × ROI
# ROI contrast separately within ASD and NT
diffs_ROI <- contrast(
  em_alpha_pow_ROI,
  method = "revpairwise",
  by     = "group",
  adjust = "none"
)
# Difference between the group-specific ROI contrasts
interaction_contrasts_ROI <- contrast(
  diffs_ROI,
  method = "pairwise",
  by     = NULL,
  adjust = "none"
)
# Unadjusted 95% CI and raw p-value
df_interaction_contrasts_ROI <- summary(
  interaction_contrasts_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Attach the FDR-adjusted p-value from the group × ROI omnibus test
interaction_ROI_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "alpha",
    hemisphere == "ipsi",
    effect == "group:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_ROI_p_FDR) == 1)
df_interaction_contrasts_ROI <- df_interaction_contrasts_ROI %>%
  mutate(
    p_FDR = interaction_ROI_p_FDR
  )
df_interaction_contrasts_ROI


## Contrateral alpha power: ROI × validity
# Estimated marginal means
em_alpha_pow_valid_ROI <- emmeans(
  m_alpha_contra,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)
# ROI differences within each validity condition and validity effects within each ROI
contr_alpha_pow_valid_ROI <- contrast(
  em_alpha_pow_valid_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs + FDR-adjusted p-values
df_alpha_pow_valid_ROI <- summary(
  contr_alpha_pow_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_alpha_pow_valid_ROI

## Contrast of contrasts: ROI × validity
diffs <- contrast(
  em_alpha_pow_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Difference between the two validity contrasts
interaction_contrasts <- contrast(
  diffs,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# CI and raw p-value
df_interaction_contrasts <- summary(
  interaction_contrasts,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the omnibus interaction
interaction_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "alpha",
    hemisphere == "contra",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_p_FDR) == 1)
df_interaction_contrasts <- df_interaction_contrasts %>%
  mutate(
    p_FDR = interaction_p_FDR
  )
df_interaction_contrasts  


###### Beta power post-hoc ####
## ipsiteral beta power: ROI × validity
# Estimated marginal means
em_beta_pow_valid_ROI <- emmeans(
  m_beta_ipsi,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)
# ROI differences within each validity condition
# and validity effects within each ROI
contr_beta_pow_valid_ROI <- contrast(
  em_beta_pow_valid_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs + FDR-adjusted p-values
df_beta_pow_valid_ROI <- summary(
  contr_beta_pow_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_beta_pow_valid_ROI

## contrast of contrast: ROI × validity
diffs <- contrast(
  em_beta_pow_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Difference between the two validity ipsists
interaction_contrast <- contrast(
  diffs,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# CI and raw p-value
df_interaction_contrast <- summary(
  interaction_contrast,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the omnibus interaction
interaction_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "beta",
    hemisphere == "ipsi",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_p_FDR) == 1)
df_interaction_contrast <- df_interaction_contrast %>%
  mutate(
    p_FDR = interaction_p_FDR
  )
df_interaction_contrast  

## Ipsilateral beta power: group × ROI
em_beta_pow_ipsi_ROI <- emmeans(
  m_beta_ipsi,
  ~ group | ROI,
  lmer.df = "asymptotic"
)
# Obtain group differences within each ROI and ROI differences within each group
contr_beta_pow_ipsi_ROI <- contrast(
  em_beta_pow_ipsi_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted asymptotic 95% CIs and FDR-adjusted p-values
df_beta_pow_ipsi_ROI <- summary(
  contr_beta_pow_ipsi_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_beta_pow_ipsi_ROI

## Contrast of contrasts: group × ROI 
# ROI contrast separately within ASD and NT
diffs_beta_ipsi_ROI <- contrast(
  em_beta_pow_ipsi_ROI,
  method = "revpairwise",
  by     = "group",
  adjust = "none"
)
# Difference between the two group-specific ROI contrasts
interaction_contrasts_beta_ipsi_ROI <- contrast(
  diffs_beta_ipsi_ROI,
  method = "pairwise",
  by     = NULL,
  adjust = "none"
)
# Unadjusted asymptotic 95% CI and raw p-value
df_interaction_beta_ipsi_ROI <- summary(
  interaction_contrasts_beta_ipsi_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the corresponding omnibus interaction
interaction_beta_ipsi_ROI_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "beta",
    hemisphere == "ipsi",
    effect == "group:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_beta_ipsi_ROI_p_FDR) == 1)
df_interaction_beta_ipsi_ROI <- df_interaction_beta_ipsi_ROI %>%
  mutate(
    p_FDR = interaction_beta_ipsi_ROI_p_FDR
  )

df_interaction_beta_ipsi_ROI


## Contralateral beta power: ROI × validity
# Estimated marginal means
em_beta_pow_valid_ROI <- emmeans(
  m_beta_contra,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)
# ROI differences within each validity condition
# and validity effects within each ROI
contr_beta_pow_valid_ROI <- contrast(
  em_beta_pow_valid_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted 95% CIs + FDR-adjusted p-values
df_beta_pow_valid_ROI <- summary(
  contr_beta_pow_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(p.value, method = "BH")
  )
df_beta_pow_valid_ROI

## contrast of contrast: ROI × validity
diffs <- contrast(
  em_beta_pow_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Difference between the two validity contrasts
interaction_contrast <- contrast(
  diffs,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# CI and raw p-value
df_interaction_contrast <- summary(
  interaction_contrast,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the omnibus interaction
interaction_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "beta",
    hemisphere == "contra",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_p_FDR) == 1)
df_interaction_contrast <- df_interaction_contrast %>%
  mutate(
    p_FDR = interaction_p_FDR
  )
df_interaction_contrast  


## Contralateral beta power: group × ROI
em_beta_pow_contra_ROI <- emmeans(
  m_beta_contra,
  ~ group | ROI,
  lmer.df = "asymptotic"
)
# Obtain group differences within each ROI and ROI differences within each group
contr_beta_pow_contra_ROI <- contrast(
  em_beta_pow_contra_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted asymptotic 95% CIs and FDR-adjusted p-values
df_beta_pow_contra_ROI <- summary(
  contr_beta_pow_contra_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_beta_pow_contra_ROI
## Contrast of contrasts: group × ROI
# ROI contrast separately within ASD and NT
diffs_beta_contra_ROI <- contrast(
  em_beta_pow_contra_ROI,
  method = "revpairwise",
  by     = "group",
  adjust = "none"
)
# Difference between the two group-specific ROI contrasts
interaction_contrasts_beta_contra_ROI <- contrast(
  diffs_beta_contra_ROI,
  method = "pairwise",
  by     = NULL,
  adjust = "none"
)
# Unadjusted asymptotic 95% CI and raw p-value
df_interaction_beta_contra_ROI <- summary(
  interaction_contrasts_beta_contra_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the corresponding omnibus interaction
interaction_beta_contra_ROI_p_FDR <- all_power_tests_FDR %>%
  filter(
    band == "beta",
    hemisphere == "contra",
    effect == "group:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_beta_contra_ROI_p_FDR) == 1)
df_interaction_beta_contra_ROI <- df_interaction_beta_contra_ROI %>%
  mutate(
    p_FDR = interaction_beta_contra_ROI_p_FDR
  )
df_interaction_beta_contra_ROI



###### Sensitivity analysis: Power after exclusion of participants with above-cutoff score in FBB-ADHD (increased ADHD scores) ####
# exclude participants from data frame
setwd("...\\input_files")
load("adhd_id")
trial_full_df_adhd<-trial_full_df[!trial_full_df$id %in% adhd_id,]

# rerun models with significant differences
## Prepare data
power_df_adhd <- trial_full_df_adhd %>%
  filter(
    band %in% c("alpha", "beta"),
    hemi %in% c("contra", "ipsi")
  ) %>%
  mutate(
    # Explicit orders determine the directions of post-hoc contrasts
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI  = factor(ROI),
    sex  = factor(sex),
    id   = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    hemi = factor(
      hemi,
      levels = c("contra", "ipsi")
    )
  ) %>%
  droplevels()
# check levels
levels(power_df_adhd$group)
levels(power_df_adhd$validity)
levels(power_df_adhd$ROI)

###### Models ####
# Alpha contralateral
m_alpha_contra <- lmer(scale(dB_total) ~ group * validity * ROI +
                         scale(age) + scale(IQ) + sex + (1 | id), data = power_df_adhd %>% filter(band == "alpha", hemi == "contra"))
# Alpha ipsilateral
m_alpha_ipsi <- lmer(scale(dB_total) ~ group * validity * ROI +
                       scale(age) + scale(IQ) + sex + (1 | id), data = power_df_adhd %>% filter(band == "alpha", hemi == "ipsi"))
# Beta contralateral
m_beta_contra <- lmer(scale(dB_total) ~ group * validity * ROI +
                        scale(age) + scale(IQ) + sex + (1 | id), data = power_df_adhd %>% filter(band == "beta", hemi == "contra"))
# Beta ipsilateral
m_beta_ipsi <- lmer(scale(dB_total) ~ group * validity * ROI +
                      scale(age) + scale(IQ) + sex + (1 | id), data = power_df_adhd %>% filter(band == "beta", hemi == "ipsi"))

## Define & extract omnibus tests 
effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)
# Function to extract Type III omnibus tests 
extract_power_tests <- function(models) {
  
  imap_dfr(
    models,
    function(model, hemisphere_name) {
      
      anova(
        model,
        type = 3,
        ddf  = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          hemisphere = hemisphere_name,
          effect      = effect,
          df_num      = NumDF,
          df_den      = DenDF,
          F_value     = `F value`,
          p_raw       = `Pr(>F)`
        ) %>%
        filter(
          effect %in% effects_of_interest
        )
    }
  )
}
# Extract alpha omnibus tests
alpha_power_models <- list(
  contra = m_alpha_contra,
  ipsi   = m_alpha_ipsi
)
alpha_power_tests <- extract_power_tests(
  alpha_power_models
)
alpha_power_tests
# Extract beta omnibus tests
beta_power_models <- list(
  contra = m_beta_contra,
  ipsi   = m_beta_ipsi
)

beta_power_tests <- extract_power_tests(
  beta_power_models
)
beta_power_tests

## Run omnibus test
# merge tests
all_power_tests <- bind_rows(
  alpha_power_tests %>%
    mutate(
      band = "alpha",
      .before = hemisphere
    ),
  beta_power_tests %>%
    mutate(
      band = "beta",
      .before = hemisphere
    )
)
# For each omnibus effect, apply BH-FDR across all four models
all_power_tests_FDR <- all_power_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    hemisphere = factor(
      hemisphere,
      levels = c("contra", "ipsi")
    ),
    effect = factor(
      effect,
      levels = effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    effect,
    band,
    hemisphere
  )

all_power_tests_FDR
# Format p-values
fmt_p <- function(x) {
  ifelse(x < .001, "< .001", sub("^0", "", sprintf("%.3f", x)))
}
# Arrange rows by model, then by omnibus effect
power_table <- all_power_tests_FDR %>%
  arrange(
    match(
      paste(band, hemisphere),
      c(
        "alpha ipsi",
        "alpha contra",
        "beta ipsi",
        "beta contra"
      )
    ),
    effect
  ) %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Hemi = recode(
      as.character(hemisphere),
      ipsi   = "Ipsi",
      contra = "Contra"
    ),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value, df_num, df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(as.data.frame(power_table), row.names = FALSE, right = TRUE)

# save as word table
output_file <- paste0(
  ".../Tables and Figures/",
  "Sensitivity_ADHD_Table_Power_Omnibus_Tests_FDR.docx"
)
ft <- flextable(power_table) %>% autofit()
save_as_docx(
  "Power omnibus tests" = ft,
  path = output_file
)



#### 2.2 LATERALIZATION  --------------------------------------------------------------------------------------------------------------
# calculated from ipsilateral minus contralateral activity

# plot
ggplot(targ_li_df, aes(x=LI_ipsi_minus_contra, fill=group)) + geom_density(alpha=0.2) + facet_grid(ROI ~ band)

## Prepare data 
li_df <- targ_li_df %>%
  filter(
    band %in% c("alpha", "beta")
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(ROI),
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    sex = factor(sex),
    id = factor(id)
  ) %>%
  droplevels()
# Check factor levels and observations
levels(li_df$group)
levels(li_df$validity)
levels(li_df$ROI)
levels(li_df$band)
li_df %>%
  count(band, ROI, validity, group)
stopifnot(
  nlevels(li_df$group) == 2,
  nlevels(li_df$validity) == 2,
  nlevels(li_df$ROI) == 2,
  nlevels(li_df$band) == 2
)
## Fit models
# Alpha lateralization
m_li_alpha <- lmer(
  scale(LI_ipsi_minus_contra) ~
    group * validity * ROI +
    scale(age) + scale(IQ) + sex +
    (1 | id),
  data = li_df %>%
    filter(band == "alpha")
)
# Beta lateralization
m_li_beta <- lmer(
  scale(LI_ipsi_minus_contra) ~
    group * validity * ROI +
    scale(age) + scale(IQ) + sex +
    (1 | id),
  data = li_df %>%
    filter(band == "beta")
)
## Inspect original model tables
anova(
  m_li_alpha,
  type = 3,
  ddf = "Satterthwaite"
)
anova(
  m_li_beta,
  type = 3,
  ddf = "Satterthwaite"
)
## Define omnibus effects
li_effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)
## Extract Type III omnibus tests
extract_li_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% li_effects_of_interest
        )
    }
  )
}
## Combine alpha and beta results
li_models <- list(
  alpha = m_li_alpha,
  beta = m_li_beta
)
all_li_tests <- extract_li_tests(
  li_models
)
all_li_tests
## BH-FDR correction
# For each omnibus effect, correct across alpha and beta
all_li_tests_FDR <- all_li_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = li_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_li_tests_FDR
## Table
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
li_table <- all_li_tests_FDR %>%
  arrange(
    band,
    effect
  ) %>%
  transmute(
    Band = tools::toTitleCase(
      as.character(band)
    ),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(
  as.data.frame(li_table),
  row.names = FALSE,
  right = TRUE
)
## Save as Word document
output_file_li <- paste0(
  ".../Tables and Figures/",
  "Table_Lateralization_Omnibus_Tests_FDR.docx"
)
ft_li <- flextable(li_table) %>% autofit()
save_as_docx(
  "Lateralization omnibus tests" = ft_li,
  path = output_file_li
)


###### Alpha lateralization post-hoc ####
## Validity × ROI
em_li_alpha_valid_ROI <- emmeans(
  m_li_alpha,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)

contr_li_alpha_valid_ROI <- contrast(
  em_li_alpha_valid_ROI,
  method = "consec",
  simple = "each",
  combine = TRUE,
  adjust = "none"
)

df_li_alpha_valid_ROI <- summary(
  contr_li_alpha_valid_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )

df_li_alpha_valid_ROI
# Valid minus invalid separately within each ROI
diffs_li_alpha_valid_ROI <- contrast(
  em_li_alpha_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Frontal validity effect minus posterior validity effect
interaction_li_alpha_valid_ROI <- contrast(
  diffs_li_alpha_valid_ROI,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
df_interaction_li_alpha_valid_ROI <- summary(
  interaction_li_alpha_valid_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Retrieve corrected omnibus validity × ROI p-value
interaction_li_alpha_valid_ROI_p_FDR <- all_li_tests_FDR %>%
  filter(
    band == "alpha",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(
  length(interaction_li_alpha_valid_ROI_p_FDR) == 1
)
df_interaction_li_alpha_valid_ROI <-
  df_interaction_li_alpha_valid_ROI %>%
  mutate(
    p_FDR = interaction_li_alpha_valid_ROI_p_FDR
  )
df_interaction_li_alpha_valid_ROI


## Alpha: Group × ROI post-hoc
em_li_alpha_group_ROI <- emmeans(
  m_li_alpha,
  ~ group | ROI,
  lmer.df = "asymptotic"
)
contr_li_alpha_group_ROI <- contrast(
  em_li_alpha_group_ROI,
  method = "consec",
  simple = "each",
  combine = TRUE,
  adjust = "none"
)
df_li_alpha_group_ROI <- summary(
  contr_li_alpha_group_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_li_alpha_group_ROI
# Posterior minus frontal separately within ASD and NT
diffs_li_alpha_group_ROI <- contrast(
  em_li_alpha_group_ROI,
  method = "revpairwise",
  by = "group",
  adjust = "none"
)
# Difference between group-specific regional contrasts
interaction_li_alpha_group_ROI <- contrast(
  diffs_li_alpha_group_ROI,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
df_interaction_li_alpha_group_ROI <- summary(
  interaction_li_alpha_group_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
interaction_li_alpha_group_ROI_p_FDR <- all_li_tests_FDR %>%
  filter(
    band == "alpha",
    effect == "group:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(
  length(interaction_li_alpha_group_ROI_p_FDR) == 1
)
df_interaction_li_alpha_group_ROI <-
  df_interaction_li_alpha_group_ROI %>%
  mutate(
    p_FDR = interaction_li_alpha_group_ROI_p_FDR
  )

df_interaction_li_alpha_group_ROI


###### Beta lateralization post hoc ####
## Validity × ROI Post-hoc
em_li_beta_valid_ROI <- emmeans(
  m_li_beta,
  ~ ROI | validity,
  lmer.df = "asymptotic"
)
contr_li_beta_valid_ROI <- contrast(
  em_li_beta_valid_ROI,
  method = "consec",
  simple = "each",
  combine = TRUE,
  adjust = "none"
)
df_li_beta_valid_ROI <- summary(
  contr_li_beta_valid_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_li_beta_valid_ROI
# Valid minus invalid separately within each ROI
diffs_li_beta_valid_ROI <- contrast(
  em_li_beta_valid_ROI,
  method = "revpairwise",
  by = "ROI",
  adjust = "none"
)
# Frontal validity effect minus posterior validity effect
interaction_li_beta_valid_ROI <- contrast(
  diffs_li_beta_valid_ROI,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
df_interaction_li_beta_valid_ROI <- summary(
  interaction_li_beta_valid_ROI,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
interaction_li_beta_valid_ROI_p_FDR <- all_li_tests_FDR %>%
  filter(
    band == "beta",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(
  length(interaction_li_beta_valid_ROI_p_FDR) == 1
)
df_interaction_li_beta_valid_ROI <-
  df_interaction_li_beta_valid_ROI %>%
  mutate(
    p_FDR = interaction_li_beta_valid_ROI_p_FDR
  )

df_interaction_li_beta_valid_ROI


###### Sensitivity analysis: Lateralization after exclusion of participants with above cut-off ADHD symptoms (parent-rated: FBB-ADHD) ####
setwd("...\\input_files")
load("adhd_id")
targ_li_df_adhd<-targ_li_df[!targ_li_df$id %in% adhd_id,]

## Prepare data 
li_df_adhd <- targ_li_df_adhd %>%
  filter(
    band %in% c("alpha", "beta")
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(ROI),
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    sex = factor(sex),
    id = factor(id)
  ) %>%
  droplevels()
# Check factor levels and observations
levels(li_df_adhd$group)
levels(li_df_adhd$validity)
levels(li_df_adhd$ROI)
levels(li_df_adhd$band)
li_df_adhd %>%
  count(band, ROI, validity, group)
stopifnot(
  nlevels(li_df_adhd$group) == 2,
  nlevels(li_df_adhd$validity) == 2,
  nlevels(li_df_adhd$ROI) == 2,
  nlevels(li_df_adhd$band) == 2
)
## Fit models
# Alpha lateralization
m_li_alpha <- lmer(
  scale(LI_ipsi_minus_contra) ~
    group * validity * ROI +
    scale(age) + scale(IQ) + sex +
    (1 | id),
  data = li_df_adhd %>%
    filter(band == "alpha")
)
# Beta lateralization
m_li_beta <- lmer(
  scale(LI_ipsi_minus_contra) ~
    group * validity * ROI +
    scale(age) + scale(IQ) + sex +
    (1 | id),
  data = li_df_adhd %>%
    filter(band == "beta")
)
## Inspect original model tables
anova(
  m_li_alpha,
  type = 3,
  ddf = "Satterthwaite"
)
anova(
  m_li_beta,
  type = 3,
  ddf = "Satterthwaite"
)
## Define omnibus effects
li_effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)
## Extract Type III omnibus tests
extract_li_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% li_effects_of_interest
        )
    }
  )
}
## Combine alpha and beta results
li_models <- list(
  alpha = m_li_alpha,
  beta = m_li_beta
)
all_li_tests <- extract_li_tests(
  li_models
)
all_li_tests
## BH-FDR correction
# For each omnibus effect, correct across alpha and beta
all_li_tests_FDR <- all_li_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = li_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_li_tests_FDR
## Table
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
li_table <- all_li_tests_FDR %>%
  arrange(
    band,
    effect
  ) %>%
  transmute(
    Band = tools::toTitleCase(
      as.character(band)
    ),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(
  as.data.frame(li_table),
  row.names = FALSE,
  right = TRUE
)
## Save as Word document
output_file_li <- paste0(
  ".../Tables and Figures/",
  "Sensitivity_ADHD_Table_Lateralization_Omnibus_Tests_FDR.docx"
)
ft_li <- flextable(li_table) %>% autofit()
save_as_docx(
  "Lateralization omnibus tests" = ft_li,
  path = output_file_li
)

#### 2.3 INTER-TRIAL PHASE COHERENCE (ITPC) -----------------------------------------------------------------------------------------------------------------------------
# ITPC matched for n-trials (n=20 within validity, SOA, target side conditions) -> with iterative resampling

# plot
ggplot(targ_itc_matchedN[targ_itc_matchedN$hemi=="contra",], aes(x=ITC_mean, fill=group)) + geom_density(alpha=0.2) + facet_grid(ROI ~ band)

## FDR correction across the two contralateral ITPC models
## prepare data 
itpc_df <- targ_itc_matchedN %>%
  filter(
    ROI != "pari",
    band %in% c("alpha", "beta"),
    hemi == "contra"
  ) %>%
  mutate(
    # Explicit orders determine post-hoc contrast directions
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(
      ROI,
      levels = c("front", "occi")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()
# Verify that the expected two ROIs are included
stopifnot(
  nlevels(itpc_df$group) == 2,
  nlevels(itpc_df$validity) == 2,
  nlevels(itpc_df$ROI) == 2
)
levels(itpc_df$ROI)

###### Models ####
# Contralateral alpha ITPC
m_itpc_alpha <- lmer(scale(ITC_mean) ~ group * validity * ROI +
                       scale(age) + scale(IQ) + sex + (1 | id), data = itpc_df %>% filter(band == "alpha"))
# Contralateral beta ITPC
m_itpc_beta <- lmer(scale(ITC_mean) ~ group * validity * ROI +
                      scale(age) + scale(IQ) + sex + (1 | id), data = itpc_df %>% filter(band == "beta"))

## Omnibus test
# Define omnibus effects of interest
itpc_effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)

# Extract Type III omnibus tests
extract_itpc_tests <- function(models) {
  
  imap_dfr(
    models,
    function(model, band_name) {
      
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% itpc_effects_of_interest
        )
    }
  )
}

# Combine alpha and beta results
itpc_models <- list(
  alpha = m_itpc_alpha,
  beta  = m_itpc_beta
)
all_itpc_tests <- extract_itpc_tests(
  itpc_models
)

## BH-FDR correction
# For each omnibus effect, correct across the alpha and beta ITPC models
all_itpc_tests_FDR <- all_itpc_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = itpc_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_itpc_tests_FDR


# Prepare table
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
itpc_table <- all_itpc_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(
  as.data.frame(itpc_table),
  row.names = FALSE,
  right = TRUE
)

# Save as Word table 
output_file_itpc <- paste0(
  ".../Tables and Figures/",
  "Table_ITPC_Omnibus_Tests_FDR.docx"
)
ft_itpc <- flextable(itpc_table) %>% autofit()
save_as_docx(
  "ITPC omnibus tests" = ft_itpc,
  path = output_file_itpc
)

###### ITPC Post hoc ####
## Contralateral beta ITPC: group × ROI
em_beta_itpc_contra_ROI <- emmeans(
  m_itpc_beta,
  ~ group | ROI,
  lmer.df = "asymptotic"
)
# Obtain group differences within each ROI and ROI differences within each group
contr_beta_itpc_contra_ROI <- contrast(
  em_beta_itpc_contra_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted asymptotic 95% CIs and FDR-adjusted p-values
df_beta_itpc_contra_ROI <- summary(
  contr_beta_itpc_contra_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_beta_itpc_contra_ROI
## Contrast of contrasts: group × ROI
# ROI contrast separately within ASD and NT
diffs_beta_contra_ROI <- contrast(
  em_beta_itpc_contra_ROI,
  method = "revpairwise",
  by     = "group",
  adjust = "none"
)
# Difference between the two group-specific ROI contrasts
interaction_contrasts_beta_contra_ROI <- contrast(
  diffs_beta_contra_ROI,
  method = "pairwise",
  by     = NULL,
  adjust = "none"
)
# Unadjusted asymptotic 95% CI and raw p-value
df_interaction_beta_contra_ROI <- summary(
  interaction_contrasts_beta_contra_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the corresponding omnibus interaction
interaction_beta_contra_ROI_p_FDR <- all_itpc_tests_FDR %>%
  filter(
    band == "beta",
    effect == "group:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_beta_contra_ROI_p_FDR) == 1)
df_interaction_beta_contra_ROI <- df_interaction_beta_contra_ROI %>%
  mutate(
    p_FDR = interaction_beta_contra_ROI_p_FDR
  )
df_interaction_beta_contra_ROI


## Contralateral beta ITPC: validity × ROI
em_beta_itpc_contra_valid_ROI <- emmeans(
  m_itpc_beta,
  ~ validity | ROI,
  lmer.df = "asymptotic"
)
# Obtain group differences within each ROI and ROI differences within each group
contr_beta_itpc_contra_valid_ROI <- contrast(
  em_beta_itpc_contra_valid_ROI,
  method  = "consec",
  simple  = "each",
  combine = TRUE,
  adjust  = "none"
)
# Unadjusted asymptotic 95% CIs and FDR-adjusted p-values
df_beta_itpc_contra_valid_ROI <- summary(
  contr_beta_itpc_contra_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_beta_itpc_contra_valid_ROI
## Contrast of contrasts: validity × ROI
diffs_beta_contra_valid_ROI <- contrast(
  em_beta_itpc_contra_valid_ROI,
  method = "revpairwise",
  by     = "validity",
  adjust = "none"
)
# Difference between the two group-specific ROI contrasts
interaction_contrasts_beta_contra_valid_ROI <- contrast(
  diffs_beta_contra_valid_ROI,
  method = "pairwise",
  by     = NULL,
  adjust = "none"
)
# Unadjusted asymptotic 95% CI and raw p-value
df_interaction_beta_contra_valid_ROI <- summary(
  interaction_contrasts_beta_contra_valid_ROI,
  infer  = c(TRUE, TRUE),
  level  = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add FDR-adjusted p-value from the corresponding omnibus interaction
interaction_beta_contra_valid_ROI_p_FDR <- all_itpc_tests_FDR %>%
  filter(
    band == "beta",
    effect == "validity:ROI"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_beta_contra_valid_ROI_p_FDR) == 1)
df_interaction_beta_contra_valid_ROI <- df_interaction_beta_contra_valid_ROI %>%
  mutate(
    p_FDR = interaction_beta_contra_valid_ROI_p_FDR
  )
df_interaction_beta_contra_valid_ROI


###### Sensitivity analysis: Exclusion of participants with above-cutoff score in FBB-ADHD (increased ADHD scores) ####
# exclude participants from data frame
setwd("...\\input_files")
load("adhd_id")
targ_itc_matchedN_adhd<-targ_itc_matchedN[!targ_itc_matchedN$id %in% adhd_id,]

## FDR correction across the two contralateral ITPC models
## prepare data 
itpc_df_adhd <- targ_itc_matchedN_adhd %>%
  filter(
    ROI != "pari",
    band %in% c("alpha", "beta"),
    hemi == "contra"
  ) %>%
  mutate(
    # Explicit orders determine post-hoc contrast directions
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    ROI = factor(
      ROI,
      levels = c("front", "occi")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()
# Verify that the expected two ROIs are included
stopifnot(
  nlevels(itpc_df_adhd$group) == 2,
  nlevels(itpc_df_adhd$validity) == 2,
  nlevels(itpc_df_adhd$ROI) == 2
)
levels(itpc_df_adhd$ROI)

###### Models ####
# Contralateral alpha ITPC
m_itpc_alpha <- lmer(scale(ITC_mean) ~ group * validity * ROI +
                       scale(age) + scale(IQ) + sex + (1 | id), data = itpc_df_adhd %>% filter(band == "alpha"))
# Contralateral beta ITPC
m_itpc_beta <- lmer(scale(ITC_mean) ~ group * validity * ROI +
                      scale(age) + scale(IQ) + sex + (1 | id), data = itpc_df_adhd %>% filter(band == "beta"))

## Omnibus test
# Define omnibus effects of interest
itpc_effects_of_interest <- c(
  "group", "validity", "ROI",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity", "group:ROI", "validity:ROI",
  "group:validity:ROI"
)

# Extract Type III omnibus tests
extract_itpc_tests <- function(models) {
  
  imap_dfr(
    models,
    function(model, band_name) {
      
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% itpc_effects_of_interest
        )
    }
  )
}

# Combine alpha and beta results
itpc_models <- list(
  alpha = m_itpc_alpha,
  beta  = m_itpc_beta
)
all_itpc_tests <- extract_itpc_tests(
  itpc_models
)

## BH-FDR correction
# For each omnibus effect, correct across the alpha and beta ITPC models
all_itpc_tests_FDR <- all_itpc_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = itpc_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_itpc_tests_FDR


# Prepare table
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
itpc_table <- all_itpc_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      ROI = "ROI",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity",
      `group:ROI` = "Group × ROI",
      `validity:ROI` = "Validity × ROI",
      `group:validity:ROI` = "Group × validity × ROI"
    ),
    
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(
  as.data.frame(itpc_table),
  row.names = FALSE,
  right = TRUE
)

# Save as Word table 
output_file_itpc <- paste0(
  ".../Tables and Figures/",
  "Sensitivity_ADHD_ITPC omnibus tests.docx"
)
ft_itpc <- flextable(itpc_table) %>% autofit()
save_as_docx(
  "Sensitivity_ADHD_ITPC omnibus tests" = ft_itpc,
  path = output_file_itpc
)



#### 2.4 POWER CONNECTIVITY (Long-range) -----------------------------------------------------------------------------------------------------------------------
# Occipital–frontal, contralateral
# FDR correction across alpha and beta within each omnibus-effect family

# plot distribution
ggplot(df_powconnect[df_powconnect$hemi=="contra",], aes(x=pow_cor, fill=group)) + geom_density(alpha=0.2) + facet_grid(ROI_interaction ~ band)

## Omnibus test
## Prepare data
powcon_df <- df_powconnect %>%
  filter(
    band %in% c("alpha", "beta"),
    ROI_interaction == "front.occi",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()

stopifnot(
  nlevels(powcon_df$group) == 2,
  nlevels(powcon_df$validity) == 2
)

###### Models ####
# Alpha power connectivity
m_powcon_alpha <- lmer(scale(pow_cor) ~ group * validity +
                         scale(age) + scale(IQ) + sex + (1 | id), data = powcon_df %>% filter(band == "alpha"))
# Beta power connectivity
m_powcon_beta <- lmer(scale(pow_cor) ~ group * validity +
                        scale(age) + scale(IQ) + sex + (1 | id), data = powcon_df %>% filter(band == "beta"))


## Define and extract omnibus tests
powcon_effects_of_interest <- c(
  "group", "validity",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity"
)
extract_powcon_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% powcon_effects_of_interest
        )
    }
  )
}
powcon_models <- list(
  alpha = m_powcon_alpha,
  beta = m_powcon_beta
)
all_powcon_tests <- extract_powcon_tests(
  powcon_models
)
## BH-FDR correction
# Correct each omnibus effect across the alpha and beta models
all_powcon_tests_FDR <- all_powcon_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = powcon_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_powcon_tests_FDR

# Format p-values
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
powcon_table <- all_powcon_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity"
    ),
    
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(
  as.data.frame(powcon_table),
  row.names = FALSE,
  right = TRUE
)
# save as word table
output_file_powcon <- paste0(
  ".../Tables and Figures/",
  "Table_PowerConnectivity_Omnibus_Tests_FDR.docx"
)
ft_powcon <- flextable(powcon_table) %>% autofit()
save_as_docx(
  "Power-connectivity omnibus tests" = ft_powcon,
  path = output_file_powcon
)


###### Sensitivity analysis: Exclusion of participants with above-cutoff score in FBB-ADHD (increased ADHD scores) ####
# exclude participants from data frame
setwd("...\\input_files")
load("adhd_id")
df_powconnect_adhd<-df_powconnect[!df_powconnect$id %in% adhd_id,]

## Omnibus test
## Prepare data
powcon_df_adhd <- df_powconnect_adhd %>%
  filter(
    band %in% c("alpha", "beta"),
    ROI_interaction == "front.occi",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()

stopifnot(
  nlevels(powcon_df_adhd$group) == 2,
  nlevels(powcon_df_adhd$validity) == 2
)

###### Models ####
# Alpha power connectivity
m_powcon_alpha <- lmer(scale(pow_cor) ~ group * validity +
                         scale(age) + scale(IQ) + sex + (1 | id), data = powcon_df_adhd %>% filter(band == "alpha"))
# Beta power connectivity
m_powcon_beta <- lmer(scale(pow_cor) ~ group * validity +
                        scale(age) + scale(IQ) + sex + (1 | id), data = powcon_df_adhd %>% filter(band == "beta"))


## Define and extract omnibus tests
powcon_effects_of_interest <- c(
  "group", "validity",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity"
)
extract_powcon_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% powcon_effects_of_interest
        )
    }
  )
}
powcon_models <- list(
  alpha = m_powcon_alpha,
  beta = m_powcon_beta
)
all_powcon_tests <- extract_powcon_tests(
  powcon_models
)
## BH-FDR correction
# Correct each omnibus effect across the alpha and beta models
all_powcon_tests_FDR <- all_powcon_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = powcon_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_powcon_tests_FDR

# Format p-values
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
powcon_table <- all_powcon_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity"
    ),
    
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(
  as.data.frame(powcon_table),
  row.names = FALSE,
  right = TRUE
)
# save as word table
output_file_powcon <- paste0(
  ".../Tables and Figures/",
  "Sensitivity_ADHD_Table_PowerConnectivity_Omnibus_Tests_FDR.docx"
)
ft_powcon <- flextable(powcon_table) %>% autofit()
save_as_docx(
  "Power-connectivity omnibus tests" = ft_powcon,
  path = output_file_powcon
)


# 2.5 PHASE CONNECTIVITY (long range) -----------------------------------------------------------------------------------
# Contralateral Occipital–frontal Inter-site-phase-coherence (ISPC) ~ Phase-locking value (PLV)
# FDR correction across alpha and beta within each omnibus-effect family

# plot
ggplot(targ_ispc_matchedN[targ_ispc_matchedN$hemi=="contra",], aes(x=ISPC_mean, fill=group)) + geom_density(alpha=0.2) + facet_grid(roi_pair ~ band)

## Omnibus test
## Prepare data
phasecon_df <- targ_ispc_matchedN %>%
  filter(
    band %in% c("alpha", "beta"),
    roi_pair == "occi_front",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()

stopifnot(
  nlevels(phasecon_df$group) == 2,
  nlevels(phasecon_df$validity) == 2
)
###### Mdels ####
# Alpha phase connectivity
m_phasecon_alpha <- lmer(scale(ISPC_mean) ~ group * validity +
                           scale(age) + scale(IQ) + sex + (1 | id), data = phasecon_df %>% filter(band == "alpha"))
# Beta phase connectivity
m_phasecon_beta <- lmer(scale(ISPC_mean) ~ group * validity +
                          scale(age) + scale(IQ) + sex + (1 | id), data = phasecon_df %>% filter(band == "beta"))

## Define and extract omnibus tests
phasecon_effects_of_interest <- c(
  "group", "validity",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity"
)
extract_phasecon_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% phasecon_effects_of_interest
        )
    }
  )
}
phasecon_models <- list(
  alpha = m_phasecon_alpha,
  beta  = m_phasecon_beta
)
all_phasecon_tests <- extract_phasecon_tests(
  phasecon_models
)

## BH-FDR correction 
# Correct each omnibus effect across alpha and beta
all_phasecon_tests_FDR <- all_phasecon_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = phasecon_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_phasecon_tests_FDR

# Format p-values
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
phasecon_table <- all_phasecon_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(
  as.data.frame(phasecon_table),
  row.names = FALSE,
  right = TRUE
)
# save as word table
output_file_phasecon <- paste0(
  ".../Tables and Figures/",
  "Table_PhaseConnectivity_Omnibus_Tests_FDR.docx"
)
ft_phasecon <- flextable(phasecon_table) %>% autofit()
save_as_docx(
  "Phase-connectivity omnibus tests" = ft_phasecon,
  path = output_file_phasecon
)


###### ISPC post-hoc ####
## Alpha ISPC: general validity effect
# Estimated marginal means averaged across groups
em_alpha_ispc_validity <- emmeans(
  m_phasecon_alpha,
  ~ validity,
  lmer.df = "asymptotic"
)
em_alpha_ispc_validity
# Explicit contrast: valid minus invalid
contr_alpha_ispc_validity <- contrast(
  em_alpha_ispc_validity,
  method = list(
    "valid - invalid" = c(-1, 1)
  ),
  adjust = "none"
)
df_alpha_ispc_validity <- summary(
  contr_alpha_ispc_validity,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Add the FDR-adjusted p-value from the omnibus validity test
alpha_ispc_validity_p_FDR <- all_phasecon_tests_FDR %>%
  filter(
    band == "alpha",
    effect == "validity"
  ) %>%
  pull(p_FDR)
stopifnot(length(alpha_ispc_validity_p_FDR) == 1)
df_alpha_ispc_validity <- df_alpha_ispc_validity %>%
  mutate(
    p_FDR = alpha_ispc_validity_p_FDR
  )
df_alpha_ispc_validity


## Beta ISPC: group × validity
em_phasecon_beta_valid <- emmeans(m_phasecon_beta, ~ group | validity, lmer.df = "asymptotic")
# Group differences within validity and validity effects within group
contr_phasecon_beta_valid <- contrast(
  em_phasecon_beta_valid,
  method = "consec",
  simple = "each",
  combine = TRUE,
  adjust = "none"
)
# Unadjusted asymptotic 95% CIs and jointly FDR-adjusted p-values
df_phasecon_beta_valid <- summary(
  contr_phasecon_beta_valid,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = p.adjust(
      p.value,
      method = "BH"
    )
  )
df_phasecon_beta_valid
# Validity effect separately within ASD and NT
diffs_phasecon_beta_valid <- contrast(
  em_phasecon_beta_valid,
  method = "revpairwise",
  by = "group",
  adjust = "none"
)
# Difference between the group-specific validity effects
interaction_phasecon_beta_valid <- contrast(
  diffs_phasecon_beta_valid,
  method = "pairwise",
  by = NULL,
  adjust = "none"
)
# Unadjusted asymptotic 95% CI and raw p-value
df_interaction_phasecon_beta_valid <- summary(
  interaction_phasecon_beta_valid,
  infer = c(TRUE, TRUE),
  level = 0.95,
  adjust = "none"
) %>%
  as.data.frame()
# Retrieve the corrected beta group × validity omnibus p-value
interaction_phasecon_beta_p_FDR <- all_phasecon_tests_FDR %>%
  filter(
    band == "beta",
    effect == "group:validity"
  ) %>%
  pull(p_FDR)
stopifnot(length(interaction_phasecon_beta_p_FDR) == 1)
df_interaction_phasecon_beta_valid <-
  df_interaction_phasecon_beta_valid %>%
  mutate(
    p_FDR = interaction_phasecon_beta_p_FDR
  )

df_interaction_phasecon_beta_valid


###### Sensitivity analysis: Exclusion of participants with above-cutoff score in FBB-ADHD (increased ADHD scores) ####
# exclude participants from data frame
setwd("...\\input_files")
load("adhd_id")
targ_ispc_matchedN_adhd<-targ_ispc_matchedN[!targ_ispc_matchedN$id %in% adhd_id,]


## Omnibus test
## Prepare data
phasecon_df_adhd <- targ_ispc_matchedN_adhd %>%
  filter(
    band %in% c("alpha", "beta"),
    roi_pair == "occi_front",
    hemi == "contra"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("ASD", "NT")
    ),
    validity = factor(
      validity,
      levels = c("invalid", "valid")
    ),
    sex = factor(sex),
    id = factor(id),
    band = factor(
      band,
      levels = c("alpha", "beta")
    )
  ) %>%
  droplevels()

stopifnot(
  nlevels(phasecon_df_adhd$group) == 2,
  nlevels(phasecon_df_adhd$validity) == 2
)

###### Models ####
# Alpha phase connectivity
m_phasecon_alpha <- lmer(scale(ISPC_mean) ~ group * validity +
                           scale(age) + scale(IQ) + sex + (1 | id), data = phasecon_df_adhd %>% filter(band == "alpha"))
# Beta phase connectivity
m_phasecon_beta <- lmer(scale(ISPC_mean) ~ group * validity +
                          scale(age) + scale(IQ) + sex + (1 | id), data = phasecon_df_adhd %>% filter(band == "beta"))

## Define and extract omnibus tests
phasecon_effects_of_interest <- c(
  "group", "validity",
  "scale(age)", "scale(IQ)", "sex",
  "group:validity"
)
extract_phasecon_tests <- function(models) {
  imap_dfr(
    models,
    function(model, band_name) {
      anova(
        model,
        type = 3,
        ddf = "Satterthwaite"
      ) %>%
        as.data.frame() %>%
        rownames_to_column("effect") %>%
        transmute(
          band = band_name,
          effect = effect,
          df_num = NumDF,
          df_den = DenDF,
          F_value = `F value`,
          p_raw = `Pr(>F)`
        ) %>%
        filter(
          effect %in% phasecon_effects_of_interest
        )
    }
  )
}
phasecon_models <- list(
  alpha = m_phasecon_alpha,
  beta  = m_phasecon_beta
)
all_phasecon_tests <- extract_phasecon_tests(
  phasecon_models
)

## BH-FDR correction 
# Correct each omnibus effect across alpha and beta
all_phasecon_tests_FDR <- all_phasecon_tests %>%
  group_by(effect) %>%
  mutate(
    p_FDR = p.adjust(
      p_raw,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    band = factor(
      band,
      levels = c("alpha", "beta")
    ),
    effect = factor(
      effect,
      levels = phasecon_effects_of_interest
    ),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(
    band,
    effect
  )
all_phasecon_tests_FDR

# Format p-values
fmt_p <- function(x) {
  ifelse(
    x < .001,
    "< .001",
    sub("^0", "", sprintf("%.3f", x))
  )
}
phasecon_table <- all_phasecon_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      group = "Group",
      validity = "Validity",
      `scale(age)` = "Age",
      `scale(IQ)` = "IQ",
      sex = "Sex",
      `group:validity` = "Group × validity"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value,
      df_num,
      df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )
print(
  as.data.frame(phasecon_table),
  row.names = FALSE,
  right = TRUE
)
# save as word table
output_file_phasecon <- paste0(
  ".../Tables and Figures/",
  "Sensitivity_ADHD_Table_PhaseConnectivity_Omnibus_Tests_FDR.docx"
)
ft_phasecon <- flextable(phasecon_table) %>% autofit()
save_as_docx(
  "Phase-connectivity omnibus tests" = ft_phasecon,
  path = output_file_phasecon
)

#### 2.6 Figure 2 - Group differences ------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(emmeans)
library(patchwork)
library(grid)

## style
base_family <- "sans"

col_nt  <- "#1F7A83"   
col_asd <- "#C46B45"   
cols <- c("NT" = col_nt, "ASD" = col_asd)

theme_fig <- theme_minimal(base_family = base_family, base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "grey35", linewidth = .4),
    axis.ticks = element_line(colour = "grey35", linewidth = .4),
    axis.title = element_text(size = 10, colour = "grey15"),
    axis.text = element_text(size = 10, colour = "grey15"),
    plot.title = element_text(face = "bold", size = 11.5, colour = "grey10"),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    TRUE     ~ ""
  )
}

fmt_p_fig <- function(p) {
  ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))
}

pull_pfdr <- function(df, eff, band = NULL, hemi = NULL) {
  out <- df
  if (!is.null(band)) out <- out[out$band == band, , drop = FALSE]
  if (!is.null(hemi)) out <- out[out$hemisphere == hemi, , drop = FALSE]
  out <- out[out$effect == eff, , drop = FALSE]
  stopifnot(nrow(out) == 1)
  out$p_FDR
}

section_header <- function(label, fill = "#F5F3F0") {
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = fill, colour = NA) +
    annotate("text", x = 0, y = 0, label = label,
             family = base_family, fontface = "bold",
             size = 5, colour = "grey15") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") +
    theme_void()
}

## helpers: emmeans / contrasts 
get_plot_emm <- function(model, xvar) {
  out <- emmeans(
    model,
    as.formula(paste("~ group *", xvar)),
    lmer.df = "asymptotic"
  ) %>%
    summary(infer = c(TRUE, TRUE), level = .95) %>%
    as.data.frame()
  
  if ("asymp.LCL" %in% names(out)) {
    out <- out %>%
      rename(
        lower.CL = asymp.LCL,
        upper.CL = asymp.UCL
      )
  }
  
  out %>%
    mutate(
      group = factor(as.character(group), levels = c("NT", "ASD"))
    )
}

get_panel_tests <- function(model, xvar) {
  
  if (xvar == "validity") {
    within_lab <- "valid - invalid"
  } else if (xvar == "ROI") {
    within_lab <- "posterior - frontal"
  } else {
    stop("xvar must be 'validity' or 'ROI'")
  }
  
  em_within <- emmeans(
    model,
    as.formula(paste("~", xvar, "| group")),
    lmer.df = "asymptotic"
  )
  
  df_within <- summary(
    contrast(
      em_within,
      method = setNames(list(c(-1, 1)), within_lab),
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    level = .95,
    adjust = "none"
  ) %>%
    as.data.frame()
  
  em_between <- emmeans(
    model,
    as.formula(paste("~ group |", xvar)),
    lmer.df = "asymptotic"
  )
  
  df_between <- summary(
    contrast(
      em_between,
      method = list("NT - ASD" = c(-1, 1)),
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    level = .95,
    adjust = "none"
  ) %>%
    as.data.frame()
  
  if ("asymp.LCL" %in% names(df_within)) {
    df_within <- df_within %>%
      rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  if ("asymp.LCL" %in% names(df_between)) {
    df_between <- df_between %>%
      rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  
  p_all <- p.adjust(c(df_within$p.value, df_between$p.value), method = "BH")
  df_within$p_FDR  <- p_all[seq_len(nrow(df_within))]
  df_between$p_FDR <- p_all[nrow(df_within) + seq_len(nrow(df_between))]
  
  list(
    within  = df_within  %>% mutate(stars = p_to_stars(p_FDR)),
    between = df_between %>% mutate(stars = p_to_stars(p_FDR))
  )
}

get_interaction_df <- function(model, xvar, omnibus_p_fdr) {
  
  if (xvar == "validity") {
    x1 <- "invalid"
    x2 <- "valid"
  } else if (xvar == "ROI") {
    x1 <- "front"
    x2 <- "occi"
  } else {
    stop("xvar must be 'validity' or 'ROI'")
  }
  
  em <- emmeans(
    model,
    as.formula(paste("~ group *", xvar)),
    lmer.df = "asymptotic"
  )
  
  em_df <- as.data.frame(em)
  coef <- rep(0, nrow(em_df))
  
  coef[em_df$group == "ASD" & em_df[[xvar]] == x1] <- -1
  coef[em_df$group == "ASD" & em_df[[xvar]] == x2] <-  1
  coef[em_df$group == "NT"  & em_df[[xvar]] == x1] <-  1
  coef[em_df$group == "NT"  & em_df[[xvar]] == x2] <- -1
  
  out <- summary(
    contrast(
      em,
      method = list("interaction" = coef),
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    level = .95,
    adjust = "none"
  ) %>%
    as.data.frame()
  
  if ("asymp.LCL" %in% names(out)) {
    out <- out %>%
      rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  
  out %>%
    mutate(
      p_FDR = omnibus_p_fdr,
      stars = p_to_stars(p_FDR)
    )
}

## EMM dot-and-whiskers panels
make_emm_panel <- function(
    df_emm, tests,
    xvar, x_levels, x_labels,
    title, subtitle, ylab,
    data_range = NULL) {
  
  # NT and ASD slightly separated horizontally
  dodge  <- .085
  x_base <- seq_along(x_levels)
  
  df_plot <- df_emm %>%
    filter(.data[[xvar]] %in% x_levels) %>%
    mutate(
      group = factor(as.character(group), levels = c("NT", "ASD")),
      x_index = match(.data[[xvar]], x_levels),
      x = x_base[x_index] + ifelse(group == "NT", -dodge, dodge)
    )
  
  # Use shared range when supplied, so panels within a row have the same y-scale
  if (is.null(data_range)) {
    data_range <- range(
      c(df_plot$lower.CL, df_plot$upper.CL),
      na.rm = TRUE
    )
  }
  
  lo <- data_range[1]
  hi <- data_range[2]
  rr <- hi - lo
  if (rr == 0) rr <- 1
  
  # Fixes annotation lanes
  y_nt  <- hi + .31 * rr   # Teal bracket slightly higher
  y_asd <- hi + .13 * rr   # Orange bracket slightly lower
  y_between <- lo - .16 * rr
  
  tick <- .040 * rr
  bracket_star_gap_col <- .010 * rr
  bracket_star_gap_blk <- .020 * rr
  
  # Keep the ORIGINAL y-limits to preserve alignment
  y_lim <- c(
    lo - .30 * rr,
    hi + .39 * rr
  )
  
  p <- ggplot(
    df_plot,
    aes(x = x, y = emmean, colour = group)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey82",
      linewidth = .45
    ) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = .035,
      linewidth = .52
    ) +
    geom_point(
      size = 1.8
    ) +
    scale_colour_manual(
      values = cols,
      breaks = c("NT", "ASD"),
      labels = c("Neurotypical (NT)", "Autistic (ASD)")
    ) +
    scale_x_continuous(
      breaks = x_base,
      labels = x_labels,
      limits = c(.68, length(x_levels) + .32),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = y_lim,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = ylab
    ) +
    theme_fig +
    theme(
      legend.position = "none",
      plot.title = element_text(
        size = 9.7,
        face = "bold",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = 8.0,
        colour = "grey38",
        margin = margin(b = 3)
      ),
      axis.title = element_text(size = 8.7),
      axis.text = element_text(size = 8.0),
      plot.margin = margin(5, 7, 4, 4)
    )
  
  # Within-group contrasts
  # Long brackets, always in fixed top lanes
  within_sig <- tests$within %>%
    filter(stars != "")
  
  if (nrow(within_sig) > 0) {
    
    for (i in seq_len(nrow(within_sig))) {
      
      g <- as.character(within_sig$group[i])
      
      this_col <- if (g == "NT") cols["NT"] else cols["ASD"]
      offset   <- if (g == "NT") -dodge else dodge
      
      y <- if (g == "NT") y_nt else y_asd
      
      x1 <- x_base[1] + offset
      x2 <- x_base[2] + offset
      
      p <- p +
        annotate(
          "segment",
          x = x1, xend = x2,
          y = y, yend = y,
          colour = this_col,
          linewidth = .52
        ) +
        annotate(
          "segment",
          x = x1, xend = x1,
          y = y, yend = y - tick,
          colour = this_col,
          linewidth = .52
        ) +
        annotate(
          "segment",
          x = x2, xend = x2,
          y = y, yend = y - tick,
          colour = this_col,
          linewidth = .52
        ) +
        annotate(
          "text",
          x = mean(c(x1, x2)),
          y = y + bracket_star_gap_col,
          label = within_sig$stars[i],
          colour = this_col,
          family = base_family,
          fontface = "bold",
          size = 3.45,
          vjust = 0
        )
    }
  }
  
  # Between-group simple effects
  between_sig <- tests$between %>%
    filter(
      stars != "",
      .data[[xvar]] %in% x_levels
    )
  
  if (nrow(between_sig) > 0) {
    
    for (i in seq_len(nrow(between_sig))) {
      
      lev <- as.character(between_sig[[xvar]][i])
      xi  <- match(lev, x_levels)
      
      x1 <- x_base[xi] - dodge
      x2 <- x_base[xi] + dodge
      
      p <- p +
        annotate(
          "segment",
          x = x1, xend = x2,
          y = y_between, yend = y_between,
          colour = "black",
          linewidth = .48
        ) +
        annotate(
          "segment",
          x = x1, xend = x1,
          y = y_between, yend = y_between + tick,
          colour = "black",
          linewidth = .48
        ) +
        annotate(
          "segment",
          x = x2, xend = x2,
          y = y_between, yend = y_between + tick,
          colour = "black",
          linewidth = .48
        ) +
        annotate(
          "text",
          x = mean(c(x1, x2)),
          y = y_between - bracket_star_gap_blk,
          label = between_sig$stars[i],
          family = base_family,
          fontface = "bold",
          size = 3.35,
          colour = "black",
          vjust = 1
        )
    }
  }
  
  p
}


## Header / Divider helpers
make_section_header <- function(label, fill, text_size = 4.0) {
  
  ggplot() +
    annotate(
      "rect",
      xmin = 0, xmax = 1,
      ymin = 0, ymax = 1,
      fill = fill,
      colour = NA
    ) +
    annotate(
      "text",
      x = .5, y = .5,
      label = label,
      family = base_family,
      fontface = "bold",
      size = text_size,
      colour = "grey18"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    theme_void()
}


make_divider <- function() {
  
  ggplot() +
    annotate(
      "segment",
      x = .5, xend = .5,
      y = 0, yend = 1,
      colour = "grey88",
      linewidth = .45
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    theme_void()
}


make_text_cell <- function(label, hjust = .5, size = 2.8, face = "bold") {
  
  ggplot() +
    annotate(
      "text",
      x = if (hjust == 0) .02 else .5,
      y = .5,
      label = label,
      hjust = hjust,
      family = base_family,
      fontface = face,
      size = size,
      colour = "grey20"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    theme_void()
}

## Extract all EMMs / post-hoc tests / interactions
## Alpha validity
alpha_valid_emm <- get_plot_emm(
  m_alpha_ipsi,
  "validity"
)

alpha_valid_test <- get_panel_tests(
  m_alpha_ipsi,
  "validity"
)

alpha_valid_int <- get_interaction_df(
  m_alpha_ipsi,
  "validity",
  pull_pfdr(
    all_power_tests_FDR,
    eff = "group:validity",
    band = "alpha",
    hemi = "ipsi"
  )
)

## Alpha spatial
alpha_roi_emm <- get_plot_emm(
  m_alpha_ipsi,
  "ROI"
)

alpha_roi_test <- get_panel_tests(
  m_alpha_ipsi,
  "ROI"
)

alpha_roi_int <- get_interaction_df(
  m_alpha_ipsi,
  "ROI",
  pull_pfdr(
    all_power_tests_FDR,
    eff = "group:ROI",
    band = "alpha",
    hemi = "ipsi"
  )
)


## Beta ipsilateral
beta_ipsi_emm <- get_plot_emm(
  m_beta_ipsi,
  "ROI"
)

beta_ipsi_test <- get_panel_tests(
  m_beta_ipsi,
  "ROI"
)

beta_ipsi_int <- get_interaction_df(
  m_beta_ipsi,
  "ROI",
  pull_pfdr(
    all_power_tests_FDR,
    eff = "group:ROI",
    band = "beta",
    hemi = "ipsi"
  )
)


## Beta contralateral
beta_contra_emm <- get_plot_emm(
  m_beta_contra,
  "ROI"
)

beta_contra_test <- get_panel_tests(
  m_beta_contra,
  "ROI"
)

beta_contra_int <- get_interaction_df(
  m_beta_contra,
  "ROI",
  pull_pfdr(
    all_power_tests_FDR,
    eff = "group:ROI",
    band = "beta",
    hemi = "contra"
  )
)


## Beta ITPC
beta_itpc_emm <- get_plot_emm(
  m_itpc_beta,
  "ROI"
)

beta_itpc_test <- get_panel_tests(
  m_itpc_beta,
  "ROI"
)

beta_itpc_int <- get_interaction_df(
  m_itpc_beta,
  "ROI",
  pull_pfdr(
    all_itpc_tests_FDR,
    eff = "group:ROI",
    band = "beta"
  )
)


## Beta ISPC
beta_ispc_emm <- get_plot_emm(
  m_phasecon_beta,
  "validity"
)

beta_ispc_test <- get_panel_tests(
  m_phasecon_beta,
  "validity"
)

beta_ispc_int <- get_interaction_df(
  m_phasecon_beta,
  "validity",
  pull_pfdr(
    all_phasecon_tests_FDR,
    eff = "group:validity",
    band = "beta"
  )
)


## Shared y-scales
alpha_power_range <- range(
  c(
    alpha_valid_emm$lower.CL,
    alpha_valid_emm$upper.CL,
    alpha_roi_emm$lower.CL,
    alpha_roi_emm$upper.CL
  ),
  na.rm = TRUE
)


beta_power_range <- range(
  c(
    beta_ipsi_emm$lower.CL,
    beta_ipsi_emm$upper.CL,
    beta_contra_emm$lower.CL,
    beta_contra_emm$upper.CL
  ),
  na.rm = TRUE
)


phase_range <- range(
  c(
    beta_itpc_emm$lower.CL,
    beta_itpc_emm$upper.CL,
    beta_ispc_emm$lower.CL,
    beta_ispc_emm$upper.CL
  ),
  na.rm = TRUE
)


## Main EMM panels
panel_final_theme <- theme(
  plot.title = element_text(
    size = 8.5,
    face = "bold",
    margin = margin(b = .5)
  ),
  plot.subtitle = element_text(
    size = 8.0,
    colour = "grey38",
    margin = margin(b = 2)
  ),
  axis.title = element_text(size = 8.0),
  axis.text = element_text(size = 8.0),
  plot.margin = margin(3, 5, 3, 3)
)


# A - Alpha power
p_A <- make_emm_panel(
  alpha_valid_emm,
  alpha_valid_test,
  xvar = "validity",
  x_levels = c("invalid", "valid"),
  x_labels = c("Invalid", "Valid"),
  title = "A  Alpha power",
  subtitle = "Ipsilateral · validity modulation",
  ylab = "Power (z)",
  data_range = alpha_power_range
) + panel_final_theme


# second alpha panel
p_alpha_spatial <- make_emm_panel(
  alpha_roi_emm,
  alpha_roi_test,
  xvar = "ROI",
  x_levels = c("front", "occi"),
  x_labels = c("Frontal", "Posterior"),
  title = " ",
  subtitle = "Ipsilateral · spatial profile",
  ylab = NULL,
  data_range = alpha_power_range
) + panel_final_theme


# B — Beta power
p_B <- make_emm_panel(
  beta_ipsi_emm,
  beta_ipsi_test,
  xvar = "ROI",
  x_levels = c("front", "occi"),
  x_labels = c("Frontal", "Posterior"),
  title = "B  Beta power",
  subtitle = "Ipsilateral · spatial profile",
  ylab = "Power (z)",
  data_range = beta_power_range
) + panel_final_theme


# second beta-power panel
p_beta_contra <- make_emm_panel(
  beta_contra_emm,
  beta_contra_test,
  xvar = "ROI",
  x_levels = c("front", "occi"),
  x_labels = c("Frontal", "Posterior"),
  title = " ",
  subtitle = "Contralateral · spatial profile",
  ylab = NULL,
  data_range = beta_power_range
) + panel_final_theme


# C — Beta ITPC
p_C <- make_emm_panel(
  beta_itpc_emm,
  beta_itpc_test,
  xvar = "ROI",
  x_levels = c("front", "occi"),
  x_labels = c("Frontal", "Posterior"),
  title = "C  Beta ITPC",
  subtitle = "Spatial profile",
  ylab = "ITPC (z)",
  data_range = phase_range
) + panel_final_theme


# D — Beta ISPC
p_D <- make_emm_panel(
  beta_ispc_emm,
  beta_ispc_test,
  xvar = "validity",
  x_levels = c("invalid", "valid"),
  x_labels = c("Invalid", "Valid"),
  title = "D  Beta ISPC",
  subtitle = "Validity modulation",
  ylab = "ISPC (z)",
  data_range = phase_range
) + panel_final_theme


## Conceptual header grobs
make_single_header_grob <- function(label, fill) {
  
  grobTree(
    rectGrob(
      x = .5,
      y = .5,
      width = 1,
      height = 1,
      gp = gpar(
        fill = fill,
        col = NA
      )
    ),
    textGrob(
      label,
      x = .5,
      y = .5,
      gp = gpar(
        fontsize = 8.2,
        fontface = "bold",
        col = "grey18",
        fontfamily = base_family
      )
    )
  )
}


# Two separate rectangles with a real white gap.
make_split_header_grob <- function(
    label_left,
    fill_left,
    label_right,
    fill_right,
    gap = .035) {
  
  w <- (1 - gap) / 2
  
  grobTree(
    
    # left rectangle
    rectGrob(
      x = w / 2,
      y = .5,
      width = w,
      height = 1,
      gp = gpar(
        fill = fill_left,
        col = NA
      )
    ),
    
    textGrob(
      label_left,
      x = w / 2,
      y = .5,
      gp = gpar(
        fontsize = 8.2,
        fontface = "bold",
        col = "grey18",
        fontfamily = base_family
      )
    ),
    
    # right rectangle
    rectGrob(
      x = 1 - w / 2,
      y = .5,
      width = w,
      height = 1,
      gp = gpar(
        fill = fill_right,
        col = NA
      )
    ),
    
    textGrob(
      label_right,
      x = 1 - w / 2,
      y = .5,
      gp = gpar(
        fontsize = 8.2,
        fontface = "bold",
        col = "grey18",
        fontfamily = base_family
      )
    )
  )
}


hdr_amplitude <- make_single_header_grob(
  "Local amplitude dynamics",
  "#F5F1EB"
)


hdr_phase <- make_split_header_grob(
  label_left  = "Local phase dynamics",
  fill_left   = "#EEF5F2",
  label_right = "Network phase dynamics",
  fill_right  = "#EEF3F7",
  gap = .040
)

make_forest_title <- function(label, fill) {
  grobTree(
    rectGrob(
      gp = gpar(fill = fill, col = NA)
    ),
    textGrob(
      label,
      x = unit(0.025, "npc"),
      y = unit(0.5, "npc"),
      just = "left",
      gp = gpar(
        fontsize = 8.2,
        fontface = "bold",
        fontfamily = base_family
      )
    ),
    textGrob(
      expression(p[adj]),
      x = unit(0.975, "npc"),
      y = unit(0.5, "npc"),
      just = "right",
      gp = gpar(
        fontsize = 8.2,
        fontface = "bold",
        fontfamily = base_family
      )
    )
  )
}

forest_title_alpha <- make_forest_title("Alpha power", "#F5F1EB")
forest_title_beta  <- make_forest_title("Beta power",  "#F5F1EB")
forest_title_itpc  <- make_forest_title("Beta ITPC",   "#EEF5F2")
forest_title_ispc  <- make_forest_title("Beta ISPC",   "#EEF3F7")


## left-side rows
row_alpha <- (
  p_A |
    p_alpha_spatial
) +
  plot_layout(
    widths = c(1, 1)
  )

row_beta <- (
  p_B |
    p_beta_contra
) +
  plot_layout(
    widths = c(1, 1)
  )

row_phase <- (
  p_C |
    p_D
) +
  plot_layout(
    widths = c(1, 1)
  )


## forest data helpers
forest_subheader_row <- function(label, section) {
  
  tibble(
    section = section,
    row_type = "subheader",
    effect_kind = NA_character_,
    label = label,
    estimate_plot = NA_real_,
    lower_plot = NA_real_,
    upper_plot = NA_real_,
    p_plot = NA_real_
  )
}


forest_between_rows <- function(
    test_object,
    xvar,
    x_levels,
    x_labels,
    section) {
  
  test_object$between %>%
    filter(
      .data[[xvar]] %in% x_levels
    ) %>%
    mutate(
      ordering = match(
        as.character(.data[[xvar]]),
        x_levels
      )
    ) %>%
    arrange(ordering) %>%
    transmute(
      section = section,
      row_type = "effect",
      effect_kind = "simple",
      label = x_labels[ordering],
      
      # original contrast = NT − ASD
      # forest direction = ASD − NT
      estimate_plot = -estimate,
      lower_plot = -upper.CL,
      upper_plot = -lower.CL,
      
      p_plot = p_FDR
    )
}


forest_interaction_row <- function(
    interaction_df,
    label,
    section) {
  
  interaction_df %>%
    slice(1) %>%
    transmute(
      section = section,
      row_type = "effect",
      effect_kind = "interaction",
      label = label,
      estimate_plot = estimate,
      lower_plot = lower.CL,
      upper_plot = upper.CL,
      p_plot = p_FDR
    )
}


prepare_forest <- function(df) {
  df %>%
    mutate(
      y = rev(seq_len(n()))
    )
}


## common forest spacing
forest_row_step     <- 0.82
forest_title_height <- 0.070

forest_y_pad_power <- 0.35

prepare_forest <- function(
    df,
    row_step = forest_row_step) {
  
  df %>%
    mutate(
      # Start at zero so padding below/above the rows is truly symmetric
      y = rev(seq(
        0,
        by = row_step,
        length.out = n()
      ))
    )
}


forest_center_span <- function(df) {
  diff(range(df$y, na.rm = TRUE))
}

## alpha forest
forest_alpha <- bind_rows(
  
  forest_subheader_row(
    "Validity modulation",
    "alpha"
  ),
  
  forest_between_rows(
    alpha_valid_test,
    xvar = "validity",
    x_levels = c("invalid", "valid"),
    x_labels = c("Invalid", "Valid"),
    section = "alpha"
  ),
  
  forest_interaction_row(
    alpha_valid_int,
    label = "\u0394 valid\u2212invalid",
    section = "alpha"
  ),
  
  forest_subheader_row(
    "Spatial profile",
    "alpha"
  ),
  
  forest_between_rows(
    alpha_roi_test,
    xvar = "ROI",
    x_levels = c("front", "occi"),
    x_labels = c("Frontal", "Posterior"),
    section = "alpha"
  ),
  
  forest_interaction_row(
    alpha_roi_int,
    label = "\u0394 posterior\u2212frontal",
    section = "alpha"
  )
  
) %>%
  prepare_forest()

# beta power forest
forest_beta <- bind_rows(
  
  forest_subheader_row(
    "Ipsilateral spatial profile",
    "beta"
  ),
  
  forest_between_rows(
    beta_ipsi_test,
    xvar = "ROI",
    x_levels = c("front", "occi"),
    x_labels = c("Frontal", "Posterior"),
    section = "beta"
  ),
  
  forest_interaction_row(
    beta_ipsi_int,
    label = "\u0394 posterior\u2212frontal",
    section = "beta"
  ),
  
  forest_subheader_row(
    "Contralateral spatial profile",
    "beta"
  ),
  
  forest_between_rows(
    beta_contra_test,
    xvar = "ROI",
    x_levels = c("front", "occi"),
    x_labels = c("Frontal", "Posterior"),
    section = "beta"
  ),
  
  forest_interaction_row(
    beta_contra_int,
    label = "\u0394 posterior\u2212frontal",
    section = "beta"
  )
  
) %>%
  prepare_forest()


# Beta ITPC forest
forest_itpc <- bind_rows(
  
  forest_between_rows(
    beta_itpc_test,
    xvar = "ROI",
    x_levels = c("front", "occi"),
    x_labels = c("Frontal", "Posterior"),
    section = "itpc"
  ),
  
  forest_interaction_row(
    beta_itpc_int,
    label = "\u0394 posterior\u2212frontal",
    section = "itpc"
  )
  
) %>%
  prepare_forest()


# Beta ISPC forest
forest_ispc <- bind_rows(
  
  forest_between_rows(
    beta_ispc_test,
    xvar = "validity",
    x_levels = c("invalid", "valid"),
    x_labels = c("Invalid", "Valid"),
    section = "ispc"
  ),
  
  forest_interaction_row(
    beta_ispc_int,
    label = "\u0394 valid\u2212invalid",
    section = "ispc"
  )
  
) %>%
  prepare_forest()

# Common forest scale + p format
forest_all_effects <- bind_rows(
  forest_alpha,
  forest_beta,
  forest_itpc,
  forest_ispc
) %>%
  filter(row_type == "effect")

forest_max <- max(
  abs(
    c(
      forest_all_effects$lower_plot,
      forest_all_effects$upper_plot
    )
  ),
  na.rm = TRUE
)

forest_xlim <- c(-1, 1) *
  max(
    .25,
    ceiling(forest_max * 10) / 10
  )

format_forest <- function(df) {
  df %>%
    mutate(
      p_text = ifelse(
        row_type == "effect",
        ifelse(
          p_plot < .001,
          "<.001",
          sub("^0", "", sprintf("%.3f", p_plot))
        ),
        ""
      )
    )
}

forest_alpha <- format_forest(forest_alpha)
forest_beta  <- format_forest(forest_beta)
forest_itpc  <- format_forest(forest_itpc)
forest_ispc  <- format_forest(forest_ispc)

# Forest body
forest_col_widths <- c(
  1.68,
  1.30,
  .42
)

make_forest_body <- function(df, y_pad) {
  
  ymin <- min(df$y, na.rm = TRUE) - y_pad
  ymax <- max(df$y, na.rm = TRUE) + y_pad
  
  effects <- df %>%
    filter(row_type == "effect")
  
  # label column
  p_labels <- ggplot() +
    
    geom_text(
      data = df %>% filter(row_type == "subheader"),
      aes(x = .02, y = y, label = label),
      hjust = 0,
      family = base_family,
      fontface = "italic",
      size = 2.85,
      colour = "grey42"
    ) +
    
    geom_text(
      data = df %>%
        filter(row_type == "effect",
               effect_kind == "simple"),
      aes(x = .045, y = y, label = label),
      hjust = 0,
      family = base_family,
      size = 2.85,
      colour = "grey18"
    ) +
    
    geom_text(
      data = df %>%
        filter(row_type == "effect",
               effect_kind == "interaction"),
      aes(x = .045, y = y, label = label),
      hjust = 0,
      family = base_family,
      fontface = "bold",
      size = 2.85,
      colour = "grey18"
    ) +
    
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(ymin, ymax),
      clip = "off"
    ) +
    
    theme_void(base_family = base_family) +
    theme(
      plot.margin = margin(0, 4, 0, 1)
    )
  
  # forest column
  p_forest <- ggplot() +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey76",
      linewidth = .30
    ) +
    
    geom_segment(
      data = effects,
      aes(
        x = lower_plot,
        xend = upper_plot,
        y = y,
        yend = y
      ),
      linewidth = .40,
      colour = "grey15"
    ) +
    
    geom_segment(
      data = effects,
      aes(
        x = lower_plot,
        xend = lower_plot,
        y = y - .08,
        yend = y + .08
      ),
      linewidth = .29,
      colour = "grey15"
    ) +
    
    geom_segment(
      data = effects,
      aes(
        x = upper_plot,
        xend = upper_plot,
        y = y - .08,
        yend = y + .08
      ),
      linewidth = .29,
      colour = "grey15"
    ) +
    
    geom_point(
      data = effects %>% filter(effect_kind == "simple"),
      aes(x = estimate_plot, y = y),
      shape = 16,
      size = 1.15,
      colour = "black"
    ) +
    
    geom_point(
      data = effects %>% filter(effect_kind == "interaction"),
      aes(x = estimate_plot, y = y),
      shape = 18,
      size = 1.35,
      colour = "black"
    ) +
    
    scale_x_continuous(
      limits = forest_xlim,
      breaks = pretty(forest_xlim, n = 3)
    ) +
    
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    
    labs(x = NULL, y = NULL) +
    theme_fig +
    
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 8.0),
      axis.ticks.x = element_line(
        colour = "grey40",
        linewidth = .27
      ),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 2, 0, 2)
    )
  
  # pFDR column
  p_p <- ggplot() +
    
    geom_text(
      data = effects,
      aes(
        x = .5,
        y = y,
        label = p_text
      ),
      hjust = .5,
      family = base_family,
      size = 2.85,
      colour = "grey38"
    ) +
    
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(ymin, ymax)
    ) +
    
    theme_void(base_family = base_family) +
    theme(
      plot.margin = margin(0, 1, 0, 1)
    )
  
  
  (
    p_labels |
      p_forest |
      p_p
  ) +
    plot_layout(
      widths = forest_col_widths
    )
}

# forest vertical alignment + shared axis label
+# Fraction of one complete result row occupied by the forest bodies
power_body_fraction <- 1 - forest_title_height

# Two title strips + two equally tall phase forest bodies
phase_body_fraction <- (
  1 - 2 * forest_title_height
) / 2

# Physical row spacing in the power forests
power_span <- mean(c(
  forest_center_span(forest_alpha),
  forest_center_span(forest_beta)
)) + 2 * forest_y_pad_power

common_row_gap <- (
  power_body_fraction *
    forest_row_step /
    power_span
)

# Calculate the padding needed in the short ITPC/ISPC forests
# to produce the same physical row spacing
phase_center_span <- mean(c(
  forest_center_span(forest_itpc),
  forest_center_span(forest_ispc)
))

phase_required_span <- (
  phase_body_fraction *
    forest_row_step /
    common_row_gap
)

forest_y_pad_phase <- max(
  0,
  (phase_required_span - phase_center_span) / 2
)


# Shared forest-axis label
forest_axis_label_plot <- ggplot() +
  annotate(
    "text",
    x = .68,
    y = .92,
    label = "Standardized effect (ASD \u2212 NT)",
    family = base_family,
    size = 2.85,
    colour = "grey25",
    hjust = .5
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 2, 0, 2)
  )

# complete forest plot
alpha_forest_plot <- (
  wrap_elements(full = forest_title_alpha) /
    wrap_elements(
      full = make_forest_body(
        forest_alpha,
        y_pad = forest_y_pad_power
      )
    )
) +
  plot_layout(
    heights = c(
      forest_title_height,
      power_body_fraction
    )
  )

beta_forest_plot <- (
  wrap_elements(full = forest_title_beta) /
    wrap_elements(
      full = make_forest_body(
        forest_beta,
        y_pad = forest_y_pad_power
      )
    )
) +
  plot_layout(
    heights = c(
      forest_title_height,
      power_body_fraction
    )
  )

phase_forest_plot <- (
  wrap_elements(full = forest_title_itpc) /
    wrap_elements(
      full = make_forest_body(
        forest_itpc,
        y_pad = forest_y_pad_phase
      )
    ) /
    wrap_elements(full = forest_title_ispc) /
    wrap_elements(
      full = make_forest_body(
        forest_ispc,
        y_pad = forest_y_pad_phase
      )
    )
) +
  plot_layout(
    heights = c(
      forest_title_height,
      phase_body_fraction,
      forest_title_height,
      phase_body_fraction
    )
  )

# Divider
divider <- grobTree(
  segmentsGrob(
    x0 = .5,
    x1 = .5,
    y0 = 0,
    y1 = 1,
    gp = gpar(
      col = "grey88",
      lwd = .6
    )
  )
)

# final layout
main_left_width <- 1.43
divider_width   <- .016
forest_width    <- 1.10

concept_header_height <- .070

amplitude_header_complete <- (
  wrap_elements(full = hdr_amplitude) |
    wrap_elements(full = divider) |
    plot_spacer()
) +
  plot_layout(
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

alpha_row_complete <- (
  wrap_elements(full = row_alpha) |
    wrap_elements(full = divider) |
    wrap_elements(full = alpha_forest_plot)
) +
  plot_layout(
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

beta_row_complete <- (
  wrap_elements(full = row_beta) |
    wrap_elements(full = divider) |
    wrap_elements(full = beta_forest_plot)
) +
  plot_layout(
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

phase_header_complete <- (
  wrap_elements(full = hdr_phase) |
    wrap_elements(full = divider) |
    plot_spacer()
) +
  plot_layout(
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

phase_row_complete <- (
  wrap_elements(full = row_phase) |
    wrap_elements(full = divider) |
    wrap_elements(full = phase_forest_plot)
) +
  plot_layout(
    nrow = 1,
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

# legend
legend_df <- tibble(
  group = factor(
    c("NT", "ASD"),
    levels = c("NT", "ASD")
  ),
  x = c(1, 2),
  y = 1
)

legend_plot <- ggplot(
  legend_df,
  aes(
    x = x,
    y = y,
    colour = group
  )
) +
  
  geom_point(
    size = 2.15,
    alpha = 0
  ) +
  
  scale_colour_manual(
    values = cols,
    labels = c(
      "Neurotypical (NT)",
      "Autistic (ASD)"
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = 1,
        size = 2.15
      )
    )
  ) +
  
  theme_void(
    base_family = base_family
  ) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 8.0
    )
  )

legend_row <- (
  wrap_elements(full = legend_plot) |
    plot_spacer() |
    wrap_elements(full = forest_axis_label_plot)
) +
  plot_layout(
    widths = c(
      main_left_width,
      divider_width,
      forest_width
    )
  )

## complete figure
fig_group_dashboard <- (
  amplitude_header_complete /
    alpha_row_complete /
    beta_row_complete /
    phase_header_complete /
    phase_row_complete /
    legend_row
) +
  plot_layout(
    heights = c(
      concept_header_height,
      1,
      1,
      concept_header_height,
      1,
      .075
    )
  )

fig_group_dashboard

# 22) SAVE
output_dir <- paste0(
  "...",
  ".../Tables and Figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


pdf_file <- file.path(
  output_dir,
  "Figure_group_differences_dashboard.pdf"
)


tiff_file <- file.path(
  output_dir,
  "Figure_group_differences_dashboard.tiff"
)

tmp_pdf <- tempfile(
  fileext = ".pdf"
)

tmp_tiff <- tempfile(
  fileext = ".tiff"
)

fig_width  <- 190 / 25.4

fig_height <- 7.65

ggsave(
  filename = tmp_pdf,
  plot = fig_group_dashboard,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)


ggsave(
  filename = tmp_tiff,
  plot = fig_group_dashboard,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 1000,
  compression = "lzw",
  bg = "white"
)


ok_pdf <- file.copy(
  tmp_pdf,
  pdf_file,
  overwrite = TRUE
)


ok_tiff <- file.copy(
  tmp_tiff,
  tiff_file,
  overwrite = TRUE
)


if (!ok_pdf) {
  warning(
    "PDF could not be copied to the network folder."
  )
}

if (!ok_tiff) {
  warning(
    "TIFF could not be copied to the network folder."
  )
}

if (ok_pdf) {
  
  shell.exec(
    normalizePath(
      pdf_file,
      winslash = "\\",
      mustWork = TRUE
    )
  )
}

# Export PNG 
png_file <- file.path(
  output_dir,
  "Figure_group_differences_dashboard.png"
)

ggsave(
  filename = png_file,
  plot = fig_group_dashboard,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 600,
  bg = "white"
)


# 3 ASSOCIATION OSCILLATIONS - RT  ######################################################################################################################################################

#### 3.1 PREDICTION: RT - POWER -----------------------------------------------------------------------------------------------------------------------
# Trial-level RT – power 

## Prepare trial-level data
power_trial <- trial_full_df %>%
  filter(
    band %in% c("alpha", "beta"),
    hemi %in% c("contra", "ipsi"),
    ROI %in% c("front", "occi")
  ) %>%
  select(
    id, group, sex, age, IQ, validity, band, hemi,
    epoch_idx, ROI, dB_total, RT_ms
  ) %>%
  distinct() %>%
  mutate(ROI = recode(as.character(ROI), occi = "post")) %>%
  pivot_wider(
    names_from = ROI,
    values_from = dB_total,
    names_prefix = "dB_"
  ) %>%
  mutate(
    group = factor(group, levels = c("ASD", "NT")),
    validity = factor(validity, levels = c("invalid", "valid")),
    band = factor(band, levels = c("alpha", "beta")),
    hemi = factor(hemi, levels = c("contra", "ipsi")),
    sex = factor(sex),
    id = factor(id)
  ) %>%
  drop_na(RT_ms, dB_post, dB_front, age, IQ, sex) %>%
  group_by(band, hemi) %>%
  mutate(
    RT_z = as.numeric(scale(RT_ms)),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ)),
    dB_front_z = as.numeric(scale(dB_front)),
    dB_post_z = as.numeric(scale(dB_post))
  ) %>%
  ungroup() %>%
  droplevels()

###### Models ####
m1_alpha_contra_trial <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "alpha" &
      power_trial$hemi == "contra", ]
)

m1_alpha_ipsi_trial <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "alpha" &
      power_trial$hemi == "ipsi", ]
)

m1_beta_contra_trial <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "beta" &
      power_trial$hemi == "contra", ]
)

m1_beta_ipsi_trial <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "beta" &
      power_trial$hemi == "ipsi", ]
)

## Extract omnibus tests
power_RT_trial_effects <- c(
  "validity", "dB_post_z", "dB_front_z",
  "age_z", "IQ_z", "sex",
  "validity:dB_post_z",
  "validity:dB_front_z"
)

power_RT_trial_models <- list(
  alpha_contra = m1_alpha_contra_trial,
  alpha_ipsi   = m1_alpha_ipsi_trial,
  beta_contra  = m1_beta_contra_trial,
  beta_ipsi    = m1_beta_ipsi_trial
)

all_power_RT_trial_tests <- imap_dfr(
  power_RT_trial_models,
  function(model, model_name) {
    anova(model, type = 3, ddf = "Satterthwaite") %>%
      as.data.frame() %>%
      rownames_to_column("effect") %>%
      transmute(
        model = model_name,
        band = sub("_.*", "", model_name),
        hemi = sub(".*_", "", model_name),
        effect,
        df_num = NumDF,
        df_den = DenDF,
        F_value = `F value`,
        p_raw = `Pr(>F)`
      ) %>%
      filter(effect %in% power_RT_trial_effects)
  }
)

## FDR correction across the four models
all_power_RT_trial_tests_FDR <- all_power_RT_trial_tests %>%
  group_by(effect) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    band = factor(band, levels = c("alpha", "beta")),
    hemi = factor(hemi, levels = c("contra", "ipsi")),
    effect = factor(effect, levels = power_RT_trial_effects),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(band, hemi, effect)

## Format table
fmt_p <- function(x) {
  ifelse(x < .001, "< .001", sub("^0", "", sprintf("%.3f", x)))
}

power_RT_trial_table <- all_power_RT_trial_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Hemisphere = recode(
      as.character(hemi),
      contra = "Contralateral",
      ipsi = "Ipsilateral"
    ),
    Effect = recode(
      as.character(effect),
      validity = "Validity",
      dB_post_z = "Posterior power",
      dB_front_z = "Frontal power",
      age_z = "Age",
      IQ_z = "IQ",
      sex = "Sex",
      `validity:dB_post_z` = "Validity × posterior power",
      `validity:dB_front_z` = "Validity × frontal power"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value, df_num, df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(
  as.data.frame(power_RT_trial_table),
  row.names = FALSE,
  right = TRUE
)

## Export to Word
output_file_power_RT_trial <- paste0(
  ".../Tables and Figures/",
  "Table_Power_RT_TrialLevel_Omnibus_Tests_FDR.docx"
)

save_as_docx(
  "Trial-level power–RT omnibus tests" =
    flextable(power_RT_trial_table) %>% autofit(),
  path = output_file_power_RT_trial
)


###### Post-hoc RT - Power ####
## Alpha ipsilateral posterior power: slopes within validity
alpha_ipsi_post <- emtrends(
  m1_alpha_ipsi_trial,
  ~ validity,
  var = "dB_post_z",
  lmer.df = "asymptotic"
)

alpha_ipsi_post_df <- summary(
  alpha_ipsi_post,
  infer = c(TRUE, TRUE),
  adjust = "none"
) %>%
  as.data.frame() %>%
  rename(estimate = dB_post_z.trend) %>%
  mutate(p_FDR = p.adjust(p.value, method = "BH"))

alpha_ipsi_post_df

## Difference between valid and invalid slopes
alpha_ipsi_post_interaction <- contrast(
  alpha_ipsi_post,
  method = "revpairwise",
  adjust = "none"
)

summary(
  alpha_ipsi_post_interaction,
  infer = c(TRUE, TRUE)
)


## Beta contralateral posterior power: overall slope
beta_contra_post_df <- emtrends(
  m1_beta_contra_trial,
  ~ 1,
  var = "dB_post_z",
  weights = "equal",
  lmer.df = "asymptotic"
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as.data.frame() %>%
  rename(estimate = dB_post_z.trend) %>%
  mutate(p_FDR = all_power_RT_trial_tests_FDR %>%
           filter(band == "beta", hemi == "contra",
                  as.character(effect) == "dB_post_z") %>%
           pull(p_FDR))

## Beta ipsilateral posterior power: overall slope
beta_ipsi_post_df <- emtrends(
  m1_beta_ipsi_trial,
  ~ 1,
  var = "dB_post_z",
  weights = "equal",
  lmer.df = "asymptotic"
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as.data.frame() %>%
  rename(estimate = dB_post_z.trend) %>%
  mutate(p_FDR = all_power_RT_trial_tests_FDR %>%
           filter(band == "beta", hemi == "ipsi",
                  as.character(effect) == "dB_post_z") %>%
           pull(p_FDR))

beta_contra_post_df
beta_ipsi_post_df

###### exploratory: interaction effect with group added only for neural effects that were significant in across-group models ####
## Alpha ipsilateral posterior power:
## test whether the validity-dependent association differs by group
m_alpha_ipsi_group <- lmer(
  RT_z ~ group * validity * dB_post_z +
    validity * dB_front_z + age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "alpha" & power_trial$hemi == "ipsi", ]
)

anova(m_alpha_ipsi_group, type = 3, ddf = "Satterthwaite")

## Beta contralateral posterior power:
## test whether the overall association differs by group
m_beta_contra_group <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    group * dB_post_z + age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "beta" & power_trial$hemi == "contra", ]
)

anova(m_beta_contra_group, type = 3, ddf = "Satterthwaite")

## Beta ipsilateral posterior power:
## test whether the overall association differs by group
m_beta_ipsi_group <- lmer(
  RT_z ~ validity * (dB_post_z + dB_front_z) +
    group * dB_post_z + age_z + IQ_z + sex + (1 | id),
  data = power_trial[
    power_trial$band == "beta" & power_trial$hemi == "ipsi", ]
)

anova(m_beta_ipsi_group, type = 3, ddf = "Satterthwaite")



#### 3.2 Prediction: RT - Lateralization (LI) -----------------------------------------------------------------------------
# Trial-level RT – lateralization

## Prepare trial-level data
li_trial <- targ_li_df %>%
  filter(band %in% c("alpha", "beta"), ROI %in% c("front", "occi")) %>%
  select(id, group, sex, age, IQ, validity, band, epoch_idx,
         ROI, LI_ipsi_minus_contra, RT_ms) %>%
  distinct() %>%
  mutate(ROI = recode(as.character(ROI), occi = "post")) %>%
  pivot_wider(names_from = ROI, values_from = LI_ipsi_minus_contra,
              names_prefix = "LI_") %>%
  mutate(
    validity = factor(validity, levels = c("invalid", "valid")),
    band = factor(band, levels = c("alpha", "beta")),
    group = factor(group, levels = c("ASD", "NT")),
    sex = factor(sex), id = factor(id)
  ) %>%
  drop_na(RT_ms, LI_post, LI_front, age, IQ, sex) %>%
  group_by(band) %>%
  mutate(
    RT_z = as.numeric(scale(RT_ms)),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ)),
    LI_post_z = as.numeric(scale(LI_post)),
    LI_front_z = as.numeric(scale(LI_front))
  ) %>%
  ungroup() %>%
  droplevels()

###### Models ####
m1_LI_alpha_trial <- lmer(
  RT_z ~ validity * (LI_post_z + LI_front_z) + age_z + IQ_z + sex + (1 | id),
  data = li_trial[li_trial$band == "alpha", ]
)

m1_LI_beta_trial <- lmer(
  RT_z ~ validity * (LI_post_z + LI_front_z) + age_z + IQ_z + sex + (1 | id),
  data = li_trial[li_trial$band == "beta", ]
)

## Extract omnibus tests
LI_RT_trial_effects <- c(
  "validity", "LI_post_z", "LI_front_z", "age_z", "IQ_z", "sex",
  "validity:LI_post_z", "validity:LI_front_z"
)

LI_RT_trial_models <- list(
  alpha = m1_LI_alpha_trial,
  beta = m1_LI_beta_trial
)

all_LI_RT_trial_tests <- imap_dfr(
  LI_RT_trial_models,
  function(model, band_name) {
    anova(model, type = 3, ddf = "Satterthwaite") %>%
      as.data.frame() %>%
      rownames_to_column("effect") %>%
      transmute(
        band = band_name, effect,
        df_num = NumDF, df_den = DenDF,
        F_value = `F value`, p_raw = `Pr(>F)`
      ) %>%
      filter(effect %in% LI_RT_trial_effects)
  }
)

## FDR correction across alpha and beta
all_LI_RT_trial_tests_FDR <- all_LI_RT_trial_tests %>%
  group_by(effect) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    band = factor(band, levels = c("alpha", "beta")),
    effect = factor(effect, levels = LI_RT_trial_effects),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(band, effect)

## Format table
LI_RT_trial_table <- all_LI_RT_trial_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      validity = "Validity",
      LI_post_z = "Posterior lateralization",
      LI_front_z = "Frontal lateralization",
      age_z = "Age", IQ_z = "IQ", sex = "Sex",
      `validity:LI_post_z` = "Validity × posterior lateralization",
      `validity:LI_front_z` = "Validity × frontal lateralization"
    ),
    `F(df1, df2)` = sprintf("%.2f (%.0f, %.2f)", F_value, df_num, df_den),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(as.data.frame(LI_RT_trial_table), row.names = FALSE, right = TRUE)

## Export
output_file_LI_RT_trial <- paste0(
  ".../Tables and Figures/",
  "Table_Lateralization_RT_TrialLevel_Omnibus_Tests_FDR.docx"
)

save_as_docx(
  "Trial-level lateralization–RT omnibus tests" =
    flextable(LI_RT_trial_table) %>% autofit(),
  path = output_file_LI_RT_trial
)

###### Post-hoc RT ~ lateralization ####
## Alpha posterior lateralization
LI_alpha_post <- emtrends(
  m1_LI_alpha_trial,
  ~ 1,
  var = "LI_post_z",
  weights = "equal",
  lmer.df = "asymptotic"
)

LI_alpha_post_df <- summary(
  LI_alpha_post,
  infer = c(TRUE, TRUE),
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    Band = "Alpha",
    Effect = "Posterior lateralization",
    p_FDR = all_LI_RT_trial_tests_FDR %>%
      filter(band == "alpha", effect == "LI_post_z") %>%
      pull(p_FDR)
  )

## Beta frontal lateralization
LI_beta_front <- emtrends(
  m1_LI_beta_trial,
  ~ 1,
  var = "LI_front_z",
  weights = "equal",
  lmer.df = "asymptotic"
)

LI_beta_front_df <- summary(
  LI_beta_front,
  infer = c(TRUE, TRUE),
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    Band = "Beta",
    Effect = "Frontal lateralization",
    p_FDR = all_LI_RT_trial_tests_FDR %>%
      filter(band == "beta", effect == "LI_front_z") %>%
      pull(p_FDR)
  )

LI_significant_slopes <- bind_rows(
  LI_alpha_post_df,
  LI_beta_front_df
)

LI_significant_slopes


###### Exploratory group moderation: RT ~ lateralization ####

## Alpha posterior lateralization
m_LI_alpha_group <- lmer(
  RT_z ~ validity * (LI_post_z + LI_front_z) +
    group * LI_post_z + age_z + IQ_z + sex + (1 | id),
  data = li_trial[li_trial$band == "alpha", ]
)

anova(m_LI_alpha_group, type = 3, ddf = "Satterthwaite")

## Beta frontal lateralization
m_LI_beta_group <- lmer(
  RT_z ~ validity * (LI_post_z + LI_front_z) +
    group * LI_front_z + age_z + IQ_z + sex + (1 | id),
  data = li_trial[li_trial$band == "beta", ]
)

anova(m_LI_beta_group, type = 3, ddf = "Satterthwaite")


#### 3.3 RT - ITPC ----------------------------------------------------------------------------
## ITPC preparation
ITPC_n <- targ_itc_matchedN %>%
  filter(hemi == "contra", ROI %in% c("front", "occi")) %>%
  select(id, validity, band, ROI, ITC_mean, RT_mean, age, IQ, sex, group) %>%
  mutate(ROI = recode(ROI, occi = "post")) %>%
  pivot_wider(names_from = ROI, values_from = ITC_mean, names_prefix = "ITPC_") %>%
  group_by(band) %>%
  mutate(
    RT_mean_z = as.numeric(scale(RT_mean)),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ)),
    ITPC_front_z = as.numeric(scale(ITPC_front)),
    ITPC_post_z = as.numeric(scale(ITPC_post))
  ) %>%
  ungroup()

###### Models ####
m1_ITPC_alpha_contra <- lmer(
  RT_mean_z ~ validity * (ITPC_post_z + ITPC_front_z) + age_z + IQ_z + sex + (1 | id),
  data = ITPC_n[ITPC_n$band == "alpha", ]
)

m1_ITPC_beta_contra <- lmer(
  RT_mean_z ~ validity * (ITPC_post_z + ITPC_front_z) + age_z + IQ_z + sex + (1 | id),
  data = ITPC_n[ITPC_n$band == "beta", ]
)

## Omnibus tests and FDR
ITPC_RT_effects <- c(
  "validity", "ITPC_post_z", "ITPC_front_z",
  "age_z", "IQ_z", "sex",
  "validity:ITPC_post_z", "validity:ITPC_front_z"
)

ITPC_RT_models <- list(
  alpha = m1_ITPC_alpha_contra,
  beta  = m1_ITPC_beta_contra
)

all_ITPC_RT_tests <- imap_dfr(ITPC_RT_models, function(model, band_name) {
  anova(model, type = 3, ddf = "Satterthwaite") %>%
    as.data.frame() %>%
    rownames_to_column("effect") %>%
    transmute(
      band = band_name, effect,
      df_num = NumDF, df_den = DenDF,
      F_value = `F value`, p_raw = `Pr(>F)`
    ) %>%
    filter(effect %in% ITPC_RT_effects)
})

all_ITPC_RT_tests_FDR <- all_ITPC_RT_tests %>%
  group_by(effect) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    band = factor(band, levels = c("alpha", "beta")),
    effect = factor(effect, levels = ITPC_RT_effects)
  ) %>%
  arrange(band, effect)

## Format and export
ITPC_RT_table <- all_ITPC_RT_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      validity = "Validity",
      ITPC_post_z = "Posterior ITPC",
      ITPC_front_z = "Frontal ITPC",
      age_z = "Age", IQ_z = "IQ", sex = "Sex",
      `validity:ITPC_post_z` = "Validity × posterior ITPC",
      `validity:ITPC_front_z` = "Validity × frontal ITPC"
    ),
    `F(df1, df2)` = sprintf("%.2f (%.0f, %.2f)", F_value, df_num, df_den),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(as.data.frame(ITPC_RT_table), row.names = FALSE, right = TRUE)

output_file_ITPC_RT <- paste0(
  ".../Tables and Figures/",
  "Table_ITPC_RT_Omnibus_Tests_FDR.docx"
)

save_as_docx(
  "ITPC–RT omnibus tests" = flextable(ITPC_RT_table) %>% autofit(),
  path = output_file_ITPC_RT
)


###### Post-hoc RT - ITPC ####
front_alpha_trend <- emtrends(
  m1_ITPC_alpha_contra,
  ~ 1,
  var = "ITPC_front_z",
  weights = "equal",
  lmer.df = "asymptotic"
)

front_alpha_df <- summary(
  front_alpha_trend,
  infer = c(TRUE, TRUE),
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    p_FDR = all_ITPC_RT_tests_FDR %>%
      filter(band == "alpha", effect == "ITPC_front_z") %>%
      pull(p_FDR)
  )

front_alpha_df

###### Exploratory group moderation ####
m_ITPC_alpha_group <- lmer(
  RT_mean_z ~ group * ITPC_front_z +
    validity * (ITPC_post_z + ITPC_front_z) +
    age_z + IQ_z + sex + (1 | id),
  data = ITPC_n[ITPC_n$band == "alpha", ]
)

anova(m_ITPC_alpha_group, type = 3, ddf = "Satterthwaite")


#### 3.4 RT - power connectivity ----------------------------------------------------------------------------

## Prepare data
df_powconnect_z <- df_powconnect %>%
  filter(
    band %in% c("alpha", "beta"),
    ROI_interaction == "front.occi",
    hemi == "contra"
  ) %>%
  mutate(
    validity = factor(validity, levels = c("invalid", "valid")),
    band = factor(band, levels = c("alpha", "beta")),
    sex = factor(sex),
    id = factor(id)
  ) %>%
  group_by(band) %>%
  mutate(
    pow_cor_z = as.numeric(scale(pow_cor)),
    RT_mean_z = as.numeric(scale(RT_mean)),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ))
  ) %>%
  ungroup() %>%
  droplevels()

###### Models ####
m1_powcon_alpha <- lmer(
  RT_mean_z ~ validity * pow_cor_z + age_z + IQ_z + sex + (1 | id),
  data = df_powconnect_z[df_powconnect_z$band == "alpha", ]
)

m1_powcon_beta <- lmer(
  RT_mean_z ~ validity * pow_cor_z + age_z + IQ_z + sex + (1 | id),
  data = df_powconnect_z[df_powconnect_z$band == "beta", ]
)

## Extract omnibus tests
powcon_RT_effects <- c(
  "validity", "pow_cor_z", "age_z", "IQ_z", "sex",
  "validity:pow_cor_z"
)

powcon_RT_models <- list(
  alpha = m1_powcon_alpha,
  beta = m1_powcon_beta
)

all_powcon_RT_tests <- imap_dfr(powcon_RT_models, function(model, band_name) {
  anova(model, type = 3, ddf = "Satterthwaite") %>%
    as.data.frame() %>%
    rownames_to_column("effect") %>%
    transmute(
      band = band_name,
      effect,
      df_num = NumDF,
      df_den = DenDF,
      F_value = `F value`,
      p_raw = `Pr(>F)`
    ) %>%
    filter(effect %in% powcon_RT_effects)
})

## FDR correction across alpha and beta
all_powcon_RT_tests_FDR <- all_powcon_RT_tests %>%
  group_by(effect) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    band = factor(band, levels = c("alpha", "beta")),
    effect = factor(effect, levels = powcon_RT_effects),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(band, effect)

## Format table
fmt_p <- function(x) {
  ifelse(x < .001, "< .001", sub("^0", "", sprintf("%.3f", x)))
}

powcon_RT_table <- all_powcon_RT_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      validity = "Validity",
      pow_cor_z = "Power connectivity",
      age_z = "Age",
      IQ_z = "IQ",
      sex = "Sex",
      `validity:pow_cor_z` = "Validity × power connectivity"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value, df_num, df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(as.data.frame(powcon_RT_table), row.names = FALSE, right = TRUE)

## Export to Word
output_file_powcon_RT <- paste0(
  ".../Tables and Figures/",
  "Table_PowerConnectivity_RT_Omnibus_Tests_FDR.docx"
)

save_as_docx(
  "Power-connectivity–RT omnibus tests" =
    flextable(powcon_RT_table) %>% autofit(),
  path = output_file_powcon_RT
)  



#### 3.5 RT - phase connectivity ----------------------------------------------------------------------------

## Prepare data
ISPC_matchedN_z <- targ_ispc_matchedN %>%
  filter(
    band %in% c("alpha", "beta"),
    roi_pair == "occi_front",
    hemi == "contra"
  ) %>%
  mutate(
    validity = factor(validity, levels = c("invalid", "valid")),
    band = factor(band, levels = c("alpha", "beta")),
    sex = factor(sex),
    id = factor(id)
  ) %>%
  group_by(band) %>%
  mutate(
    ISPC_mean_z = as.numeric(scale(ISPC_mean)),
    RT_mean_z = as.numeric(scale(RT_mean)),
    age_z = as.numeric(scale(age)),
    IQ_z = as.numeric(scale(IQ))
  ) %>%
  ungroup() %>%
  droplevels()

###### Models ####
m1_ISPC_alpha <- lmer(
  RT_mean_z ~ validity * ISPC_mean_z + age_z + IQ_z + sex + (1 | id),
  data = ISPC_matchedN_z[ISPC_matchedN_z$band == "alpha", ]
)

m1_ISPC_beta <- lmer(
  RT_mean_z ~ validity * ISPC_mean_z + age_z + IQ_z + sex + (1 | id),
  data = ISPC_matchedN_z[ISPC_matchedN_z$band == "beta", ]
)

## Extract omnibus tests
ISPC_RT_effects <- c(
  "validity", "ISPC_mean_z", "age_z", "IQ_z", "sex",
  "validity:ISPC_mean_z"
)

ISPC_RT_models <- list(
  alpha = m1_ISPC_alpha,
  beta = m1_ISPC_beta
)

all_ISPC_RT_tests <- imap_dfr(ISPC_RT_models, function(model, band_name) {
  anova(model, type = 3, ddf = "Satterthwaite") %>%
    as.data.frame() %>%
    rownames_to_column("effect") %>%
    transmute(
      band = band_name,
      effect,
      df_num = NumDF,
      df_den = DenDF,
      F_value = `F value`,
      p_raw = `Pr(>F)`
    ) %>%
    filter(effect %in% ISPC_RT_effects)
})

## FDR correction across alpha and beta
all_ISPC_RT_tests_FDR <- all_ISPC_RT_tests %>%
  group_by(effect) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    band = factor(band, levels = c("alpha", "beta")),
    effect = factor(effect, levels = ISPC_RT_effects),
    significant_FDR = p_FDR < .05
  ) %>%
  arrange(band, effect)

## Format table
ISPC_RT_table <- all_ISPC_RT_tests_FDR %>%
  transmute(
    Band = tools::toTitleCase(as.character(band)),
    Effect = recode(
      as.character(effect),
      validity = "Validity",
      ISPC_mean_z = "Phase connectivity",
      age_z = "Age",
      IQ_z = "IQ",
      sex = "Sex",
      `validity:ISPC_mean_z` = "Validity × phase connectivity"
    ),
    `F(df1, df2)` = sprintf(
      "%.2f (%.0f, %.2f)",
      F_value, df_num, df_den
    ),
    p = fmt_p(p_raw),
    `p FDR` = fmt_p(p_FDR)
  )

print(as.data.frame(ISPC_RT_table), row.names = FALSE, right = TRUE)

## Export to Word
output_file_ISPC_RT <- paste0(
  ".../Tables and Figures/",
  "Table_PhaseConnectivity_RT_Omnibus_Tests_FDR.docx"
)

save_as_docx(
  "Phase-connectivity–RT omnibus tests" =
    flextable(ISPC_RT_table) %>% autofit(),
  path = output_file_ISPC_RT
)


# 4 ASSOCIATIONS: Neural markers ###########################################################################################################################################
## Test correlations between neural markers that were examined in group comparisons between ASD and NT 
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)

#### 4.1 General settings ---------------------------------------------------------------------

n_permutations_global <- 50000
n_permutations_direct <- 50000
permutation_seed <- 1234

#### 4.2 Define neural measures ----------------------------------------
power_base <- expand.grid(
  band = c("alpha", "beta"),
  ROI = c("front", "occi"),
  hemi = c("ipsi", "contra"),
  stringsAsFactors = FALSE
) %>%
  transmute(base = paste("power", band, ROI, hemi, sep = "_")) %>%
  pull(base)

# ITPC
itpc_base <- expand.grid(
  band = c("alpha", "beta"),
  ROI = c("front", "occi"),
  stringsAsFactors = FALSE
) %>%
  transmute(base = paste("ITPC", band, ROI, "contra", sep = "_")) %>%
  pull(base)

# posterior–frontal contralateral connectivity
find_connectivity_base <- function(df, prefix, band) {
  candidates <- names(df)[
    grepl(paste0("^", prefix, "_", band, "_"), names(df)) &
      grepl("front", names(df)) &
      grepl("occi", names(df)) &
      grepl("_contra_valid$", names(df))
  ]
  
  if (length(candidates) != 1) {
    stop(
      "Expected exactly one ", prefix, " ", band,
      " posterior–frontal variable, but found: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  sub("_valid$", "", candidates)
}

powconn_base <- vapply(
  c("alpha", "beta"),
  \(band) find_connectivity_base(shared_neural_df, "powconn", band),
  character(1)
)

ispc_base <- vapply(
  c("alpha", "beta"),
  \(band) find_connectivity_base(shared_neural_df, "ispc", band),
  character(1)
)

feature_base <- unname(c(power_base, itpc_base, powconn_base, ispc_base))

# Use the original condition-specific variables
vars_valid <- paste0(feature_base, "_valid")
vars_invalid <- paste0(feature_base, "_invalid")

missing_vars <- setdiff(c(vars_valid, vars_invalid), names(shared_neural_df))

if (length(missing_vars) > 0) {
  stop("The following required variables are missing:\n", paste(missing_vars, collapse = "\n"))
}

length(vars_valid)    
length(vars_invalid)  

#### 4.3 Correlation matrices, p-values, sample sizes and within-matrix FDR ------------------
# Pearson correlations are computed using pairwise-complete observations.
# Benjamini–Hochberg FDR correction is applied across the 120 unique correlations within each
# group and condition
correlation_results <- function(df, vars) {
  x <- df[, vars, drop = FALSE]
  
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("All selected neural variables must be numeric.")
  }
  
  k <- length(vars)
  r_mat <- cor(x, use = "pairwise.complete.obs", method = "pearson")
  p_mat <- matrix(NA_real_, k, k, dimnames = list(vars, vars))
  n_mat <- matrix(NA_integer_, k, k, dimnames = list(vars, vars))
  
  for (i in seq_len(k)) {
    for (j in i:k) {
      keep <- complete.cases(x[[i]], x[[j]])
      n_pair <- sum(keep)
      
      n_mat[i, j] <- n_mat[j, i] <- n_pair
      
      if (i == j) {
        p_mat[i, j] <- 0
      } else if (n_pair >= 3 && sd(x[[i]][keep]) > 0 && sd(x[[j]][keep]) > 0) {
        test <- suppressWarnings(cor.test(x[[i]][keep], x[[j]][keep], method = "pearson"))
        p_mat[i, j] <- p_mat[j, i] <- test$p.value
      }
    }
  }
  
  upper <- upper.tri(p_mat)
  p_fdr <- matrix(NA_real_, k, k, dimnames = list(vars, vars))
  p_fdr[upper] <- p.adjust(p_mat[upper], method = "BH")
  p_fdr[lower.tri(p_fdr)] <- t(p_fdr)[lower.tri(p_fdr)]
  diag(p_fdr) <- 0
  
  stars <- matrix("", k, k, dimnames = list(vars, vars))
  stars[!is.na(p_fdr) & p_fdr < .05] <- "*"
  stars[!is.na(p_fdr) & p_fdr < .01] <- "**"
  stars[!is.na(p_fdr) & p_fdr < .001] <- "***"
  diag(stars) <- ""
  
  list(r = r_mat, p = p_mat, p_fdr = p_fdr, n = n_mat, stars = stars)
}


#### 4.4 Pooled dendrogram within each condition ----------------------------------------------------
# For each condition, one dendrogram is calculated from the pooled ASD and NT correlation matrix.
# That same order is then applied to both group matrices
make_pooled_dendrogram <- function(df, vars) {
  pooled_r <- cor(df[, vars, drop = FALSE], use = "pairwise.complete.obs", method = "pearson")
  
  if (anyNA(pooled_r)) {
    stop("The pooled correlation matrix contains missing correlations.")
  }
  
  distance <- as.dist(sqrt(2 * (1 - pooled_r)))
  hclust(distance, method = "average")
}


#### 4.5 valid and invalid matrices -------------------------------------------------
prepare_condition_results <- function(df, vars) {
  dat <- df %>% filter(group %in% c("ASD", "NT"))
  
  list(
    NT = correlation_results(dat %>% filter(group == "NT"), vars),
    ASD = correlation_results(dat %>% filter(group == "ASD"), vars),
    hc = make_pooled_dendrogram(dat, vars)
  )
}

res_valid <- prepare_condition_results(shared_neural_df, vars_valid)
res_invalid <- prepare_condition_results(shared_neural_df, vars_invalid)

#### 4.6 Descriptive valid-invalid comaprison within each group -----------------------------
# describe how closely each group's valid correlation architecture resembles its
# invalid architecture. Formal group comparison of reconfiguration is performed later

compare_valid_invalid_matrices <- function(valid_result, invalid_result, group_name) {
  valid_r <- valid_result$r[upper.tri(valid_result$r)]
  invalid_r <- invalid_result$r[upper.tri(invalid_result$r)]
  keep <- complete.cases(valid_r, invalid_r)
  
  tibble(
    group = group_name,
    matrix_similarity = cor(valid_r[keep], invalid_r[keep], method = "pearson"),
    RMSE_valid_vs_invalid = sqrt(mean((valid_r[keep] - invalid_r[keep])^2)),
    mean_abs_r_valid = mean(abs(valid_r[keep])),
    mean_abs_r_invalid = mean(abs(invalid_r[keep])),
    mean_signed_r_valid = mean(valid_r[keep]),
    mean_signed_r_invalid = mean(invalid_r[keep]),
    n_FDR_pairs_valid = sum(valid_result$p_fdr[upper.tri(valid_result$p_fdr)] < .05, na.rm = TRUE),
    n_FDR_pairs_invalid = sum(invalid_result$p_fdr[upper.tri(invalid_result$p_fdr)] < .05, na.rm = TRUE)
  )
}

condition_comparison <- bind_rows(
  compare_valid_invalid_matrices(res_valid$NT, res_invalid$NT, "Neurotypical"),
  compare_valid_invalid_matrices(res_valid$ASD, res_invalid$ASD, "Autistic")
)

print(condition_comparison, n = Inf)


#### 4.7 ASD–NT matrix correspondence -----------------------------------------------
# Matrix similarity is Pearson correlation between corresponding ASD and NT coefficients
# Global RMSE quantifies the overall magnitude of the discrepancy between the two matrices

upper_values <- function(mat) mat[upper.tri(mat)]

global_matrix_test <- function(df, vars, analysis_name, n_permutations = 5000, seed = 1234) {
  dat <- df %>% filter(group %in% c("ASD", "NT"))
  get_cor <- function(rows) cor(dat[rows, vars, drop = FALSE], use = "pairwise.complete.obs", method = "pearson")
  
  rows_ASD <- which(dat$group == "ASD")
  rows_NT <- which(dat$group == "NT")
  
  r_ASD <- get_cor(rows_ASD)
  r_NT <- get_cor(rows_NT)
  
  vec_ASD <- upper_values(r_ASD)
  vec_NT <- upper_values(r_NT)
  keep <- complete.cases(vec_ASD, vec_NT)
  
  matrix_similarity <- cor(vec_ASD[keep], vec_NT[keep], method = "pearson")
  observed_rmse <- sqrt(mean((vec_ASD[keep] - vec_NT[keep])^2))
  
  set.seed(seed)
  group_original <- dat$group
  permutation_rmse <- numeric(n_permutations)
  
  for (b in seq_len(n_permutations)) {
    group_permuted <- sample(group_original)
    
    r_perm_ASD <- get_cor(which(group_permuted == "ASD"))
    r_perm_NT <- get_cor(which(group_permuted == "NT"))
    
    vec_perm_ASD <- upper_values(r_perm_ASD)
    vec_perm_NT <- upper_values(r_perm_NT)
    keep_perm <- complete.cases(vec_perm_ASD, vec_perm_NT)
    
    permutation_rmse[b] <- sqrt(mean((vec_perm_ASD[keep_perm] - vec_perm_NT[keep_perm])^2))
  }
  
  permutation_p <- (sum(permutation_rmse >= observed_rmse, na.rm = TRUE) + 1) /
    (sum(!is.na(permutation_rmse)) + 1)
  
  tibble(
    analysis = analysis_name,
    matrix_similarity = matrix_similarity,
    global_RMSE = observed_rmse,
    permutation_p = permutation_p,
    n_permutations = n_permutations
  )
}

matrix_results_condition <- bind_rows(
  global_matrix_test(shared_neural_df, vars_valid, "Valid", n_permutations_global, permutation_seed),
  global_matrix_test(shared_neural_df, vars_invalid, "Invalid", n_permutations_global, permutation_seed)
)

print(matrix_results_condition, n = Inf)


#### 4.8 ASD–NT comparison of overall correlation strength --------------------------------------------------
# Mean absolute correlation strength: irrespective of direction
# Mean signed correlation:  whether the matrix is generally shifted toward more positive/negative relationships

compare_correlation_strength <- function(df, vars, analysis_name, n_permutations = 5000, seed = 1234) {
  dat <- df %>% filter(group %in% c("ASD", "NT"))
  
  get_summary <- function(group_labels) {
    r_ASD <- cor(dat[group_labels == "ASD", vars, drop = FALSE], use = "pairwise.complete.obs")
    r_NT <- cor(dat[group_labels == "NT", vars, drop = FALSE], use = "pairwise.complete.obs")
    
    v_ASD <- r_ASD[upper.tri(r_ASD)]
    v_NT <- r_NT[upper.tri(r_NT)]
    
    c(
      mean_abs_ASD = mean(abs(v_ASD), na.rm = TRUE),
      mean_abs_NT = mean(abs(v_NT), na.rm = TRUE),
      difference_abs = mean(abs(v_ASD), na.rm = TRUE) - mean(abs(v_NT), na.rm = TRUE),
      mean_signed_ASD = mean(v_ASD, na.rm = TRUE),
      mean_signed_NT = mean(v_NT, na.rm = TRUE),
      difference_signed = mean(v_ASD, na.rm = TRUE) - mean(v_NT, na.rm = TRUE)
    )
  }
  
  observed <- get_summary(dat$group)
  
  set.seed(seed)
  permutation_results <- replicate(n_permutations, get_summary(sample(dat$group))) %>% t()
  
  p_abs <- (sum(abs(permutation_results[, "difference_abs"]) >= abs(observed["difference_abs"])) + 1) /
    (n_permutations + 1)
  
  p_signed <- (sum(abs(permutation_results[, "difference_signed"]) >= abs(observed["difference_signed"])) + 1) /
    (n_permutations + 1)
  
  tibble(
    analysis = analysis_name,
    mean_abs_ASD = observed["mean_abs_ASD"],
    mean_abs_NT = observed["mean_abs_NT"],
    difference_abs = observed["difference_abs"],
    p_abs = p_abs,
    mean_signed_ASD = observed["mean_signed_ASD"],
    mean_signed_NT = observed["mean_signed_NT"],
    difference_signed = observed["difference_signed"],
    p_signed = p_signed,
    n_permutations = n_permutations
  )
}

strength_results_condition <- bind_rows(
  compare_correlation_strength(shared_neural_df, vars_valid, "Valid", n_permutations_global, permutation_seed),
  compare_correlation_strength(shared_neural_df, vars_invalid, "Invalid", n_permutations_global, permutation_seed)
)

print(strength_results_condition, n = Inf)


#### 4.9 Direct tests of matrix architecture ---------------------------------------------------
# 1) Is the ASD–NT matrix discrepancy larger in valid than invalid trials?
# 2) Is valid–invalid neural reconfiguration larger in NT than ASD?

matrix_rmse_between_groups <- function(dat, vars, group_labels) {
  r_ASD <- cor(dat[group_labels == "ASD", vars, drop = FALSE], use = "pairwise.complete.obs")
  r_NT <- cor(dat[group_labels == "NT", vars, drop = FALSE], use = "pairwise.complete.obs")
  
  v_ASD <- r_ASD[upper.tri(r_ASD)]
  v_NT <- r_NT[upper.tri(r_NT)]
  keep <- complete.cases(v_ASD, v_NT)
  
  sqrt(mean((v_ASD[keep] - v_NT[keep])^2))
}

matrix_rmse_between_conditions <- function(dat, vars_valid, vars_invalid, group_labels, group_name) {
  r_valid <- cor(dat[group_labels == group_name, vars_valid, drop = FALSE], use = "pairwise.complete.obs")
  r_invalid <- cor(dat[group_labels == group_name, vars_invalid, drop = FALSE], use = "pairwise.complete.obs")
  
  v_valid <- r_valid[upper.tri(r_valid)]
  v_invalid <- r_invalid[upper.tri(r_invalid)]
  keep <- complete.cases(v_valid, v_invalid)
  
  sqrt(mean((v_valid[keep] - v_invalid[keep])^2))
}

test_condition_dependent_architecture <- function(df, vars_valid, vars_invalid,
                                                  n_permutations = 50000, seed = 1234) {
  dat <- df %>% filter(group %in% c("ASD", "NT"))
  group_original <- dat$group
  
  rmse_group_valid <- matrix_rmse_between_groups(dat, vars_valid, group_original)
  rmse_group_invalid <- matrix_rmse_between_groups(dat, vars_invalid, group_original)
  delta_group_discrepancy <- rmse_group_valid - rmse_group_invalid
  
  rmse_condition_NT <- matrix_rmse_between_conditions(dat, vars_valid, vars_invalid, group_original, "NT")
  rmse_condition_ASD <- matrix_rmse_between_conditions(dat, vars_valid, vars_invalid, group_original, "ASD")
  delta_reconfiguration <- rmse_condition_NT - rmse_condition_ASD
  
  set.seed(seed)
  perm_delta_group <- numeric(n_permutations)
  perm_delta_reconfiguration <- numeric(n_permutations)
  
  for (b in seq_len(n_permutations)) {
    group_permuted <- sample(group_original)
    
    perm_valid <- matrix_rmse_between_groups(dat, vars_valid, group_permuted)
    perm_invalid <- matrix_rmse_between_groups(dat, vars_invalid, group_permuted)
    perm_delta_group[b] <- perm_valid - perm_invalid
    
    perm_NT <- matrix_rmse_between_conditions(dat, vars_valid, vars_invalid, group_permuted, "NT")
    perm_ASD <- matrix_rmse_between_conditions(dat, vars_valid, vars_invalid, group_permuted, "ASD")
    perm_delta_reconfiguration[b] <- perm_NT - perm_ASD
  }
  
  tibble(
    test = c(
      "ASD–NT discrepancy: valid minus invalid",
      "Valid–invalid reconfiguration: NT minus ASD"
    ),
    observed_difference = c(delta_group_discrepancy, delta_reconfiguration),
    p_permutation = c(
      (sum(abs(perm_delta_group) >= abs(delta_group_discrepancy)) + 1) / (n_permutations + 1),
      (sum(abs(perm_delta_reconfiguration) >= abs(delta_reconfiguration)) + 1) / (n_permutations + 1)
    ),
    n_permutations = n_permutations
  )
}

condition_architecture_test <- test_condition_dependent_architecture(
  shared_neural_df,
  vars_valid,
  vars_invalid,
  n_permutations_direct,
  permutation_seed
)

print(condition_architecture_test, n = Inf)

#### 4.10 Figure 3: Neural association Figure -------------------------------------------------------------------
base_family <- "sans"

## labels
pub_strip_suffix <- function(x) sub("_(valid|invalid)$", "", x)

pub_make_labels <- function(base_vars) {
  tibble(
    base = base_vars,
    label = vapply(base_vars, function(x) {
      band <- if (grepl("_alpha_", x)) "α" else "β"
      
      if (grepl("^power_", x)) {
        ROI <- if (grepl("_occi_", x)) "Post" else "Front"
        hemi <- if (grepl("_ipsi$", x)) "i" else "c"
        return(paste(band, "pow", ROI, hemi))
      }
      
      if (grepl("^ITPC_", x)) {
        ROI <- if (grepl("_occi_", x)) "Post" else "Front"
        return(paste(band, "ITPC", ROI))
      }
      
      if (grepl("^powconn_", x)) return(paste(band, "conn", "P–F"))
      if (grepl("^ispc_", x)) return(paste(band, "ISPC", "P–F"))
      
      x
    }, character(1))
  )
}

pub_key <- pub_make_labels(feature_base)


## style
pub_col_fun <- circlize::colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#274C9B", "#8FB6DE", "#FAFAFA", "#E3A18D", "#A51C30")
)

pub_upper_stars <- function(stars, hc) {
  pos <- integer(nrow(stars))
  pos[hc$order] <- seq_len(nrow(stars))
  keep <- outer(pos, pos, FUN = "<")
  
  out <- stars
  out[!keep] <- ""
  out
}

pub_spread_dendrogram <- function(hc) {
  dend <- as.dendrogram(hc)
  
  get_heights <- function(x) {
    if (is.leaf(x)) return(numeric(0))
    c(attr(x, "height"), unlist(lapply(x, get_heights), use.names = FALSE))
  }
  
  original_heights <- sort(unique(get_heights(dend)))
  display_heights <- if (length(original_heights) == 1) 1 else seq(0.08, 1, length.out = length(original_heights))
  
  replace_heights <- function(x) {
    if (is.leaf(x)) {
      attr(x, "height") <- 0
      return(x)
    }
    
    original_height <- attr(x, "height")
    for (i in seq_along(x)) x[[i]] <- replace_heights(x[[i]])
    attr(x, "height") <- display_heights[match(original_height, original_heights)]
    x
  }
  
  replace_heights(dend)
}


## heatmap panel
pub_heatmap_grob <- function(result, hc, key, panel_title,
                             show_row_labels = TRUE, show_column_labels = TRUE,
                             row_names_side = "right") {
  mat <- result$r
  stars_display <- pub_upper_stars(result$stars, hc)
  
  row_map <- setNames(key$label, key$base)
  current_base <- pub_strip_suffix(rownames(mat))
  current_labels <- unname(row_map[current_base])
  
  if (anyNA(current_labels)) {
    stop("Some publication-figure labels could not be matched.")
  }
  
  dend_display <- pub_spread_dendrogram(hc)
  
  ht <- Heatmap(
    mat,
    name = "Pearson r",
    col = pub_col_fun,
    
    cluster_rows = dend_display,
    cluster_columns = dend_display,
    row_dend_reorder = FALSE,
    column_dend_reorder = FALSE,
    
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_side = "left",
    column_dend_side = "top",
    
    row_dend_width = unit(10, "mm"),
    column_dend_height = unit(12, "mm"),
    row_dend_gp = gpar(col = "#555555", lwd = 0.5),
    column_dend_gp = gpar(col = "#555555", lwd = 0.5),
    
    show_row_names = show_row_labels,
    row_names_side = row_names_side,
    row_labels = current_labels,
    row_names_gp = gpar(
      fontsize = 8,
      fontfamily = base_family
    ),
    row_names_max_width = unit(24, "mm"),
    
    show_column_names = show_column_labels,
    column_labels = current_labels,
    column_names_gp = gpar(
      fontsize = 8,
      fontfamily = base_family
    ),
    column_names_rot = 55,
    column_names_centered = FALSE,
    column_names_max_height = unit(24, "mm"),
    
    heatmap_width = unit(84, "mm"),
    heatmap_height = unit(86, "mm"),
    
    rect_gp = gpar(col = NA, lwd = 0),
    border = FALSE,
    use_raster = FALSE,
    
    column_title = panel_title,
    column_title_gp = gpar(fontsize = 10.5, fontface = "bold", fontfamily = base_family),
    
    show_heatmap_legend = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- stars_display[i, j]
      
      if (!is.na(lab) && lab != "") {
        grid.text(
          lab, x, y,
          gp = gpar(
            fontsize = 8,
            fontface = "plain",
            fontfamily = base_family,
            col = "#555555"
          )
        )
      }
    }
  )
  
  grid.grabExpr(
    draw(
      ht,
      newpage = FALSE,
      padding = unit(c(2, 2, 3, 2), "mm")
    )
  )
}


## create four panels
pub_A <- pub_heatmap_grob(
  result = res_valid$NT,
  hc = res_valid$hc,
  key = pub_key,
  panel_title = "A  Neurotypical"
)

pub_B <- pub_heatmap_grob(
  result = res_valid$ASD,
  hc = res_valid$hc,
  key = pub_key,
  panel_title = "B  Autistic"
)

pub_C <- pub_heatmap_grob(
  result = res_invalid$NT,
  hc = res_invalid$hc,
  key = pub_key,
  panel_title = "C  Neurotypical"
)

pub_D <- pub_heatmap_grob(
  result = res_invalid$ASD,
  hc = res_invalid$hc,
  key = pub_key,
  panel_title = "D  Autistic"
)

## headers and legend
pub_format_p <- function(p) {
  if (p < .001) "< .001" else paste0("= ", sprintf("%.3f", p))
}

pub_make_header <- function(title, stats_row) {
  
  grobTree(
    
    rectGrob(
      x = .5,
      y = .5,
      width = 1,
      height = 1,
      gp = gpar(
        fill = "#F2F2F2",
        col = NA
      )
    ),
    
    textGrob(
      title,
      x = .5,
      y = .72,
      gp = gpar(
        fontsize = 10.5,
        fontface = "bold",
        fontfamily = base_family,
        col = "grey18"
      )
    ),
    
    textGrob(
      paste0(
        "Matrix correspondence r = ",
        sprintf("%.3f", stats_row$matrix_similarity),
        "   |   RMSE = ",
        sprintf("%.3f", stats_row$global_RMSE),
        "   |   permutation p ",
        pub_format_p(stats_row$permutation_p)
      ),
      x = .5,
      y = .38,
      gp = gpar(
        fontsize = 8.5,
        fontfamily = base_family,
        col = "grey35"
      )
    )
  )
}

pub_stat_valid <- matrix_results_condition %>% filter(analysis == "Valid") %>% slice(1)
pub_stat_invalid <- matrix_results_condition %>% filter(analysis == "Invalid") %>% slice(1)

pub_header_valid <- pub_make_header("Valid cue condition", pub_stat_valid)
pub_header_invalid <- pub_make_header("Invalid cue condition", pub_stat_invalid)

pub_legend <- grid.grabExpr(
  draw(
    Legend(
      title = "Pearson r",
      col_fun = pub_col_fun,
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("−1.0", "−0.5", "0", "0.5", "1.0"),
      title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = base_family),
      labels_gp = gpar(fontsize = 8,   fontfamily = base_family),
      legend_height = unit(25, "mm")
    )
  )
)


## assemble figure
pub_top_row <- arrangeGrob(
  pub_A,
  nullGrob(),
  pub_B,
  ncol = 3,
  widths = unit.c(unit(1, "null"), unit(2, "mm"), unit(1, "null")),
  padding = unit(0, "mm")
)

pub_bottom_row <- arrangeGrob(
  pub_C,
  nullGrob(),
  pub_D,
  ncol = 3,
  widths = unit.c(unit(1, "null"), unit(2, "mm"), unit(1, "null")),
  padding = unit(0, "mm")
)

pub_legend_block <- arrangeGrob(
  pub_legend,
  nullGrob(),
  ncol = 2,
  widths = c(0.98, 0.02)
)

pub_main <- arrangeGrob(
  pub_header_valid,
  pub_top_row,
  nullGrob(),
  pub_header_invalid,
  pub_bottom_row,
  ncol = 1,
  heights = unit.c(
    unit(11, "mm"),
    unit(87, "mm"),
    unit(4, "mm"),
    unit(11, "mm"),
    unit(87, "mm")
  ),
  padding = unit(0, "mm")
)

pub_figure <- arrangeGrob(
  pub_main,
  pub_legend_block,
  ncol = 2,
  widths = unit.c(unit(1, "null"), unit(16, "mm")),
  padding = unit(0, "mm")
)


## open figure in R
graphics.off()

fig_width  <- 190 / 25.4
fig_height <- 210 / 25.4

windows(
  width = fig_width,
  height = fig_height
)

grid.newpage()
grid.draw(pub_figure)


## export
output_dir <- paste0(
  ".../Tables and Figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(output_dir)) {
  stop("The output directory could not be accessed:\n", output_dir)
}

pdf_file <- file.path(
  output_dir,
  "Figure_neural_correlations_valid_invalid.pdf"
)

png_file <- file.path(
  output_dir,
  "Figure_neural_correlations_valid_invalid.png"
)

# PDF export
cairo_pdf(
  filename = pdf_file,
  width = fig_width,
  height = fig_height,
  family = base_family
)

grid.newpage()
grid.draw(pub_figure)
dev.off()

# PNG export
png(
  filename = png_file,
  width = fig_width,
  height = fig_height,
  units = "in",
  res = 600,
  type = "cairo",
  bg = "white"
)

grid.newpage()
grid.draw(pub_figure)
dev.off()

shell.exec(normalizePath(png_file, winslash = "\\", mustWork = TRUE))


# 5 PHENOTYPIC ASSOCIATIONS - NEURAL MARKERS  ###########################################################################################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(lmtest)
library(sandwich)
library(car)

#### 5.1 Construct the four neural markers --------------------------------------------------------------------------
# One marker represents each neural domain identified in the group analyses:
# 1) Ipsilateral alpha power, averaged across frontal and posterior ROIs
# 2) Frontal beta power, averaged across hemispheres
# 3) Posterior contralateral beta ITPC
# 4) Posterior–frontal contralateral beta ISPC
#
# Each marker is represented both as:
# - the equally weighted mean across valid and invalid conditions
# - valid-minus-invalid modulation
phen_neural_df <- shared_neural_df %>%
  mutate(
    alpha_ipsi_mean = (power_alpha_front_ipsi_mean_valid_invalid + power_alpha_occi_ipsi_mean_valid_invalid) / 2,
    alpha_ipsi_diff = (power_alpha_front_ipsi_diff_valid_minus_invalid + power_alpha_occi_ipsi_diff_valid_minus_invalid) / 2,
    
    beta_front_mean = (power_beta_front_ipsi_mean_valid_invalid + power_beta_front_contra_mean_valid_invalid) / 2,
    beta_front_diff = (power_beta_front_ipsi_diff_valid_minus_invalid + power_beta_front_contra_diff_valid_minus_invalid) / 2,
    
    beta_ITPC_post_mean = ITPC_beta_occi_contra_mean_valid_invalid,
    beta_ITPC_post_diff = ITPC_beta_occi_contra_diff_valid_minus_invalid,
    
    beta_ISPC_mean = ispc_beta_occi_front_contra_mean_valid_invalid,
    beta_ISPC_diff = ispc_beta_occi_front_contra_diff_valid_minus_invalid
  ) %>%
  filter(group == "ASD") %>%
  mutate(sex = factor(sex))


#### 5.2. Define the two neural-marker sets -------------------------------------------------
mean_vars <- c("alpha_ipsi_mean", "beta_front_mean", "beta_ITPC_post_mean", "beta_ISPC_mean")
diff_vars <- c("alpha_ipsi_diff", "beta_front_diff", "beta_ITPC_post_diff", "beta_ISPC_diff")

## 5.3 Prepare complete-case and standardized model data ------------------------------------
prepare_model_data <- function(data, outcome, neural_vars) {
  continuous_vars <- c(outcome, neural_vars, "age", "IQ")
  
  data %>%
    select(all_of(c(continuous_vars, "sex"))) %>%
    drop_na() %>%
    mutate(across(all_of(continuous_vars), ~as.numeric(scale(.x))))
}


#### 5.4 Prepare data for models --------------------------------------
dat_adhd_mean <- prepare_model_data(phen_neural_df, "ADHS_aufm", mean_vars)
dat_adhd_diff <- prepare_model_data(phen_neural_df, "ADHS_aufm", diff_vars)

dat_rrb_mean <- prepare_model_data(phen_neural_df, "Severity_RRB", mean_vars)
dat_rrb_diff <- prepare_model_data(phen_neural_df, "Severity_RRB", diff_vars)

dat_sa_mean <- prepare_model_data(phen_neural_df, "Severity_SA", mean_vars)
dat_sa_diff <- prepare_model_data(phen_neural_df, "Severity_SA", diff_vars)


#### 5.5 Covariate-adjusted models ----------------------------------------------------------

## ADHD inattention
m_adhd_mean <- lm(
  ADHS_aufm ~ alpha_ipsi_mean + beta_front_mean + beta_ITPC_post_mean + beta_ISPC_mean + age + IQ + sex,
  data = dat_adhd_mean
)

m_adhd_diff <- lm(
  ADHS_aufm ~ alpha_ipsi_diff + beta_front_diff + beta_ITPC_post_diff + beta_ISPC_diff + age + IQ + sex,
  data = dat_adhd_diff
)

## ADOS restricted and repetitive behaviours
m_rrb_mean <- lm(
  Severity_RRB ~ alpha_ipsi_mean + beta_front_mean + beta_ITPC_post_mean + beta_ISPC_mean + age + IQ + sex,
  data = dat_rrb_mean
)

m_rrb_diff <- lm(
  Severity_RRB ~ alpha_ipsi_diff + beta_front_diff + beta_ITPC_post_diff + beta_ISPC_diff + age + IQ + sex,
  data = dat_rrb_diff
)

## ADOS social affect 
m_sa_mean <- lm(
  Severity_SA ~ alpha_ipsi_mean + beta_front_mean + beta_ITPC_post_mean + beta_ISPC_mean + age + IQ + sex,
  data = dat_sa_mean
)

m_sa_diff <- lm(
  Severity_SA ~ alpha_ipsi_diff + beta_front_diff + beta_ITPC_post_diff + beta_ISPC_diff + age + IQ + sex,
  data = dat_sa_diff
)


#### 5.6 Store models -------------------------------------------------------------------
phenotype_models <- tibble(
  phenotype = c("ADHD inattention", "ADHD inattention",
                "ADOS RRB", "ADOS RRB",
                "ADOS social affect", "ADOS social affect"),
  
  score_set = c("Condition average", "Valid-minus-invalid",
                "Condition average", "Valid-minus-invalid",
                "Condition average", "Valid-minus-invalid"),
  
  predictors = list(mean_vars, diff_vars, mean_vars, diff_vars, mean_vars, diff_vars),
  
  model = list(m_adhd_mean, m_adhd_diff,
               m_rrb_mean, m_rrb_diff,
               m_sa_mean, m_sa_diff)
)


#### 5.7 Extract neural coefficients  ------------------------------------------
# using HC3 standard errors

extract_HC3_coefficients <- function(model, neural_vars) {
  robust_covariance <- sandwich::vcovHC(model, type = "HC3")
  robust_test <- lmtest::coeftest(model, vcov. = robust_covariance)
  
  estimate <- robust_test[, 1]
  robust_SE <- robust_test[, 2]
  critical_t <- qt(.975, df.residual(model))
  
  tibble(
    predictor = rownames(robust_test),
    beta = estimate,
    SE_HC3 = robust_SE,
    CI_low_HC3 = estimate - critical_t * robust_SE,
    CI_high_HC3 = estimate + critical_t * robust_SE,
    t_HC3 = robust_test[, 3],
    p_HC3 = robust_test[, 4]
  ) %>%
    filter(predictor %in% neural_vars)
}

HC3_coefficients <- phenotype_models %>%
  mutate(
    n = map_int(model, nobs),
    R2 = map_dbl(model, ~broom::glance(.x)$r.squared),
    adj_R2 = map_dbl(model, ~broom::glance(.x)$adj.r.squared),
    coefficients = map2(model, predictors, extract_HC3_coefficients)
  ) %>%
  select(phenotype, score_set, n, R2, adj_R2, coefficients) %>%
  unnest(coefficients) %>%
  mutate(marker = sub("_(mean|diff)$", "", predictor)) %>%
  
  # BH-FDR across the three phenotypes, separately for each neural marker
  # and separately for condition-average and valid-minus-invalid scores
  group_by(score_set, marker) %>%
  mutate(p_FDR = p.adjust(p_HC3, method = "BH")) %>%
  ungroup() %>%
  
  arrange(phenotype, score_set, p_FDR)

print(HC3_coefficients, n = Inf)


#### 5.8 HC3-robust omnibus tests of neural markers ----------------------

extract_HC3_omnibus <- function(model, neural_vars) {
  hypotheses <- paste0(neural_vars, " = 0")
  
  robust_test <- car::linearHypothesis(
    model,
    hypothesis.matrix = hypotheses,
    vcov. = sandwich::vcovHC(model, type = "HC3"),
    test = "F"
  )
  
  robust_table <- as.data.frame(robust_test)
  result_row <- robust_table[nrow(robust_table), , drop = FALSE]
  
  tibble(
    df_num = result_row[["Df"]],
    df_den = result_row[["Res.Df"]],
    F_HC3 = result_row[["F"]],
    p_raw = result_row[["Pr(>F)"]]
  )
}

HC3_omnibus <- phenotype_models %>%
  mutate(
    n = map_int(model, nobs),
    R2 = map_dbl(model, ~broom::glance(.x)$r.squared),
    adj_R2 = map_dbl(model, ~broom::glance(.x)$adj.r.squared),
    omnibus = map2(model, predictors, extract_HC3_omnibus)
  ) %>%
  select(phenotype, score_set, n, R2, adj_R2, omnibus) %>%
  unnest(omnibus) %>%
  
  # Two separate BH-FDR families:
  # 1) the three condition-average omnibus tests
  # 2) the three valid-minus-invalid omnibus tests
  group_by(score_set) %>%
  mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  arrange(score_set, p_FDR)

print(HC3_omnibus, n = Inf)


#### 5.9 Plot phenotype distributions --------------------------------------------
phen_neural_df %>%
  select(ADHS_aufm, Severity_RRB, Severity_SA) %>%
  pivot_longer(everything(), names_to = "phenotype", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10, boundary = 0, colour = "white") +
  facet_wrap(~phenotype, scales = "free") +
  labs(x = "Score", y = "Count") +
  theme_classic()


#### 5.10 Diagnostic plots for covariate-adjusted models --------------------------

diagnostic_models <- list(
  "ADHD inattention - condition average" = m_adhd_mean,
  "ADHD inattention - valid minus invalid" = m_adhd_diff,
  "ADOS RRB - condition average" = m_rrb_mean,
  "ADOS RRB - valid minus invalid" = m_rrb_diff,
  "ADOS social affect - condition average" = m_sa_mean,
  "ADOS social affect - valid minus invalid" = m_sa_diff
)

pdf("phenotype_model_diagnostics_HC3.pdf", width = 10, height = 8)

for (model_name in names(diagnostic_models)) {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.5, 1), oma = c(0, 0, 3, 0))
  plot(diagnostic_models[[model_name]], which = c(1, 2, 3, 5))
  mtext(model_name, outer = TRUE, line = 1, font = 2, cex = 1.3)
}

dev.off()


#### 5.11 Numerical diagnostic summary -----------------------------------------
diagnostic_summary <- phenotype_models %>%
  mutate(
    n = map_int(model, nobs),
    shapiro_p = map_dbl(model, ~shapiro.test(residuals(.x))$p.value),
    breusch_pagan_p = map_dbl(model, ~lmtest::bptest(.x)$p.value),
    max_Cooks_D = map_dbl(model, ~max(cooks.distance(.x))),
    n_Cooks_D_above_4_over_n = map_int(model, ~sum(cooks.distance(.x) > 4 / nobs(.x)))
  ) %>%
  select(phenotype, score_set, n, shapiro_p, breusch_pagan_p,
         max_Cooks_D, n_Cooks_D_above_4_over_n)

print(diagnostic_summary, n = Inf)


#### 5.12 inspect complete model #####
lmtest::coeftest(m_adhd_mean, vcov. = sandwich::vcovHC(m_adhd_mean, type = "HC3"))
lmtest::coeftest(m_adhd_diff, vcov. = sandwich::vcovHC(m_adhd_diff, type = "HC3"))

lmtest::coeftest(m_rrb_mean, vcov. = sandwich::vcovHC(m_rrb_mean, type = "HC3"))
lmtest::coeftest(m_rrb_diff, vcov. = sandwich::vcovHC(m_rrb_diff, type = "HC3"))

lmtest::coeftest(m_sa_mean, vcov. = sandwich::vcovHC(m_sa_mean, type = "HC3"))
lmtest::coeftest(m_sa_diff, vcov. = sandwich::vcovHC(m_sa_diff, type = "HC3"))


#### 5.13 Figure 4: Phenotype association figure ----------------------------------------------
library(dplyr)
library(ggplot2)
library(grid)

base_family <- "sans"
## Prepare plotting data
forest_df <- HC3_coefficients %>%
  mutate(
    phenotype = factor(
      phenotype,
      levels = c("ADOS social affect", "ADOS RRB", "ADHD inattention"),
      labels = c(
        labels = c(
          "ADOS social affect",
          "ADOS RRB",
          "FBB-ADHD inattention"
        )
      )
    ),
    
    score_set = factor(
      score_set,
      levels = c("Condition average", "Valid-minus-invalid"),
      labels = c(
        "Condition average",
        "Validity modulation"
      )
    ),
    
    # Factor order is defined from bottom to top
    marker = factor(
      marker,
      levels = c(
        "beta_ISPC",
        "beta_ITPC_post",
        "beta_front",
        "alpha_ipsi"
      ),
      labels = c(
        "β ISPC, posterior–frontal",
        "β ITPC, posterior",
        "β power, frontal",
        "α power, ipsilateral"
      )
    ),
    
    sig_status = factor(
      case_when(
        p_FDR < .05 ~ "adj",
        p_HC3 < .05 ~ "raw",
        TRUE ~ "ns"
      ),
      levels = c("ns", "raw", "adj")
    )
  )


## Forest plot
legend_df <- tibble(
  sig_status = factor(
    c("ns", "raw", "adj"),
    levels = c("ns", "raw", "adj")
  )
)


p_forest <- ggplot(forest_df, aes(y = marker)) +
  
  # Zero-effect reference
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey70",
    linewidth = 0.35
  ) +
  
  # HC3 confidence intervals
  geom_segment(
    aes(
      x = CI_low_HC3,
      xend = CI_high_HC3,
      yend = marker
    ),
    colour = "#555555",
    linewidth = 0.50,
    lineend = "round"
  ) +
  
  # Invisible points outside the plotting range, used only to construct
  # a complete three-category legend
  geom_point(
    data = legend_df,
    aes(x = -99, y = 1, fill = sig_status),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.5,
    colour = "#252525",
    stroke = 0.6,
    show.legend = TRUE
  ) +
  
  # Actual model estimates
  geom_point(
    aes(x = beta, fill = sig_status),
    shape = 21,
    size = 2.8,
    colour = "#252525",
    stroke = 0.65,
    show.legend = FALSE
  ) +
  
  facet_grid(
    rows = vars(score_set),
    cols = vars(phenotype)
  ) +
  
  scale_fill_manual(
    values = c(
      "ns" = "white",
      "raw" = "#8F8F8F",
      "adj" = "black"
    ),
    limits = c("ns", "raw", "adj"),
    breaks = c("ns", "raw", "adj"),
    labels = expression(
      p >= .05,
      p < .05,
      p[adj] < .05
    ),
    drop = FALSE,
    guide = guide_legend(
      direction = "horizontal",
      override.aes = list(
        shape = 21,
        size = 2.5,
        colour = "#252525",
        stroke = 0.6,
        fill = c("white", "#8F8F8F", "black")
      )
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(-0.8, 0.8, by = 0.4),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  
  coord_cartesian(
    xlim = c(-0.85, 0.85),
    clip = "off"
  ) +
  
  labs(
    x = "Standardized regression coefficient (β)",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 10,
    base_family = base_family
  ) +
  
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.y = element_text(
      colour = "#252525",
      size = 8.5
    ),
    axis.text.x = element_text(
      colour = "#252525",
      size = 8.5
    ),
    axis.title.x = element_text(
      margin = margin(t = 7),
      size = 9
    ),
    
    panel.grid.major.y = element_line(
      colour = "#EEEEEE",
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    
    strip.background.x = element_rect(
      fill = "#F4F4F4",
      colour = NA
    ),
    strip.background.y = element_blank(),
    
    strip.text.x = element_text(
      face = "bold",
      size = 9,
      lineheight = .95,
      margin = margin(3, 4, 3, 4)
    ),
    
    strip.text.y.right = element_text(
      face = "plain",
      size = 8.5,
      colour = "grey35",
      angle = 90,
      margin = margin(4, 3, 4, 3)
    ),
    
    panel.spacing.x = unit(6, "pt"),
    panel.spacing.y = unit(8, "pt"),
    
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = unit(12, "pt"),
    legend.spacing.x = unit(4, "pt"),
    legend.margin = margin(2, 0, 0, 0),
    legend.text = element_text(
      size = 8,
      colour = "#252525"
    ),
    
    plot.margin = margin(7, 8, 7, 8)
  )

p_forest

## Export 
output_dir <- paste0(
  ".../Tables and Figures"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(output_dir)) {
  stop("Output directory could not be accessed:\n", output_dir)
}

pdf_file <- file.path(
  output_dir,
  "Figure_phenotype_associations_forest_plot.pdf"
)

png_file <- file.path(
  output_dir,
  "Figure_phenotype_associations_forest_plot.png"
)

fig_width  <- 190 / 25.4
fig_height <- 120 / 25.4

ggsave(
  filename = pdf_file,
  plot = p_forest,
  device = cairo_pdf,
  width = fig_width,
  height = fig_height,
  units = "in",
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = p_forest,
  device = "png",
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 600,
  bg = "white"
)

shell.exec(normalizePath(png_file, winslash = "\\", mustWork = TRUE))