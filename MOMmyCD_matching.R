library(MatchIt)
library(readxl)
library(dplyr)
library(tableone)
# Methods: https://kosukeimai.github.io/MatchIt/articles/assessing-balance.html#recommendations-for-balance-assessment-1
# Data cleaning and combine---- 
MOMmyCDdb <- read_excel("Data/hui_control_selection.xlsx")[,c(1,2,6,15,16)]
MOMmyCDdb$Sub_ID <- substr(MOMmyCDdb$Sub_ID,0,9)

# Change to readRDS("MOMmyCD_control_db_01112021.rds") to read MOMmyCDcon db directly
# MOMmyCDcon <- PNdb %>% filter(!Sub_ID %in% MOMmyCDdb$Sub_ID)
# MOMmyCDcon <- MOMmyCDcon[,-c(11,12,20,21,23:60,89:128,180)]
# MOMmyCDcon <- MOMmyCDcon %>% filter(Region=="HK"&Any_disease=="No"&complete.cases(BMI)&
#                                       complete.cases(Del_mode)&complete.cases(Age))
MOMmyCD_match <- MOMmyCDcon %>% select(2,15,20,62)
MOMmyCD_match <- data.frame(MOMmyCD_match[,c(1:3)],
                            SRD = "No",
                            MOD = MOMmyCD_match[,4])

MOMmyCD_match$MOD <- ifelse(MOMmyCD_match$MOD=="陰道分娩","Vag",
                            ifelse(MOMmyCD_match$MOD=="剖腹產","CS","VOT"))
colnames(MOMmyCDdb)[c(2,4,5)] <- c("Age","SRD","MOD")
MOMmyCD_match <- data.frame(rbind(MOMmyCDdb,
                       MOMmyCD_match))

MOMmyCD_match$SRD_bin <- ifelse(MOMmyCD_match$SRD=="IBD",1,0)

# Matching----
# 1. Match using nearest, 1:8, set seed == 100
set.seed(100) # Must run this before matching!
MOMmyCD_matched <- matchit(SRD_bin ~ Age + BMI + MOD,
                         data = MOMmyCD_match, ratio = 2.2,
                         caliper=0.2, discard= "both",
                         method = "nearest", distance = "logit")
# summary(MOMmyCD_matched)
MOMmyCD_matched <- match.data(MOMmyCD_matched, group = "all")


# 2. Check data before and after matching
Tabb4 <-CreateTableOne(vars = c("Age","BMI","MOD"), data = MOMmyCD_match, 
                       strata = "SRD")
Tabb4 <- data.frame(print(Tabb4, nonnormal = c("Age","BMI"), exact = "stage", smd = T)) 
Tabaft <-CreateTableOne(vars = c("Age","BMI","MOD"), data = MOMmyCD_matched, 
                       strata = "SRD")
Tabaft <- data.frame(print(Tabaft, nonnormal = c("Age","BMI"), exact = "stage", smd = T)) 
# All SMD of co-variate <=0.1, indicates matching satisfied

writexl::write_xlsx(MOMmyCD_matched,"MOMmyCD_matched_01112021.xlsx")
saveRDS(MOMmyCDcon, "MOMmyCD_control_db_01112021.rds")
demodb_pre <-  demodb_new %>% filter (substr(Sub_ID,6,9) %in% (substr(preterm_ID,6,9)))

# Matching 2nd round (20211117)----
MOMmyCDdb <- read_excel("Data/MOMmy_CD/wing_control_selection_20211117.xlsx",
                        sheet = "Baseline_CD")[,c(1,3,7:9)]
MOMmyCDdb$Sub_ID <- substr(MOMmyCDdb$Sub_ID,0,9)

MOMmyCDcon <- readRDS("Data/MOMmy_CD/MOMmyCD_control_db_01112021.rds") 
MOMmyCD_match <- MOMmyCDcon %>% select(2,15,20,62)
MOMmyCD_match <- data.frame(MOMmyCD_match[,c(1:3)],
                            SRD = "No",
                            MOD = MOMmyCD_match[,4])

MOMmyCD_match$MOD <- ifelse(MOMmyCD_match$MOD=="陰道分娩","Vag",
                            ifelse(MOMmyCD_match$MOD=="剖腹產","CS","VOT"))
colnames(MOMmyCDdb)[c(2,4,5)] <- c("Age","SRD","MOD")



# Include the "super healthy" families only
Healthy_fams <- readRDS("Data/MOMmy_CD/Healthy_fams_18112021.rds")
MOMmyCD_match <- MOMmyCD_match %>% filter(substr(Sub_ID,6,9) %in% Healthy_fams)
MOMmyCD_match <- data.frame(rbind(MOMmyCDdb,
                                  MOMmyCD_match))
MOMmyCD_match$SRD_bin <- ifelse(MOMmyCD_match$SRD=="IBD",1,0)


# Matching
set.seed(10000) # Must run this before matching!
MOMmyCD_matched <- matchit(SRD_bin ~ Age + BMI + MOD,
                           data = MOMmyCD_match, ratio = 4,
                           caliper = 0.4, discard= "both",
                           method = "nearest", distance = "logit")
# summary(MOMmyCD_matched)
MOMmyCD_matched <- match.data(MOMmyCD_matched, group = "all")
writexl::write_xlsx(MOMmyCD_matched,"Data/MOMmy_CD/MOMmyCD_CD_matched_18112021.xlsx")


# New set for IBD matching
MOMmyCDdb <- data.frame(read_excel("Data/MOMmy_CD/wing_control_selection_03122021.xlsx",
                        sheet = "Baseline_IBD")[,c(1,3,7,16,17)])
colnames(MOMmyCDdb)[c(2,4,5)] <- c("Age","SRD","MOD")

MOMmyCDcon <- readRDS("Data/MOMmy_CD/MOMmyCD_control_db_01112021.rds") 
MOMmyCD_match <- MOMmyCDcon %>% select(2,15,20,62)
MOMmyCD_match <- data.frame(MOMmyCD_match[,c(1:3)],
                            SRD = "No",
                            MOD = MOMmyCD_match[,4])

MOMmyCD_match$MOD <- ifelse(MOMmyCD_match$MOD=="陰道分娩","Vag",
                            ifelse(MOMmyCD_match$MOD=="剖腹產","CS","VOT"))

Healthy_fams <- readRDS("Data/MOMmy_CD/Healthy_fams_18112021.rds")
MOMmyCD_match <- MOMmyCD_match %>% filter(substr(Sub_ID,6,9) %in% Healthy_fams)
MOMmyCD_match <- data.frame(rbind(MOMmyCDdb,
                                  MOMmyCD_match))

MOMmyCD_match$SRD_bin <- ifelse(MOMmyCD_match$SRD=="IBD",1,0)
set.seed(100) # Must run this before matching!
MOMmyCD_matched <- matchit(SRD_bin ~ Age + BMI + MOD,
                           data = MOMmyCD_match, ratio = 2,
                           caliper = 0.4, discard= "both",
                           method = "nearest", distance = "logit")
# summary(MOMmyCD_matched)
MOMmyCD_matched <- match.data(MOMmyCD_matched, group = "all")
writexl::write_xlsx(MOMmyCD_matched,"Data/MOMmy_CD/MOMmyCD_IBD_matched_03122021.xlsx")







