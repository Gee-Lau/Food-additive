library(vegan)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggvegan)
library(reshape2)
library(ggpubr)
library(readxl)
library(stringr)
library(ggplot2)
library(ape)
library(hagis)
library(tableone)
library(MicrobiotaProcess)
library(phyloseq)

# Food additive nMDS----

HK_FS_FA <- MomFA[,c(1:10)]

FS_FA <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_MOMmy_Annual_Food_Additives_2022.03.21.xlsx"))
FS_subID <- data.frame(rbind(read.table("Data/FOSHAN data20220128/Mom_subID.csv",
                       sep = ",", header = TRUE),
                       rbind(read.table("Data/FOSHAN data20220128/Dad_subID.csv",
                                        sep = ",", header = TRUE))))
FS_subID <- FS_subID %>% filter(startsWith(Sub_ID,"F"))
FS_subID$Record.ID[1:12] <- paste0("21-",str_pad(1:12, 2, pad = "0"))
FS_FA$Sub_ID <- FS_subID$Sub_ID[match(FS_FA$Additives..mg.,FS_subID$Record.ID)]
FS_FA <- FS_FA[,c(12,2:4,6:11)]
HK_FS_FA <- rbind(HK_FS_FA,FS_FA)
HK_FS_FA$Area <- ifelse(startsWith(HK_FS_FA$Sub_ID,"P"),"HK","FS")
HK_FS_FA$Sex <- ifelse(substr(HK_FS_FA$Sub_ID,4,4)=="1","Female","Male")

gender.ds <- vegdist(HK_FS_FA[,c(2:10)],na.rm = T)

gender.md <- monoMDS(gender.ds)
gender.md <- data.frame(cbind(Sub_ID = HK_FS_FA$Sub_ID,
                              Sex = HK_FS_FA$Sex,
                              Area = HK_FS_FA$Area,
                              gender.md[["points"]]))
gender.md$MDS1 <- as.numeric(gender.md$MDS1)
gender.md$MDS2 <- as.numeric(gender.md$MDS2)
gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
filter(abs(MDS2)<0.1)
# saveRDS(gender.md,"HK_FS_FA_NMDS.rds")

ggplot(gender.md, aes(MDS1, MDS2, col = Sex)) +
  geom_point() +
  scale_color_manual(name = "Sex",
                     label = c("Female","Male"),
                     values=c("#018571","#0571b0")) +
  stat_ellipse(linetype = 2, size = 1.5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("MDS of weight adjusted annual estimated food additive") 
invisible(dev.off())

HK_FS_FA_Mom <- HK_FS_FA %>% filter(Sex=="Female") %>% 
  filter(Sub_ID %in% gender.md$Sub_ID)

gender.anosim <- anosim(as.matrix(HK_FS_FA_Mom[,c(2:10)]),  distance = "bray",
                        HK_FS_FA_Mom$Area,
                        permutations = 999)

# Baseline demographic----
MOMdemo_FS <- data.frame(read_excel("Data/FOSHAN data20220128/Mother baseline3.18.xlsx"))
MOM_HK_FS_demo <- MOMcomb_sub %>% select(colnames(MOMcomb_sub)[colnames(MOMcomb_sub) %in% colnames(MOMdemo_FS)]) %>% 
  filter(substr(Sub_ID,4,4)=="1")
MOMdemo_FS$Weight <- str_replace(MOMdemo_FS$Weight,"（孕早期）","")
MOMdemo_FS$Weight <- as.numeric(MOMdemo_FS$Weight)
MOMdemo_FS <- MOMdemo_FS %>% select(colnames(MOMdemo_FS)[colnames(MOMdemo_FS) %in% colnames(MOM_HK_FS_demo)]) %>% 
  mutate(Area = "FS")
MOM_HK_FS_demo$Area <- "HK"

MOMdemo_FS$Income_fam_new[MOMdemo_FS$Income_fam_new=="< 2,000"|MOMdemo_FS$Income_fam_new=="2,000 - 4,999"|MOMdemo_FS$Income_fam_new=="5,000 - 9,999"|MOMdemo_FS$Income_fam_new=="10,000 - 14,999"|MOMdemo_FS$Income_fam_new=="15,000 - 29,999"] <- "<30,000"
MOMdemo_FS$Income_fam_new[MOMdemo_FS$Income_fam_new=="30,000 - 49,999"] <- "30,000 - 49,999"

MOMdemo_FS$Income_fam_new <- factor(as.factor(MOMdemo_FS$Income_fam_new),
                                levels(as.factor(MOMdemo_FS$Income_fam_new))[c(1,3,4,2)])
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="High school or bellow"] <- "Highschool_bellow"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Bachelor or equivalent"] <- "Bachelor/degree"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Master or above"] <- "Master_above"

MOMdemo_FS$Edu <- factor(as.factor(MOMdemo_FS$Edu),
                     levels(as.factor(MOMdemo_FS$Edu))[c(2,1,3)])
substr(MOMdemo_FS$Work_stat,1,1) <- toupper(substr(MOMdemo_FS$Work_stat,1,1) )
MOMdemo_FS$Work_stat <- factor(as.factor(MOMdemo_FS$Work_stat),
                           levels(as.factor(MOMdemo_FS$Work_stat))[c(1,4,2,3)])
MOMdemo_FS$PE_Group <- ifelse(MOMdemo_FS$PE>2,"≥3",
                                  paste(MOMdemo_FS$PE))
MOMdemo_FS$PE_Group <- factor(as.factor(MOMdemo_FS$PE_Group),
                          levels(as.factor(MOMdemo_FS$PE_Group))[c(2,3,1)])
MOMdemo_FS[,c(13:30)] <- apply(MOMdemo_FS[,c(13:30)], 2,
                               function(x) ifelse(x=="NA",NA,
                                                  ifelse(x=="No","1No","2Yes")))
MOM_HK_FS_demo$Any_disease <- ifelse(MOM_HK_FS_demo$Any_disease=="No","1No","2Yes")
MOM_HK_FS_demo <- rbind(MOM_HK_FS_demo,MOMdemo_FS)

MOM_HK_FS_demo$BMI <- round(MOM_HK_FS_demo$Weight/(MOM_HK_FS_demo$Height/100)^2,1)
MOM_HK_FS_demo$BMI_Group <- ifelse(is.na(MOM_HK_FS_demo$BMI),NA,
                            ifelse(MOM_HK_FS_demo$BMI<18.5,"aUnderweight",
                                   ifelse(MOM_HK_FS_demo$BMI<23,"bHealthy weight",
                                          ifelse(MOM_HK_FS_demo$BMI<25,"cOverweight","Obesity"))))
tabcat <- colnames(MOM_HK_FS_demo)[c(5,6,8,10,11,13:25,27:30,34)]
tabcon <- colnames(MOM_HK_FS_demo)[c(9,33)]
demotab <- CreateTableOne(vars = c(tabcat,tabcon), 
                          data = MOM_HK_FS_demo, strata = "Area")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"MOM_HK_FS_demo_mom_03042022.xlsx")


# Foshan BMI----
FS_demo_F <- data.frame(read_excel("Data/FOSHAN data20220128/Father_baseline_20220409.xls"))
FS_demo_F$BMI <- round(FS_demo_F$Weight/(FS_demo_F$Height/100)^2,1)
FS_demo_F <- FS_demo_F[,-c(10:19,73:80)]
FS_demo_F <- FS_demo_F %>% filter(substr(Sub_ID,4,4)==3)
FS_demo_F <- left_join(FS_demo_F,FS_FA,"Sub_ID")
FS_demo_F$Emu <- as.numeric(apply(FS_demo_F[,64:66],1,sum))
FS_demo_F$Sweet <- as.numeric(apply(FS_demo_F[,70:72],1,sum))
FS_demo_F$Other_FA <- as.numeric(apply(FS_demo_F[,67:69],1,sum))
FS_demo_F$FA_all <- as.numeric(apply(FS_demo_F[,64:72],1,sum))
colnames(FS_demo_F)[c(59:61)] <- c("MoverV", "VoverM", "VMeven")
FS_demo_F$VMgroup <- as.character(apply(FS_demo_F[,59:61],1,
                                 function(x) colnames(FS_demo_F)[c(59:61)][x=="Yes"]))
FS_demo_F$VMgroup[grepl("c",FS_demo_F$VMgroup)] <- NA

# FS_demo_F[,c(8:10,12:58)] <- apply(FS_demo_F[,c(8:10,12:58)],2,
#                                    function(x) ifelse(x=="Yes","2Yes","2No"))


FS_demo_F$Income[FS_demo_F$Income=="< 2,000"|FS_demo_F$Income=="2,000 - 4,999"|FS_demo_F$Income_fam_new=="5,000 - 9,999"|FS_demo_F$Income_fam_new=="10,000 - 14,999"|FS_demo_F$Income_fam_new=="15,000 - 29,999"] <- "<30,000"
FS_demo_F$Income[FS_demo_F$Income=="30,000 - 49,999"] <- "30,000 - 49,999"

FS_demo_F$Income <- factor(as.factor(FS_demo_F$Income),
                                   levels(as.factor(FS_demo_F$Income))[c(1,3,4,2)])
FS_demo_F$Edu[FS_demo_F$Edu=="High school or bellow"] <- "Highschool_bellow"
FS_demo_F$Edu[FS_demo_F$Edu=="Bachelor or equivalent"] <- "Bachelor/degree"
FS_demo_F$Edu[FS_demo_F$Edu=="Master or above"] <- "Master_above"
FS_demo_F$Edu <- factor(as.factor(FS_demo_F$Edu),
                        levels(as.factor(FS_demo_F$Edu))[c(2,1,3)])
substr(FS_demo_F$Work_stat,1,1) <- toupper(substr(FS_demo_F$Work_stat,1,1) )
FS_demo_F$Work_stat <- factor(as.factor(FS_demo_F$Work_stat),
                              levels(as.factor(FS_demo_F$Work_stat))[c(1,4,2,3)])
FS_demo_F$PE_Group <- ifelse(FS_demo_F$PE>2,"≥3",
                             paste(FS_demo_F$PE))
FS_demo_F$PE_Group <- factor(as.factor(FS_demo_F$PE_Group),
                             levels(as.factor(FS_demo_F$PE_Group))[c(2,3,1)])
FS_demo_F$BMI_Group <- ifelse(is.na(FS_demo_F$BMI),NA,
                                   ifelse(FS_demo_F$BMI<18.5,"aUnderweight",
                                          ifelse(FS_demo_F$BMI<23,"bHealthy weight",
                                                 ifelse(FS_demo_F$BMI<25,"cOverweight","Obesity"))))



uni_c <- data.frame()
for (i in colnames(FS_demo_F)[c(5:11,77,12:26,28,29,40,64:76)]){
  
  if (all(FS_demo_F[,i]=="No")){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~",i)),
                  FS_demo_F,family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                    exp(cbind(OR = coef(temp_a), confint(temp_a,level = 0.95))))
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
uni_c$Pr...t.. <- ifelse(uni_c$Pr...t..<0.001, paste("<0.001"),paste(uni_c$Pr...t..))
# uni_c$Pr...z.. <- ifelse(uni_c$Pr...z..<0.001, paste("<0.001"),paste(uni_c$Pr...z..))

uni_c <- uni_c %>% filter(Pr...t..=="<0.001"|Pr...t..<0.05)
# uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.05)

uni_c$Var <- rownames(uni_c)
colnames(uni_c) <- c("Est","p.value","OR","Lower","Upper","Var")
writexl::write_xlsx(uni_c[,c(6,1:5)],"FS_Father_BMI_Uni.xlsx")



MOMdemo_FS <- data.frame( read_excel("Data/FOSHAN data20220128/Mother_baseline_20220409.xls"))
MOMdemo_FS$Weight <- as.numeric(MOMdemo_FS$Weight)
MOMdemo_FS$BMI <- round(MOMdemo_FS$Weight/(MOMdemo_FS$Height/100)^2,1)
MOMdemo_FS <- MOMdemo_FS[,-c(10:19,73:80)]
MOMdemo_FS <- MOMdemo_FS %>% filter(substr(Sub_ID,4,4)==1)
MOMdemo_FS <- left_join(MOMdemo_FS,FS_FA,"Sub_ID")
MOMdemo_FS$Emu <- as.numeric(apply(MOMdemo_FS[,65:67],1,sum))
MOMdemo_FS$Sweet <- as.numeric(apply(MOMdemo_FS[,71:73],1,sum))
MOMdemo_FS$Other_FA <- as.numeric(apply(MOMdemo_FS[,68:70],1,sum))
MOMdemo_FS$FA_all <- as.numeric(apply(MOMdemo_FS[,65:73],1,sum))
colnames(MOMdemo_FS)[c(59:61)] <- c("MoverV", "VoverM", "VMeven")
MOMdemo_FS$VMgroup <- as.character(apply(MOMdemo_FS[,59:61],1,
                                        function(x) colnames(MOMdemo_FS)[c(59:61)][x=="Yes"]))
# MOMdemo_FS$VMgroup[grepl("c",MOMdemo_FS$VMgroup)] <- NA

# FS_demo_F[,c(8:10,12:58)] <- apply(FS_demo_F[,c(8:10,12:58)],2,
#                                    function(x) ifelse(x=="Yes","2Yes","2No"))


MOMdemo_FS$Income[MOMdemo_FS$Income=="< 2,000"|MOMdemo_FS$Income=="2,000 - 4,999"|MOMdemo_FS$Income_fam_new=="5,000 - 9,999"|MOMdemo_FS$Income_fam_new=="10,000 - 14,999"|MOMdemo_FS$Income_fam_new=="15,000 - 29,999"] <- "<30,000"
MOMdemo_FS$Income[MOMdemo_FS$Income=="30,000 - 49,999"] <- "30,000 - 49,999"

MOMdemo_FS$Income <- factor(as.factor(MOMdemo_FS$Income),
                           levels(as.factor(MOMdemo_FS$Income))[c(1,3,4,2)])
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="High school or bellow"] <- "Highschool_bellow"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Bachelor or equivalent"] <- "Bachelor/degree"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Master or above"] <- "Master_above"
MOMdemo_FS$Edu <- factor(as.factor(MOMdemo_FS$Edu),
                        levels(as.factor(MOMdemo_FS$Edu))[c(2,1,3)])
substr(MOMdemo_FS$Work_stat,1,1) <- toupper(substr(MOMdemo_FS$Work_stat,1,1) )
MOMdemo_FS$Work_stat <- factor(as.factor(MOMdemo_FS$Work_stat),
                              levels(as.factor(MOMdemo_FS$Work_stat))[c(1,4,2,3)])
MOMdemo_FS$PE_Group <- ifelse(MOMdemo_FS$PE>2,"≥3",
                             paste(MOMdemo_FS$PE))
MOMdemo_FS$PE_Group <- factor(as.factor(MOMdemo_FS$PE_Group),
                             levels(as.factor(MOMdemo_FS$PE_Group))[c(2,3,1)])
MOMdemo_FS$BMI_Group <- ifelse(is.na(MOMdemo_FS$BMI),NA,
                              ifelse(MOMdemo_FS$BMI<18.5,"aUnderweight",
                                     ifelse(MOMdemo_FS$BMI<23,"bHealthy weight",
                                            ifelse(MOMdemo_FS$BMI<25,"cOverweight","Obesity"))))





uni_c <- data.frame()
for (i in colnames(MOMdemo_FS)[c(5:10,78,11:25,27,28,39,65:78)]){
  
  if (all(MOMdemo_FS[,i]=="No")){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~",i)),
                  MOMdemo_FS,family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                    exp(cbind(OR = coef(temp_a), confint(temp_a,level = 0.95))))
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
uni_c$Pr...t.. <- ifelse(uni_c$Pr...t..<0.001, paste("<0.001"),paste(uni_c$Pr...t..))
# uni_c$Pr...z.. <- ifelse(uni_c$Pr...z..<0.001, paste("<0.001"),paste(uni_c$Pr...z..))

uni_c <- uni_c %>% filter(Pr...t..=="<0.001"|Pr...t..<0.05)
# uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.05)

uni_c$Var <- rownames(uni_c)
colnames(uni_c) <- c("Est","p.value","OR","Lower","Upper","Var")
writexl::write_xlsx(uni_c[,c(6,1:5)],"FS_Mother_BMI_Uni.xlsx")



