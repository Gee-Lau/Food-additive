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
library(ComplexHeatmap)


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
HK_FS_FA$Area <- ifelse(substr(HK_FS_FA$Sub_ID,1,1)=="P","HK","FS")
HK_FS_FA$Sex <- ifelse(substr(HK_FS_FA$Sub_ID,4,4)=="1","Female","Male")
HK_FS_FA_Mom <- HK_FS_FA %>% filter(Sex=="Female") %>% select(-Sex) 
HK_FS_FA_Mom <- HK_FS_FA_Mom %>% 
  mutate(SumFA = as.numeric(apply(HK_FS_FA_Mom[,c(2:10)],1,sum))) %>% 
  filter(SumFA != 0) %>% select(-SumFA)


GBAdis <- vegdist(HK_FS_FA_Mom[,c(2:10)],na.rm = T, method = "bray") 

# library(ade4)
# GBApcoa <- dudi.pco(GBAdis, scan = FALSE,nf=3)
# GBApcoa_eig <- (GBApcoa$eig)[1:2] / sum(GBApcoa$eig)
# # GBAsite$FA <- rownames(GBApcoa)
# GBAsite <- data.frame({GBApcoa$li})[1:2]
# sample_site$names <- rownames(sample_site)
# names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

GBAnmds <- monoMDS(GBAdis)

HKFS.md <- data.frame(cbind(Sub_ID = HK_FS_FA_Mom$Sub_ID,
                              Area = HK_FS_FA_Mom$Area,
                            GBAnmds[["points"]]))
HKFS.md$MDS1 <- as.numeric(HKFS.md$MDS1)
HKFS.md$MDS2 <- as.numeric(HKFS.md$MDS2)



HK.FA.md <- vegdist(t(HK_FS_FA_Mom[HK_FS_FA_Mom$Area=="HK",c(2:10)])
                    ,na.rm = T, method = "bray")
HK.FA.md <- monoMDS(HK.FA.md)
HK.FA.md <- data.frame(cbind(FA = colnames(HK_FS_FA_Mom)[2:10],
                            Area = "HK",
                            HK.FA.md[["points"]]))

FS.FA.md <- vegdist(t(HK_FS_FA_Mom[HK_FS_FA_Mom$Area=="FS",c(2:10)])
                    ,na.rm = T, method = "bray")
FS.FA.md <- monoMDS(FS.FA.md)
FS.FA.md <- data.frame(cbind(FA = colnames(HK_FS_FA_Mom)[2:10],
                             Area = "FS",
                             FS.FA.md[["points"]]))
GBA_FA_md <- rbind(HK.FA.md,FS.FA.md)
colnames(GBA_FA_md)[3:4] <- c("FA_MDS1","FA_MDS2")
GBA_FA_md$FA_MDS1 <- as.numeric(GBA_FA_md$FA_MDS1)
GBA_FA_md$FA_MDS2 <- as.numeric(GBA_FA_md$FA_MDS2)
cent <- aggregate(cbind(MDS1,MDS2) ~ Area, data = HKFS.md, FUN = mean)
colnames(cent)[2:3] = c('cent1','cent2')
GBA_FA_md <- merge(GBA_FA_md, cent, by = 'Area', all.x = TRUE)

GBAper <- adonis2(HK_FS_FA_Mom[,c(2:10)] ~ Area, 
                  data = HK_FS_FA_Mom, 
                  permutation=999)

GBA_FA_md_plt <- ggplot(HKFS.md, aes(MDS1, MDS2, col = Area)) +
  geom_point() +
  scale_color_manual(name = "Region",
                     label = c("Foshan","Hong Kong"),
                     values=c("#00A087FF","#8491B4FF")) +
  scale_fill_manual(name = "",
                     label = c("Foshan","Hong Kong"),
                     values=c("#00A087FF","#8491B4FF")) +
  scale_x_continuous(breaks = seq(-2,2,1)) +
  scale_y_continuous(breaks = seq(-1,1,0.5)) +
  ylim(-1,1) +
  xlim(-2,2) +
  # geom_segment(mapping = aes(xend = cent1, yend = cent2),
  #              alpha = 0.05,
  #              show.legend = FALSE) +
  geom_point(data = cent, aes(cent1, cent2),
             size = 5, show.legend = FALSE) +
  geom_segment(data = GBA_FA_md,
               aes(x = cent1, y = cent2, 
                   xend = FA_MDS1, yend = FA_MDS2,
                   col = Area), 
               size = 0.8, alpha = 0.4) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  geom_point(data = GBA_FA_md, shape = 24, 
             aes(FA_MDS1, FA_MDS2, 
                 fill = Area),
             size = 4, color = "black", show.legend = FALSE) +
  geom_text(data = GBA_FA_md,
            aes(x = FA_MDS1-0.05, y = FA_MDS2+0.08),
            col = "gray30",
            label = GBA_FA_md$FA, show.legend = FALSE) +
  theme_classic() +
  theme(legend.position = c(0.87, 0.95),
        legend.background = element_rect(color ="black"),
        axis.text = element_text(color = "black", size =14),
        axis.title = element_text(size = 16, face = "bold")) +
  ggtitle("Permanova p<0.001, R = 0.11") 





# Baseline demographic----
MOMdemo_FS <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/MotherBaseline2022-04-16-372.xls"))
MOM_HK_FS_demo <- MOMcomb_sub %>% select(colnames(MOMcomb_sub)[colnames(MOMcomb_sub) %in% colnames(MOMdemo_FS)]) %>%
  filter(substr(Sub_ID,4,4)=="1")
# MOMdemo_FS$Weight <- str_replace(MOMdemo_FS$Weight,"（孕早期）","")
MOMdemo_FS$Weight <- as.numeric(MOMdemo_FS$Weight)
MOMdemo_FS <- MOMdemo_FS %>% select(colnames(MOMdemo_FS)[colnames(MOMdemo_FS) %in% colnames(MOM_HK_FS_demo)]) %>% 
  mutate(Area = "FS")
MOM_HK_FS_demo$Area <- "HK"
MOMdemo_FS[MOMdemo_FS=="NA"] <- NA
a$Income_fam_new[a$Income=="< 2,000"|a$Income=="2,000 - 4,999"] <- "<5,000"
a$Income_fam_new[a$Income=="5,000 - 9,999"] <- "5,000 - 9,999"
a$Income_fam_new[a$Income=="10,000 - 14,999"|a$Income=="15,000 - 29,999"|a$Income=="30,000 - 49,999"|a$Income=="50,000 - 99,999"] <- ">10,000"
MOMdemo_FS$Income_fam_new <- a$Income_fam_new[match(MOMdemo_FS$Sub_ID,a$Sub_ID)]
MOMdemo_FS$Income_fam_new <- factor(as.factor(MOMdemo_FS$Income_fam_new),
                                levels(as.factor(MOMdemo_FS$Income_fam_new))[c(1,3,2)])
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="High school or bellow"] <- "Highschool_bellow"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Bachelor or equivalent"] <- "Bachelor_degree"
MOMdemo_FS$Edu[MOMdemo_FS$Edu=="Master or above"] <- "Master_above"

MOMdemo_FS$Edu <- factor(as.factor(MOMdemo_FS$Edu),
                     levels(as.factor(MOMdemo_FS$Edu))[c(2,1,3)])
substr(MOMdemo_FS$Work_stat,1,1) <- toupper(substr(MOMdemo_FS$Work_stat,1,1) )
MOMdemo_FS$Work_stat <- factor(as.factor(MOMdemo_FS$Work_stat),
                           levels(as.factor(MOMdemo_FS$Work_stat))[c(1,4,2,3)])
MOMdemo_FS$PE_Group <- ifelse(MOMdemo_FS$PE==0,"1No",
                                  "2Yes")
MOMdemo_FS$PE_Group <- factor(as.factor(MOMdemo_FS$PE_Group),
                          levels(as.factor(MOMdemo_FS$PE_Group))[c(2,3,1)])
# MOMdemo_FS[,c(9:26)] <- apply(MOMdemo_FS[,c(9:26)], 2,
#                                function(x) ifelse(x=="NA",NA,
#                                                   ifelse(x=="No","1No","2Yes")))
# MOM_HK_FS_demo <- rbind(MOM_HK_FS_demo,MOMdemo_FS)

MOMdemo_FS$BMI <- round(MOMdemo_FS$Weight/(MOMdemo_FS$Height/100)^2,1)
MOMdemo_FS$BMI_Group <- ifelse(is.na(MOMdemo_FS$BMI),NA,
                            ifelse(MOMdemo_FS$BMI<18.5,"aUnderweight",
                                   ifelse(MOMdemo_FS$BMI<23,"bHealthy weight",
                                          ifelse(MOMdemo_FS$BMI<25,"cOverweight","Obesity"))))



tabcat <- colnames(MOM_HK_FS_demo)[c(5,6,8,10,11,13:25,27:30,34)]
tabcon <- colnames(MOM_HK_FS_demo)[c(9,33)]
demotab <- CreateTableOne(vars = c(tabcat,tabcon), 
                          data = MOM_HK_FS_demo, strata = "Area")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"MOM_HK_FS_demo_mom_03042022.xlsx")


# Foshan BMI----
FS_demo_F <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/FatherBaseline2022-04-16-353.xls"))
FS_demo_F$BMI <- round(FS_demo_F$Weight/(FS_demo_F$Height/100)^2,1)
FS_demo_F <- FS_demo_F[,c(1:8,80,21:23,66:79,25,26,37,54:56)]
FS_demo_F <- FS_demo_F %>% filter(substr(Sub_ID,4,4)==3)
FS_demo_F <- left_join(FS_demo_F,FS_FA,"Sub_ID")
FS_demo_F[FS_demo_F=="NA"] <- NA
FS_demo_F$Emu <- as.numeric(apply(FS_demo_F[,33:35],1,sum))
FS_demo_F$Sweet <- as.numeric(apply(FS_demo_F[,39:41],1,sum))
FS_demo_F$Other_FA <- as.numeric(apply(FS_demo_F[,36:38],1,sum))
FS_demo_F$FA_all <- as.numeric(apply(FS_demo_F[,33:41],1,sum))
colnames(FS_demo_F)[c(30:32)] <- c("MoverV", "VoverM", "VMeven")
FS_demo_F$VMgroup <- as.character(apply(FS_demo_F[,30:32],1,
                                 function(x) colnames(FS_demo_F)[c(30:32)][x=="Yes"]))
FS_demo_F$VMgroup[grepl("c",FS_demo_F$VMgroup)] <- NA

# FS_demo_F[,c(8:10,12:58)] <- apply(FS_demo_F[,c(8:10,12:58)],2,
#                                    function(x) ifelse(x=="Yes","2Yes","2No"))
FS_demo_F$Income_fam_new[FS_demo_F$Income=="< 2,000"|FS_demo_F$Income=="2,000 - 4,999"] <- "<5,000"
FS_demo_F$Income_fam_new[FS_demo_F$Income=="5,000 - 9,999"] <- "5,000 - 9,999"
FS_demo_F$Income_fam_new[FS_demo_F$Income=="10,000 - 14,999"|FS_demo_F$Income=="15,000 - 29,999"|FS_demo_F$Income=="30,000 - 49,999"|FS_demo_F$Income=="50,000 - 99,999"] <- ">10,000"
FS_demo_F$Income_fam_new <- factor(as.factor(FS_demo_F$Income_fam_new),
                                    levels(as.factor(FS_demo_F$Income_fam_new))[c(1,3,2)])
FS_demo_F$Edu[FS_demo_F$Edu=="High school or bellow"] <- "Highschool_bellow"
FS_demo_F$Edu[FS_demo_F$Edu=="Bachelor or equivalent"] <- "Bachelor/degree"
FS_demo_F$Edu[FS_demo_F$Edu=="Master or above"] <- "Master_above"
FS_demo_F$Edu <- factor(as.factor(FS_demo_F$Edu),
                        levels(as.factor(FS_demo_F$Edu))[c(2,1,3)])
substr(FS_demo_F$Work_stat,1,1) <- toupper(substr(FS_demo_F$Work_stat,1,1) )
FS_demo_F$Work_stat <- factor(as.factor(FS_demo_F$Work_stat),
                              levels(as.factor(FS_demo_F$Work_stat))[c(1,4,2,3)])
FS_demo_F$PE_Group <- ifelse(FS_demo_F$PE>0,"2Yes",
                             "1No")

FS_demo_F$BMI_Group <- ifelse(is.na(FS_demo_F$BMI),NA,
                                   ifelse(FS_demo_F$BMI<18.5,"aUnderweight",
                                          ifelse(FS_demo_F$BMI<23,"bHealthy weight",
                                                 ifelse(FS_demo_F$BMI<25,"cOverweight","Obesity"))))
FS_FA_past_F <-  data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/FatherDAT1_2022-04-16.xls"))
colnames(FS_FA_past_F)[1] <- "Sub_ID"
FS_demo_F[,12:29] <- apply(FS_demo_F[,12:29],2,function(x)
  ifelse(x=="Yes","2Yes","1No"))
FS_demo_F <- left_join(FS_demo_F,FS_FA_past_F,"Sub_ID")

FS_demo_F$Alco <- FS_Alco$Alcohol.Unit.Per.Week[match(FS_demo_F$Sub_ID,
                                                      FS_Alco$Sub_ID)]

uni_c <- data.frame()
for (i in colnames(FS_demo_F)[c(4:6,8,10,11,47,88,46,49,12:29,33:45,
                                50,55,51:54,56:87)]){
  
  if (all(FS_demo_F[,i]=="1No"|FS_demo_F[,i]=="No"|
          FS_demo_F[,i][complete.cases(FS_demo_F[,i])]=="Fresh")){
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

# uni_c <- uni_c %>% filter(Pr...t..=="<0.001"|Pr...t..<0.05)
# uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.05)

uni_c$Var <- rownames(uni_c)
colnames(uni_c) <- c("Est","p.value","OR","Lower","Upper","Var")
writexl::write_xlsx(uni_c[,c(6,1:5)],"Results/FA_FS/FS_Father_BMI_Uni_16052022.xlsx")


mglmtb <- data.frame()
for(i in colnames(FS_demo_F)[c(33,35,41,42,45)]){
  mglm <- glm(as.formula(paste0("BMI ~ Age + Alco + Digest...66 + HG + HT +
                                ProcVeg_10_18y + Fastfood_5_10y +
                                SoftDrink_5_10y +", i)),
              FS_demo_F,family = "gaussian")
  temp.db <- cbind(summary(mglm)$coefficients[,c(1)],confint(mglm),
                   summary(mglm)$coefficients[,c(4)])
  temp.db <- apply(temp.db, 2, function(x) round(x,4))
  temp.db <- data.frame(temp.db)
  temp.db <- temp.db[rownames(temp.db)==i,]
  mglmtb <- rbind(mglmtb,temp.db)
}


colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))

mglmtb$Var <- rownames(mglmtb)

writexl::write_xlsx(mglmtb,"Results/FA_FS/FS_Father_BMI_multi_16052022.xlsx")


# MOMdemo_FS <- data.frame( read_excel("Data/FOSHAN data20220128/Mother_baseline_20220409.xls"))
# MOMdemo_FS$Weight <- as.numeric(MOMdemo_FS$Weight)
# MOMdemo_FS$BMI <- round(MOMdemo_FS$Weight/(MOMdemo_FS$Height/100)^2,1)
# MOMdemo_FS <- MOMdemo_FS[,-c(11:20,40:48,51:69,74:81)]
MOMdemo_FS <- left_join(MOMdemo_FS,FS_FA,"Sub_ID")
MOMdemo_FS$Emu <- as.numeric(apply(MOMdemo_FS[,31:33],1,sum))
MOMdemo_FS$Sweet <- as.numeric(apply(MOMdemo_FS[,37:39],1,sum))
MOMdemo_FS$Other_FA <- as.numeric(apply(MOMdemo_FS[,34:36],1,sum))
MOMdemo_FS$FA_all <- as.numeric(apply(MOMdemo_FS[,31:39],1,sum))
a <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/MotherBaseline2022-04-16-372.xls"))

colnames(a)[c(70:72)] <- c("MoverV", "VoverM", "VMeven")
a$VMgroup <- as.character(apply(a[,70:72],1,
                                        function(x) colnames(a)[c(70:72)][x=="Yes"]))
# MOMdemo_FS$VMgroup[grepl("c",MOMdemo_FS$VMgroup)] <- NA

# FS_demo_F[,c(8:10,12:58)] <- apply(FS_demo_F[,c(8:10,12:58)],2,
#                                    function(x) ifelse(x=="Yes","2Yes","2No"))
MOMdemo_FS$VMgroup <- a$VMgroup[match(MOMdemo_FS$Sub_ID,a$Sub_ID)]
MOMdemo_FS$BF[MOMdemo_FS$BF=="NA"] <- NA
MOMdemo_FS$BF <- ifelse(MOMdemo_FS$BF=="Yes",
                        "2Yes","1No")

FS_FA_past <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/MotherDAT12022-04-16.xls"))
colnames(FS_FA_past)[1] <- "Sub_ID"
MOMdemo_FS <- left_join(MOMdemo_FS,FS_FA_past,"Sub_ID")




FS_Alco <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/FoShan Alcohol Unit Calculation (1)(1).xlsx"))
MOMdemo_FS$Alco <- FS_Alco$Alcohol.Unit.Per.Week[match(MOMdemo_FS$Sub_ID,
                                                       FS_Alco$Sub_ID)]

uni_c <- data.frame()
for (i in colnames(MOMdemo_FS)[c(4:8,84,85,45,44,9:22,24:26,
                                 31:43,46,51,47:50,52:83)]){
  
  if (all(MOMdemo_FS[,i]=="No"|MOMdemo_FS[,i]=="1No"|
          MOMdemo_FS[,i]=="Fresh")){
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


uni_c$Var <- rownames(uni_c)
colnames(uni_c) <- c("Est","p.value","OR","Lower","Upper","Var")
writexl::write_xlsx(uni_c[,c(6,1:5)],"Results/FA_FS/FS_Mother_BMI_Uni_16052022.xlsx")

mglmtb <- data.frame()
for(i in colnames(MOMdemo_FS)[c(31:33,37:41,43)]){
  mglm <- glm(as.formula(paste0("BMI ~ Age + Alco + Work_stat + Dairy_5_10y +
                                Fastfood_5_10y + ", i)),
              MOMdemo_FS,family = "gaussian")
  temp.db <- cbind(summary(mglm)$coefficients[,c(1)],confint(mglm),
                  summary(mglm)$coefficients[,c(4)])
  temp.db <- apply(temp.db, 2, function(x) round(x,4))
  temp.db <- data.frame(temp.db)
  temp.db <- temp.db[rownames(temp.db)==i,]
  mglmtb <- rbind(mglmtb,temp.db)
}


colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))

mglmtb$Var <- rownames(mglmtb)

writexl::write_xlsx(mglmtb,"Results/FA_FS/FS_Mother_BMI_multi_16052022.xlsx")




FS_FA_raw <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/Mother-DAT2_2022-04-16.xls"))
FS_FA_raw <- rbind(FS_FA_raw,
                   data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/FatherDAT2_2022-04-16.xls")))

FS_FA_raw <- FS_FA_raw %>% filter(FS_FA_raw$研究纳入编号. %!in% FS_FA$Sub_ID)
writexl::write_xlsx(FS_FA_raw,"Data/FOSHAN data20220128/Foshan_raw_16042022/FS_FA_Raw_04052022.xlsx")

# FS maternal outcomes
FS_FA_Del <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/delivery2(1).xlsx")) %>% 
  filter(substr(Sub_ID,4,4)==2)
FS_FA_Del <- FS_FA_Del[,c(1,580:639)]
FS_FA_Del <- FS_FA_Del[,-c(8,9,18:28,35,38,54)]
FS_FA_BL <- data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/T123-0621.xlsx"))[,c(1:9)]
colnames(FS_FA_BL) <- c("Samp_ID","Height","Pre_Weight","EDOB","Tri",
                        "GA","Date","Weight","Weight_Unit")
colnames(FS_FA_Del) <- c(colnames(PNdb)[c(1,103:109)],"Preterm",
                         colnames(PNdb)[c(110:145)])
FS_FA_Del$Pre_Weight <- FS_FA_BL$Pre_Weight[match(FS_FA_Del$Sub_ID,FS_FA_BL$Samp_ID)]


FS_FA_T2 <-  data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/T123-0621.xlsx"))[,c(1,12:16)]
FS_FA_T3 <-  data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/T123-0621.xlsx"))[,c(1,19:23)]
colnames(FS_FA_T2) <-  c("Samp_ID","Tri",
                         "GA","Date","Weight","Weight_Unit")
FS_FA_DM <-  data.frame(read_excel("Data/FOSHAN data20220128/Foshan_raw_16042022/T2 and T3.xlsx"))

FS_FA_Del$T2_weight <- FS_FA_T2$Weight[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_T2$Samp_ID,5,9)
)]

FS_FA_Del$T3_weight <- FS_FA_T3$Weight[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_T3$Samp_ID,5,9)
)]

FS_FA_Del$T2_GA <- FS_FA_T2$GA[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_T2$Samp_ID,5,9)
)]
FS_FA_Del$T3_GA <- FS_FA_T3$GA[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_T3$Samp_ID,5,9)
)]
FS_FA_Del[,c(154,155)][FS_FA_Del[,c(154,155)]==0] <- NA

FS_FA_Del$T2_GA[FS_FA_Del$T2_GA>FS_FA_Del$T3_GA|FS_FA_Del$T2_GA==FS_FA_Del$T3_GA] <- NA

FS_FA_DM$DM...11 <- ifelse(is.na(FS_FA_DM$DM...11),"No","Yes")
FS_FA_DM$DM...5 <- ifelse(is.na(FS_FA_DM$DM...5),"No","Yes")

FS_FA_Del$T2_GDM <- FS_FA_DM$DM...5[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_DM$Study.ID...1,5,9)
)]
FS_FA_Del$T3_GDM <- FS_FA_DM$DM...5[match(
  substr(FS_FA_Del$ID,5,9),
  substr(FS_FA_DM$Study.ID...1,5,9)
)]

FS_FA_Del$GDM_all <- as.character(apply(FS_FA_Del[,c(5,50,51)], 1,
                                        function(x) 
                                          ifelse(any(x=="Yes"),"Yes","No")))
FS_FA_Del$Sub_ID <- paste0("FFP1",substr(FS_FA_Del$ID,5,9))
FS_FA_Del <- left_join(MOMdemo_FS,FS_FA_Del,"Sub_ID")



FS_FA_Del$GDM_Bi <- ifelse(FS_FA_Del$GDM_all=="Yes",1,0)
FS_FA_Del$Preterm_bi <- ifelse(is.na(FS_FA_Del$Preterm),0,1)


for(i in 31:43){
  FS_FA_Del[,i+108] <- ifelse(ntile(FS_FA_Del[,i],4)==4,
                              "2High","1Normal")
  colnames(FS_FA_Del)[i+108] <- paste0(colnames(FS_FA_Del)[i],"_cut")
}

FS_FA_Del <- FS_FA_Del %>% mutate(
  T2_WC = round((as.numeric(T2_weight)-as.numeric(Pre_Weight))/T2_GA,1),
  T3_WC = round((as.numeric(T3_weight)-as.numeric(Pre_Weight))/T3_GA,1)
)


FS_FA_Del$T2_WC[FS_FA_Del$T2_WC<0|FS_FA_Del$T2_WC>1] <- NA
FS_FA_Del$T3_WC[FS_FA_Del$T3_WC<0|FS_FA_Del$T3_WC>1] <- NA
set.seed(100)
FS_FA_Del$GBS_Bi <- ifelse(FS_FA_Del$GBS=="NA",NA,
                           ifelse(FS_FA_Del$GBS=="Yes",1,0))

FS_FA_Del[,c(113:121)] <- apply(FS_FA_Del[,c(113:121)],2,
                                function(x) ifelse(x=="No",0,1))
uni_c <- data.frame()
for(j in colnames(FS_FA_Del)[c(107)]){
  if (length(unique(FS_FA_Del[,j][complete.cases(FS_FA_Del[,j])]))==1) {
    next
  }else{
    
    for (i in colnames(FS_FA_Del)[c(4:8,84,85,45,44,9:22,24:26,
                                    139:151)]){
      if (all(FS_FA_Del[,i]=="No"|FS_FA_Del[,i]=="1No"|
              FS_FA_Del[,i]=="Fresh")){
        print(i)
      }
      else {
        print(i)
        temp_a <- glm(as.formula(paste0(j, "~",i)),
                      FS_FA_Del,family ="gaussian")
        temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                        confint(temp_a,level = 0.95),
                        summary(temp_a)$coefficients[,c(4)],
                        Outcome = j,
                        Var = i)
        uni_c <- rbind(uni_c,temp_b)
      }
    }
  }
}

for(j in colnames(FS_FA_Del)[c(113:121,138)]){
  if (length(unique(FS_FA_Del[,j][complete.cases(FS_FA_Del[,j])]))==1){
    next
  }else{
    for (i in colnames(FS_FA_Del)[c(4:8,84,85,45,44,9:22,24:26,
                                    139:151)]){
      if (all(FS_FA_Del[,i]=="No"|FS_FA_Del[,i]=="1No"|
              FS_FA_Del[,i]=="Fresh")){
        print(i)
      }
      else {
        print(i)
        temp_a <- glm(as.formula(paste0(j, " ~",i)),
                      FS_FA_Del,family ="binomial")
        temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                        exp(cbind(OR = coef(temp_a), 
                                  confint(temp_a))),
                        Outcome = j,
                        Var = i)
        temp_b <- temp_b[,c(3:5,2,6,7)]
        colnames(temp_b) <- colnames(uni_c)
        uni_c <- rbind(uni_c,temp_b)
      }
    }
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c[,c(1:4)] <- apply(uni_c[,c(1:4)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)
# uni_c$Pr...z.. <- ifelse(uni_c$Pr...z..<0.001, paste("<0.001"),paste(uni_c$Pr...z..))

colnames(uni_c) <- c("Est/OR","Lower","Upper","p.value","Outcome","Var")
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))
uni_c$Var_strat <- rownames(uni_c)
# writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_Gcomp_Uni_21062022_all.xlsx")
uni_c <- uni_c %>% filter(p.value<0.06)
writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_PN_Uni_22072022_sig.xlsx")


mglmtb <- data.frame()

for(j in colnames(FS_FA_Del)[c(107)]){
  for (i in colnames(FS_FA_Del)[31:43]){
    if (paste0(i,"_cut") %in% uni_c$Var[uni_c$Outcome==j]){
      print(i)
      temp_a <- glm(as.formula(paste0(j, "~",i, "_cut + Age  +", 
                                      paste0(uni_c$Var[!grepl("cut",uni_c$Var)&uni_c$Outcome==j],
                                             collapse = "+"))),
                    FS_FA_Del,family ="gaussian")
      temp_b <- data.frame(cbind(summary(temp_a)$coefficients[,c(1)],
                      confint(temp_a,level = 0.95),
                      summary(temp_a)$coefficients[,c(4)],
                      Outcome = j,
                      Var = i))
      mglmtb <- rbind(mglmtb,temp_b[2,])
    }
    else {
      next
    }
  }
}

for(j in colnames(FS_FA_Del)[c(113:121,138)]){
  for (i in colnames(FS_FA_Del)[31:43]){
    if(paste0(i,"_cut") %in% uni_c$Var[uni_c$Outcome==j]&
       length(uni_c$Var[!grepl("cut",uni_c$Var)&uni_c$Outcome==j])==0){
      temp_a <- glm(as.formula(paste0(j, "~",i, "_cut + Age  + SexOB")),
                    FS_FA_Del,family ="binomial")
      temp_b <- data.frame(cbind(summary(temp_a)$coefficients[,c(1,4)],
                                 exp(cbind(OR = coef(temp_a), 
                                           confint(temp_a))),
                                 Outcome = j,
                                 Var = i))
      temp_b <- temp_b[,c(3:5,2,6,7)]
      colnames(temp_b) <- colnames(mglmtb)
      mglmtb <- rbind(mglmtb,temp_b[2,])
    }else{
      if (paste0(i,"_cut") %in% uni_c$Var[uni_c$Outcome==j]){
        print(i)
        temp_a <- glm(as.formula(paste0(j, "~",i, "_cut + Age + SexOB +", 
                                        paste0(uni_c$Var[!grepl("cut",uni_c$Var)&uni_c$Outcome==j],
                                               collapse = "+"))),
                      FS_FA_Del,family ="binomial")
        temp_b <- data.frame(cbind(summary(temp_a)$coefficients[,c(1,4)],
                                   exp(cbind(OR = coef(temp_a), 
                                             confint(temp_a))),
                                   Outcome = j,
                                   Var = i))
        temp_b <- temp_b[,c(3:5,2,6,7)]
        colnames(temp_b) <- colnames(mglmtb)
        mglmtb <- rbind(mglmtb,temp_b[2,])
      }
      else {
        next
      }
    }
  }
}





colnames(mglmtb) <-  c("Est/OR","Lower","Upper","p.value","Outcome","Var")

mglmtb[,1:4] <- apply(mglmtb[,1:4], 2, function(x) round(as.numeric(x),4))


mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))

mglmtb$Var_Strat <- rownames(mglmtb)
writexl::write_xlsx(mglmtb,"Results/FA_FS/FS_Mother_PN_Mul_22072022_all.xlsx")

# FS VS HK -----
# 1. Early Processed food intake

PDHsum$Var <- str_replace(PDHsum$Var,"\\-","_")

PDH_FS_sum <- FS_FA_Del[,c(45:83)]
PDH_FS_sum$BF <- ifelse(PDH_FS_sum$BF=="2Yes","Yes","No")
PDH_FS_sum$Infancy_food_source <- ifelse(PDH_FS_sum$Infancy_food_source=="Homemade","Yes","No")
PDH_FS_sum$Adolo_food_source <- ifelse(PDH_FS_sum$Adolo_food_source=="Fresh","Yes","No")

PDH_FS_sum <- data.frame(Var = colnames(PDH_FS_sum),
                         Perc = as.numeric(apply(PDH_FS_sum,2,function(x) 
                           round(sum(x=="Yes",na.rm = T)/sum(complete.cases(x))*100,1))))


PDH_GBA <- PDHsum 
colnames(PDH_GBA)[2] <- "HK"
PDH_GBA$FS <- PDH_FS_sum$Perc[match(PDH_GBA$Var,PDH_FS_sum$Var)]
# PDH_FS_sum$Var[1] <- "Breast feeding"
# PDH_FS_sum$Var[2] <- "Infancy_Homemade"
# PDH_FS_sum$Var[7] <- "Adolescent_Fresh"
PDH_GBA$Var <- str_replace(PDH_GBA$Var,"\\_",".")
PDH_GBA$Timepoint <- as.character(lapply(PDH_GBA$Var,
                                         function(x) strsplit(x,"\\.")[[1]][2]))
PDH_GBA$Timepoint <- str_replace(PDH_GBA$Timepoint,"y"," Years")
PDH_GBA$Timepoint <- str_replace(PDH_GBA$Timepoint,"m"," Months")
PDH_GBA <- PDH_GBA[c(1,2,7,3:6,8:nrow(PDH_GBA)),]
PDH_GBA$Timepoint[1:2] <- "Infancy"
PDH_GBA$Timepoint[3] <- "Adolescent"
PDH_GBA$Item <- NA
PDH_GBA$Item[1:3] <- c("Breast fed","Homemade food","Fresh food")
PDH_GBA$Item[4:7] <- rep(c("Home-grown vegetables"),4)
PDH_GBA$Item[8:11] <- rep(c("Processed dairy"),4)
PDH_GBA$Item[12:15] <- rep(c("Processed seafood"),4)
PDH_GBA$Item[16:19] <- rep(c("Processed carbohydrates"),4)
PDH_GBA$Item[20:23] <- rep(c("Processed fruit"),4)
PDH_GBA$Item[24:27] <- rep(c("Processed vegetables"),4)
PDH_GBA$Item[28:31] <- rep(c("Fast food"),4)
PDH_GBA$Item[32:35] <- rep(c("Soft drinks"),4)
PDH_GBA$Item[36:39] <- rep(c("Snacks"),4)
PDH_GBA$Order <- 1:nrow(PDH_GBA)
PDH_GBA$Item <- reorder(PDH_GBA$Item,PDH_GBA$Order)
PDH_GBA$Timepoint <- str_replace(PDH_GBA$Timepoint,"\\_","\\-")
col_fun = circlize::colorRamp2(c(0,20,100), c("white", "#fcbba1", "#a50f15"))
col_fun2 = circlize::colorRamp2(c(0,50,100), c("white", "#78c679", "#006837"))
lgd1 = Legend(col_fun = col_fun, direction = "horizontal",
             at = c(0, 30,60,100),
             title = "Early-life dietary\nhabits (%)",
             legend_height = unit(1.8, "cm"), border = "black")
lgd2 = Legend(col_fun = col_fun2, direction = "horizontal",
              at = c(0, 30,60,100), border = "black",
              legend_height = unit(1.8, "cm"))
lgd <- packLegend(lgd1,lgd2)
PDH_HM1 <- Heatmap(PDH_GBA[c(8:nrow(PDH_GBA)),c(2,3)], 
        na_col = "white",
        rect_gp = gpar(col = "gray30"),
        border = T,
        col = col_fun,
        row_split = PDH_GBA[8:nrow(PDH_GBA), 5],
        row_title_rot = 0,
        row_labels = PDH_GBA$Timepoint[8:nrow(PDH_GBA)],
        row_names_gp = gpar(fontsize = 10,
                            col = c(rep(c("#74a9cf","#3690c0","#0570b0","#0570b0"), 9))),
        show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        height = unit(16, "cm"),
        width = unit(5, "cm"),
        show_heatmap_legend = F) 

PDH_HM2 <- Heatmap(PDH_GBA[1:7,c(2,3)], 
                   na_col = "white",
                   rect_gp = gpar(col = "gray30"),
                   border = T,
                   col = col_fun2,
                   row_split = PDH_GBA[1:7, 5],
                   row_title_rot = 0,
                   row_labels = PDH_GBA$Timepoint[1:7],
                   row_names_gp = gpar(fontsize = 10,
                                       col = c("#74a9cf","#74a9cf","#045a8d",
                                               "#74a9cf","#3690c0","#0570b0","#0570b0")),
                   column_labels = c("Hong Kong\n","Foshan\n"),
                   column_names_rot = 0,
                   column_names_gp = gpar(fontsize = 14, hjust = 0),
                   column_names_side = "top",
                   column_gap = unit(1, "points"),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   height = unit(4, "cm"),
                   width = unit(5, "cm"),
                   show_heatmap_legend = F) 

ht_list = PDH_HM2 %v%PDH_HM1 

pdf("Results/FA_FS/PDH_GBA_04072022.pdf",width = 10,height = 9)
draw(ht_list,
     heatmap_legend_list = list(lgd))
invisible(dev.off())

# 2. The impact of food additives to clinical outcomes

FS_FA_Del <- data.frame(FS_FA_Del)
# Plot CMC and SUC
plot_list = list()
my.data <- data.frame(x = FS_FA_Del[,140][complete.cases(FS_FA_Del[,140])],
                      y = FS_FA_Del$T2_WC[complete.cases(FS_FA_Del[140])])
j <- colnames(FS_FA_Del)[140]
p1 <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
  # geom_violin() +
  geom_boxplot(width = 0.6, size = 1.2,
               position = position_dodge(0.9)) +
  scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="1Normal"&complete.cases(my.data$y))),
                             paste0("High \nN = ",sum(my.data$x=="2High"&complete.cases(my.data$y))))) +
  geom_point(shape=16, position=position_jitter(0.3),
             size = 2, alpha = 0.3) +
  scale_color_manual(values=c("#00A087FF","#E64B35FF")) +
  theme_classic()+
  stat_compare_means(method = "wilcox", size = 5,
                     label.x = 0.7,label.y = 0.9) +
  # stat_compare_means(comparisons = statecomp,label = "p.signif",
  #                    method = "wilcox",hide.ns = T) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 15)) +
  xlab(paste0("Annual CMC intake group")) +
  ylab("Trimester 2 weight gain\n(kg/gestational week)") 

plot_list <- c(plot_list, list(p1))
set.seed(100)
my.data <- data.frame(x = FS_FA_Del[,146][complete.cases(FS_FA_Del[,146])],
                      y = FS_FA_Del$T2_WC[complete.cases(FS_FA_Del[146])])
j <- colnames(FS_FA_Del)[146]
p2 <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
  # geom_violin() +
  geom_boxplot(width = 0.6, size = 1.2,
               position = position_dodge(0.9)) +
  scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="1Normal"&complete.cases(my.data$y))),
                             paste0("High \nN = ",sum(my.data$x=="2High"&complete.cases(my.data$y))))) +
  geom_point(shape=16, position=position_jitter(0.3),
             size = 2, alpha = 0.3) +
  scale_color_manual(values=c("#00A087FF","#E64B35FF")) +
  theme_classic()+
  stat_compare_means(method = "wilcox", size = 5,
                     label.x = 0.7,label.y = 0.9) +
  # stat_compare_means(comparisons = statecomp,label = "p.signif",
  #                    method = "wilcox",hide.ns = T) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 12),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 15)) +
  xlab(paste0("Annual SUC intake group")) +
  ylab("Trimester 2 weight gain\n(kg/gestational week)") 
plot_list <- c(plot_list, list(p2))


pdf(file="Results/T2_WC_FA_FS_05072022.pdf", width = 11.2, height = 7)
(p1+p2)/
  FS_mglmtb + plot_layout(heights = unit(c(4,1), c('inch', 'null')))
invisible(dev.off())


FS_mglmtb <- mglmtb[1:2,-c(5:7)]
rownames(FS_mglmtb) <- c("Carboxymethyl cellulose", 
                         "Sucralose")
colnames(FS_mglmtb) <- c("Estimate","95% CI (lower)",
                         "95%CI (upper)","p-value")

FS_mglmtb <- ggtexttable(FS_mglmtb, rows = rownames(FS_mglmtb),
                         theme = ttheme("blank", base_size = 18,
                                        padding = unit(c(1,1), "cm"),
                                        rownames.style = rownames_style(face = "bold",
                                                                        hjust=0, x=0.1,
                                                                        size = 18))) %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 5, linetype = 1) %>%
  tab_add_hline(at.row = c(3), row.side = "bottom", linewidth = 5, linetype = 1) 

# FS PN ----
FS_PN <- data.frame(read_xlsx("Data/FOSHAN data20220128/20220718-PN.xlsx"))

