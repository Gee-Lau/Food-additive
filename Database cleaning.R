library("reshape2")
library("Hmisc")
library("car")
library("dplyr")
library("tableone")
library("readxl")
library("tidyverse")
library("lubridate")
library("export")
library("rgl")
library("ggpubr")
library("ggplot2")
# Find extreme values for all food additives----

checksub <- matrix(nrow=30)
checksub <- data.frame(checksub)
for (i in c(3:12)){
  j <- colnames(MomFA)[i]
  checksub[,j] <- head(MomFA[,i][order(MomFA[,i],decreasing = T)],30)
}

for (i in 2:11){
  j <- colnames(checksub)[i]
  checksub[,i+10] <- head(MomFA$Sub_ID[order(MomFA[,j],decreasing = T)],30)
  colnames(checksub)[i+10] <- paste0(colnames(checksub)[i],"_ID")
}
checksub <- checksub[,-1]
checksub <- checksub[c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]
writexl::write_xlsx(checksub,"FA_Extreme values_03072021.xlsx")

# Sum up 
MomFAsum <- data.frame()

for (i in c(2:10)) {
  temp.db <- summary(MomFA[,i])
  MomFAsum <- rbind(MomFAsum,temp.db)
}
colnames(MomFAsum) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
MomFAsum$FA_cat <- colnames(MomFA)[2:10]
MomFAsum <- gather(MomFAsum, key = "FA_cat")
colnames(MomFAsum)[1] <- "Var"
MomFAsum$FA_cat <- rep(colnames(MomFA)[2:10],6)

# Clean up database (Dietary habit)----
clinicdb <- data.frame(read_xlsx("Data/0010_MOM_WL20210423A1 (Data only)_20210428 Clinical Team.xlsx",skip = 1))

# Past Dietary Habits
# Rename and select variables
PDHdb <- clinicdb[,c(1:50)]
colnames(PDHdb)[c(5:43,48)] <- c("BF","Infancy_food_source","HG_4-12m","HG_1-5y","HG_5-10y",
                                 "HG_10-18y","Adolo_food_source","Dairy_4-12m","Dairy_1-5y","Dairy_5-10y",
                                 "Dairy_10-18y", "Seafood_4-12m","Seafood_1-5y","Seafood_5-10y", "Seafood_10-18y","CarboHfood_4-12m",
                                 "CarboHfood_1-5y","CarboHfood_5-10y", "CarboHfood_10-18y", "ProcFruit_4-12m", "ProcFruit_1-5y",
                                 "ProcFruit_5-10y", "ProcFruit_10-18y", "ProcVeg_4-12m","ProcVeg_1-5y","ProcVeg_5-10y", "ProcVeg_10-18y",
                                 "Fastfood_4-12m","Fastfood_1-5y","Fastfood_5-10y", "Fastfood_10-18y", "SoftDrink_4-12m","SoftDrink_1-5y",
                                 "SoftDrink_5-10y", "SoftDrink_10-18y", "Snacks_4-12m","Snacks_1-5y","Snacks_5-10y", "Snacks_10-18y", "Complete")
PDHdb <- PDHdb %>% select(-c(1,3,4,44:47,49,50))
colnames(PDHdb)[1] <- "Sub_ID"

# Check if all binary data and convert to Eng
for (i in c(2,4:7,9:40)) {
  print(colnames(PDHdb)[i])
  print(table(PDHdb[,i]))
}

for (i in c(2,4:7,9:40)) {
  PDHdb[,i] <- ifelse(PDHdb[,i]=="否","No","Yes")
}

PDHdb$Infancy_food_source <- ifelse(grepl("商店", PDHdb$Infancy_food_source),"From shops","Homemade")
PDHdb$Adolo_food_source <- ifelse(grepl("便利店", PDHdb$Adolo_food_source),"Preserved", "Fresh")

# Build up data frame with summary dietary history overtime
PDHdb$Infancy_food_source <- ifelse(PDHdb$Infancy_food_source=="Homemade","Yes","No")
PDHdb$Adolo_food_source <- ifelse(PDHdb$Adolo_food_source=="Fresh","Yes","No")
PDHsum <- data.frame(Var = colnames(PDHdb)[-c(1,41)],
                     Perc = as.numeric(apply(PDHdb[,-c(1,41)],2,function(x) 
                       round(sum(x=="Yes",na.rm = T)/sum(complete.cases(x))*100,1))))

for (i in c(2,4:7,9:40)) {
  temp.db <- data.frame(rbind(cbind(Var = paste0(colnames(PDHdb[i]),"_M"), 
                                    t(table(PDHdb[,i][PDHdb$Sex=="M"]))),
                         cbind(Var = paste0(colnames(PDHdb[i]),"_F"), 
                               t(table(PDHdb[,i][PDHdb$Sex=="F"])))))
  PDHsum <- rbind(PDHsum,temp.db)
}

PDHsum <- rbind(PDHsum, cbind(Var = "AInfancy_food_source_M",
            No = as.numeric(table(PDHdb$Infancy_food_source[PDHdb$Sex=="M"])[1]),
            Yes = as.numeric(table(PDHdb$Infancy_food_source[PDHdb$Sex=="M"])[2])),
      cbind(Var =  "AInfancy_food_source_F",
            No = as.numeric(table(PDHdb$Infancy_food_source[PDHdb$Sex=="F"])[1]),
            Yes = as.numeric(table(PDHdb$Infancy_food_source[PDHdb$Sex=="F"])[2])),
      cbind(Var =  "Adolo_food_source_M",
            No = as.numeric(table(PDHdb$Adolo_food_source[PDHdb$Sex=="F"])[2]),
            Yes = as.numeric(table(PDHdb$Adolo_food_source[PDHdb$Sex=="F"])[1])),
      cbind(Var =  "Adolo_food_source_F",
            No = as.numeric(table(PDHdb$Adolo_food_source[PDHdb$Sex=="M"])[2]),
            Yes = as.numeric(table(PDHdb$Adolo_food_source[PDHdb$Sex=="M"])[1])))


PDHsum$Yes <- as.numeric(PDHsum$Yes)
PDHsum$No <- as.numeric(PDHsum$No)

PDHsum$Perc <- round(PDHsum$Yes/(PDHsum$Yes+PDHsum$No)*100,1)
PDHsum$Item <- as.character(lapply(as.list(PDHsum$Var),function(x) strsplit(x,"_")[[1]][1]))
PDHsum$Sex <- as.character(lapply(as.list(PDHsum$Var),function(x) strsplit(x,"_")[[1]][3]))
PDHsum$Sex[c(1,75,77)] <- "M"
PDHsum$Sex[c(2,76,78)] <- "F"

PDHsum$Var <- paste0(as.character(lapply(as.list(PDHsum$Var),function(x) strsplit(x,"_")[[1]][1])),"_",
                as.character(lapply(as.list(PDHsum$Var),function(x) strsplit(x,"_")[[1]][2])))
PDHsum$Timepoint <- as.character(lapply(as.list(PDHsum$Var),function(x) strsplit(x,"_")[[1]][2]))
PDHsum$Timepoint <- factor(as.factor(PDHsum$Timepoint),
                           levels(as.factor(PDHsum$Timepoint))[c(3,1,4,2,5:7)])
PDHsum <- PDHsum[order(PDHsum$Item,PDHsum$Timepoint),]
PDHsum <- PDHsum[order(PDHsum$Item,PDHsum$Timepoint),]
PDHsum$Var[PDHsum$Var=="Adolo_food"] <- "Adolescent_Fresh" 
PDHsum$Var[PDHsum$Var=="AInfancy_food"] <- "Infancy_Homemade" 
PDHsum$Var[grepl("BF",PDHsum$Var)] <- "Breast feeding" 

PDHsum$ID <- 1:length(PDHsum$Var)
PDHsum <- PDHsum %>% mutate(Item_od = fct_reorder(Var, desc(ID))) 




ggplot(data = PDHsum, aes( x = Timepoint,
                           y = Perc, color = Item)) +
  # scale_color_manual(values=c("steelblue2","red4")) +
  geom_point(size = 2) +
  geom_path(aes(group = Item)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))+
  ylab("Yes %") +
  xlab("Life stage")  


dadPDHdb <- data.frame(read_excel("Data/Baseline DAT 1 + Alcohol Unit (send to Wing on Jul 7).xlsx"))
dadPDHdb <- dadPDHdb %>% select(-record_id, -redcap_event_name) %>% 
  filter(substr(study_no,4,4)==3)
colnames(dadPDHdb)[c(1:40)] <- colnames(PDHdb)[c(1:40)] 
dadPDHdb <- dadPDHdb %>% filter(Sub_ID %in% dadcombdb$Sub_ID) %>% 
  filter(complete.cases(BF))
dadPDHdb <- dadPDHdb[,-c(41:45)]
dadPDHdb$Early_FA <- as.character(apply(dadPDHdb[,c(9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38)], 1,
                                     function(x) if (any(x[complete.cases(x)]==1)) {1} else {0}))



# Clean up database (Demographic)----
demodb <- clinicdb[,c(52:54,59,201,202,62,63,77:79,83,84,89,104:106,125,126,199,200,60,61,250,251)]
colnames(demodb) <- c("Sub_ID","Date","Region","Height","Weight","Weight_unit","DOB","Race","Edu",
                      "Work_stat","Work_pos","Income_ind","Income_fam","Pet","Smoke_ind","Smoke_fam",
                      "Alco","Drug","PE","FullTerm","Del_mode","Weight_T1","Weight_T1_Unit",
                      "Weight_T3","Weight_T3_Unit")

demodb$Region <- ifelse(demodb$Region=="香港","HK",
                        ifelse(demodb$Region=="廣州","GZ","Other"))

demodb$Weight[grepl("kg",tolower(demodb$Weight))] <- substr(demodb$Weight[grepl("kg",tolower(demodb$Weight))],1,
                                                            nchar(demodb$Weight[grepl("kg",tolower(demodb$Weight))])-2)
demodb$Weight <- as.numeric(demodb$Weight)
demodb$Weight[complete.cases(demodb$Weight)&demodb$Weight_unit=="市斤"] <- round(demodb$Weight[complete.cases(demodb$Weight)&demodb$Weight_unit=="市斤"]/2,1)
demodb$Weight[complete.cases(demodb$Weight)&demodb$Weight_unit=="磅 (lb)"] <- round(demodb$Weight[complete.cases(demodb$Weight)&demodb$Weight_unit=="磅 (lb)"]*0.453,1)
demodb$Height[grepl("cm",tolower(demodb$Height))] <- substr(demodb$Height[grepl("cm",tolower(demodb$Height))],1,
                                                            nchar(demodb$Height[grepl("cm",tolower(demodb$Height))])-2)

demodb$Weight_T1[grepl("kg",tolower(demodb$Weight_T1))] <- substr(demodb$Weight_T1[grepl("kg",tolower(demodb$Weight_T1))],1,
                                                            nchar(demodb$Weight_T1[grepl("kg",tolower(demodb$Weight_T1))])-2)
demodb$Weight_T1 <- as.numeric(demodb$Weight_T1)
demodb$Weight_T1[complete.cases(demodb$Weight_T1)&demodb$Weight_T1_Unit=="市斤"] <- round(demodb$Weight_T1[complete.cases(demodb$Weight_T1)&demodb$Weight_T1_Unit=="市斤"]/2,1)
demodb$Weight_T1[complete.cases(demodb$Weight_T1)&demodb$Weight_T1_Unit=="磅 (lb)"] <- round(demodb$Weight_T1[complete.cases(demodb$Weight_T1)&demodb$Weight_T1_Unit=="磅 (lb)"]*0.453,1)

demodb$Weight_T3[grepl("kg",tolower(demodb$Weight_T3))] <- substr(demodb$Weight_T3[grepl("kg",tolower(demodb$Weight_T3))],1,
                                                                  nchar(demodb$Weight_T3[grepl("kg",tolower(demodb$Weight_T3))])-2)
demodb$Weight_T3 <- as.numeric(demodb$Weight_T3)
demodb$Weight_T3[complete.cases(demodb$Weight_T3)&demodb$Weight_T3_Unit=="市斤"] <- round(demodb$Weight_T3[complete.cases(demodb$Weight_T3)&demodb$Weight_T3_Unit=="市斤"]/2,1)
demodb$Weight_T3[complete.cases(demodb$Weight_T3)&demodb$Weight_T3_Unit=="磅 (lb)"] <- round(demodb$Weight_T3[complete.cases(demodb$Weight_T3)&demodb$Weight_T3_Unit=="磅 (lb)"]*0.453,1)

demodb$Weight[demodb$Weight<20|demodb$Weight>200] <- NA


demodb$Height <- as.numeric(demodb$Height)
demodb$Height[demodb$Height<10&complete.cases(demodb$Height)] <-  demodb$Height[demodb$Height<10&complete.cases(demodb$Height)]*100
demodb <- cbind(demodb[,c(1:5,7)],
                BMI = round(demodb$Weight/(demodb$Height/100)^2,1),
                Age = round(as.numeric(ymd(demodb$Date)-ymd(demodb$DOB))/365),
                Group = NA,
                demodb[,c(8:21,22,24)])

demodb$Age[demodb$Age==0] <- NA # Remove uncertain age
demodb$BMI[demodb$BMI>50] <- NA
demodb$BMI[demodb$BMI<10] <- NA
demodb$Edu <- ifelse(grepl("大學",demodb$Edu),"Bachelor/degree",
                     ifelse(grepl("碩士",demodb$Edu),"Master_above","Highschool_bellow"))
demodb$Work_stat <- ifelse(grepl("全職家庭",demodb$Work_stat),"Housewife",
                           ifelse(grepl("全職",demodb$Work_stat),"Fulltime",
                                  ifelse(grepl("兼職",demodb$Work_stat),"Parttime","Others")))
demodb$Income_fam_new <- demodb$Income_fam
demodb$Income_fam_new[grepl("不方便|幣",demodb$Income_fam_new)] <- NA
demodb$Income_fam_new[demodb$Income_fam_new=="< 2,000"|demodb$Income_fam_new=="2,000 - 4,999"|demodb$Income_fam_new=="5,000 - 9,999"|demodb$Income_fam_new=="10,000 - 14,999"|demodb$Income_fam_new=="15,000 - 29,999"] <- "<30,000"
demodb$Income_fam_new[demodb$Income_fam_new=="30,000 - 49,999"] <- "30,000 - 49,999"
demodb$DOB <- as.Date(as.numeric(demodb$DOB),origin = "1900-01-01")
demodb$Age <- round(as.numeric(ymd(demodb$Date)-ymd(demodb$DOB))/365)
demodb$Age[demodb$Age<10] <- NA
demodb$Height[demodb$Height>1000|demodb$Height<140] <- NA
demodb$BMI = round(demodb$Weight/(demodb$Height/100)^2,1)
demodb$Sex <- ifelse(substr(demodb$Sub_ID,4,4)==1,"F",
                     "M")
for (i in c(16,18:20)) {
  demodb[,i] <- ifelse(demodb[,i]=="有", "Yes", "No")
}
demodb$Smoke_indhis <- ifelse(demodb$Smoke_ind=="否","No","Yes")
demodb$PE_Group <- ifelse(demodb$PE>2,"≥3",paste(demodb$PE))
demodb$EDOC <- clinicdb$X3..預產期.[match(demodb$Sub_ID,clinicdb$Study.No....4)]
demodb$TTB <- as.numeric(ymd(demodb$EDOC)-ymd(demodb$Date))
demodb$FullTerm <- ifelse(demodb$FullTerm=="是", "Yes",
                          ifelse(demodb$FullTerm=="否","No","Unknown"))
demodb$Del_mode <- ifelse(demodb$Del_mode=="順產", "Vag",
                          ifelse(demodb$Del_mode=="剖腹","CS","Unknown"))
demodb$Weight_T1BL <- demodb$Weight_T1-demodb$Weight
demodb$Weight_T3BL <- demodb$Weight_T3-demodb$Weight


tabcat <- colnames(demodb)[c(11,12,16,18,19,22:26)]
tabcon <- colnames(demodb)[c(7,8,28)]
demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = demodb[demodb$Group!="Normal",], strata = "Group")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 




demodb$Group <- NA
demodb$Group[demodb$Sub_ID %in% paste0("PWH100",mhighlist)] <- "High"
demodb$Group[demodb$Sub_ID %in% paste0("PWH100",mlowlist)] <- "Low"
demodb$Group[is.na(demodb$Group)] <- "Normal"
# demodb$Group <- NA
# demodb$Group[demodb$Sub_ID %in% bothhigh$Sub_ID] <- "High"
# demodb$Group[demodb$Sub_ID %in% bothlow$Sub_ID] <- "Low"
# demodb$Group[is.na(demodb$Group)] <- "Normal"

setEPS(width = 6, height = 5)
postscript("BMI_in_FAgroups_both_03042021.eps",onefile = FALSE)
ggstripchart(demodb[demodb$Group!="Normal",],
             x = "Group", y = "TTB", size = 1.5,
             color = "Group", add = c("boxplot"),
             add.params = list(size = 1,
                               linetype = "solid", shape="0.5")) +
  geom_hline(aes(yintercept = 21.8), 
             color="black", size = 0.5, linetype = "dashed") +
  scale_color_manual(values=c("#993404","#7fcdbb")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  stat_compare_means() +
  scale_x_discrete(labels=c("High food additive group \n(N = 47)",
                            "Low food additive group \n(N = 36)")) #H:111 Low:105
dev.off()


daddemodb <- clinicdb[,c(306,307,54,310:314,319,320,321,325,333,334,353,354,426,427)]
colnames(daddemodb) <- c("Sub_ID","Date","Region","Height","Weight","Weight_unit","DOB","Race","Edu",
                         "Work_stat","Work_pos","Income_ind","Smoke_ind",
                         "Alco","Drug","PE","FullTerm","Del_mode")

daddemodb$Region <- ifelse(daddemodb$Region=="香港","HK",
                        ifelse(daddemodb$Region=="廣州","GZ","Other"))

daddemodb$Weight[grepl("kg",tolower(daddemodb$Weight))] <- substr(daddemodb$Weight[grepl("kg",tolower(daddemodb$Weight))],1,
                                                                  nchar(daddemodb$Weight[grepl("kg",tolower(daddemodb$Weight))])-2)
daddemodb$Weight <- as.numeric(daddemodb$Weight)
daddemodb$Weight[complete.cases(daddemodb$Weight)&daddemodb$Weight_unit=="市斤"] <- round(daddemodb$Weight[complete.cases(daddemodb$Weight)&daddemodb$Weight_unit=="市斤"]/2,1)
daddemodb$Weight[complete.cases(daddemodb$Weight)&daddemodb$Weight_unit=="磅 (lb)"] <- round(daddemodb$Weight[complete.cases(daddemodb$Weight)&daddemodb$Weight_unit=="磅 (lb)"]*0.453,1)
daddemodb$Height[grepl("cm",tolower(daddemodb$Height))] <- substr(daddemodb$Height[grepl("cm",tolower(daddemodb$Height))],1,
                                                                  nchar(daddemodb$Height[grepl("cm",tolower(daddemodb$Height))])-2)
daddemodb$Height <- as.numeric(daddemodb$Height)
daddemodb$Height[daddemodb$Height<10] <-  daddemodb$Height[daddemodb$Height<10]*100

daddemodb$Age[daddemodb$Age==0] <- NA # Remove uncertain age
daddemodb$Weight[daddemodb$Weight>100] <- NA # Remove uncertain weight
daddemodb$BMI[daddemodb$BMI>50] <- NA
daddemodb$Edu <- ifelse(grepl("大學",daddemodb$Edu),"Bachelor/degree",
                        ifelse(grepl("碩士",daddemodb$Edu),"Master_above","Highschool_bellow"))
daddemodb$Work_stat <- ifelse(grepl("全職家庭",daddemodb$Work_stat),"Housewife",
                              ifelse(grepl("全職",daddemodb$Work_stat),"Fulltime",
                                     ifelse(grepl("兼職",daddemodb$Work_stat),"Parttime","Others")))

for (i in c(13:15)) {
  daddemodb[,i] <- ifelse(daddemodb[,i]=="有", "Yes", "No")
}
daddemodb$Smoke_indhis <- ifelse(daddemodb$Smoke_ind=="否","No","Yes")
daddemodb$PE_Group <- ifelse(daddemodb$PE>2,"≥3",paste(daddemodb$PE))
daddemodb$EDOC <- clinicdb$X3..預產期.[match(daddemodb$Sub_ID,clinicdb$Study.No....4)]
daddemodb$TTB <- as.numeric(ymd(daddemodb$EDOC)-ymd(daddemodb$Date))
daddemodb$FullTerm <- ifelse(daddemodb$FullTerm=="是", "Yes",
                             ifelse(daddemodb$FullTerm=="否","No","Unknown"))
daddemodb$Del_mode <- ifelse(daddemodb$Del_mode=="順產", "Vag",
                             ifelse(daddemodb$Del_mode=="剖腹","CS","Unknown"))

daddemodb <- cbind(daddemodb[,c(1:5,7)],
                   BMI = round(daddemodb$Weight/(as.numeric(daddemodb$Height)/100)^2,1),
                   Age = round(as.numeric(ymd(daddemodb$Date)-ymd(as.Date(as.numeric(daddemodb$DOB),origin = "1900-01-01")))/365),
                   Group = NA,
                   daddemodb[,c(8:18)])



# Clean up database (Clinical)----
phdb <- clinicdb[,c(52,127,128,142:152,185:196,129,153:157,130,158:163,131,164:166,132,167:171,
                    133,172:175,134,176:178,135,179,180,136:138,181,139,182,140,183,141,273:275)]
colnames(phdb) <- c("Sub_ID","Any_dis","Digest","GI","IBD","IBS","LI","Polyps","Piles","Gastritis_GU","Digest_Other",
                    "IBD_type","Gastritis_GU_type","Digest_Other_type","Gast_OP1","Gast_OP_Rea1","Gast_OP_Year1",
                    "Gast_OP2","Gast_OP_Rea2","Gast_OP_Year2","Gast_OP3","Gast_OP_Rea3","Gast_OP_Year3",
                    "Gast_OP4","Gast_OP_Rea4","Gast_OP_Year4","Immune","RA","SLE","Psoriasis","Immune_Other",
                    "Immune_Other_type","Allergy","Asthma","ADE","FoodA","AR","Allergy_Other","Allergy_Other_type",
                    "HD","CAD","HD_Other","HD_Other_type","Liver","HepB","Other_Hep","FL","Liver_Other",
                    "Liver_Other_type","Renal","Glom","DN","Renal_Other","Renal_Other_type","Repro","POS","Repro_Other",
                    "Repro_Other_type","HG","HG_type","HG_Age","HT","HL","Thyroid","Thyroid_type","Cancer","Cancer_type",
                    "Psycho","Psycho_type","Other_Dis","New_dis","New_dis_type","New_dis_other")


dadphdb <- clinicdb[,c(306,355,356,370:377,379,380,378,412:423,357,381:385,358,386:391,359,392:394,
                       360,395:399,361,400:403,362,404,405,363,406,407,364:366,408,367,409,368,
                       410,411)]
colnames(dadphdb) <- c("Sub_ID","Any_dis","Digest","GI","IBD","IBS","LI","Polyps","Piles","Gastritis_GU","Digest_Other",
                    "IBD_type","Gastritis_GU_type","Digest_Other_type","Gast_OP1","Gast_OP_Rea1","Gast_OP_Year1",
                    "Gast_OP2","Gast_OP_Rea2","Gast_OP_Year2","Gast_OP3","Gast_OP_Rea3","Gast_OP_Year3",
                    "Gast_OP4","Gast_OP_Rea4","Gast_OP_Year4","Immune","RA","SLE","Psoriasis","Immune_Other",
                    "Immune_Other_type","Allergy","Asthma","ADE","FoodA","AR","Allergy_Other","Allergy_Other_type",
                    "HD","CAD","HD_Other","HD_Other_type","Liver","HepB","Other_Hep","FL","Liver_Other",
                    "Liver_Other_type","Renal","Glom","DN","Renal_Other","Renal_Other_type","Repro", "Repro_Other",
                    "Repro_Other_type","HG","HG_type","HG_Age","HT","HL","Thyroid","Thyroid_type","Cancer","Cancer_type",
                    "Psycho","Psycho_type","Other_Dis")



# Check answers
for(i in c(3:11,27:31,33:38,40:42,44:48,50:53,55,58,61:63,65,67)){
  print(colnames(dadphdb)[i])
  print(table(dadphdb[,i]))
}

for(i in c(3:11,27:31,33:38,40:42,44:48,50:53,55,58,61:63,65,67)){
  dadphdb[,i] <- ifelse(dadphdb[,i]=="Checked",1,0)
}

# Main type columns: 3,27,33,40,44,50,55,59,62:64,66,68,70

dadphdb$Any_dis_sum <- as.numeric(apply(dadphdb[,c(3:11,27:31,33:38,40:42,44:48,50:53,55,58,61:63,65,67)],1,
                                     function (x) sum(as.numeric(x))))
dadphdb$Any_dis_sum[is.na(dadphdb$Any_dis)] <- NA


for (i in c(18:26)){
  print(colnames(dadphdb)[i])
  if (all(is.na(dadphdb[,i]))){
    print("All NA")
  }else{
    print("With record")
  }
} # Check all NA columns for digestive sys operation

dadphdb <- dadphdb[,-c(24:26)]
dadphdb[is.na(dadphdb$Any_dis),c(3:64)] <- NA

# Check description in other type and recode
for (i in colnames(dadphdb)[grepl("_type",colnames(dadphdb))]){
  print(i)
  print(table(dadphdb[,i]))
}


phdb$IBD_type[grepl("Croh",phdb$IBD_type)] <- "Crohn's Disease"
phdb$IBD_type[grepl("克隆氏",phdb$Digest_Other_type)] <- "Crohn's Disease"
phdb$Digest_Other_type[grepl("克隆氏",phdb$Digest_Other_type)] <- NA
phdb$Digest_Other_type[grepl("幽門螺旋菌",phdb$Digest_Other_type)] <- "HP"
phdb$Digest_Other_type[grepl("膽石",phdb$Digest_Other_type)] <- "Gallstones"
phdb$Immune_Other_type[complete.cases(phdb$Immune_Other_type)] <- "White matter lesions"
phdb$Allergy_Other_type[complete.cases(phdb$Allergy_Other_type)] <- "Rubella"
phdb$Repro_Other_type[grepl("線|腺",phdb$Repro_Other_type)] <- "Adenomyosis"
phdb$Repro_Other_type[grepl("子宮肌瘤",phdb$Repro_Other_type)] <- "Uterine fibroids"
phdb$HG_type[grepl("一型",phdb$HG_type)] <- "Type1"
phdb$HG_type[grepl("二型",phdb$HG_type)] <- "Type2"
phdb$HG_type[grepl("GDM",phdb$HG_type)] <- "GDM"
phdb$Thyroid_type[grepl("切除",phdb$Thyroid_type)] <- "Thyroidectomy"
phdb$Thyroid_type[grepl("切除",phdb$Thyroid_type)] <- "Thyroidectomy"
phdb$Thyroid_type[grepl("低|減|ypot",phdb$Thyroid_type)] <- "Hypothyroidism"
phdb$Thyroid_type[grepl("高|亢",phdb$Thyroid_type)] <- "Hyperthyroidism"
phdb$Thyroid_type[grepl("thyrotoxcosis",phdb$Thyroid_type)] <- "Thyrotoxcosis"
phdb$Thyroid_type[grepl("no med|未",phdb$Thyroid_type)] <- "No medication_Not firm"
phdb$Psycho_type[grepl("depression|抑鬱",phdb$Psycho_type)] <- "Depression"
phdb$Psycho_type[grepl("焦慮",phdb$Psycho_type)] <- "Anxiety"
phdb$New_dis_type[grepl("其他",phdb$New_dis_type)] <- "Others"
phdb$New_dis_type[grepl("糖尿病",phdb$New_dis_type)] <- "GDM"
phdb$New_dis_type[grepl("高血壓",phdb$New_dis_type)] <- "GH"
phdb$New_dis <- ifelse(phdb$New_dis=="有","Yes",
                       ifelse(phdb$New_dis=="沒有","No",NA))

phdb$Any_dis <- ifelse(phdb$Any_dis=="是","Yes",
                       ifelse(phdb$Any_dis=="否","No",NA))
phdb$GDM_DM <- as.character(apply(phdb[,c(54,66)],1,function(x) if (any(grepl("GDM",x))|complete.cases(x[1])) {return(1)} else {return(0)}))


# Check answers
for(i in c(3:11,27:31,33:38,40:42,44:48,50:53,55:57,59,62:64,66,68,70)){
  print(colnames(dadphdb)[i])
  print(table(dadphdb[,i]))
}

for(i in c(3:11,27:31,33:38,40:42,44:48,50:53,55:57,59,62:64,66,68,70)){
  phdb[,i] <- ifelse(phdb[,i]=="Checked",1,0)
}

# Main type columns: 3,27,33,40,44,50,55,59,62:64,66,68,70

phdb$Any_dis_sum <- as.numeric(apply(phdb[,c(3,27,33,40,44,50,55,59,62:64,66,68,70)],1,
                                     function (x) sum(as.numeric(x))))
phdb$Any_dis_sum[is.na(phdb$Any_dis)] <- NA


for (i in c(18:26)){
  print(colnames(phdb)[i])
  if (all(is.na(phdb[,i]))){
    print("All NA")
  }else{
    print("With record")
  }
} # Check all NA columns for digestive sys operation

phdb <- phdb[,-c(21:26)]
phdb[is.na(phdb$Any_dis),c(3:64)] <- NA

# Check description in other type and recode
for (i in colnames(phdb)[grepl("_type",colnames(phdb))]){
  print(i)
  print(table(phdb[,i]))
}


dadphdb$Gastritis_GU_type[grepl("Gastric",dadphdb$IBD_type)] <- "Gastric ulcer"
dadphdb$Digest_Other_type[grepl("幽門螺旋",dadphdb$Digest_Other_type)] <- "HP"
dadphdb$Digest_Other_type[grepl("胃抽筋",dadphdb$Digest_Other_type)] <- "Stomach cramps"
dadphdb$Digest_Other_type[grepl("腸化生",dadphdb$Digest_Other_type)] <- "Intestinal metaplasia"
dadphdb$Immune_Other_type[complete.cases(dadphdb$Immune_Other_type)] <- "Rubella"
dadphdb$Allergy_Other_type[grepl("allergy|過敏",dadphdb$Allergy_Other_type)] <- "Drug/alch allergy"
dadphdb$RA[grepl("鼻敏感",dadphdb$Digest_Other_type)] <- 1
dadphdb$HD_Other_type[grepl("心漏",dadphdb$HD_Other_type)] <- "Congenital mitral valve anomalies"
dadphdb$HD_Other_type[grepl("不清楚|輕微",dadphdb$HD_Other_type)] <- "Not calrified"
dadphdb$Liver_Other_type[grepl("輕微",dadphdb$Liver_Other_type)] <- "Not calrified"
dadphdb$Repro_Other_type[grepl("疣",dadphdb$Repro_Other_type)] <- "Wart (Low risk)"
dadphdb$HG_type[grepl("二型",dadphdb$HG_type)] <- "Type2"
dadphdb$HG_type[grepl("GDM",dadphdb$HG_type)] <- "GDM"
dadphdb$Thyroid_type[grepl("切除",dadphdb$Thyroid_type)] <- "Thyroidectomy"
dadphdb$Thyroid_type[grepl("抗",dadphdb$Thyroid_type)] <- "Hyperthyroidism"
dadphdb$Psycho_type[grepl("|抑鬱",dadphdb$Psycho_type)] <- "Depression"
dadphdb$Cancer_type[grepl("腸",dadphdb$Cancer_type)] <- "Bowel cancer"


dadphdb$Any_dis <- ifelse(dadphdb$Any_dis=="是","Yes",
                       ifelse(dadphdb$Any_dis=="否","No",NA))

# Clean up database (Current dietary habit)----
CDHdb <- clinicdb[,c(52,203:238)]
colnames(CDHdb) <- c("Sub_ID","FoodA","FoddA_type","MoverV","VoverM","VMeven","Vegan","Vegeta","GlutenFree","LactFree",
                     "Lowcarb","LowNa","Lowcal","HighPro","Other","Other_type","Avg_meal","Avg_takeout","Any_Syn",
                     "Bif","Lact","Cheese","Yogurt","Yakult","Yogurt_drink","Miso","Natto","Kimchi","Suancai","Othey_Syn","Enzyme","DieFrib","Aojiru",
                     "Other_Pro","Other_Syn_type","Other_Pro_type","Stool_type")
dadCDHdb <- clinicdb[,c(306,428:443)]

CDHdb$Avg_meal <- ifelse(is.na(CDHdb$Avg_meal),NA,
                         ifelse(CDHdb$Avg_meal=="1-2 餐","1-2",
                                ifelse(CDHdb$Avg_meal=="3餐","3",">3")))
CDHdb$Avg_takeout <- ifelse(is.na(CDHdb$Avg_takeout),NA,
                            ifelse(CDHdb$Avg_takeout=="1-2 餐","1-2",
                                   ifelse(CDHdb$Avg_takeout=="3餐","3",
                                          ifelse(CDHdb$Avg_takeout=="4餐或以上",">3","None"))))
CDHdb$Any_Syn <- ifelse(CDHdb$Any_Syn=="有","Yes","No")

for(i in c(20:34)){
  CDHdb[,i] <- ifelse(is.na(CDHdb[,i]),NA,
                      ifelse(CDHdb[,i]=="從不","0",
                             ifelse(CDHdb[,i]=="每星期幾次","4",
                                    ifelse(CDHdb[,i]=="每月幾次","3",
                                           ifelse(CDHdb[,i]=="每月少過或等於一次","2",
                                                  ifelse(CDHdb[,i]=="差不多每天","5","1"))))))
}

CDHdb$Dairy <- as.character(apply(CDHdb[,c(22:25)],1,function(x) if(all(is.na(x))) {return(NA)} else {x[which.max(as.numeric(x))]}))
CDHdb$Presfood <-  as.character(apply(CDHdb[,c(26:29)],1,function(x) if(all(is.na(x))) {return(NA)} else {x[which.max(as.numeric(x))]}))


for(i in c(20:34,38,39)){
  CDHdb[,i] <- ifelse(is.na(CDHdb[,i]),NA,
                      ifelse(CDHdb[,i]=="0","Never",
                             ifelse(CDHdb[,i]=="4","Weekly",
                                    ifelse(CDHdb[,i]=="3","Monthly",
                                           ifelse(CDHdb[,i]=="2","<Monthly",
                                                  ifelse(CDHdb[,i]=="5","Everyday","Irregular"))))))
}


CDHdb$Dairy <- as.character(apply(CDHdb[,c(22:25)],1,function(x) if(grepl("Everyday",x)) {return("Everyday")} else {
  if(grepl("Weekly",x)) {return("Weekly")} else {
    if(any(x[complete.cases(x)]=="Monthly")) {return("Monthly")} else {
      if(grepl("<",x)) {return("<Monthly")} else {
        if(grepl("Irregular",x)) {return("Irregular")} else {
          if(grepl("Never",x)) {return("Never")} else {return(NA)}
        }
      }
    }
  }
}))

CDHdb$Stool_type_Group <- ifelse(grepl("一|二",CDHdb$Stool_type),"Constipation",
                                 ifelse(grepl("三|四",CDHdb$Stool_type),"Normal","Lacking fibre/Diarrhea"))

for(i in c(4:6)){
  CDHdb[,i] <- ifelse(CDHdb[,i]=="Checked","Yes","No")
}

CDHdb$Veg_group <- as.character(apply(CDHdb[,c(7,8)],1,function(x) if (any(grepl("Checked",x))) {return("Yes")} else {return("No")}))
CDHdb$SpeDie_group <- as.character(apply(CDHdb[,c(9:14)],1,function(x) if (any(grepl("Checked",x))) {return("Yes")} else {return("No")}))
CDHdb$Dairy_Group <- ifelse(is.na(CDHdb$Dairy),NA,
                            ifelse(grepl("<|Irre|Never",CDHdb$Dairy),"Rarely","Frequently"))
CDHdb$Presfood_Group <- ifelse(is.na(CDHdb$Presfood),NA,
                            ifelse(grepl("<|Irre|Never",CDHdb$Presfood),"Rarely","Frequently"))


dadCDHdb <- clinicdb[,c(306,428:443)]
colnames(dadCDHdb) <- c("Sub_ID","FoodA","FoddA_type","MoverV","VoverM","VMeven","Vegan","Vegeta","GlutenFree","LactFree",
                     "Lowcarb","LowNa","Lowcal","HighPro","Other","Other_type","Stool_type")

for(i in c(4:6)){
  dadCDHdb[,i] <- ifelse(dadCDHdb[,i]=="Checked","Yes","No")
}

dadCDHdb$Veg_group <- as.character(apply(dadCDHdb[,c(7,8)],1,function(x) if (any(grepl("Checked",x))) {return("Yes")} else {return("No")}))
dadCDHdb$SpeDie_group <- as.character(apply(dadCDHdb[,c(9:14)],1,function(x) if (any(grepl("Checked",x))) {return("Yes")} else {return("No")}))
dadCDHdb$Stool_type_Group <- ifelse(grepl("一|二",dadCDHdb$Stool_type),"Constipation",
                                    ifelse(grepl("三|四",dadCDHdb$Stool_type),"Normal","Lacking fibre/Diarrhea"))
dadCDHdb$Any_Syn <- ifelse(dadCDHdb$Any_Syn=="有","Yes","No")



# Merge database----
combdb <- left_join(demodb[,c(1,2,6,27,4,5,7,8,11,12,24,25,18,19,26,28,31,32)],
                    phdb[,c(1,2,68,3,21,27,34,38,44,49,53,54,69,56:58,60,62,65,67)],
                    by = "Sub_ID")
combdb <- left_join(combdb,
                    CDHdb[,c(1,2,4:6,41,42,17:19,38,39,40)],
                    by = "Sub_ID")
combdb <- left_join(combdb,
                    momdb[,c(1,23:26)],
                    by = "Sub_ID")

# check factor levels
for(i in 2:ncol(combdb)){
  if (is.numeric(combdb[,i])) {
    print(colnames(combdb)[i])
  }else{
    print(colnames(combdb)[i])
    print(levels(as.factor(combdb[,i])))
  }
}

combdb$Edu <- factor(as.factor(combdb$Edu),
                           levels(as.factor(combdb$Edu))[c(2,1,3)])
combdb$Work_stat <- factor(as.factor(combdb$Work_stat),
                     levels(as.factor(combdb$Work_stat))[c(1,4,2,3)])
combdb$Income_fam_new <- factor(as.factor(combdb$Income_fam_new),
                     levels(as.factor(combdb$Income_fam_new))[c(1,3,4,2)])
combdb$PE_Group <- factor(as.factor(combdb$PE_Group),
                     levels(as.factor(combdb$PE_Group))[c(2,3,1)])
combdb$Avg_meal <- factor(as.factor(combdb$Avg_meal),
                     levels(as.factor(combdb$Avg_meal))[c(2,3,1)])
combdb$Avg_takeout <- factor(as.factor(combdb$Avg_takeout),
                     levels(as.factor(combdb$Avg_takeout))[c(4,2,3,1)])
combdb$Dairy <- factor(as.factor(combdb$Dairy),
                     levels(as.factor(combdb$Dairy))[c(5,3,1,4,6,2)])
combdb$Presfood <- factor(as.factor(combdb$Presfood),
                     levels(as.factor(combdb$Presfood))[c(4,2,1,3,5)])
combdb$Stool_type_Group <- factor(as.factor(combdb$Stool_type_Group),
                     levels(as.factor(combdb$Stool_type_Group))[c(1,3,2)])

combdb$Group <- ifelse(substr(combdb$Sub_ID,7,9) %in% mhighlist,"High",
                       ifelse(substr(combdb$Sub_ID,7,9) %in% mlowlist,"Low","Normal"))

tabcat <- colnames(combdb)[c(9:10,13:16,19:49)]
tabcon <- colnames(combdb)[c(7,8,17,18,50:53)]

for (i in c(9:10,13:16,19:49)) {
  combdb[,i] <- as.character(combdb[,i])
}

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = combdb[combdb$Group!="Normal",], strata = "Group")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 

combdb$BMI <- as.numeric(combdb$BMI)
combdb$BMI_Group <- "NA"
combdb$BMI_Group[complete.cases(combdb$BMI)&combdb$BMI<18.5] <- "Underweight"
combdb$BMI_Group[complete.cases(combdb$BMI)&combdb$BMI>18.4&combdb$BMI<23] <- "Normal"
combdb$BMI_Group[complete.cases(combdb$BMI)&combdb$BMI>22.9&combdb$BMI<27.1] <- "Overweight"
combdb$BMI_Group[complete.cases(combdb$BMI)&combdb$BMI>27] <- "Obese"


 # Dad
dadcombdb <- left_join(daddemodb,dadphdb[,-33],"Sub_ID")
dadcombdb <- left_join(dadcombdb,dadCDHdb,"Sub_ID")

tabcat <- colnames(dadcombdb)[c(3,10:12,14:21,23:25,37:40)]
tabcon <- colnames(dadcombdb)[c(7,8)]
dadcombdb$PE_group <- ifelse(dadcombdb$PE>2,"≥3",paste(dadcombdb$PE))

dadcombdb$Group[dadcombdb$Sub_ID %in% paste0("PWH300",dhighlist)] <- "High"
dadcombdb$Group[dadcombdb$Sub_ID %in% paste0("PWH300",dlowlist)] <- "Low"
dadcombdb$Group[is.na(dadcombdb$Group)] <- "Normal"

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = dadcombdb[dadcombdb$Group!="Normal",], strata = "Group")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 


dadcombdb$Avg_meal <- combdb$Avg_meal[match(substr(dadcombdb$Sub_ID,7,9),substr(combdb$Sub_ID,7,9))]
dadcombdb$Avg_takeout <- combdb$Avg_takeout[match(substr(dadcombdb$Sub_ID,7,9),substr(combdb$Sub_ID,7,9))]

dadcombdb <- dadcombdb %>% select(colnames(dadcombdb)[colnames(dadcombdb) %in% colnames(combdb)])
# Plot overall demo


plot_list = list()
for (i in c(9:11,12,14,15,42,43,47)) {
  temp.db <- data.frame(Var = combdb[,i])
  p <- ggplot(temp.db,aes(x = Var)) + 
    geom_histogram( stat = "count", binwidth = 1, bins = 5,
                    color = "steelblue1", fill = "steelblue1") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 14),
          axis.title = element_text(size = 18)) +
    ylab("Count") +
    xlab(paste(colnames(combdb)[i]))
  plot_list <- c(plot_list, list(p))
}

wrap_plots(plot_list,ncol = 3)

demofig <- combdb[,c(17:34)]
colnames(demofig)[c(6,10,11,13,14)] <- c("Heart disease","Hyperglycaemia","Hyperglycaemia_type","Hypertension","Hyperlipidemia")
demofig[,c(3:10,12:17)] <- apply(demofig[,c(3:10,12:17)],2,function(x) ifelse(x==1,"Yes","No"))

plot_list = list()
for (i in c(1:18)) {
  temp.db <- data.frame(Var = demofig[,i])
  p <- ggplot(temp.db,aes(x = Var)) + 
    geom_histogram( stat = "count", binwidth = 1, bins = 5,
                    color = "tomato1", fill = "tomato1") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 14),
          axis.title = element_text(size = 18)) +
    ylab("Count") +
    xlab(paste(colnames(demofig)[i]))
  plot_list <- c(plot_list, list(p))
}

wrap_plots(plot_list,ncol = 4)


# Add up latest database with quantitative alcohol drinking----
dadalco <-  data.frame(read_excel("Data/Baseline DAT 1 + Alcohol Unit (send to Wing on Jul 7).xlsx", 
                                  sheet = "Father"))
dadalco <- dadalco %>% select(Study.No.,Alcohol.Unit.Per.Week)
momalco <- data.frame(read_excel("Data/Baseline DAT 1 + Alcohol Unit (send to Wing on Jul 7).xlsx", 
                                 sheet = "Mother"))
momalco <- momalco %>% select(Record.ID,Alcohol.Unit.Per.Week)
HLMdb$Alco_con <- momalco$Alcohol.Unit.Per.Week[match(HLMdb$Sub_ID,momalco$Record.ID)]
HLMdb$Alco_con[is.na(HLMdb$Alco_con)] <-dadalco$Alcohol.Unit.Per.Week[match(HLMdb$Sub_ID[is.na(HLMdb$Alco_con)],
                                                                            dadalco$Study.No.)]


# Updated database----
PDHdb <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                               sheet = "Combined Baseline"))[,c(1,3:41)]
colnames(PDHdb) <- c("Sub_ID","BF","Infancy_food_source","HG_4-12m","HG_1-5y","HG_5-10y",
  "HG_10-18y","Adolo_food_source","Dairy_4-12m","Dairy_1-5y","Dairy_5-10y",
  "Dairy_10-18y", "Seafood_4-12m","Seafood_1-5y","Seafood_5-10y", "Seafood_10-18y","CarboHfood_4-12m",
  "CarboHfood_1-5y","CarboHfood_5-10y", "CarboHfood_10-18y", "ProcFruit_4-12m", "ProcFruit_1-5y",
  "ProcFruit_5-10y", "ProcFruit_10-18y", "ProcVeg_4-12m","ProcVeg_1-5y","ProcVeg_5-10y", "ProcVeg_10-18y",
  "Fastfood_4-12m","Fastfood_1-5y","Fastfood_5-10y", "Fastfood_10-18y", "SoftDrink_4-12m","SoftDrink_1-5y",
  "SoftDrink_5-10y", "SoftDrink_10-18y", "Snacks_4-12m","Snacks_1-5y","Snacks_5-10y", "Snacks_10-18y")
PDHdb$Sex <- ifelse(substr(PDHdb$Sub_ID,4,4)==1,"F","M")
PDHdb <- PDHdb %>% filter(complete.cases(BF))

MomFA <- data.frame(read_excel("Data/MOMmy Diet Related REDCap Data 2021.08.05.xlsx", 
                                sheet = "Combined Baseline"))[,c(1,43:52,66)]
MomFA$Comp.case <- as.numeric(apply(MomFA[,c(2:11)],1,function(x) ifelse(all(is.na(x)),1,0)))
MomFA <- MomFA %>% 
  filter(Comp.case==0) %>% 
  select(-Comp.case,-MDX)
colnames(MomFA)[1] <- "Sub_ID"
MomFA$All_FA <- as.numeric(apply(MomFA[,c(2:10)],1,function(x) sum(as.numeric(x),
                                                                   na.rm = T)))
colnames(MomFA)[11] <- "Weight"
MomFA$Weight <- round(as.numeric(MomFA$Weight),1)

for (i in 2:10){
  j <- colnames(MomFA)[i]
  MomFA[,i+11] <- round(MomFA[,i]/MomFA$Weight,1)
  colnames(MomFA)[i+11] <- paste0(colnames(MomFA)[i],"_W")
}

MomFA$All_FA_W <- round(MomFA$All_FA/MomFA$Weight,1)
MomFA$Sex <- ifelse(substr(MomFA$Sub_ID,4,4)==1,"F","M")

a <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                sheet = "MOMmy-FaBaseline_DATA_LABELS_20"))
demodb <- a[,c(1,4,5,10:14,28:30,35,55,57,76,77,150,151)]
daddemodb <- a[,c(1,5,6,8:12,17:19,21,31,32,51,52,124,125)]
colnames(daddemodb) <- c("Sub_ID","Date","Region","Height","Weight","Weight_unit","DOB","Race","Edu",
                         "Work_stat","Work_pos","Income_fam","Smoke_ind",
                         "Alco","Drug","PE","FullTerm","Del_mode")
daddemodb$Income_fam <- demodb$Income_fam[match(substr(daddemodb$Sub_ID,7,9),
                                                substr(demodb$Sub_ID,7,9))]
daddemodb$Region <- demodb$Region[match(substr(daddemodb$Sub_ID,7,9),
                                                substr(demodb$Sub_ID,7,9))]
a <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                           sheet = "Combined Baseline"))
# demodb <- rbind(demodb,daddemodb)
demodb$Alco <- a$Alcohol.Unit.Per.Week[match(demodb$Sub_ID,a$Study.No.)]
colnames(a)[1] <- "Sub_ID"
a$SpeDie_group <-  as.character(apply(a[,c(56:63)],1,
                                      function(x) if (any(x==1))
                                        {return("Yes")} else {return("No")}))

a$Diet_type <- as.character(apply(a[,c(53:55)],1,function(x) paste(colnames(a)[c(53:55)][x==1],collapse = ",")))
a$VM_group <- as.character(lapply(as.list(a$Diet_type),function(x) strsplit(x,",")[[1]][1]))
a$VM_group <- ifelse(is.na(a$VM_group),NA,
                     ifelse(a$VM_group=="肉類.蔬菜水果","MoverV",
                            ifelse(a$VM_group=="蔬菜水果.肉類","VoverM","VMeven")))
demodb <- left_join(demodb,a[,c(1,88,90)], by = "Sub_ID")
demodb$Sex <- ifelse(substr(demodb$Sub_ID,4,4)==1,"F","M")


tabcat <- colnames(demodb)[c(2,8:10,17:22,25)]
tabcon <- colnames(demodb)[c(14,26,27)]

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = demodb, strata = "Sex")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"demo_29072021.xlsx")

for (i in c(2:10)){
  MomFA <- cbind(MomFA,i=log2(MomFA[,i]))
  
}

colnames(MomFA)[c(24:32)] <- paste0("log2",colnames(MomFA)[c(2:10)])
MOMcomb <- left_join(MOMcomb,MomFA,"Sub_ID")
MOMcomb <- MOMcomb[,-c(109:130)]
colnames(MOMcomb)[c(14,60:70)] <- str_replace_all(colnames(MOMcomb)[c(14,60:70)], ".x","")
MOMcomb[,c(109:117,)][MOMcomb[,c(109:117)]=="-Inf"] <- NA

# Build up combined database----
demodb <- demodb %>% select(-Weight,-Weight_unit,-BMI,-PE,-Work_pos,-Smoke_ind,
                            -Drug,-PE,-Income_fam)
PDHdb <- PDHdb %>% select(-Sex)
MOMcomb <- left_join(demodb,PDHdb,by = "Sub_ID")
MOMcomb <- left_join(MOMcomb,MomFA[,c(1,11,13:22)], by = "Sub_ID")
MOMcomb$BMI <- round(MOMcomb$Weight/((MOMcomb$Height/100)^2))
MOMcomb$BMI_Group <- ifelse(is.na(MOMcomb$BMI),NA,
                            ifelse(MOMcomb$BMI<18.5,"aUnderweight",
                                   ifelse(MOMcomb$BMI<23,"bHealthy weight",
                                          ifelse(MOMcomb$BMI<25,"cOverweight","Obesity"))))
MOMcomb$BMI_Congroup <- ifelse(is.na(MOMcomb$BMI),NA,
                               ifelse(MOMcomb$BMI<18.5,1,
                                      ifelse(MOMcomb$BMI<23,2,
                                             ifelse(MOMcomb$BMI<25,3,4))))


for (i in 61:70){
  j <- colnames(MOMcomb)[i]
  MOMcomb[,i+12] <- log2(MOMcomb[,i])
  colnames(MOMcomb)[i+12] <- paste0("log2",colnames(MOMcomb)[i])
}


MOMcomb$Emu_W <- MOMcomb$P80_W+MOMcomb$CMC_W+MOMcomb$CRN_W

MOMcomb$Sweet_W <- MOMcomb$SUC_W+MOMcomb$SAC_W+MOMcomb$ASP_W
MOMcomb$OtherFA_W <- as.numeric(apply(MOMcomb[,c(64:66)],1,
                                      function(x) sum(x,na.rm = T)))
for (i in 83:85){
  j <- colnames(MOMcomb)[i]
  MOMcomb[,i+3] <- log2(MOMcomb[,i])
  colnames(MOMcomb)[i+3] <- paste0("log2",colnames(MOMcomb)[i])
}


MOMcomb$Edu <- factor(as.factor(MOMcomb$Edu),
                      levels(as.factor(MOMcomb$Edu))[c(2,1,3)])
MOMcomb$Work_stat <- factor(as.factor(MOMcomb$Work_stat),
                            levels(as.factor(MOMcomb$Work_stat))[c(1,4,2,3)])
MOMcomb$Income_fam_new <- factor(as.factor(MOMcomb$Income_fam_new),
                                 levels(as.factor(MOMcomb$Income_fam_new))[c(1,3,4,2)])
MOMcomb$PE_Group <- factor(as.factor(MOMcomb$PE_Group),
                           levels(as.factor(MOMcomb$PE_Group))[c(2:4,1)])

MOMcomb[,c(71:88)][MOMcomb[,c(71:88)]=="-Inf"] <- NA
MOMcomb$VM_group <- factor(as.factor(MOMcomb$VM_group),
                           levels(as.factor(MOMcomb$VM_group))[c(3,2,1)])

MOMcomb$EDOC <- a$X3..預產期.[match(MOMcomb$Sub_ID,a$Record.ID)]
MOMcomb$TTB <- as.numeric(ymd(MOMcomb$EDOC)-ymd(MOMcomb$Date))

# Add self report disease

dadSRD <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Father Baseline Q_20210716.xlsx"))[,c(5,55:68,70,71,86)]
colnames(dadSRD) <- c("Sub_ID","Digest","Immune","Allergy","HD","Liver","Renal","Repro","HG","HT","HL","Thyroid",
                      "Cancer","Psycho","Other_Dis","IBD","IBS","Eczema")
momSRD <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Mother Baseline Q_20210716.xlsx"))[,c(4,80:93,95,96,111)]
colnames(momSRD) <- c("Sub_ID","Digest","Immune","Allergy","HD","Liver","Renal","Repro","HG","HT","HL","Thyroid",
                      "Cancer","Psycho","Other_Dis","IBD","IBS","Eczema")
RSD <- rbind(momSRD,dadSRD)


MOMcomb <- left_join(MOMcomb,RSD,by = "Sub_ID")
MOMcomb <- left_join(MOMcomb,RSD[,c(1,16:18)],by = "Sub_ID")

MOMcomb[,c(89:108)][MOMcomb[,c(89:108)]=="Checked"] <- "2Yes"
MOMcomb[,c(89:108)][MOMcomb[,c(89:108)]=="Unchecked"] <- "1No"
MOMcomb$Any_disease <- as.character(apply(MOMcomb[,c(90:108)],1,function(x) ifelse(any(x=="2Yes"),"Yes","No")))

MOMcomb$OVW <- ifelse(MOMcomb$BMI_Group=="cOverweight"|MOMcomb$BMI_Group=="Obesity",
                      "Overweight","Normal")

tabcat <- colnames(MOMcomb)[c(3,6:8,12,13,15,17,18,90:108,72)]
tabcon <- colnames(MOMcomb)[c(16,9,20,71,70)]

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = MOMcomb, strata = "Sex")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"Table1_update1.xlsx")

MOMcomb$OVW <- as.character(MOMcomb$OVW)

tabcat <- colnames(MOMcomb)[c(3,6:8,118,13,15,17,18,90:108,72)]
tabcon <- colnames(MOMcomb)[c(16,9,20,71,132:135)]

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = MOMcomb, strata = "Sex")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"Table_all_update_17092021.xlsx")

# Add individual income
a <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                           sheet = "MOMmy-FaBaseline_DATA_LABELS_20"))
indi_income <- a[,c(4,23)]
colnames(indi_income) <- c("Sub_ID","Indi_income")
a <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                           sheet = "MOMmy-MoBaseline_DATA_LABELS_20"))

a <- a[,c(3,34)]
colnames(a) <- c("Sub_ID","Indi_income")
indi_income <- rbind(a,indi_income)
indi_income$Indi_income <- ifelse(is.na(indi_income$Indi_income)|indi_income$Indi_income=="不方便作答","NA",
                                  ifelse(grepl("2,000|- 9,999|10,000 -",indi_income$Indi_income),"<15,000",
                                         ifelse(grepl("50,000 -|≥",indi_income$Indi_income),"≥50,000",paste0(indi_income$Indi_income))))
MOMcomb$Indi_income <- indi_income$Indi_income[match(MOMcomb$Sub_ID,indi_income$Sub_ID)]
rm(indi_income)
# Add log2 overall FA and 3 subgroups
MomFA$log2Emu <- log2(MomFA$P80+MomFA$CMC+MomFA$CRN)
MomFA$log2Sweet <- log2(MomFA$ASP+MomFA$SAC+MomFA$SUC)
MomFA$log2Other <- log2(MomFA$AlSiO+MomFA$TiO2+MomFA$X.SO3.2)
MomFA$log2_FA <- log2(MomFA$All_FA)

MomFA$Emu <- MomFA$P80+MomFA$CMC+MomFA$CRN
MomFA$Sweet <- MomFA$ASP+MomFA$SAC+MomFA$SUC
MomFA$Other <- MomFA$AlSiO+MomFA$TiO2+MomFA$X.SO3.2



MOMcomb <- left_join(MOMcomb,MomFA[,c(1,33:36)],"Sub_ID")
MOMcomb[,c(119:122)][MOMcomb[,c(119:122)]=="-Inf"] <- NA
MOMcomb[,c(119:122)][MOMcomb[,c(119:122)]=="Inf"] <- NA
MOMcomb <- left_join(MOMcomb,MomFA[,c(1:10,12,37:39)])

MOMcomb <- left_join(MOMcomb,MomFA[,c(1,40:51)],"Sub_ID")



# Add pre-preganancy weight
a <- a[,c(3,152,153)]
colnames(a) <- c("Sub_ID","Pre_weight","Pre_weight_unit")
a$Pre_weight <- str_replace_all(a$Pre_weight,"kg","")
a$Pre_weight <- as.numeric(a$Pre_weight)
a <- a %>% filter(complete.cases(Pre_weight_unit)&complete.cases(Pre_weight))
a$Pre_weight[a$Pre_weight_unit=="市斤"] <- round(a$Pre_weight[a$Pre_weight_unit=="市斤"]/2,1)
a$Pre_weight[a$Pre_weight_unit=="磅 (lb)"] <- round(a$Pre_weight[a$Pre_weight_unit=="磅 (lb)"]*0.453,1)
a$Pre_weight[a$Pre_weight>200] <- NA
a$Pre_weight[a$Pre_weight<40] <- NA
MOMcomb$Pre_Weight <- a$Pre_weight[match(MOMcomb$Sub_ID,a$Sub_ID)]
MOMcomb$Pre_BMI <- round(MOMcomb$Pre_Weight/(MOMcomb$Height/100)^2,1) 
MOMcomb$Pre_OVW <- ifelse(is.na(MOMcomb$Pre_BMI),NA,
                          ifelse(MOMcomb$Pre_BMI<23,"0","1"))


# Update demo for table one-----
a <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Mother Baseline Q_20210716_Wing_20210827_Checked.xlsx"))[,-1]
b <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Father Baseline Q_20210716_Wing_20200819_Checked.xlsx"))[,-1]

a <- a %>% filter(Complete.=="Complete")
a <- a[,c(1,4,5,10:14,28:30,35,55,57,76,77,150,151)]
b <- b %>% filter(Complete.=="Complete")
b <- b[,c(1,5,6,8:12,17:19,21,31,32,51,52,124,125)]
colnames(b) <- c("Sub_ID","Date","Region","Height","Weight","Weight_unit","DOB","Race","Edu",
                         "Work_stat","Work_pos","Income_fam","Smoke_ind",
                         "Alco","Drug","PE","FullTerm","Del_mode")
b$Region <- a$Region[match(substr(b$Sub_ID,7,9),substr(a$Sub_ID,7,9))]
demodb_new <- rbind(a,b)
demodb_new$DOB <- as.Date(as.numeric(demodb_new$DOB),origin="1900-01-01")

a <- data.frame(read_excel("Data/MOMmy_Diet_Related_REDCap_Data_2021.07.27.xlsx", 
                           sheet = "Combined Baseline"))
colnames(a)[1] <- "Sub_ID"
a$SpeDie_group <-  as.character(apply(a[,c(56:63)],1,
                                      function(x) if (any(x==1))
                                      {return("Yes")} else {return("No")}))

a$Diet_type <- as.character(apply(a[,c(53:55)],1,function(x) paste(colnames(a)[c(53:55)][x==1],collapse = ",")))
a$VM_group <- as.character(lapply(as.list(a$Diet_type),function(x) strsplit(x,",")[[1]][1]))

a$VM_group <- as.character(lapply(as.list(a$Diet_type),function(x) strsplit(x,",")[[1]][1]))

a$VM_group <- ifelse(is.na(a$VM_group),NA,
                     ifelse(a$VM_group=="肉類.蔬菜水果","MoverV",
                            ifelse(a$VM_group=="蔬菜水果.肉類","VoverM","VMeven")))
demodb_new <- left_join(demodb_new,a[,c(1,87,90)], by = "Sub_ID")

demodb_new$Income_indi <- indi_income$Indi_income[match(demodb_new$Sub_ID,indi_income$Sub_ID)]

demodb_new$Region <- ifelse(demodb_new$Region=="香港","HK",
                            ifelse(demodb_new$Region=="廣州","GZ","Other"))

demodb_new$Weight <- str_replace_all(demodb_new$Weight,"kg","")
demodb_new$Weight <- as.numeric(demodb_new$Weight)
demodb_new$Weight_unit[demodb_new$Sub_ID=="PWH100004"] <- "公斤 (kg)"
demodb_new$Weight[complete.cases(demodb_new$Weight)&demodb_new$Weight_unit=="市斤"] <- round(demodb_new$Weight[complete.cases(demodb_new$Weight)&demodb_new$Weight_unit=="市斤"]/2,1)
demodb_new$Weight[complete.cases(demodb_new$Weight)&demodb_new$Weight_unit=="磅 (lb)"] <- round(demodb_new$Weight[complete.cases(demodb_new$Weight)&demodb_new$Weight_unit=="磅 (lb)"]*0.453,1)
demodb_new$Height <- str_replace_all(demodb_new$Height,"cm","")
demodb_new$Height <- as.numeric(demodb_new$Height)
demodb_new$Height[demodb_new$Height<10] <-  demodb_new$Height[demodb_new$Height<10]*100
demodb_new$Height[demodb_new$Height==1700] <-  170
demodb_new$Weight[demodb_new$Weight>100] <- NA
demodb_new$Weight[demodb_new$Weight<40] <- NA # Remove uncertain weight
demodb_new <- cbind(demodb_new[,c(1:5,7)],
                    BMI = round(demodb_new$Weight/(as.numeric(demodb_new$Height)/100)^2,1),
                    Age = round(as.numeric(ymd(demodb_new$Date)-ymd(demodb_new$DOB))/365),
                    Group = NA,
                    demodb_new[,c(8:21)])
demodb_new$Age[demodb_new$Age<10] <- NA # Remove uncertain age
demodb_new$BMI[demodb_new$BMI>50] <- NA
demodb_new$Edu <- ifelse(grepl("大學",demodb_new$Edu),"Bachelor/degree",
                         ifelse(grepl("碩士",demodb_new$Edu),"Master_above","Highschool_bellow"))
demodb_new$Work_stat <- ifelse(grepl("全職家庭",demodb_new$Work_stat),"Housewife",
                               ifelse(grepl("全職",demodb_new$Work_stat),"Fulltime",
                                      ifelse(grepl("兼職",demodb_new$Work_stat),"Parttime","Others")))

demodb_new$Smoke_indhis <- ifelse(demodb_new$Smoke_ind=="否","No","Yes")
demodb_new$PE_Group <- ifelse(demodb_new$PE>2,"≥3",paste(demodb_new$PE))
demodb_new$EDOC <- a$X3..預產期.[match(demodb_new$Sub_ID,a$Study.No.)]
demodb_new$TTB <- as.numeric(ymd(demodb_new$EDOC)-ymd(demodb_new$Date))
demodb_new$FullTerm <- ifelse(demodb_new$FullTerm=="是", "Yes",
                              ifelse(demodb_new$FullTerm=="否","No","Unknown"))
demodb_new$Del_mode <- ifelse(demodb_new$Del_mode=="順產", "Vag",
                              ifelse(demodb_new$Del_mode=="剖腹","CS","Unknown"))

demodb_new <- left_join(demodb_new,RSD,"Sub_ID")


demodb_new[,c(28:44)][demodb_new[,c(28:44)]=="Checked"] <- "2Yes"
demodb_new[,c(28:44)][demodb_new[,c(28:44)]=="Unchecked"] <- "1No"
demodb_new$Any_disease <- as.character(apply(demodb_new[,c(28:41)],1,function(x) ifelse(any(x=="2Yes"),"2Yes","1No")))
demodb_new$Mar_status <- a$X5..你現時的婚姻狀況.[match(demodb_new$Sub_ID,a$Study.No.)]
demodb_new$BMI_Group <- ifelse(is.na(demodb_new$BMI),NA,
                            ifelse(demodb_new$BMI<18.5,"aUnderweight",
                                   ifelse(demodb_new$BMI<23,"bHealthy weight",
                                          ifelse(demodb_new$BMI<25,"cOverweight","Obesity"))))

demodb_new$OVW <- ifelse(is.na(demodb_new$BMI),NA,
                         ifelse(demodb_new$BMI<23,"1No","2Yes"))

demodb_new$Sex <- ifelse(substr(demodb_new$Sub_ID,4,4)==1,"F","M")
demodb_new$Smoke_ind <- ifelse(grepl("仍",demodb_new$Smoke_ind),"Actively",
                               ifelse(grepl("戒掉",demodb_new$Smoke_ind),"Has quit","Never"))


tabcat <- colnames(demodb_new)[c(3,10:12,23,16,25,15,22,26,28:45,47,48,46)]
tabcon <- colnames(demodb_new)[c(7,8,21,27)]

demotab <- CreateTableOne(vars = c(tabcat,tabcon), data = wsamplist, strata = "Sex")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$var <- rownames(a)
writexl::write_xlsx(a,"Table1_wsl_19092021.xlsx")

# Extract subjects with stool samples
wsamplist <- data.frame(read_xlsx("Data/Meta_wsl_with sample.xlsx"))
wsamplist <- demodb_new %>% filter(Sub_ID %in% wsamplist$Sub_ID)


# Add weight change for mothers-----
T2 <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Mother T2 T3 Q_20210716.xlsx"))[,-1]
T2 <- T2 %>% filter(grepl("二",X2..孕期.)) %>%
  filter(Complete.=="Complete") %>% 
  select(Record.ID,X1..今天的日期.,
         X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.,X3a..體重單位.)
T3 <- data.frame(read_excel("~/Documents/MOMmy_REDCap/MOMmy-Mother T2 T3 Q_20210716.xlsx"))[,-1]
T3 <- T3 %>% filter(grepl("三",X2..孕期.)) %>% 
  filter(Complete.=="Complete") %>% 
  select(Record.ID,X1..今天的日期.,
         X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.,X3a..體重單位.)
T2$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤. <- as.numeric(T2$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.)
T2$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T2$X3a..體重單位.=="磅 (lb)"] <- round(T2$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T2$X3a..體重單位.=="磅 (lb)"]*0.453,1)
T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤. <- as.numeric(T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.)
T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T3$X3a..體重單位.=="磅 (lb)"] <- round(T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T3$X3a..體重單位.=="磅 (lb)"]*0.453,1)
T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T3$X3a..體重單位.=="市斤"] <- round(T3$X3..體重..公斤.kg....nbsp.磅.lb....nbsp.市斤.[T3$X3a..體重單位.=="市斤"]/2,1)

colnames(T2) <- c("Sub_ID","Date_T2","Weight_T2","Unit")
colnames(T3) <- c("Sub_ID","Date_T3","Weight_T3","Unit")
T2$Weight_T2[T2$Weight_T2>100|T2$Weight_T2<40] <- NA
T3$Weight_T3[T3$Weight_T3>100|T3$Weight_T3<40] <- NA

T2$EDOC <- MOMcomb$EDOC[match(T2$Sub_ID,MOMcomb$Sub_ID)]
T3$EDOC <- MOMcomb$EDOC[match(T3$Sub_ID,MOMcomb$Sub_ID)]
T2$EDOC[T2$Sub_ID=="PWH100945"] <- NA
T2$GA_T2 <- 40-round(as.numeric(ymd(T2$EDOC)-ymd(T2$Date_T2))/7)
T2$GA_T2[T2$GA_T2<1|T2$GA_T2>35] <- NA
T3$GA_T3 <- 40-round(as.numeric(ymd(T3$EDOC)-ymd(T3$Date_T3))/7)
T3$GA_T3[T3$GA_T3<1|T3$GA_T3>42] <- NA
T2$Pre_Weiht <- MOMcomb$Pre_Weight[match(T2$Sub_ID,MOMcomb$Sub_ID)]
T2$Weight_change <- round((T2$Weight_T2-T2$Pre_Weiht)/T2$GA_T2,2)
T2$Weight_change[T2$Weight_change<0] <- NA
T3$Pre_Weiht <- MOMcomb$Pre_Weight[match(T3$Sub_ID,MOMcomb$Sub_ID)]
T3$Weight_change <- round((T3$Weight_T3-T3$Pre_Weiht)/T3$GA_T3,2)
T3$Weight_change[T3$Weight_change<0] <- NA

MOMcomb$T2_WC <- T2$Weight_change[match(MOMcomb$Sub_ID,T2$Sub_ID)]
MOMcomb$T3_WC <- T3$Weight_change[match(MOMcomb$Sub_ID,T3$Sub_ID)]

T3$Height <- MOMcomb$Height[match(T3$Sub_ID,MOMcomb$Sub_ID)]
T3$Pre_BMI <- round(T3$Pre_Weiht/(T3$Height/100)^2,1)
T3$Weight_change_raw <- T3$Weight_T3-T3$Pre_Weiht 
T3$IOM_GWG <- ifelse((T3$Pre_BMI<18.5&T3$Weight_change_raw>18)|
                       (T3$Pre_BMI>18.5&T3$Pre_BMI<23&T3$Weight_change_raw>16)|
                       (T3$Pre_BMI>22.99&T3$Pre_BMI<25&T3$Weight_change_raw>11.5)|
                       (T3$Pre_BMI>24.99&T3$Weight_change_raw>9),"Over",
                     ifelse((T3$Pre_BMI<18.5&T3$Weight_change_raw<12.5)|
                              (T3$Pre_BMI>18.5&T3$Pre_BMI<23&T3$Weight_change_raw<11.5)|
                              (T3$Pre_BMI>22.99&T3$Pre_BMI<25&T3$Weight_change_raw<7)|
                              (T3$Pre_BMI>24.99&T3$Weight_change_raw<5),"Under","Normal"))

T3$BMI_Group <- MOMcomb$BMI_Group[match(T3$Sub_ID,MOMcomb$Sub_ID)]
MOMcomb$IOM_GWG <- T3$IOM_GWG[match(MOMcomb$Sub_ID,T3$Sub_ID)]





# Select urine samples for Emu test----

FAsamp_db <- data.frame(
  read_excel("Data/0015_MOM_LZ_20210421B1_20210714_revised_20210811B3_20211118_WING.xlsx",
             skip = 1))
FAsamp_db_fil <- FAsamp_db %>% filter(complete.cases(STYPE...2)) %>% 
  filter(complete.cases(STYPE...21)) %>% 
  filter(complete.cases(STYPE...40))
FAsamp_db_fil$GA_Urine <- 40-(round(as.numeric(ymd(MOMcomb_sub$EDOC[match(FAsamp_db_fil$STUDY_NO...1,MOMcomb_sub$Sub_ID)])-
  ymd(FAsamp_db_fil$CDATE...8)),1)/7)
FAsamp_db_fil$GA_Urine[FAsamp_db_fil$GA_Urine<0] <- NA

FAsamp_db_fil$OVA_Emu <- MomFA$Emu[match(FAsamp_db_fil$STUDY_NO...1,MomFA$Sub_ID)]
FAsamp_db_fil <- FAsamp_db_fil %>% filter(complete.cases(OVA_Emu))
FAsamp_db_fil <- FAsamp_db_fil[order(FAsamp_db_fil$OVA_Emu,decreasing = TRUE),]
FAsamp_db_fil$Order <- 1:nrow(FAsamp_db_fil)
FAsamp_db_fil <- FAsamp_db_fil %>% filter(Order %in% c(1:50,
                                                      (nrow(FAsamp_db_fil)-49):(nrow(FAsamp_db_fil)),
                                                       (round(nrow(FAsamp_db_fil)/2)-25):
                                                         (round(nrow(FAsamp_db_fil)/2)+24))) %>% 
  select(STUDY_NO...1,OVA_Emu,Order)
FAsamp_db_fil$Group[1:50] <- "High"
FAsamp_db_fil$Group[51:100] <- "Middle"
FAsamp_db_fil$Group[101:150] <- "Low"

writexl::write_xlsx(FAsamp_db_fil,"Data/FA_Samptest_pilot_list.xlsx")
FAsamp_db$Sample.SNO...Selected.by.applicants....15 <- ifelse(FAsamp_db$STUDY_NO...1 %in% FAsamp_db_fil$STUDY_NO...1,
                                                              "Yes",NA)
writexl::write_xlsx(FAsamp_db,
                    "Data/0015_MOM_LZ_20210421B1_20210714_revised_20210811B3_20211118_WING_20211208.xlsx")
FAsamp_db <- data.frame(read_excel("Data/FA_Samptest_pilot_list.xlsx"))
stoollistz <- read_excel("Data/2. sample list summary_20211130.xlsx")
FAsamp_db_fil$Stool <- stoollistz$`MOTHER T3`[match(FAsamp_db_fil$STUDY_NO...1,stoollistz$`MOTHER T3`)]
FAsamp_db_fil$Stool2 <- stoollistz$T2[match(FAsamp_db_fil$STUDY_NO...1,stoollistz$T3)]
FAsamp_db$stool <- stoollistz$`MOTHER T3`[match(FAsamp_db$STUDY_NO,stoollistz$`MOTHER T3`)]
FAsamp_db$stool[is.na(FAsamp_db$stool)] <- stoollistz$`MOTHER T2`[match(FAsamp_db$STUDY_NO[is.na(FAsamp_db$stool)],stoollistz$`MOTHER T2`)]
FAsamp_db$stool[is.na(FAsamp_db$stool)] <- stoollistz$`Mother T1`[match(FAsamp_db$STUDY_NO[is.na(FAsamp_db$stool)],stoollistz$`Mother T1`)]
writexl::write_xlsx(FAsamp_db,"Stool_MTG_21022022.xlsx")
  
  
  
  
# Clinical data 4 Y1 report----
Y1db <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx")[,-c(1,5:7)])

# Delivery mode/ other demo
PdDel <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                               sheet = "REDCap Qs (D01-M12)")) %>% 
  filter(Event.Name == "Date of Born (Arm 2: Infant)")

PdDel <- PdDel[,c(1,6:26)]

colnames(PdDel) <- c("MOM_SNO","Del_mode","CS_mode","Weight_child","Weight_unit_child",
                     "NComp","NTrauma","NHypoG","NInfec","NDeformity",
                     "NMAS","NJaundice","NAsphyxia","NRDS","NComp_Other",
                     "NComp_Other_specific","NComp_Treatment",
                     "Light","Medication","NComp_Treatment_Other","NComp_Treatment_Other_spe",
                     "Del_feeding_mode")
PdDel$Weight_child[complete.cases(PdDel$Weight_child)&PdDel$Weight_unit_child=="磅 (lb)"] <- round(PdDel$Weight_child[complete.cases(PdDel$Weight_child)&PdDel$Weight_unit_child=="磅 (lb)"]*453.6,1)
PdDel <- PdDel %>% select(-Weight_unit_child)
PdDel$Del_mode <- ifelse(PdDel$Del_mode=="陰道分娩", "Vag",
                         ifelse(PdDel$Del_mode=="剖腹產","CS","Vac_or_forcep"))
PdDel$CS_mode <- ifelse(is.na(PdDel$CS_mode),NA,
                        ifelse(PdDel$CS_mode=="選擇性","Elective","Emergency"))

PdDel[,c(5,16)] <- apply(PdDel[,c(5,16)],2,function(x)
  ifelse(is.na(x),NA,ifelse(x=="否",0,1)))

PdDel[,c(6:14,17:19)] <- apply(PdDel[,c(6:14,17:19)],2,function(x)
  ifelse(is.na(x),NA,ifelse(x=="Unchecked",0,1)))
PdDel$Del_feeding_mode <- ifelse(PdDel$Del_feeding_mode=="混合","Mix",
                                 ifelse(grepl("母乳",PdDel$Del_feeding_mode),"BF","Formula"))
Y1db <- left_join(Y1db,PdDel,"MOM_SNO")

# BF hist
PdBF <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                              sheet = "REDCap Qs (D01-M12)")) %>% 
  filter(Event.Name != "Date of Born (Arm 2: Infant)")
PdBF <- PdBF[,-c(2:28,34:89)]
colnames(PdBF) <- c("MOM_SNO","Timepoint","Current_feeding_mode",
                    "BF_mode","Feeding_mode","Weekly_BF_perc")
PdBF <- PdBF %>% select(-Current_feeding_mode) %>% 
  filter(complete.cases(Timepoint))
PdBF$Timepoint <- ifelse(PdBF$Timepoint=="1個月","M1","M2")
PdBF$BF_mode <- ifelse(PdBF$BF_mode=="混合 <貼身/埋身餵哺> 和 <使用奶瓶餵哺>","Mix",
                       ifelse(PdBF$BF_mode=="幾乎全使用奶瓶餵哺 (≥90%時間使用奶瓶)","Pump","Direct"))
PdBF$Feeding_mode <- ifelse(PdBF$Feeding_mode=="混合","Mix",
                            ifelse(grepl("母乳",PdBF$Feeding_mode),"BF","Formula"))
PdBF$Weekly_BF_perc <- str_replace(PdBF$Weekly_BF_perc,"多於",">")
PdBF <- reshape(PdBF, idvar = "MOM_SNO", 
                timevar = "Timepoint", direction = "wide")
Y1db <- left_join(Y1db,PdBF,"MOM_SNO")


PdBF_M6Y1 <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                                sheet = "REDCap Qs (D01-M12)")) %>% 
  filter(Event.Name != "Date of Born (Arm 2: Infant)")
PdBF_M6Y1 <- PdBF_M6Y1[,c(1,54,55)]
colnames(PdBF_M6Y1) <- c("MOM_SNO","Timepoint","Feeding_mode")
PdBF_M6Y1 <- PdBF_M6Y1 %>% filter(complete.cases(Timepoint))
PdBF_M6Y1$Timepoint <- ifelse(PdBF_M6Y1$Timepoint=="6個月","M6","M12")
PdBF_M6Y1$Feeding_mode <- ifelse(PdBF_M6Y1$Feeding_mode=="混合","Mix",
                            ifelse(grepl("母乳",PdBF_M6Y1$Feeding_mode),"BF","Formula"))
PdBF_M6Y1 <- reshape(PdBF_M6Y1, idvar = "MOM_SNO", 
                timevar = "Timepoint", direction = "wide")
Y1db <- left_join(Y1db,PdBF_M6Y1,"MOM_SNO")
# Disease summary (clinical notes)
Pdnotes <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                                 sheet = "Paed Notes")[,-c(1:2)]) %>% 
  filter(ATTEND != "NO") %>% 
  select(-VISIT_S, -ATTEND, -Anthropometric.Measurement_B,
         -OVERALL_R, -Disease.Status.Changed1_B_NO,
         -Disease.Status.Changed2_B_NO, -Disease.Status.Changed3_B_NO)
colnames(Pdnotes)[11:13] <- c("Head.Circumference","Body.Weight",
                              "Body.Height")
Pdnotes <- reshape(Pdnotes, idvar = "MOM_SNO", 
                   timevar = "VISIT", direction = "wide")
Pdnotes <- Pdnotes[,c(1:12,24:45,13:23)]
Y1db <- left_join(Y1db,Pdnotes,"MOM_SNO")

# Antibiotic use
PdATB <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                                sheet = "ATB Record")[,-c(1,2)])
PdATB$STOPD[is.na(PdATB$STOPD)] <- ymd(substr(PdATB$STARTD[is.na(PdATB$STOPD)] +
  PdATB$DURATION[is.na(PdATB$STOPD)],0,10))
PdATB <- PdATB %>% mutate(ATB_Period = paste0(STARTD," - ",STOPD),
                          ATB_Dur = as.numeric(ymd(STOPD)-ymd(STARTD))+1) %>% 
  group_by(MOM_SNO) %>% 
  summarise(ATB_used = paste0(ATB_NAME, collapse = ";"),
            ATB_indication = paste0(IND_OTH, collapse = ";"),
            ATB_period = paste0(ATB_Period, collapse = ";"),
            ATB_dur = paste0(ATB_Dur, collapse = ";"),
            ATB_dursum = sum(ATB_Dur,na.rm = T),
            ATB_route = paste0(ROUTE[!duplicated(ROUTE)],collapse = ";")) %>% 
  filter(ATB_route!="TOPICAL"|ATB_route!="OTHER")
Y1db <- left_join(Y1db,PdATB,"MOM_SNO")

# Other medication
PdMed <- data.frame(read_excel("Data/Clincal data_M12/0015_MOM_LZ_20210421B5F5 (Clinical data).xlsx", 
                               sheet = "Med Record")[,-c(1:2)])
PdMed$STOPD[is.na(PdMed$STOPD)] <- ymd(substr(PdMed$STARTD[is.na(PdMed$STOPD)] +
                                                PdMed$DURATION[is.na(PdMed$STOPD)],0,10))
PdMed <- PdMed %>% mutate(Med_Period = paste0(STARTD," - ",STOPD),
                          Med_Dur = as.numeric(ymd(STOPD)-ymd(STARTD))+1) %>% 
  group_by(MOM_SNO) %>% 
  summarise(Med_used = paste0(MED_NAME, collapse = ";"),
            Med_indication = paste0(IND_OTH, collapse = ";"),
            Med_period = paste0(Med_Period, collapse = ";"),
            Med_dur = paste0(Med_Dur, collapse = ";"),
            Med_dursum = sum(Med_Dur,na.rm = T),
            Med_dursum = sum(Med_Dur,na.rm = T),
            Med_route = paste0(ROUTE[!duplicated(ROUTE)],collapse = ";")) %>% 
  filter(Med_route!="TOPICAL"|ATB_route!="OTHER")
Y1db <- left_join(Y1db,PdMed,"MOM_SNO")
Y1db$Any_ATB_Y1 <- ifelse(is.na(Y1db$ATB_used),0,1)
Y1db$Any_Med_Y1 <- ifelse(is.na(Y1db$Med_used),0,1)

Y1metadata <- Y1db[,c(1:6,25,28,30,31,32:39,43:50,54:61,
                      65:72,88,89)]
Y1metadata <- Y1metadata %>% select(colnames(Y1metadata)[!grepl("_B_R.",colnames(Y1metadata))])

Y1metadata$Samp_name <- colnames(OTU_list)[match(substr(Y1metadata$MOM_SNO,7,9),
                                                 substr(colnames(OTU_list),2,4))]
Y1metadata$Samp_name <- str_replace(Y1metadata$Samp_name,"X","")
Y1metadata <- Y1metadata[,c(33,1:32)]
Y1metadata$Samp_name <- str_replace(Y1metadata$Samp_name,"\\.","-")
# 
# Y1metadata$Any_ATB_Y1 <- Y1db$Any_ATB_Y1[match(Y1metadata$MOM_SNO,Y1db$MOM_SNO)]
# Y1metadata$Any_Med_Y1 <- Y1db$Any_Med_Y1[match(Y1metadata$MOM_SNO,Y1db$MOM_SNO)]
Y1metadata <- Y1metadata %>% select(-Feeding_mode.M12,-Feeding_mode.M6)
Y1metadata$BMI_Group <- MOMcomb[MOMcomb$Sex=="F",]$BMI_Group[match(substr(Y1metadata$MOM_SNO,6,9),
                                                                 substr(MOMcomb[MOMcomb$Sex=="F",]$Sub_ID,6,9))]
Y1metadata[,c(3,10,15,20,25)] <- apply(Y1metadata[,c(3,12,17,22,27)],2,
                                       function(x) as.character(ymd(x)))
writexl::write_xlsx(Y1metadata,"Y1_metadata_22032022.xlsx")

saveRDS(Y1db,"Data/Y1_Child_dbsum.rds")

Refgroup <- data.frame(read_excel("Data/Year 1/MOMmy_M12 Microbial Report_Reference Group Clinical Data 220403_TF 220411.xlsx"))
# Refgroup <- Refgroup %>% filter(ROUTE!="TOPICAL") %>%
#   filter(!duplicated(Study.No.))
Y1metadata$Refgroup <- Refgroup$Reference[match(Y1metadata$MOM_SNO,Refgroup$MOM_SNO)]
Y1metadata$Samp_name <- Refgroup$Filename[match(Y1metadata$MOM_SNO,Refgroup$MOM_SNO)]
Y1metadata$ATB_TF <- Refgroup$Study.No.[match(Y1metadata$MOM_SNO,Refgroup$MOM_SNO)]
Y1metadata$Med_TF <- Refgroup$Study.No.[match(Y1metadata$MOM_SNO,Refgroup$MOM_SNO)]
Y1metadata[,c(33,34)] <- apply(Y1metadata[,c(33,34)],2,function(x)
  ifelse(is.na(x),0,1))
Y1metadata$Samp_name <- colnames(Y1OTU)[match(substr(Y1metadata$MOM_SNO,7,9),
                                              substr(colnames(Y1OTU),2,4))]
Y1metadata$Samp_name <- str_replace(Y1metadata$Samp_name,"X","")
Y1metadata$Samp_name <- str_replace_all(Y1metadata$Samp_name,"\\.","-")

PD_dis_group <- data.frame(read_excel("Data/Year 1/MOMmy_M12 Microbial Report_Reference Group Clinical Data 220403_TF 16042022.xlsx", 
                                       sheet = "Research Notes"))  
PD_dis_group <- PD_dis_group %>% 
  mutate(Ecz = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
    ifelse(any(grepl("ECZEMA|EZCEMA",x)),1,0))),
    Rash = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("RASH",x)),1,0))),
    Allergy = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("ALLERG",x)),1,0))),
    GI = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("DISTENSION|STOOL|CONSTIPATION|DIARRHEA|GASTROENTERITIS|
                       GASTROINTESTIONAL|NOROVIRUS INFECTION|ROTA VIRUS",x)),1,0))),
    URTI = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("COUGHING|RHINOVIRUS INFECTION|FEVER|URT",x)),1,0))),
    Vas = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("CUTIS MARMORATA|PORT-WINE STAIN",x)),1,0))),
    UTI = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("UTI",x)),1,0))),
    GE = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("GE",x)),1,0))),
    HEART = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("HEART",x)),1,0))),
    RI = as.character(apply(PD_dis_group[,c(9:16)],1, function(x) 
      ifelse(any(grepl("ROSEOLA",x)),1,0)))) 
PD_dis_group <- PD_dis_group %>% 
  filter(complete.cases(Disease.Status.Changed1_B))

PD_dis_group <- PD_dis_group %>% 
  mutate(Any_skin = as.character(apply(PD_dis_group[,c(21:23)],1, function(x) 
    ifelse(any(x==1),1,0))),
    Other_CD = ifelse(is.na(Disease_Notinclude),0,
                      ifelse(Disease_Notinclude==1,0,1)))

for(i in unique(PD_dis_group$Study.No.)){
  if(any(PD_dis_group$Any_skin[PD_dis_group$Study.No.==i]==1)){
    PD_dis_group$Any_skin[PD_dis_group$Study.No.==i] <- 1
  }
  if(any(PD_dis_group$Ecz[PD_dis_group$Study.No.==i]==1)){
    PD_dis_group$Ecz[PD_dis_group$Study.No.==i] <- 1
  }
  if(any(PD_dis_group$Rash[PD_dis_group$Study.No.==i]==1)){
    PD_dis_group$Rash[PD_dis_group$Study.No.==i] <- 1
  }
  if(any(PD_dis_group$Allergy[PD_dis_group$Study.No.==i]==1)){
    PD_dis_group$Allergy[PD_dis_group$Study.No.==i] <- 1
  }
  if(any(PD_dis_group$Other_CD[PD_dis_group$Study.No.==i]==1)){
    PD_dis_group$Other_CD[PD_dis_group$Study.No.==i] <- 1
  }
}
PD_dis_group_M6 <- PD_dis_group %>% filter(VISIT=="M6"|Any_skin==1|Other_CD==1) %>% 
  filter(is.na(Disease_Notinclude)|Other_CD==1|Any_skin==1) %>% 
  filter(!duplicated(MOM_SNO)) %>% 
  mutate(VISIT = "M6")


PD_dis_group_M12 <- PD_dis_group %>% filter(VISIT=="M12"|Any_skin==1|Other_CD==1) %>% 
  # select(MOM_SNO,Disease_Notinclude,Any_skin,Other_CD,VISIT) %>% 
  filter(is.na(Disease_Notinclude)|Other_CD==1|Any_skin==1) %>% 
  filter(!duplicated(MOM_SNO))%>% 
  mutate(VISIT = "M12")
# colnames(PD_dis_group_M12)[24:30] <- paste0(colnames(PD_dis_group_M12)[24:30],"_M12")
PD_dis_group_M12$Disease <- apply(PD_dis_group_M12[,21:30], 1, function(x)
  paste0(colnames(PD_dis_group_M12)[21:30][x==1], collapse = ", "))
PD_dis_group_M12$Disease[PD_dis_group_M12$Disease==""] <- "Others"
PD_dis_group_M12$Disease <- ifelse(grepl(",",PD_dis_group_M12$Disease),"Multiple",
                                   PD_dis_group_M12$Disease)

Y1metadata <- left_join(Y1metadata,PdBF_M6Y1,"MOM_SNO")
Y1metadata <- Y1metadata[,c(1:7,15,16,8:13)]

Y1metadata$Disease_M6 <- ifelse(is.na(PD_dis_group_M6$MOM_SNO[match(Y1metadata$MOM_SNO,
                                                       PD_dis_group_M6$MOM_SNO)]),0,1)
Y1metadata$Disease_M12 <- ifelse(is.na(PD_dis_group_M12$MOM_SNO[match(Y1metadata$MOM_SNO,
                                                                    PD_dis_group_M12$MOM_SNO)]),0,1)
Y1metadata$Disease_Spe_M12 <- ifelse(is.na(PD_dis_group_M12$Disease[match(Y1metadata$MOM_SNO,
                                                                      PD_dis_group_M12$MOM_SNO)]),"No",PD_dis_group_M12$Disease[match(Y1metadata$MOM_SNO,
                                                                                                                                      PD_dis_group_M12$MOM_SNO)])

a <- PD_dis_group[!duplicated(PD_dis_group$MOM_SNO),]
Y1metadata <- left_join(Y1metadata,a[,c(3,21:25)],"MOM_SNO")
Y1metadata[,c(18:22)] <- apply(Y1metadata[,c(18:22)], 2, function(x)
  ifelse(is.na(x),0,as.numeric(paste(x))))
Y1metadata <- left_join(Y1metadata,PdBF,"MOM_SNO")
Y1metadata <- Y1metadata[,c(1:7,23:28,8:22)]
writexl::write_xlsx(Y1metadata,"Data/Year 1/Y1_metadata_16042022.xlsx")

# Selected samples checking (11082022)----

a <- data.frame(read_xlsx("~/Dropbox/g/FA/Logistics/Sample_application/0015_MOM_LZ_20210421B8_30062022_ELM_10072022.xlsx",
                          skip = 1))

colnames(a)[grepl("STUDY_NO",colnames(a))] <- "STUDY_NO"
colnames(a)[grepl("Sample.SNO...Selected.by.applicants",colnames(a))] <- "SNO_sel"
a <- rbind(a[,c(1,9)],
           a[,c(16,24)],
           a[,c(32,40)],
           a[,c(48,56)],
           a[,c(63,71)],
           a[,c(76,84)],
           a[,c(91,99)],
           a[,c(107,115)],
           a[,c(123,131)])
b <- data.frame(read_xlsx("~/Dropbox/g/FA/Logistics/Sample_application/a.xlsx"))
a <- a %>% filter(complete.cases(SNO_sel)) %>% 
  filter(STUDY_NO %!in% b$STUDY.NO) %>% 
  filter(startsWith(SNO_sel,"ST"))

b <- data.frame(read_xlsx("~/Dropbox/g/FA/Logistics/Sample_application/0015_MOM_LZB8_05082022_status.xlsx",
                          sheet = "Sheet1"))

b <- b %>% filter(complete.cases(Check))

a$Check <- b$Check[match(a$SNO_sel,b$Check)]
a$Check[is.na(a$Check)] <- b$...2[match(a$SNO_sel[is.na(a$Check)],b$...2)]
b$REMARKS <- str_remove(b$REMARKS,"SUBSTITUTE ")
b$REMARKS[b$...5=="MISMATCH"&!startsWith(b$REMARKS,"ST")] <- paste0("ST",b$REMARKS[b$...5=="MISMATCH"&!startsWith(b$REMARKS,"ST")])
a$Check[is.na(a$Check)] <- b$REMARKS[match(a$SNO_sel[is.na(a$Check)],b$REMARKS)]
b$PWH_No <- a$STUDY_NO[match(b$Check,a$Check)]
b$PWH_No[is.na(b$PWH_No)] <- a$STUDY_NO[match(b$REMARKS[is.na(b$PWH_No)],a$Check)]
writexl::write_xlsx(b,"~/Dropbox/g/FA/Logistics/Sample_application/MOMmy_check.xlsx")

  
  