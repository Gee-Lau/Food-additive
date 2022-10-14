library(reshape2)
library(Hmisc)
library(car)
library(dplyr)
library(tableone)
library(readxl)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(stringi)
library(ggpubr)

# Data cleaning----
PT_db <- data.frame(read_excel("Data/Preterm_Control_Questionnaires/0015_MOM_LZ_20210421B4F4 preterm (Clinical data).xlsx", 
                               sheet = "Subj Info"))
PT_db <- PT_db[,c(1,5,4,6)] %>% filter(complete.cases(Preterm.Baby)) %>% 
  select(-Preterm.Baby)
colnames(PT_db) <- c("Sub_ID_Child","DOB","Sex")


PT_clinc <- data.frame(read_excel("Data/Preterm_Control_Questionnaires/0015_MOM_LZ_20210421B4F4 preterm (Clinical data).xlsx", 
                                  sheet = "Paed Notes_Preterm"))

PNdb_PT <- data.frame(read_excel("Data/Preterm_Control_Questionnaires/BabyDeliveryQ.xlsx"))[,-c(2:5,12:15,23,24,31,56:92)]
colnames(PNdb_PT) <- colnames(PNdb)[102:145]
PT_db$GA <- PNdb_PT$GA_Week[match(PT_db$Sub_ID_Child,PNdb_PT$Sub_ID_Child)]
PT_clinc$DOB <- PT_db$DOB[match(PT_clinc$Preterm.Study.No.,PT_db$Sub_ID_Child)]
PT_clinc$TTV <- as.numeric(ymd(PT_clinc$VISIT_D)-ymd(PT_clinc$DOB))
PT_clinc$TTV <- floor(PT_clinc$TTV/7)
PT_clinc$TTGA <- PT_db$GA[match(PT_clinc$Preterm.Study.No.,PT_db$Sub_ID_Child)]+
  PT_clinc$TTV
PT_db$PT_Group <- ifelse(PT_db$GA<32,"VtE","Moederate")
PT_clinc$PT_Group <- PT_db$PT_Group[match(PT_clinc$Preterm.Study.No.,PT_db$Sub_ID_Child)]
PT_HC <- PT_clinc %>% select(MOM_SNO,Head.Circumference..cm.,
                             TTGA, PT_Group) %>% 
  filter(complete.cases(Head.Circumference..cm.))

ggplot(PT_HC, aes(x = TTGA,
                  y = Head.Circumference..cm.)) +
  geom_line(aes(group = MOM_SNO, color = PT_Group), alpha = 0.2) +
  geom_point(aes(color = PT_Group), alpha = 0.3, 
             size = 1) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  scale_x_continuous(breaks = c(35,50,65,80,95)) +
  xlab("Postmenstrual age") +
  ylab("Head Circumference (cm)") +
  geom_smooth(aes(group = PT_Group,
                  color = PT_Group,
                  fill = PT_Group),
              method='lm', formula = y ~ splines::bs(x, 3),
              size = 1.5, alpha = 0.35) +
  scale_color_manual(values = c("#fb6a4a", "#cb181d"),
                     labels = c("Moderate/low", "Very/extremely")) +
  scale_fill_manual(values = c("#fb6a4a", "#cb181d"),
                    labels = c("Moderate/low", "Very/extremely"))

PTc_clinc <- data.frame(read_excel("Data/Preterm_Control_Questionnaires/0015_MOM_LZ_20210421B4F4 preterm (Clinical data).xlsx", 
                                   sheet = "Paed Notes_Control"))
PTc_clinc$DOB <- PNdb_01112021$DOB_child[match(PTc_clinc$Preterm.Study.No.,PNdb_01112021$Sub_ID_Child)]
PTc_clinc$TTV <- as.numeric(ymd(PTc_clinc$VISIT_D)-ymd(PTc_clinc$DOB))
PTc_clinc$TTV <- floor(PTc_clinc$TTV/7)
PTc_clinc$TTGA <- PNdb_01112021$GA_Week[match(PTc_clinc$Preterm.Study.No.,PNdb_01112021$Sub_ID_Child)]+
  PTc_clinc$TTV

PT_HC <- rbind(PT_HC,PTc_clinc %>% 
                 filter(TTV>0) %>% 
                 select(MOM_SNO,Head.Circumference..cm.,
                                         TTGA) %>% 
                 filter(complete.cases(Head.Circumference..cm.)) %>% 
                 filter(Head.Circumference..cm.>20))
PT_HC$Group <- NA
PT_HC$Group[PT_HC$MOM_SNO %in% PT_db$Sub_ID_Child] <- "Preterm"
PT_HC$Group[is.na(PT_HC$Group)] <- "Fullterm"






