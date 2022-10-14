library(reshape2)
library(Hmisc)
library(tableone)
library(readxl)
library(tidyverse)
library(lubridate)
library(rgl)
library(patchwork)
library(MicrobiotaProcess)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(phyloseq)
library(vegan)
library(grid)
library(stringr)
library(microbiomeMarker)
library(Maaslin2)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(file2meco)
library(ComplexHeatmap)


# Correlation between MSS results and FA ques----

Urine_SweetT1 <- data.frame(read_xlsx("Data/MS_Sweetener/MOMmy_Sweetener_data_09062022.xlsx",sheet = "Urine_T1"))
Urine_SweetT3<- data.frame(read_xlsx("Data/MS_Sweetener/MOMmy_Sweetener_data_09062022.xlsx",sheet = "Urine_T3"))

Urine_SweetT3$Calculated.Amt..ng.ml. <- as.numeric(Urine_SweetT3$Calculated.Amt..ng.ml.)
Urine_SweetT3$FA_ques <- MomFA$SAC[match(Urine_SweetT3$ID,MomFA$Sub_ID)]
Urine_SweetT3$FA_ques[Urine_SweetT3$Compound=="Aspartame"] <- MomFA$ASP[match(Urine_SweetT3$ID[Urine_SweetT3$Compound=="Aspartame"],
                                                                              MomFA$Sub_ID)]
Urine_SweetT3$FA_ques[Urine_SweetT3$Compound=="Sucralose"] <- MomFA$SUC[match(Urine_SweetT3$ID[Urine_SweetT3$Compound=="Sucralose"],
                                                                              MomFA$Sub_ID)]

Urine_SweetT3$T1_Sweet <-  Urine_SweetT1$Calculated.Amt..ng.ml.[match(paste0(Urine_SweetT3$ID,Urine_SweetT3$Compound),
                                                                      paste0(Urine_SweetT1$ID,Urine_SweetT1$Compound))]


Urine_SweetT3CB<- data.frame(read_xlsx("Data/MS_Sweetener/MOMmy_Sweetener_data_09062022.xlsx",sheet = "CB"))
Urine_SweetT3$T3_CB_Sweet <-  Urine_SweetT3CB$Calculated.Amt..ng.ml.[match(paste0(Urine_SweetT3$ID,Urine_SweetT3$Compound),
                                                                           paste0(Urine_SweetT3CB$ID,Urine_SweetT3CB$Compound))]
Urine_SweetT3MB<- data.frame(read_xlsx("Data/MS_Sweetener/MOMmy_Sweetener_data_09062022.xlsx",sheet = "MB"))
Urine_SweetT3$T3_MB_Sweet <-  Urine_SweetT3MB$Calculated.Amt..ng.ml.[match(paste0(Urine_SweetT3$ID,Urine_SweetT3$Compound),
                                                                           paste0(Urine_SweetT3MB$ID,Urine_SweetT3MB$Compound))]
colnames(Urine_SweetT3)[3:4] <- c("Sub_ID","T3_Sweet")
Urine_SweetT3[,c(7:8)] <- apply(Urine_SweetT3[,c(7:8)],2,as.numeric)
colnames(Urine_SweetT3)[6:8] <- c("T1_Urine","T3_Cord_blood","T3_Maternal_blood")

Urine_SweetT3_new <- Urine_SweetT3
Urine_SweetT3_new[is.na(Urine_SweetT3_new)] <- 0
Urine_SweetT3_new$T1_Urine <- Urine_SweetT1$Calculated.Amt..ng.ml.[match(
  paste0(Urine_SweetT3_new$Compound,Urine_SweetT3_new$Sub_ID),
  paste0(Urine_SweetT1$Compound,Urine_SweetT1$ID)
)]
Urine_SweetT3_new$T1_Urine[Urine_SweetT3_new$T1_Urine=="ND"] <- 0
Urine_SweetT3_new$T1_Urine[Urine_SweetT3_new$T1_Urine<0] <- 0
Urine_SweetT3_new$T1_Urine <- as.numeric(Urine_SweetT3_new$T1_Urine)
# Urine_SweetT3_new$Sample.ID[Urine_SweetT3_new$Sample.ID=="MU-T3"] <- 
#   paste0(Urine_SweetT3_new$Sample.ID[Urine_SweetT3_new$Sample.ID=="MU-T3"],
#          substr(Urine_SweetT3_new$Sub_ID[Urine_SweetT3_new$Sample.ID=="MU-T3"],7,9))
# Urine_SweetT3_new <- Urine_SweetT3_new %>% filter(Sample.ID!="MU-T3868")
# Urine_SweetT3_new <- Urine_SweetT3_new %>% filter(Sample.ID!="MU-T3900")


Urine_SweetT3 <- reshape(Urine_SweetT3_new[,-c(2,5)], idvar = "Sub_ID",
                         timevar = "Compound", direction = "wide")

Urine_SweetT3[,-1] <- apply(Urine_SweetT3[,-1] , 2,
                            function(x) if (sum(x==0,na.rm = T) >
                                            length(x)/2){
                              ifelse(x==0,"ND","Pos")
                            }else{
                              return(x)
                            })
for(i in colnames(Urine_SweetT3)[-1]){
  if(any(Urine_SweetT3[,i]=="ND",na.rm = T)){
    next
  }else{
    Urine_SweetT3[,i] <- as.numeric(Urine_SweetT3[,i])
  }
}

plot_list = list()
for(j in colnames(Urine_SweetT3_new)[5:8]){
  for (i in unique(Urine_SweetT3_new$Compound)) {
    my.data <- Urine_SweetT3_new %>% select(1:4,j) %>% 
      filter(Compound==i)
    colnames(my.data)[c(4,5)] <- c("x","y")
    if(all(my.data$x==0|my.data$y==0))
      next
    a <- cor.test(my.data$x,
                  my.data$y,method = "spearman")
    if(a$p.value<0.05){
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_point(size = 2, color = "#542788") +
        theme_classic() +
        xlab(paste(i," in T3 urine")) +
        ylab(paste(i," in ", j)) +
        stat_smooth(method = "glm",formula = y~x,  
                    se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
        ggtitle(paste0("NA in T3 Urine = ",sum(my.data$x==0),
                       "\nNA in ",j, " = ",
                       sum(my.data$y==0),
                       " (No. of not detected)")) +
        theme(legend.position = "none",
              axis.text = element_text(colour = "black"),
              plot.background = element_rect(fill = "white", colour = "red3", size = 1)) +
        annotate("text", x = max(my.data$x,na.rm = T)*3/4, 
                 y = max(my.data$y,na.rm = T), size = 4, 
                 label = paste0("R =", round(a$estimate,2),
                                "\np-value = ", round(a$p.value,3)))
      plot_list <- c(plot_list, list(p))
      
    }else{
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_point(size = 2, color = "#542788") +
        theme_classic() +
        xlab(paste(i," in T3 urine")) +
        ylab(paste(i," in ", j)) +
        stat_smooth(method = "glm",formula = y~x,  
                    se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
        ggtitle(paste0("NA in T3 Urine = ",sum(my.data$x==0),
                       "\nNA in ",j, " = ",
                       sum(my.data$y==0),
                       " (No. of not detected)")) +
        theme(legend.position = "none",
              axis.text = element_text(colour = "black")) + 
        annotate("text", x = max(my.data$x,na.rm = T)*3/4,
                 y = max(my.data$y,na.rm = T), size = 4, 
                 label = paste0("R =", round(a$estimate,2),
                                "\np-value = ", round(a$p.value,3)))
      
      plot_list <- c(plot_list, list(p))
      
    }
  }
}


plot_list <- plot_list[1:6]

for (i in unique(Urine_SweetT3_new$Compound)) {
  my.data <- Urine_SweetT3_new %>% select(1:3,T1_Urine, FA_ques) %>% 
    filter(Compound==i)
  colnames(my.data)[c(4,5)] <- c("x","y")
  if(all(my.data$x==0|my.data$y==0))
    next
  a <- cor.test(my.data$x,
                my.data$y,method = "spearman")
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_point(size = 2, color = "#542788") +
      theme_classic() +
      xlab(paste(i," in T1 urine")) +
      ylab(paste(i," in T1 questionnaire")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
      ggtitle(paste0("NA in T1 Urine = ",sum(my.data$x==0))) +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black"),
            plot.background = element_rect(fill = "white", colour = "red3", size = 1)) +
      annotate("text", x = max(my.data$x,na.rm = T)*3/4, 
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3)))
    plot_list <- c(plot_list, list(p))
    
  }else{
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_point(size = 2, color = "#542788") +
      theme_classic() +
      xlab(paste(i," in T1 urine")) +
      ylab(paste(i," in T1 questionnaire")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
      ggtitle(paste0("NA in T1 Urine = ",sum(my.data$x==0))) +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black")) + 
      annotate("text", x = max(my.data$x,na.rm = T)*3/4,
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3)))
    
    plot_list <- c(plot_list, list(p))
    
  }
}

pdf(file="Results/FA_MS_Sweetener/MS_Sweetener_Sum.pdf", width = 13,
    height = 14)
wrap_plots(plot_list,ncol = 3)
invisible(dev.off())

# Check the sweetener level in different emulsifier intake groups
Samp_list <- data.frame(read_xlsx("Data/FA_Samptest_pilot_list.xlsx"))
Urine_SweetT1$Emu_group <- Samp_list$Group[match(Urine_SweetT1$ID,Samp_list$STUDY_NO)]

plot_list = list()
for (i in unique(Urine_SweetT1$Compound)) {
  my.data <- Urine_SweetT1 %>% filter(Compound==i)
  colnames(my.data)[c(4,6)] <- c("y","x")
  my.data$x <- factor(as.factor(my.data$x),
                            levels(as.factor(my.data$x))[c(1,3,2)])
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
    scale_color_manual(values=c("#5e3c99", "#fdb863","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "kruskal",label.x = 1) +
    stat_compare_means(comparisons = statecomp,label = "p.signif",
                       method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 10),
          axis.title = element_text(size = 12)) +
    xlab(paste0(i)) +
    ylab(paste(i," in cord blood")) 
  
  
  plot_list <- c(plot_list, list(p))
  
}

# MS results and outcomes----
# 1. T3 MS results
# 1.1 DB preparation -- Continuous
Urine_SweetT3_new <- gather(Urine_SweetT3[,!grepl("FA|T1",colnames(Urine_SweetT3))],
                            key = "Sub_ID")
colnames(Urine_SweetT3_new) <- c("Var","value")
Urine_SweetT3_new$Sweet <- as.character(
  lapply(Urine_SweetT3_new$Var,function(x) strsplit(x,"\\.")[[1]][2])
)
Urine_SweetT3_new$Source <- str_remove(as.character(
  lapply(Urine_SweetT3_new$Var,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
Urine_SweetT3_new$Source <- str_replace(Urine_SweetT3_new$Source,"_"," ")
Urine_SweetT3_new$Source[Urine_SweetT3_new$Source=="Sweet"] <- "Urine"
Urine_SweetT3_new$value <- as.numeric(Urine_SweetT3_new$value)
Urine_SweetT3_new <- Urine_SweetT3_new %>% filter(complete.cases(value))

Urine_Sweet_T3_bi <- data.frame(Var = NA, Per = NA)

for(i in colnames(Urine_SweetT3)[!grepl("FA|T1",colnames(Urine_SweetT3))]){
  if(is.character(Urine_SweetT3[,i])){
    Urine_Sweet_T3_bi <- rbind(Urine_Sweet_T3_bi,
                               c(i,round(sum(Urine_SweetT3[,i]=="ND")/nrow(Urine_SweetT3)*100,1)))
  }
}

# 1.1 DB preparation -- Categorical
Urine_Sweet_T3_bi <- Urine_Sweet_T3_bi[-c(1,2),]
Urine_Sweet_T3_bi <- data.frame(rbind(Urine_Sweet_T3_bi,Urine_Sweet_T3_bi))
Urine_Sweet_T3_bi$Per <- as.numeric(Urine_Sweet_T3_bi$Per)
Urine_Sweet_T3_bi$Per[5:8] <- 100-Urine_Sweet_T3_bi$Per[1:4] 
Urine_Sweet_T3_bi$Sweet <- as.character(
  lapply(Urine_Sweet_T3_bi$Var,function(x) strsplit(x,"\\.")[[1]][2])
)
Urine_Sweet_T3_bi$Source <- str_remove(as.character(
  lapply(Urine_Sweet_T3_bi$Var,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
Urine_Sweet_T3_bi$Source <- str_replace(Urine_Sweet_T3_bi$Source,"_"," ")
Urine_Sweet_T3_bi$Dec <- c(rep("ND",4),rep("Detected",4))

# 1.2 Plot 
Sweet_T3_plt <- ggplot(Urine_SweetT3_new, aes(x = Source, y = value)) +
  geom_boxplot(width = 0.5, size = 2,
               fatten = 0.5,
               outlier.shape = NA,
               position = position_dodge(0.9),
               col = "#fd8d3c") +
  geom_point(shape=16, position=position_jitter(0.3),
             size = 2, alpha = 0.8,
             col = "#fd8d3c") +
  theme_classic()+
  facet_wrap(.~Sweet, scale = "free") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 12),
        axis.text.x =  element_text(color = "black",size = 12, angle = 45,
                                    hjust = 1),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="gray60",color = "black"),
        strip.text = element_text(size=15, colour="white")) +
  ylab("Trimester three sweeteners (ng/mL)") 

p1 <- ggplot(Urine_Sweet_T3_bi[Urine_Sweet_T3_bi$Sweet=="Aspartame",],
       aes(x = Source, y = Per)) +
  geom_bar(stat = "identity",
           aes(fill = Dec, alpha = 0.6), color = "gray30") +
  ggplot2::scale_fill_manual(values=c("#fd8d3c","gray85")) +
  geom_text(aes(label = c("ND: 100","ND: 92.6",NA,NA)), vjust = 1.5, color = "black",
            position = position_dodge(0.9), size=5, alpha = 0.6) +
  theme_classic()+
  ggtitle("%") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        axis.text = element_text(color = "black",size = 12),
        axis.text.x =  element_text(color = "black",size = 12, angle = 45,
                                    hjust = 1),
        axis.title = element_blank()) +
  ylab("% of not-detactable (ND)") 

p2 <- ggplot(Urine_Sweet_T3_bi[Urine_Sweet_T3_bi$Sweet=="Sucralose",],
             aes(x = Source, y = Per), 
             fill = "transparent" ) +
  geom_bar(stat = "identity",
           aes(fill = Dec, alpha = 0.6), color = "gray30") +
  ggplot2::scale_fill_manual(values=c("#fd8d3c","gray85")) +
  geom_text(aes(label = c("ND: 98.6","ND: 96.6",NA,NA)), vjust = 1.5, color = "black",
            position = position_dodge(0.9), size=5, alpha = 0.6) +
  theme_classic()+
  ggtitle("%") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        axis.text = element_text(color = "black",size = 12),
        axis.text.x =  element_text(color = "black",size = 12, angle = 45,
                                    hjust = 1),
        axis.title = element_blank()) +
  ylab("% of not-detactable (ND)") 

p3 <- ggplot(Urine_SweetT3_new[Urine_SweetT3_new$Sweet=="Saccharin"&
                           Urine_SweetT3_new$Source!="Urine",],
       aes(x = Source, y = value)) +
  geom_boxplot(width = 0.5, size = 2,
               fatten = 0.5,
               outlier.shape = NA,
               position = position_dodge(0.9),
               col = "#fd8d3c") +
  geom_point(shape=16, position=position_jitter(0.3),
             size = 2, alpha = 0.8,
             col = "#fd8d3c") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 12),
        axis.text.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        strip.background = element_rect(fill="gray60",color = "black"),
        strip.text = element_text(size=15, colour="white")) +
  ylab("Trimester three sweeteners (ng/mL)") 

pdf("Results/FA_MS_Sweetener/T3_Sweets_sum_07072022.pdf", height = 6.8, width = 15)
Sweet_T3_plt + inset_element(p1, left = 0.01, bottom = 0.5, right = 0.15, top = 0.9) +
  inset_element(p3, left = 0.36, bottom = 0.5, right = 0.55, top = 0.93) +
  inset_element(p2, left = 0.7, bottom = 0.5, right = 0.85, top = 0.9)
invisible(dev.off())



# 2. MS results and outcomes
# 2.1 DB cleaning
MB_Sweet <- data.frame(read_xlsx("Data/MOMmy_Sweetener_data.xlsx",sheet = "MB"))
MB_Sweet$Calculated.Amt..ng.ml. <- as.numeric(MB_Sweet$Calculated.Amt..ng.ml.)
MB_Sweet <- reshape(MB_Sweet[,-2], idvar = "ID", timevar = "Compound", direction = "wide")
colnames(MB_Sweet) <- c("Sub_ID","MB_Aspartame","MB_Saccharin","MB_Sucralose")


for(i in c(2:4)){
  print(sum(complete.cases(Urine_Sweet[,i])))
}

Urine_Sweet <- data.frame(read_xlsx("Data/MOMmy_Sweetener_data.xlsx",sheet = "Urine"))
Urine_Sweet$Calculated.Amt..ng.ml. <- as.numeric(Urine_Sweet$Calculated.Amt..ng.ml.)
Urine_Sweet <- reshape(Urine_Sweet[,-2], idvar = "ID", timevar = "Compound", direction = "wide")
colnames(Urine_Sweet) <- c("Sub_ID","Aspartame","Saccharin","Sucralose")

Urine_Sweet <- left_join(Urine_SweetT3, MOMcomb_sub[,c(1,13,16,64:87,100,101)],"Sub_ID")
Urine_Sweet <- left_join(Urine_Sweet, PNdb[,c(2,103:109,118:126,128:134,139,144,150:153)],"Sub_ID")
Urine_Sweet$MtGDM <- as.character(apply(Urine_Sweet[,c(43,44,48)],1,
                                        function(x) ifelse(any(grepl("Yes",x)),"Yes","No")))
Urine_Sweet$GBS <- ifelse(Urine_Sweet$GBS=="是","Yes","No")
colnames(Urine_Sweet)[69] <- "NICU"

# NInfec NHypoG Jaundice
Urine_Sweet[,c(64,63,67)] <- apply(Urine_Sweet[,c(64,63,67)], 2,
                                 function(x) ifelse(is.na(x), NA,
                                                    ifelse(grepl("1No", x), 0, 1)))

Urine_Sweet$Del_mode.y <- ifelse(is.na(Urine_Sweet$Del_mode.y), NA,
                                 ifelse(Urine_Sweet$Del_mode.y=="剖腹產",1,0))
colnames(Urine_Sweet)[104] <- "NICU"
Urine_Sweet$MtGDM_Bi <- ifelse(Urine_Sweet$MtGDM=="Yes",1,0)
Urine_Sweet$GA <- PNdb$GA_Week[match(Urine_Sweet$Sub_ID,PNdb$Sub_ID)]
colnames(Urine_Sweet) <- str_replace_all(colnames(Urine_Sweet),"T3_Sweet","T3_Urine")
Urine_Sweet$MtGDM_Bi[Urine_Sweet$HG=="2Yes"] <- NA

# 2.2 Uni/multi variable
uni_c <- data.frame()
for(j in colnames(Urine_Sweet)[c(61,63,64,67,69,71,72,75)]){
  if (length(unique(Urine_Sweet[,j]))>1) {
    for (i in colnames(Urine_Sweet)[grepl("T3_C|T3_M|T3_U",colnames(Urine_Sweet))]){
      if(length(unique(Urine_Sweet[,i]))>1&is.numeric(Urine_Sweet[,i])){
        print(i)
        temp_a <- wilcox.test(as.formula(paste0(i, "~",j)),
                              data = Urine_Sweet,
                              conf.int = T, conf.level=0.95)
        temp_b <- data.frame(OR = 1-temp_a$estimate,
                             Lower = 1-temp_a$conf.int[2],
                             Upper = 1-temp_a$conf.int[1],
                             p.value = temp_a$p.value,
                        Outcome = j,
                        Var = i)
        uni_c <- rbind(uni_c,temp_b)
      }
    } 
  }
}

for(j in colnames(Urine_Sweet)[c(61,63,67,69,71,72,75)]){
  if (length(unique(Urine_Sweet[,j]))>1) {
    for (i in colnames(Urine_Sweet)[grepl("T3_C|T3_M|T3_U",colnames(Urine_Sweet))&
                                    c(17:23,25:31,35:37,39)]){
      if(length(unique(Urine_Sweet[,i]))>1){
        temp_a <- glm(as.formula(paste0(j, " ~",i)),
                      Urine_Sweet,family ="binomial")
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

uni_c[,c(1:4)] <- apply(uni_c[,c(1:4)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)
colnames(uni_c) <- c("Est/OR","Lower","Upper","p.value","Outcome","Var")
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))
# uni_c$Var_strat <- rownames(uni_c)
# writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_Gcomp_Uni_21062022_all.xlsx")
uni_c <- uni_c %>% filter(p.value<0.06)
writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_Gcomp_Uni_21062022_sig.xlsx")


mglmtb <- data.frame()

for(j in unique(uni_c$Outcome)){
  for (i in unique(uni_c$Var)){
    if (paste0(i) %in% uni_c$Var[uni_c$Outcome==j]){
      print(i)
      temp_a <- glm(as.formula(paste0(j, "~",i, " + Age  + Pre_BMI +", 
                                      paste0(uni_c$Var[!grepl("cut",uni_c$Var)&uni_c$Outcome==j],
                                             collapse = "+"))),
                    Urine_Sweet,family ="binomial")
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

colnames(mglmtb) <-  c("Est/OR","Lower","Upper","p.value","Outcome","Var")
mglmtb[,1:4] <- apply(mglmtb[,1:4], 2, function(x) round(as.numeric(x),4))
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))
mglmtb$Var_Strat <- rownames(mglmtb)

# 2.3 Plot
T3_MS_uni <- uni_c[-c(7,8),]
T3_MS_uni$N <- NA
for(i in unique(T3_MS_uni$Outcome)){
  T3_MS_uni$N[T3_MS_uni$Outcome==i] <- paste0("N = ", sum(Urine_Sweet[,i]==1, na.rm = T),
                                              "/",sum(complete.cases(Urine_Sweet[,i])))
}
T3_MS_uni <- T3_MS_uni[-c(7:9),]

T3_MS_uni$Outcome_lab <- c("Any postnatal\ncomplication","Postnatal\ninfection",
                           rep("\nGDM",3),"\nPreterm")

T3_MS_uni$Sweet <- as.character(
  lapply(T3_MS_uni$Var,function(x) strsplit(x,"\\.")[[1]][2])
)
T3_MS_uni$Source <- str_remove(as.character(
  lapply(T3_MS_uni$Var,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
T3_MS_uni$Source <- str_replace(T3_MS_uni$Source,"_"," ")
colnames(T3_MS_uni)[1] <- "OR"
T3_MS_uni$p.value <- as.numeric(T3_MS_uni$p.value)
T3_MS_uni$Source <- factor(as.factor(T3_MS_uni$Source),
                           levels(as.factor(T3_MS_uni$Source))[c(3,2,1)])

T3_MS_uni_plt <- ggplot(T3_MS_uni, aes(x = OR-1, y = Outcome_lab)) +
  geom_bar(stat = "identity", width = 0.8, 
           position = position_dodge(),
           aes(fill = Source), color = "gray20") +
  geom_errorbar(aes(color = Source,
                    xmin = Lower-1, xmax = Upper-1), width =0.2,
                stat = "identity",
                position = position_dodge(0.8)) +
  geom_text(aes(color = Source,
                label = paste0(N,", p-value = ",
                          round(p.value,3))), 
            hjust = -0.2, vjust = c(rep(-0.5,3),-0.8,3,-0.5), 
            size = 5, fontface = "bold", show.legend = FALSE) +
  
  scale_color_manual(values = c("#fe9929","#ef3b2c","#a50f15")) +
  scale_fill_manual(values = c("#fe9929","#ef3b2c","#a50f15")) +
  scale_alpha(range = c(0.85,0.4)) +
  scale_x_continuous(labels = function(x) x + 1,
                     expand = c(0,0.05),
                     limits = c(0,16.1),
                     breaks = c(seq(0, 16, 5))) +
  theme_classic() +
  facet_wrap(. ~ Sweet,  dir = "v", scale = "free") +
  theme(legend.position = c(0.9,0.87),
        legend.box = "horizontal",
        axis.text = element_text(color = "black",size = 12),
        axis.text.y =  element_text(color = "black",size = 12, 
                                    hjust = 0.5, vjust = 1, angle = 90),
        axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        legend.background = element_rect(size=0.3, linetype="solid",
                                         colour ="black"),
        strip.background = element_rect(fill="gray60",color = "black"),
        strip.text = element_text(size=15, colour="white")) +
  xlab("Odds ratio") 

pdf("Results/FA_MS_Sweetener/T3_Sweet-Outcome_Uni_07072022.pdf", height = 8, width = 10)
T3_MS_uni_plt
invisible(dev.off())

# Urine MB sweetener correlation----
Urine_Sweet <- left_join(Urine_Sweet,MB_Sweet,"Sub_ID")
plot_list = list()
for (i in c(3)) {
  my.data <- data.frame(x = Urine_Sweet[,i],
                        y = Urine_Sweet[,i+107])
  if(all(is.na(my.data$x)))
    next
  a <- cor.test(my.data$x,
                my.data$y,method = "spearman")
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_point(size = 2, color = "#542788") +
      theme_classic() +
      xlab(paste(colnames(Urine_Sweet)[i]," in urine")) +
      ylab(paste(colnames(Urine_Sweet)[i]," in maternal blood")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black"),
            plot.background = element_rect(fill = "white", colour = "red3", size = 1)) +
      annotate("text", x = 0.07,
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3))) +
      xlim(0,max(my.data$x[complete.cases(my.data$y)]))
    plot_list <- c(plot_list, list(p))
    
  }else{
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_point(size = 2, color = "#542788") +
      theme_classic() +
      xlab(paste(colnames(Urine_Sweet)[i]," in urine")) +
      ylab(paste(colnames(Urine_Sweet)[i]," in maternal blood")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black")) + 
      stat_regline_equation(label.x = 3, label.y = 32) +
      stat_cor() +
      xlim(0,max(my.data$x[complete.cases(my.data$y)]))
    plot_list <- c(plot_list, list(p))
    
  }
}


pdf(file="MB_Urine_sweetener.pdf", width = 12, height = 6)
wrap_plots(plot_list,ncol = 1)
invisible(dev.off())


# Combine mtg data----
# Generate binary outcomes by cut offs

# Find outliers
# # get threshold values for outliers
# Tmin = med-(3*mad) 
# Tmax = med+(3*mad) 
# High_OT <- median(Urine_Sweet[,i]) +
# (1.4826 * median(abs(Urine_Sweet[,i]-median(Urine_Sweet[,i])))*3)
# Low_OT <- median(Urine_Sweet[,i]) -
#   (1.4826 * median(abs(Urine_Sweet[,i]-median(Urine_Sweet[,i])))*3)

for (i in c(2,7,10,11,12)){
  
  Urine_Sweet[,ncol(Urine_Sweet)+1] <- ifelse(Urine_Sweet[,i] %in% boxplot(Urine_Sweet[,i], plot=FALSE)$out,
                                              "OL_High", "Normal")
  colnames(Urine_Sweet)[ncol(Urine_Sweet)] <- paste0(colnames(Urine_Sweet)[i],"_cut")
}

SU_mtg <- FAmtgraw

# Below to remove the highest (extreme values for each sweeteners) 
# for(i in colnames(Urine_Sweet)[c(2,7,10,11,12)]){
#   Urine_Sweet_clean[,i][Urine_Sweet_clean[,i] %in% 
#                           head(Urine_Sweet_clean[,i][order(Urine_Sweet_clean[,i],
#                                                          decreasing = T)],3)] <- NA
# }

SU_mtg <- data.frame(left_join(data.frame(SU_mtg@sam_data)[,c(1,9:11)],
                               Urine_Sweet,"Sub_ID"))
rownames(SU_mtg) <- rownames(FAMT3Meta)[match(FAMT3Meta$Sub_ID,
                                              SU_mtg$Sub_ID)]
SU_mtg <- phyloseq(otu_table(as.matrix(FAmtgraw@otu_table), taxa_are_rows = TRUE), 
                   sample_data(SU_mtg), 
                   phyloseq::tax_table(as.matrix(FAmtgraw@tax_table)))

SU_mtg <- subset_samples(SU_mtg, complete.cases(T3_Urine.Aspartame))
SU_mtg <- prune_taxa(taxa_sums(SU_mtg) > 0, SU_mtg)
# GDM differencial species
# SU_mtg <- prune_taxa(tax_table(SU_mtg)[, "Species"] %in% GDM_Diff$f|
#                        tax_table(SU_mtg)[, "Genus"] %in% GDM_Diff$f , SU_mtg)

# 1. Alpha
SU_mtgalpha <- cbind(SU_mtg@sam_data,
                 Richness_spe = t(estimateR(ceiling(t(SU_mtg@otu_table))))[,1],
                 Richness_Chao1 = t(estimateR(ceiling(t(SU_mtg@otu_table))))[,2],
                 Evenness = diversity(t(SU_mtg@otu_table)) / log(specnumber(t(SU_mtg@otu_table))),
                 Shannon = diversity(t(SU_mtg@otu_table), index = "shannon"),
                 InvSimp = diversity(t(SU_mtg@otu_table), index = "invsimpson")) 
# FAalpha[,c(62:74)] <- apply(FAalpha[,c(62:74)], 2, function(x)
#   ifelse(x=="High","High","Normal"))
# my.comp <- list( c("Low","High"),
#                  c("Low","Normal"),
#                 c("Normal","High"))
plot_list = list()
plot_list_con = list()
for (i in colnames(SU_mtgalpha)[grepl("T3",colnames(SU_mtgalpha))&
                                grepl("\\.",colnames(SU_mtgalpha))]) {
  for(j in c(84:88)){
    if(is.numeric(SU_mtgalpha[,i])){
      my.data <- data.frame(x = SU_mtgalpha[,i], 
                            y = SU_mtgalpha[,j],
                            lab = colnames(SU_mtgalpha)[j])
      a <- cor.test(my.data$x,
                    my.data$y,method = "spearman")
      x.lab <- str_remove(i,"T3_")
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_point(size = 1,  aes(color = x)) +
        stat_smooth(method = "glm",formula = y~x,  
                    se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
        theme_bw() +      
        xlab(paste(x.lab)) +
        theme(axis.title.y = element_blank(),
              axis.text = element_text(color = "black", size =13),
              axis.title.x = element_text(size = 16, face = "bold"),
              legend.position = "none",
              panel.grid.major = element_line(colour = "gray70"),
              strip.background = element_rect(fill="gray40",color = "black"),
              strip.text = element_text(size=15, colour="white")) +
        facet_grid(. ~ lab) +
        annotate("text", x = max(my.data$x, na.rm = T)*.8,
                 y = max(my.data$y,na.rm = T), size = 4, 
                 label = paste0("R =", round(a$estimate,2),
                                "\np-value = ", round(a$p.value,3)))
      
      plot_list_con <- c(plot_list_con, list(p))
    } else{ 
      my.data <- data.frame(x = SU_mtgalpha[,i], 
                            y = SU_mtgalpha[,j],
                            lab = colnames(SU_mtgalpha)[j])
      x.lab <- str_replace(i,"_cut"," Group")
      x.lab <- str_remove(i,"T3_")
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_violin(aes(color = x), size = 1,  position = "dodge") +
        geom_point(size = 1,  aes(color = x)) +
        geom_boxplot(width = 0.15, fill="white",aes(color = x)) +
        scale_color_manual(values=c("#BC3C29FF","#00A087FF")) +
        scale_fill_manual(values=c("#BC3C29FF","#00A087FF")) +
        theme_bw() +      
        stat_compare_means(method = "wilcox") +
        stat_signif(label = "p.signif", 
                    size = 0,
                    map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2),
                    margin_top = 0.05,
                    method = "wilcox") +
        xlab(paste(x.lab)) +
        theme(axis.title.y = element_blank(),
              axis.text = element_text(color = "black", size =13),
              axis.title.x = element_text(size = 16, face = "bold"),
              legend.position = "none",
              panel.grid.major = element_line(colour = "gray70")) +
        facet_grid(. ~ lab) +
        theme(strip.background = element_rect(fill="gray40",color = "black"),
              strip.text = element_text(size=15, colour="white"))
      
      plot_list <- c(plot_list, list(p))
    }
  }
}

pdf("Results/FA_MS_Sweetener/MS_Asp_MB_alpha_sum_10072022.pdf", height = 8, width = 16)
wrap_plots(plot_list,ncol = 5) + 
  plot_annotation(title = 'T3 Aspartame (Maternal blood)',
                  theme = theme(plot.title = element_text(size = 18)))
invisible(dev.off())

pdf("Results/FA_MS_Sweetener/MS_con_alpha_sum_10072022.pdf", height = 25, width = 16)
wrap_plots(plot_list_con,ncol = 5)
invisible(dev.off())


# 2. Beta
set.seed(100)
SUnmdspre <- vegdist(t(SU_mtg@otu_table),na.rm = T)
SUnmds <- monoMDS(SUnmdspre)

SUnmds <- data.frame(cbind(data.frame(SU_mtg@sam_data),
                           SUnmds[["points"]]))

SUnmds$MDS1 <- as.numeric(SUnmds$MDS1)
SUnmds$MDS2 <- as.numeric(SUnmds$MDS2)
# gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
#   filter(abs(MDS2)<0.1)

SU_Beta_plt <- ggplot(SUnmds, aes(MDS1, MDS2)) +
  geom_point(size = 4, aes(fill = T3_Urine.Aspartame),
             shape = 21, color = "gray20") +
  # scale_color_manual(name = "CRN intake group",
  #                    values=c("#BC3C29FF","#00A087FF")) +
  scale_fill_gradient2(
    name = "T3 Aspartame (Urine)",
    low = "white",
    mid = "#fec44f",
    high = "#993404",
    midpoint = 0.55,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("") +
  annotate("text",label = paste("PERMANOVA R2 = 0.024\np = 0.011"),
           x = -1.3, y = 1.8, hjust = 0)

pdf("Results/FA_MS_Sweetener/MS_beta_U_Aspartame_11072022.pdf", height = 6, width = 7)
ggMarginal(SU_Beta_plt, colour = "#cc4c02")
invisible(dev.off())


SU_dis <- phyloseq::distance(SU_mtg, 
                             method="bray", weighted=T)
SU_Meta <- data.frame(SU_mtg@sam_data)
set.seed(1)
Uni_per <- data.frame()
for(i in colnames(SU_Meta)[grepl("T3",colnames(SU_Meta))&
                           grepl("\\.",colnames(SU_Meta))]){
  if(length(unique(SU_Meta[,i]))>1){
    temp_db <- SU_Meta %>% filter(complete.cases(SU_Meta[,i]))
    temp_db <-  subset_samples(SU_mtg, Sub_ID %in% temp_db$Sub_ID)
    temp_dis <- phyloseq::distance(temp_db, 
                                   method="bray", weighted=T)
    temp_db <-  adonis2(as.formula(paste0("temp_dis ~ ", i)), 
                       data = data.frame(temp_db@sam_data),
                       permutation=999)
    Uni_per <- rbind(Uni_per,data.frame(temp_db))
  } else {
    print(i)
  }
}


Sig_Uni_per <- Uni_per %>% filter(Pr..F.<0.1)

set.seed(1)
multiper_db <- subset_samples(FAmtgraw, Sub_ID %in% data.frame(FAMT3Meta %>% 
                                                                 select(c("Sub_ID",colnames(FAMT3Meta)[colnames(FAMT3Meta) %in% rownames(Sig_Uni_per)])) %>% 
                                                                 drop_na %>% 
                                                                 select(Sub_ID))$Sub_ID)
multiper_dis <- distance(multiper_db, method="bray", weighted=T)
multi_per <- adonis(as.formula(paste0("multiper_dis ~ ", 
                                      paste0(rownames(Sig_Uni_per),collapse = " + "))), 
                    data = data.frame(multiper_db@sam_data),
                    permutation=999)





# 3. differential species
set.seed(100)
for(i in colnames(SU_mtg@sam_data)[grepl("T3",colnames(SU_mtg@sam_data))&
                                   grepl("\\.",colnames(SU_mtg@sam_data))]){
  if(length(unique(data.frame(SU_mtg@sam_data)[,i]))==3){
    Maaslin2(
      input_data = SU_mtg@otu_table,
      input_metadata = data.frame(SU_mtg@sam_data)[data.frame(SU_mtg@sam_data)[,i]!="Normal",],
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = c(i),
      # random_effects = c(),
      output = paste0("Results/FA_MS_Sweetener/",i,"MaAsLin"),
      correction = "BH",
      reference = "Normal",
      standardize = FALSE,
      cores = 1)
  }else{
    Maaslin2(
      input_data = SU_mtg@otu_table,
      input_metadata = data.frame(SU_mtg@sam_data) ,
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = c(i),
      # random_effects = c(),
      output = paste0("Results/FA_MS_Sweetener/",i,"MaAsLin"),
      correction = "BH",
      reference = "Normal",
      standardize = FALSE,
      cores = 1)
  }
}

Maas_MS <- data.frame()
for(i in colnames(SU_mtg@sam_data)[grepl("T3",colnames(SU_mtg@sam_data))&
                                   grepl("\\.",colnames(SU_mtg@sam_data))]){
  temp_db <- read.table(paste0("Results/FA_MS_Sweetener/",i,"MaAsLin/significant_results.tsv"),
                        header = T)
  Maas_MS <- rbind(Maas_MS,temp_db)
}

Maas_MS <- Maas_MS %>% mutate(feature = as.character(lapply(strsplit(feature[grepl("\\.",feature)],"\\."),
                                                            function(x) paste0(x)[length(x)]))) 
Maas_MS$Ass_Group <- ifelse(Maas_MS$coef<0, "Neg", "Pos")
# Maas_MS$Var_fea <- paste0(Maas_MS$metadata,"-\n",Maas_MS$feature)
# Maas_Clinc$Level <- ifelse(grepl("\\_",Maas_Clinc$feature), "Species", "Genus")
Maas_MS$Sweet <- as.character(
  lapply(Maas_MS$metadata,function(x) strsplit(x,"\\.")[[1]][2])
)
Maas_MS$Sweet <- str_remove(Maas_MS$Sweet,"_cut")
Maas_MS$Source <- str_remove(as.character(
  lapply(Maas_MS$metadata,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
Maas_MS$Source <- str_replace(Maas_MS$Source,"_"," ")
Maas_MS <- Maas_MS[order(Maas_MS$coef,decreasing = T),]
Maas_MS$Order <- length(Maas_MS$feature):1
Maas_MS$feature <- reorder(Maas_MS$feature,Maas_MS$Order)
Maas_MS$GDM <- GDM_Diff$Ass_Group[match(Maas_MS$feature,GDM_Diff$f)]
Maas_MS$Coef_tx <- ifelse(Maas_MS$coef<0, Maas_MS$coef+1, Maas_MS$coef-2)

pdf("Results/FA_MS_Sweetener/MaAsLin_Suc_11072022.pdf", height = 5, width = 10)
ggplot(data = Maas_MS %>% filter(Sweet == "Sucralose"), 
       aes(x = coef, y = feature, 
           fill = Ass_Group)) +
  geom_col(width = 0.75) + 
  geom_label(aes(label = Source, 
                 x = Coef_tx,
            size=4, alpha = 0.6),
            fill = "white") +
  scale_x_continuous(expand = c(0,0.5),
                     breaks = c(seq(-8,4,4))) +
  scale_fill_manual(values=c("#74add1", "#d73027"),
                    name = "Association",
                    labels = c("Negative", "Positive")) +
  geom_vline(xintercept = 0, size = 0.8, linetype="dashed",
             color = "gray40") +
  facet_grid(. ~ Sweet, scales = "free_y", space = "free") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 14),
        axis.text.y =  element_text(size = 14,
                                    face = "italic"),
        # axis.text.y =  element_text(color = c(rep("#238b45",2),
        #                                       rep("black",2),
        #                                       rep("#0570b0",2),
        #                                       rep("black",12)),
        #                             size = 14,
        #                             face = c(rep("bold.italic",2),
        #                                      rep("italic",2),
        #                                      rep("bold.italic",2),
        #                                      rep("italic",12))),
        axis.title = element_text(size = 16),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16),
        plot.caption = element_text(size = 14),
        strip.background = element_rect(fill="gray60",color = "black"),
        strip.text = element_text(size=15, colour="white")) +
  xlab("Coef.") +
  ggtitle("") +
  labs(title = "T3 MS sucralose results and differential species",
       subtitle = "MaAsLin results",
       caption = "*No adjustment")
invisible(dev.off())

# HALLA correlation
write.table(SU_mtg@otu_table,"Data/MS_Sweetener/FAMT3_OTU_sel_alllev.txt", sep = "\t")
# Cont 5,10,13,14,18
a <- data.frame(t((SU_mtg@sam_data)[,c(5,10,13,14,15)]))
colnames(a) <- colnames(SU_mtg@otu_table)
write.table(a,"Data/MS_Sweetener/FAMT3_Meta_sel.txt", sep = "\t")

FA_Bac_CorMS <- read.table("Results/FA_MS_Sweetener/Corr_PW_alllev_pea/all_associations.txt",
                          sep = "\t", header = TRUE)

FA_Bac_CorMS <- FA_Bac_CorMS %>% #filter((association)>0.5) %>% 
  filter(q.values<0.05)

FA_Bac_CorMS$Y_features[grepl("\\|",FA_Bac_CorMS$Y_features)] <- as.character(lapply(strsplit(FA_Bac_CorMS$Y_features[grepl("\\__",FA_Bac_CorMS$Y_features)],"\\|"),
                                                                            function(x) paste0(x)[length(x)]))                                                                                 

FA_Bac_CorMS <- FA_Bac_CorMS[order(FA_Bac_CorMS$X_features,FA_Bac_CorMS$association,decreasing = T),]     
FA_Bac_CorMS <- FA_Bac_CorMS %>% filter(grepl("g__|s__",Y_features))

# CcorA correlation

test <- CCorA(data.frame(SU_mtg@sam_data)[,c(5,10,13,14,15)],
              t(data.frame(SU_mtg@otu_table)),
              permutations = 999)
library(ggcor)


# Pathway----
FA_MS_PW <- FA_PW 
colnames(FA_MS_PW) <- str_replace(colnames(FA_MS_PW),"X","")
colnames(FA_MS_PW) <- str_remove_all(colnames(FA_MS_PW),"\\_1")
FA_MS_PW <- FA_MS_PW %>% select(colnames(FA_MS_PW)[colnames(FA_MS_PW) %in%
                                         colnames(SU_mtg@otu_table)])

FA_MS_PW <- FA_MS_PW %>% mutate(sum0 = as.numeric(apply(FA_MS_PW,1,sum))) %>% 
  filter(sum0>0) %>% select(-sum0)
# write.table(FA_MS_PW,"Data/MS_Sweetener/FAMT3_PW_alllev.txt", sep = "\t")
# 
# FA_Bac_CorMSPW <- read.table("Results/FA_MS_Sweetener/Corr_PW_alllev_pea/all_associations.txt",
#                            sep = "\t", header = TRUE)
# 
# FA_Bac_CorMSPW <- FA_Bac_CorMSPW %>% #filter((association)>0.5) %>% 
#   filter(q.values<0.05)
# FA_Bac_CorMSPW$PW <- as.character(lapply(FA_Bac_CorMSPW$Y_features,
#                                          function(x)
#                                            strsplit(x,"\\|")[[1]][1]))
# FA_Bac_CorMSPW$PW_Bac <- as.character(lapply(FA_Bac_CorMSPW$Y_features,
#                                          function(x)
#                                            strsplit(x,"\\|")[[1]][2]))
# FA_Bac_CorMSPW$PW_Bac[grepl("\t",FA_Bac_CorMSPW$PW_Bac)] <-  as.character(lapply(FA_Bac_CorMSPW$PW_Bac[grepl("\t",FA_Bac_CorMSPW$PW_Bac)],
#                                                                                  function(x)
#                                                                                    strsplit(x,"\\\t0")[[1]][1]))
set.seed(100)
for(i in colnames(SU_mtg@sam_data)[grepl("T3",colnames(SU_mtg@sam_data))&
                                   grepl("\\.",colnames(SU_mtg@sam_data))]){
    Maaslin2(
      input_data = data.frame(FA_MS_PW),
      input_metadata = data.frame(SU_mtg@sam_data) ,
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = c(i),
      # random_effects = c(),
      output = paste0("Results/FA_MS_Sweetener/",i,"MaAsLin_PW"),
      correction = "BH",
      # reference = "Normal",
      standardize = FALSE,
      cores = 1)
  
}


Maas_MS_PW <- data.frame()
for(i in colnames(SU_mtg@sam_data)[grepl("T3",colnames(SU_mtg@sam_data))&
                                   grepl("\\.",colnames(SU_mtg@sam_data))]){
  temp_db <- read.table(paste0("Results/FA_MS_Sweetener/",i,"MaAsLin_PW/significant_results.tsv"),
                        header = T)
  Maas_MS_PW <- rbind(Maas_MS_PW,temp_db)
}

Maas_MS_PW <- Maas_MS_PW %>% mutate(Bac_source = as.character(lapply(strsplit(feature,"\\__"),
                                                            function(x) paste0(x)[length(x)])),
                                    PW = as.character(lapply(strsplit(feature,"\\__"),
                                                             function(x) paste0(x)[1])))

Maas_MS_PW$Bac_source[grepl("unclassified", Maas_MS_PW$Bac_source)] <- "Unclassified" 

Maas_MS_PW$Ass_Group <- ifelse(Maas_MS_PW$coef<0, "Neg", "Pos")
FA_abutest <- humann2meco(abund_table = abund_file_path, db = "MetaCyc")

a <- data.frame(environment(FA_abutest[["tidy_dataset"]])[["self"]][["tax_table"]])
FA_PW$PW <- rownames(FA_PW)
FA_PW <- data.frame(t(FA_PW))
FA_PW$Fun <- a$Superclass2[match(FA_PW$PW,rownames(a))]
# Maas_MS_PW <- Maas_MS_PW[order(Maas_MS_PW$coef,decreasing = T),]
# Maas_MS_PW$Order <- length(Maas_MS_PW$feature):1
Maas_MS_PW$Fun <- FA_PW$Fun[match(Maas_MS_PW$feature,rownames(FA_PW))]
Maas_MS_PW$PW_short <- paste0(
  lapply(Maas_MS_PW$feature,function(x) strsplit(x,"\\.")[[1]][1]),"-",
  lapply(Maas_MS_PW$feature,function(x) strsplit(x,"\\.")[[1]][2]))
a$PW_short <- as.character(lapply(rownames(a),function(x) strsplit(x,"\\:")[[1]][1]))
Maas_MS_PW$Fun[is.na(Maas_MS_PW$Fun)] <- a$Superclass2[match(Maas_MS_PW$PW_short[is.na(Maas_MS_PW$Fun)],
                                                         a$PW_short)]
Maas_MS_PW$PW <- FA_PW$PW[match(Maas_MS_PW$feature,rownames(FA_PW))]
Maas_MS_PW$PW  <- as.character(lapply(Maas_MS_PW$PW,function(x) strsplit(x,"\\:")[[1]][2]))
Maas_MS_PW$PW  <- as.character(lapply(Maas_MS_PW$PW,function(x) strsplit(x,"\\|")[[1]][1]))

Maas_MS_PW_plt <- reshape(Maas_MS_PW[,c(4,10,11)],
                idvar = "PW", timevar = "Bac_source", direction = "wide")

Maas_MS_PW_plt$Fun <- Maas_MS_PW$Fun[match(Maas_MS_PW_plt$PW,
                                           Maas_MS_PW$PW)]
rownames(Maas_MS_PW_plt) <- Maas_MS_PW_plt$PW
Maas_MS_PW_plt <- Maas_MS_PW_plt[,-1]
colnames(Maas_MS_PW_plt) <- str_remove(colnames(Maas_MS_PW_plt),"coef.")
col_fun = circlize::colorRamp2(c(-5, 0, 4), c("#4575b4", "white", "#d73027"))
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  foo = anno_block(gp = gpar(fill = c(
    "#00A087FF", "#F3987FFF")),
    labels = c("Aspartame (Urine)","Sucralose (CB)"),
    labels_gp = gpar(col = "white", fontsize = 14.5))
)
lgd <- Legend(col_fun = col_fun, direction = "horizontal",
              title = "Coef. (MaAsLin)",
              at = c(-5, 0, 4), border = "black",
              legend_height = unit(1.8, "cm"))

MaAsHM <- Heatmap(Maas_MS_PW_plt[,-c(6)], na_col = "gray95",
        rect_gp = gpar(col = "gray"),
        border = T,
        col = col_fun,
        column_split = c(rep(1, 4),2),
        top_annotation = ha, 
        row_split = Maas_MS_PW_plt[, 6],
        row_title_rot = 0,
        column_title = NULL,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = F,
        height = unit(20, "cm"))

pdf("Results/FA_MS_Sweetener/MaAs_PW_in_Bac_10072022.pdf", height = 12, width = 14)
draw(MaAsHM, heatmap_legend_list = list(lgd), heatmap_legend_side = "bottom")
invisible(dev.off())

pdf("Results/FA_MS_Sweetener/MaAs_PW_in_Bad.pdf", height = 15, width = 20)

ggplot(Maas_MS_PW_plt, aes(y = Bac, x = Pathway, 
                           fill = Coef)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-5,4),na.value = "white") +
  geom_tile(color = "black") +
  # theme_classic() +
  coord_fixed()  +
  scale_x_discrete(dup_axis(labels = Maas_MS_PW_plt$Fun)) +
  theme(title = element_text(size = 16),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        # axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.caption = element_text(size = 10, face = "italic"),
        plot.subtitle = element_text(size = 12)) 

invisible(dev.off())





# Sweetener and VT----
SU_mtg_PN <- data.frame(read_xlsx("Data/FA_Mtg/MOMmy_Abun_OTU_11072022.xlsx"))
SU_mtg_PNtaxa <- data.frame(matrix(unlist(lapply(strsplit(SU_mtg_PN$clade_name,"\\|"), function(x)
  if(length(x)<7){
    c(paste0(x),rep(NA,7-length(x)))
  }
  else{paste0(x)})), 
  ncol = 7, byrow = TRUE))
colnames(SU_mtg_PNtaxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(SU_mtg_PN) <- SU_mtg_PN$clade_name 
rownames(SU_mtg_PNtaxa) <- rownames(SU_mtg_PN)


SU_mtg_PN <- SU_mtg_PN %>% select(colnames(SU_mtg_PN)[colnames(SU_mtg_PN) %in% 
                                                        rownames(M1_Meta)])

SU_mtg_PNMeta <- M1_Meta %>% filter(Mother_Sub_ID %in% Urine_Sweet$Sub_ID) %>%
  filter(Time_point!="Mother T3")

rownames(SU_mtg_PNMeta) <- SU_mtg_PNMeta$SampleID
SU_mtg_PNMeta <- data.frame(t(SU_mtg_PNMeta))
SU_mtg_PN <- SU_mtg_PN %>% select(colnames(SU_mtg_PN)[colnames(SU_mtg_PN) %in% 
                                                        rownames(SU_mtg_PNMeta)])

SU_mtg_PN <- phyloseq(otu_table(as.matrix(SU_mtg_PN), taxa_are_rows = TRUE), 
                      sample_data(SU_mtg_PNMeta), 
                      phyloseq::tax_table(as.matrix(SU_mtg_PNtaxa)))

SU_mtg_PN <- prune_taxa(taxa_sums(SU_mtg_PN) > 0, SU_mtg_PN)
SU_mtg_PN <- prune_taxa(grepl("Bac",tax_table(SU_mtg_PN)[, "Phylum"]), SU_mtg_PN)
# 1. Alpha
SU_PN_mtgalpha <- cbind(SU_mtg_PN@sam_data,
                     Richness_spe = t(estimateR(ceiling(t(SU_mtg_PN@otu_table))))[,1],
                     Richness_Chao1 = t(estimateR(ceiling(t(SU_mtg_PN@otu_table))))[,2],
                     Evenness = diversity(t(SU_mtg_PN@otu_table)) / log(specnumber(t(SU_mtg_PN@otu_table))),
                     Shannon = diversity(t(SU_mtg_PN@otu_table), index = "shannon"),
                     InvSimp = diversity(t(SU_mtg_PN@otu_table), index = "invsimpson")) 
# FAalpha[,c(62:74)] <- apply(FAalpha[,c(62:74)], 2, function(x)
#   ifelse(x=="High","High","Normal"))
# my.comp <- list( c("Low","High"),
#                  c("Low","Normal"),
#                 c("Normal","High"))
plot_list = list()
plot_list_con = list()
for (i in colnames(SU_PN_mtgalpha)[grepl("T3",colnames(SU_PN_mtgalpha))&
                                grepl("\\.",colnames(SU_PN_mtgalpha))]) {
  for(j in c(83:87)){
    if(is.numeric(SU_PN_mtgalpha[,i])){
      my.data <- data.frame(x = SU_PN_mtgalpha[,i], 
                            y = SU_PN_mtgalpha[,j],
                            lab = colnames(SU_PN_mtgalpha)[j])
      a <- cor.test(my.data$x,
                    my.data$y,method = "spearman")
      x.lab <- str_remove(i,"T3_")
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_point(size = 1,  aes(color = x)) +
        stat_smooth(method = "glm",formula = y~x,  
                    se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
        theme_bw() +      
        xlab(paste(x.lab)) +
        theme(axis.title.y = element_blank(),
              axis.text = element_text(color = "black", size =13),
              axis.title.x = element_text(size = 16, face = "bold"),
              legend.position = "none",
              panel.grid.major = element_line(colour = "gray70"),
              strip.background = element_rect(fill="gray40",color = "black"),
              strip.text = element_text(size=15, colour="white")) +
        facet_grid(. ~ lab) +
        annotate("text", x = max(my.data$x, na.rm = T)*.8,
                 y = max(my.data$y,na.rm = T), size = 4, 
                 label = paste0("R =", round(a$estimate,2),
                                "\np-value = ", round(a$p.value,3)))
      
      plot_list_con <- c(plot_list_con, list(p))
    } else{ 
      my.data <- data.frame(x = SU_PN_mtgalpha[,i], 
                            y = SU_PN_mtgalpha[,j],
                            lab = colnames(SU_PN_mtgalpha)[j])
      x.lab <- str_replace(i,"_cut"," Group")
      x.lab <- str_remove(i,"T3_")
      p <- ggplot(my.data,  aes(x = x, y = y)) +
        geom_violin(aes(color = x), size = 1,  position = "dodge") +
        geom_point(size = 1,  aes(color = x)) +
        geom_boxplot(width = 0.15, fill="white",aes(color = x)) +
        scale_color_manual(values=c("#BC3C29FF","#00A087FF")) +
        scale_fill_manual(values=c("#BC3C29FF","#00A087FF")) +
        theme_bw() +      
        stat_compare_means(method = "wilcox") +
        stat_signif(label = "p.signif", 
                    size = 0,
                    map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2),
                    margin_top = 0.05,
                    method = "wilcox") +
        xlab(paste(x.lab)) +
        theme(axis.title = element_blank(),
              axis.text = element_text(color = "black", size =13),
              # axis.title.x = element_text(size = 16, face = "bold"),
              legend.position = "none",
              panel.grid.major = element_line(colour = "gray70")) +
        facet_grid(. ~ lab) +
        theme(strip.background = element_rect(fill="gray40",color = "black"),
              strip.text = element_text(size=15, colour="white"))
      
      plot_list <- c(plot_list, list(p))
    }
  }
}
# + 
#   plot_annotation(title = 'T3 Saccharin (Cord blood)',
#                   theme = theme(plot.title = element_text(size = 18)))
pdf("Results/FA_MS_Sweetener/MS_PN_Sac_CB_alpha_sum_11072022.pdf", height = 3, width = 13)
wrap_plots(plot_list,ncol = 5) 
invisible(dev.off())

pdf("Results/FA_MS_Sweetener/MS_PN_con_alpha_sum_11072022.pdf", height = 25, width = 16)
wrap_plots(plot_list_con,ncol = 5)
invisible(dev.off())


# 2. Beta
set.seed(100)
SU_PNnmdspre <- vegdist(t(SU_mtg_PN@otu_table),na.rm = T)
SU_PNnmds <- monoMDS(SU_PNnmdspre)

SU_PNnmds <- data.frame(cbind(data.frame(SU_mtg_PN@sam_data),
                           SU_PNnmds[["points"]]))

SU_PNnmds$MDS1 <- as.numeric(SU_PNnmds$MDS1)
SU_PNnmds$MDS2 <- as.numeric(SU_PNnmds$MDS2)
# gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
#   filter(abs(MDS2)<0.1)

SU_Beta_plt <- ggplot(SUnmds, aes(MDS1, MDS2)) +
  geom_point(size = 4, aes(fill = T3_Urine.Aspartame),
             shape = 21, color = "gray20") +
  # scale_color_manual(name = "CRN intake group",
  #                    values=c("#BC3C29FF","#00A087FF")) +
  scale_fill_gradient2(
    name = "T3 Aspartame (Urine)",
    low = "white",
    mid = "#fec44f",
    high = "#993404",
    midpoint = 0.55,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("") +
  annotate("text",label = paste("PERMANOVA R2 = 0.024\np = 0.011"),
           x = -1.3, y = 1.8, hjust = 0)

pdf("Results/FA_MS_Sweetener/MS_beta_U_Aspartame_11072022.pdf", height = 6, width = 7)
ggMarginal(SU_Beta_plt, colour = "#cc4c02")
invisible(dev.off())


set.seed(1)
Uni_per <- data.frame()
for(i in colnames(SU_mtg_PNMeta)[grepl("T3",colnames(SU_mtg_PNMeta))&
                           grepl("\\.",colnames(SU_mtg_PNMeta))]){
  if(length(unique(SU_mtg_PNMeta[,i]))>1){
    temp_db <- SU_mtg_PNMeta %>% filter(complete.cases(SU_mtg_PNMeta[,i]))
    temp_db <-  subset_samples(SU_mtg_PN, Sub_ID %in% temp_db$Sub_ID)
    temp_dis <- phyloseq::distance(temp_db, 
                                   method="bray", weighted=T)
    temp_db <-  adonis2(as.formula(paste0("temp_dis ~ ", i)), 
                        data = data.frame(temp_db@sam_data),
                        permutation=999)
    Uni_per <- rbind(Uni_per,data.frame(temp_db))
  } else {
    print(i)
  }
}


Sig_Uni_per <- Uni_per %>% filter(Pr..F.<0.1)

set.seed(1)
multiper_db <- subset_samples(FAmtgraw, Sub_ID %in% data.frame(FAMT3Meta %>% 
                                                                 select(c("Sub_ID",colnames(FAMT3Meta)[colnames(FAMT3Meta) %in% rownames(Sig_Uni_per)])) %>% 
                                                                 drop_na %>% 
                                                                 select(Sub_ID))$Sub_ID)
multiper_dis <- distance(multiper_db, method="bray", weighted=T)
multi_per <- adonis(as.formula(paste0("multiper_dis ~ ", 
                                      paste0(rownames(Sig_Uni_per),collapse = " + "))), 
                    data = data.frame(multiper_db@sam_data),
                    permutation=999)





# 3. differential species
set.seed(100)
for(i in colnames(SU_mtg_PN@sam_data)[grepl("T3",colnames(SU_mtg_PN@sam_data))&
                                   grepl("\\.",colnames(SU_mtg_PN@sam_data))]){
  if(length(unique(data.frame(SU_mtg_PN@sam_data)[,i]))==3){
    Maaslin2(
      input_data = SU_mtg_PN@otu_table,
      input_metadata = data.frame(SU_mtg_PN@sam_data)[data.frame(SU_mtg_PN@sam_data)[,i]!="Normal",],
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = c(i),
      # random_effects = c(),
      output = paste0("Results/FA_MS_Sweetener/",i,"_PN_MaAsLin"),
      correction = "BH",
      reference = "Normal",
      standardize = FALSE,
      cores = 1)
  }else{
    Maaslin2(
      input_data = SU_mtg_PN@otu_table,
      input_metadata = data.frame(SU_mtg_PN@sam_data) ,
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = c(i),
      # random_effects = c(),
      output = paste0("Results/FA_MS_Sweetener/",i,"_PN_MaAsLin"),
      correction = "BH",
      reference = "Normal",
      standardize = FALSE,
      cores = 1)
  }
}

Maas_MS_PN <- data.frame()
for(i in colnames(SU_mtg_PN@sam_data)[grepl("T3",colnames(SU_mtg_PN@sam_data))&
                                   grepl("\\.",colnames(SU_mtg_PN@sam_data))]){
  temp_db <- read.table(paste0("Results/FA_MS_Sweetener/",i,"_PN_MaAsLin/significant_results.tsv"),
                        header = T)
  Maas_MS_PN <- rbind(Maas_MS_PN,temp_db)
}

Maas_MS_PN <- Maas_MS_PN %>% mutate(feature = as.character(lapply(strsplit(feature[grepl("\\.",feature)],"\\."),
                                                            function(x) paste0(x)[length(x)]))) 
Maas_MS_PN$Ass_Group <- ifelse(Maas_MS_PN$coef<0, "Neg", "Pos")
# Maas_MS$Var_fea <- paste0(Maas_MS$metadata,"-\n",Maas_MS$feature)
# Maas_Clinc$Level <- ifelse(grepl("\\_",Maas_Clinc$feature), "Species", "Genus")
Maas_MS_PN$Sweet <- as.character(
  lapply(Maas_MS_PN$metadata,function(x) strsplit(x,"\\.")[[1]][2])
)
Maas_MS_PN$Sweet <- str_remove(Maas_MS_PN$Sweet,"_cut")
Maas_MS_PN$Source <- str_remove(as.character(
  lapply(Maas_MS_PN$metadata,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
Maas_MS_PN <- Maas_MS_PN %>% filter(coef>1)
Maas_MS_PN$Source <- str_replace(Maas_MS_PN$Source,"_"," ")
Maas_MS_PN <- Maas_MS_PN[order(Maas_MS_PN$coef,decreasing = T),]
Maas_MS_PN$Order <- length(Maas_MS_PN$feature):1
Maas_MS_PN$feature <- reorder(Maas_MS_PN$feature,Maas_MS_PN$Order)
Maas_MS_PN$Coef_tx <- ifelse(Maas_MS_PN$coef<1, Maas_MS_PN$coef+1, Maas_MS_PN$coef-1)

pdf("Results/FA_MS_Sweetener/MaAsLin_Sum_08072022.pdf", height = 7, width = 10)
ggplot(data = Maas_MS_PN, 
       aes(x = coef, y = feature, 
           fill = Ass_Group)) +
  geom_col(width = 0.75, size = 2) + 
  geom_label(aes(label = Source, 
                 x = Coef_tx,
                 size=4, alpha = 0.6),
             fill = "white") +
  scale_x_continuous(expand = c(0,0.5),
                     breaks = c(seq(-8,4,4))) +
  scale_fill_manual(values=c("#d73027"),
                    name = "Association",
                    labels = c("Negative", "Positive")) +
  geom_vline(xintercept = 0, size = 0.8, linetype="dashed",
             color = "gray40") +
  facet_grid(. ~ Sweet, scales = "free", space = "free") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 14, face = "italic"),
        axis.title = element_text(size = 16),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16),
        plot.caption = element_text(size = 14),
        strip.background = element_rect(fill="gray60",color = "black"),
        strip.text = element_text(size=15, colour="white")) +
  xlab("Coef.") +
  ggtitle("") +
  labs(title = "T3 MS sweetener results and postnatal differential species",
       subtitle = "MaAsLin results",
       caption = "*No adjustment")


invisible(dev.off())

Urine_SweetT3_pair <- left_join(Urine_Sweet[,c(1:16)],Early_paired_OTU_uni,"Sub_ID") %>% 
  filter(complete.cases(T3_Urine.Aspartame))


plot_list = list()
for(j in colnames(Urine_SweetT3_pair)[c(17:36)]){
  for (i in colnames(Urine_SweetT3_pair)[grepl("T3",colnames(Urine_SweetT3_pair))&
                                         !grepl("FA_ques",colnames(Urine_SweetT3_pair))]){
    if(is.numeric(Urine_SweetT3_pair[,i])){
      my.data <- data.frame(y = Urine_SweetT3_pair[,i], 
                            x = as.character(Urine_SweetT3_pair[,j]),
                            lab = j) %>% 
        filter(complete.cases(x))
      a <- wilcox.test(y ~ x, my.data)
      if(a$p.value<0.06){
        y.lab = str_remove(i,"T3_")
        x.lab = str_remove(j,"_1")
        p <- ggplot(my.data,  aes(x = x, y = y)) +
          geom_point(size = 1,  aes(color = x)) +
          geom_boxplot(width = 0.15, fill="white",aes(color = x)) +
          scale_color_manual(values=c("gray60","#238443")) +
          scale_fill_manual(values=c("gray60","#238443")) +
          theme_classic() +      
          stat_compare_means(method = "wilcox") +
          stat_signif(label = "p.signif", 
                      size = 0,
                      map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2),
                      margin_top = 0.05,
                      method = "wilcox") +
          ggtitle(x.lab) +
          ylab(paste0(y.lab)) +
          scale_x_discrete(labels = c(paste0("Not share","\n(N = ", sum(my.data$x==0),")"),
                                      paste0("Share","\n(N = ", sum(my.data$x==1),")"))) +
          theme(title = element_text(size = 14, face = "bold.italic"),
                axis.title.y = element_text(size = 14, face = "plain"),
                axis.text = element_text(color = "black", size =13),
                axis.title.x = element_blank(),
                legend.position = "none") 
        plot_list <- c(plot_list, list(p))
      }
        } else{
          print(paste0(paste0(c(i,j),collapse = "-"),": Not avail!"))
        }
      }
} 

+

CEdb_fisher <- data.frame()
for(j in colnames(Urine_SweetT3_pair)[c(14:33)]){
  for (i in colnames(Urine_SweetT3_pair)[grepl("T3",colnames(Urine_SweetT3_pair))]){
    if(any(Urine_SweetT3_pair[,i]=="ND")){
      if(all(Urine_SweetT3_pair[,i]=="ND")){
        next
      }else{
        
        a <- fisher.test(Urine_SweetT3_pair[,i],
                         Urine_SweetT3_pair[,j])
        if(a$p.value<0.05){
          a <- data.frame(spe = j,
                          sweet = i,
                          p_val = a$p.value,
                          OR = a$estimate,
                          OR_low = a$conf.int[1],
                          OR_up = a$conf.int[2])
          CEdb_fisher <- rbind(CEdb_fisher,a)
        } else{
          print(paste0(paste0(c(i,j),collapse = "-"),": Insig!"))
        }
      }
    } 
  }
}

# Combine GDM and Sweeteners results ----
GDM_Diff <- data.frame(read_xlsx("Results/Lefse_diff_analysis_GDM.xlsx"))
GDM_Diff$Ass_Group <- ifelse(GDM_Diff$GDM_Group=="GDM","Pos","Neg")
GDM_Diff$Conv_LDAS <- ifelse(GDM_Diff$GDM_Group=="GDM",GDM_Diff$LDAmean,
                             0-GDM_Diff$LDAmean)


Maas_MS <- data.frame()
for(i in colnames(SU_mtg@sam_data)[grepl("T3",colnames(SU_mtg@sam_data))&
                                   grepl("\\.",colnames(SU_mtg@sam_data))]){
  temp_db <- read.table(paste0("Results/FA_MS_Sweetener/",i,"MaAsLin/all_results.tsv"),
                        header = T)
  Maas_MS <- rbind(Maas_MS,temp_db)
}

Maas_MS <- Maas_MS %>% mutate(feature = as.character(lapply(strsplit(feature[grepl("\\.",feature)],"\\."),
                                                            function(x) paste0(x)[length(x)]))) 
Maas_MS$Ass_Group <- ifelse(Maas_MS$coef<0, "Neg", "Pos")
Maas_MS <- Maas_MS %>% filter(feature %in% GDM_Diff$f)
Maas_MS$Sweet <- as.character(
  lapply(Maas_MS$metadata,function(x) strsplit(x,"\\.")[[1]][2])
)
Maas_MS$Source <- str_remove(as.character(
  lapply(Maas_MS$metadata,function(x) strsplit(x,"\\.")[[1]][1])
),"T3_")
Maas_MS$Source <- str_replace(Maas_MS$Source,"_"," ")

Maas_MS_plt <- reshape(Maas_MS[,c(1,2,4)], idvar = "feature",timevar = "metadata", direction = "wide")
Maas_MS_plt_p <- reshape(Maas_MS[,c(1,2,9)], idvar = "feature",
                         timevar = "metadata", direction = "wide")
Maas_MS_plt_p[,-1][Maas_MS_plt_p[,-1]>0.05] <- NA
Maas_MS_plt_p[,-1][is.na(Maas_MS_plt_p[,-1])] <- 1
Maas_MS_plt <- Maas_MS_plt[,-c(4:7)]
Maas_MS_plt_p <- Maas_MS_plt_p[,-c(4:7)]

Maas_MS_plt$GDM <- ifelse(GDM_Diff$Ass_Group[match(Maas_MS_plt$feature,GDM_Diff$f)]=="Pos",
                          1,0)
Maas_MS_plt <- Maas_MS_plt[order(Maas_MS_plt$GDM,decreasing = T),]
colnames(Maas_MS_plt) <- str_remove(colnames(Maas_MS_plt),"coef.T3_")
Maas_MS_plt_p$GDM <- ifelse(GDM_Diff$Ass_Group[match(Maas_MS_plt_p$feature,GDM_Diff$f)]=="Pos",
                          1,0)
Maas_MS_plt_p <- Maas_MS_plt_p[order(Maas_MS_plt_p$GDM,decreasing = T),]
Maas_MS_plt <- Maas_MS_plt[,c(1,11,2,6,3,7,9,8,10,5,4)]
Maas_MS_plt_p <- Maas_MS_plt_p[,c(1,11,2,6,3,7,9,8,10,5,4)]
Maas_MS_plt$feature <- str_remove(Maas_MS_plt$feature,"\\_")

Maas_MS_colab <- c("Urine","Urine_bi","Maternal_blood","Urine_bi",
                   "Maternal_blood","Cord_blood","Urine", "Maternal_blood","Cord_blood")
col_fun = circlize::colorRamp2(c(-8, 0, 6.5), c("#2166ac", "white", "#b2182b"))
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  foo = anno_block(gp = gpar(fill = c(
    "#00A087FF", "#3C5488FF", "#F3987FFF")),
    labels = c("Aspartame","Saccharin","Sucralose"),
    labels_gp = gpar(col = "white", fontsize = 14.5))
)
lgd = Legend(col_fun = col_fun, direction = "horizontal",
             title = "Coef. (MaAsLin)",
              at = c(-8, 0, 7), border = "black",
              legend_height = unit(1.8, "cm"))

Maas_MS_HM <- Heatmap(Maas_MS_plt[,-c(1,2)], 
        na_col = "gray90",
        rect_gp = gpar(col = "gray30"),
        border = T,
        col = col_fun,
        # row_split = PDH_GBA[1:3, 2],
        row_title_rot = 0,
        row_labels = Maas_MS_plt[,1],
        row_split = Maas_MS_plt$GDM,
        right_annotation = rowAnnotation(
          foo = anno_block(gp = gpar(fill = c("#74add1", "#d73027")),
                           labels = c("GDM depleted", "GDM enriched"), 
                           labels_gp = gpar(col = "white", fontsize = 13.5))),
        row_title = NULL,
        # row_names_gp = gpar(fontsize = 10,
        #                     col = c("#74a9cf","#74a9cf","#045a8d")),
        show_column_names = FALSE, 
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(Maas_MS_colab, rot = 0,
                           gp = gpar(col = c(rep("#fe9929",2),"#ef3b2c",
                                               rep(c("#fe9929","#ef3b2c","#a50f15"),2))),
                           offset = unit(0.92, "npc"), 
                           just = "center"),
          annotation_height = max_text_width(Maas_MS_colab)
        ),
        column_gap = unit(2, "points"),
        column_split = c(rep(1:3, each = 3)),
        top_annotation = ha, 
        column_title = NULL,
        cluster_rows = FALSE, cluster_columns = FALSE,
        height = unit(16, "cm"),
        show_heatmap_legend = F,
        cell_fun =  function(j, i, x, y, w, h, fill) {
          if(Maas_MS_plt_p[,-c(1,2)][i, j] < 0.05) {
            grid.text("*", x, y)
          } else if(Maas_MS_plt_p[,-c(1,2)][i, j] < 0.1) {
            grid.text(".", x, y)
          }}) 

pdf("Results/FA_MS_Sweetener/MaAsLin_GDM_09072022.pdf", height = 9, width = 14)
draw(Maas_MS_HM,
     heatmap_legend_list = list(lgd), heatmap_legend_side = "bottom")
invisible(dev.off())

# Try media
library(mediation)
SU_Med <- cbind(SU_mtg@sam_data,
              t(data.frame(SU_mtg@otu_table)[grepl("Parabacteroides",
                                                  rownames(data.frame(SU_mtg@otu_table))),]))
colnames(SU_Med)[c(84:91)] <- as.character(lapply(colnames(SU_Med)[c(84:91)],
       function(x) str_split(x,"\\|")[[1]][length(str_split(x,"\\|")[[1]])]) )

# [81] "T3_Cord_blood.Saccharin_cut"    
# [82] "T3_Maternal_blood.Saccharin_cut"

test.med.gdm <- glm(NInfec~T3_Cord_blood.Saccharin_cut,SU_Med,
                    family = "binomial")
test.out <- glm(s__Parabacteroides_merdae~NInfec+T3_Cord_blood.Saccharin_cut,SU_Med,
                family = "gaussian")

med.test<- mediation::mediate(test.med.gdm, test.out, 
                   treat = "T3_Cord_blood.Saccharin_cut",
                   mediator ="NInfec",
                  robustSE = TRUE, sims =100)###treat填自变量，mediator填中介变量






