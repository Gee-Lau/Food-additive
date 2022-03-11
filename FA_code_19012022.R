
library(tidyverse)
library(dplyr)
library(patchwork)
library(reshape2)
library(ggpubr)
library(readxl)
library(stringr)
library(ggplot2)
library(ape)
library(hagis)
library(tableone)

# Univariable (BMI as an example)----


uni_c <- data.frame()
for (i in colnames(MOMcomb_sub)[c(3,6:8,13,15:19)]){
  if (all(MOMcomb_sub[,i]==0)){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~",i)),
                  MOMcomb_sub[MOMcomb_sub$Sex=="F",],
                  family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                    confint(temp_a,level = 0.95),
                    summary(temp_a)$coefficients[,c(4)])
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
colnames(uni_c) <- c("Est","Lower","Upper","p.value")
uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
# select p<0.05
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))
uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)

uni_c$Var <- rownames(uni_c)


# Multi-variable (BMI as an example)---- 
mglm <- glm(BMI~Age+VM_group+Any_dis+Income_fam + Avg_takeout +
              SAC+Emu_log,
            momall,family = "gaussian")
mglmtb <- cbind(summary(mglm)$coefficients[,c(1)],confint(mglm),
                summary(mglm)$coefficients[,c(4)])
mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)
colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))
mglmtb <- mglmtb %>% filter(p.value=="<0.001"|p.value<0.06)
mglmtb <- mglmtb %>% filter(!grepl("Intercep",rownames(mglmtb)))

mglmtb$Var <- rownames(mglmtb)

# Plot multi-v outcome
BMI_mul <- mglmtb
BMI_mul_plt <-  ggplot(BMI_mul, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#54278f") +  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+
  scale_y_discrete(name="", labels=c("VM_groupMoverV"="Meat > vegetable",
                                     "AlcoYes" = "Alcohole drinking",
                                     "SAC"="Saccharin")) +
  ylab("")+
  xlab("Estimate (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("BMI")

# Nutrients and emulsifier contribution----
NTdb <- left_join(NTdb,MomFA[,c(1:10,12)], by = "Sub_ID")
NTdb <- left_join(NTdb, MOMcomb_sub[,c(1,88:90,97)], by = "Sub_ID")
NTdb <- NTdb %>% filter(complete.cases(`Overall emulsifier`)) %>% 
  filter(substr(Sub_ID,4,4)==1)
Emulmtx <- NTdb %>% select(cho_1_d2,cho_2_d2,cho_3_d2,soups_d2,meats_d2,veggies_d2,
                           flavoured_milk_d2,cream_d2,baked_goods_d2,alcoholic_d2,
                           supplements_d2,egg_derived_d2,sweets_d2,creamer_d2)
Emulmtx <- apply(Emulmtx,2,as.numeric)
Emulmtx.pcoa <-  vegdist(t(Emulmtx), "jaccard", na.rm = TRUE)
Emulmtx.pcoa <- betadisper(Emulmtx.pcoa, as.factor(groups))
plot(Emulmtx.pcoa, hull = T, ellipse = T)

pcoa <- cmdscale(Emulmtx.pcoa)
efit <- envfit(pcoa, data.frame(colnames(Emulmtx)))
plot(pcoa, col = c("#BC3C2999", "#0072B599","#E1872799","#20854E99",
                   "#7876B199","#6F99AD99"), pch = c(21,24))
# groups <- c("Carbroh","Carbroh","Carbroh","Processed_fastfood","Processed_fastfood",
#             "Processed_fastfood","Processed_dairy","Processed_dairy",
#             "Processed_fastfood","Alcoholic","Supp","Egg","Supp","Supp")

plot(efit, col = "black", cex = 0.7)
biplot(prcomp(Emulmtx,scale. = T), scale = 0)
pcoa <- data.frame(pcoa)
pcoa$Group <- groups
ggplot(data.frame(pcoa), aes(x = X1, y = X2)) +
  geom_point(size = 1, alpha = .6, aes(color = Group)) +
  scale_color_manual(values = c("#BC3C2999", "#0072B599","#E1872799","#20854E99",
                                "#7876B199","#6F99AD99"),
                     name = "Food intake group") +
  stat_ellipse(aes(color = Group)) +
  ggtitle("A") +
  theme(axis.text = element_text(color = "black", size = 12),
        title = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "plain"),
        axis.title.x = element_text(size = 14, face = "plain"),
        legend.position = "bottom")
+
  geom_segment(data = fit_spp, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               col = "black")

Emulmtx.pcoa <- betadisper(Emulmtx.pcoa, as.factor(NTdb$`Overall emulsifier`))


Emulmtx <- NTdb %>% select(cho_1_d2,cho_2_d2,cho_3_d2,soups_d2,meats_d2,veggies_d2,
                           flavoured_milk_d2,cream_d2,baked_goods_d2,alcoholic_d2,
                           supplements_d2,egg_derived_d2,sweets_d2,creamer_d2,P80.x)
Emulmtx <- apply(Emulmtx, 2, as.numeric)
Emul_sum <- data.frame(Var = colnames(Emulmtx),
                       Mean = colMeans(Emulmtx))
rownames(Emul_sum) <- 1:length(Emul_sum[,1])
Emul_sum$Conv_mean <- NA
Emul_sum$Conv_mean[c(1,3:6,7,11)] <- Emul_sum$Mean[c(1,3:6,7,11)]*5000
Emul_sum$Conv_mean[c(2,7:9,13)] <- Emul_sum$Mean[c(2,7:9,13)]*3000
Emul_sum$Conv_mean[c(10)] <- Emul_sum$Mean[c(10)]*120
Emul_sum$Conv_mean[c(12)] <- Emul_sum$Mean[c(12)]*10
Emul_sum$Conv_mean[c(14)] <- Emul_sum$Mean[c(14)]*4000
Emul_sum$Per <- round(Emul_sum$Conv_mean/Emul_sum$Mean[15]*100,1)
Emul_sum <- Emul_sum[-15,]

CMCmtx <- NTdb %>% select(cho_1_d2,meats_d2,seafood_d2,veggies_d2,fruits_d2,
                          fm_d2,cream_d2,tea_d2,coffee_d2,CMC)
CMCmtx <- apply(CMCmtx, 2, as.numeric)
CMC_sum <- data.frame(Var = colnames(CMCmtx),
                      Mean = colMeans(CMCmtx))
rownames(CMC_sum) <- 1:length(CMC_sum[,1])
CMC_sum$Conv_mean <- NA
CMC_sum$Conv_mean[c(1)] <- CMC_sum$Mean[c(1)]*15000
CMC_sum$Conv_mean[c(2)] <- CMC_sum$Mean[c(2)]*5000
CMC_sum$Conv_mean[c(3)] <- CMC_sum$Mean[c(3)]*40000
CMC_sum$Conv_mean[c(4:5)] <- CMC_sum$Mean[c(4:5)]*20000
CMC_sum$Conv_mean[c(6:7)] <- CMC_sum$Mean[c(6:7)]*5000
CMC_sum$Conv_mean[c(8:9)] <- CMC_sum$Mean[c(8:9)]*4000
CMC_sum$Per <- round(CMC_sum$Conv_mean/CMC_sum$Mean[10]*100,1)
CMC_sum <- CMC_sum[-10,]


CRNmtx <- NTdb %>% select(cho_1_d2,meats_d2,seafood_d2,fruits_d2,flavoured_milk_d2,
                          fm_d2,cream_d2,sugars_d2,baked_goods_d2,egg_derived_d2,supplements_d2,CRN)
CRNmtx <- apply(CRNmtx, 2, as.numeric)
CRN_sum <- data.frame(Var = colnames(CRNmtx),
                      Mean = colMeans(CRNmtx))
rownames(CRN_sum) <- 1:length(CRN_sum[,1])
CRN_sum$Conv_mean <- NA
CRN_sum$Conv_mean[c(1,9)] <- CRN_sum$Mean[c(1,9)]*7000
CRN_sum$Conv_mean[c(2,3)] <- CRN_sum$Mean[c(2,3)]*10000
CRN_sum$Conv_mean[c(4)] <- CRN_sum$Mean[c(4)]*20000
CRN_sum$Conv_mean[c(5)] <- CRN_sum$Mean[c(5)]*2000
CRN_sum$Conv_mean[c(6:8,10)] <- CRN_sum$Mean[c(6:8,10)]*5000
CRN_sum$Conv_mean[c(11)] <- CRN_sum$Mean[c(11)]*6000
CRN_sum$Per <- round(CRN_sum$Conv_mean/CRN_sum$Mean[12]*100,1)
CRN_sum <- CRN_sum[-12,]
CRN_sum$Type <- "CRN"
CMC_sum$Type <- "CMC"
Emul_sum$Type <- "P80"
Emul_sum <- rbind(Emul_sum,CMC_sum,CRN_sum)

pdf(file="Emu_NT_12012022.pdf", width = 10, height = 10)
ggplot(Emul_sum) +
  geom_col(aes(x = Type, y = Per, fill = Var), color = "black") +
  theme_classic() +
  ylab("Food sources (%)") +
  scale_fill_manual(values = c("#d6604d", "#fddbc7","#fdae61","#fee090",
                               "#ffffbf","#8c510a","#fcc5c0","#fa9fb5",
                               "#ec7014","#f7fbff","#969696","#74c476",
                               "#a50f15","#cc4c02","#ec7014","#4292c6",
                               "#ce1256", "#e7298a","#a63603","#006d2c"),
                    name = "Food sources") +
  theme(axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "bottom")
invisible(dev.off())

# Demo between different emulsifier intake groups----
MOMcomb_sub$Pre_BMI_Group <- "NA"
MOMcomb_sub$Pre_BMI_Group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Pre_BMI<18.5] <- "Underweight"
MOMcomb_sub$Pre_BMI_Group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Pre_BMI>18.4&MOMcomb_sub$Pre_BMI<23] <- "Normal"
MOMcomb_sub$Pre_BMI_Group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Pre_BMI>22.9&MOMcomb_sub$Pre_BMI<27.1] <- "Overweight"
MOMcomb_sub$Pre_BMI_Group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Pre_BMI>27] <- "Obese"

MOMcomb_sub$Age_group <- "NA"
MOMcomb_sub$Age_group[complete.cases(MOMcomb_sub$Age)&MOMcomb_sub$Age<30] <- "<30"
MOMcomb_sub$Age_group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Age>29.9&MOMcomb_sub$Age<35] <- "30-35"
MOMcomb_sub$Age_group[complete.cases(MOMcomb_sub$Pre_BMI)&MOMcomb_sub$Age>34.9] <- ">=35"


tabcat <- colnames(MOMcomb_sub)[c(102, 7, 8, 13, 15, 17, 18, 79:80, 101)]
tabcon <- colnames(MOMcomb_sub)[c(9,16,84)]
colnames(MOMcomb_sub)[97] <- "OE_Group"
for (i in c(88:90)){
  MOMcomb_sub[,i] <- as.character(MOMcomb_sub[,i])
}
MOMcomb_sub$OE_Group <- as.character(MOMcomb_sub$OE_Group)
demotab <- CreateTableOne(vars = c(tabcat,tabcon), 
                          data = MOMcomb_sub[MOMcomb_sub$Sex=="F"&MOMcomb_sub$carrageenan!="Low",], 
                          strata = "carrageenan")
a <- data.frame(print(demotab, nonnormal = tabcon, exact = "stage", smd = T)) 
a$Var <- rownames(a)
writexl::write_xlsx(a,"Table1_cha_1.xlsx")

# Emulsifier in different maternal characteristics
MOMcomb_sub$P80_con <- as.numeric(MomFA$P80[match(MOMcomb_sub$Sub_ID,MomFA$Sub_ID)])
MOMcomb_sub$CMC_con <- MomFA$CMC[match(MOMcomb_sub$Sub_ID,MomFA$Sub_ID)]
MOMcomb_sub$CRN_con <- MomFA$CRN[match(MOMcomb_sub$Sub_ID,MomFA$Sub_ID)]
MOMcomb_sub$Emu_con <- MomFA$Emu[match(MOMcomb_sub$Sub_ID,MomFA$Sub_ID)]

ANOVA_by_group <- function(var1){
  temp_db <- data.frame(apply(MOMcomb_sub[,c(103:106)],2,function(x)
    data.frame(y = MOMcomb_sub[[var1]],x,
               Var = var1) %>%
      filter(y!="NA"&complete.cases(x)) %>% 
      mutate(p.value = round(kruskal.test(x~y)$p.value,4)) %>%
      group_by(y) %>% 
      summarise(Median = median(x),
                Q1 = summary(x)[[2]],
                Q3 = summary(x)[[5]],
                Var = rep(unique(Var),length(Q1)),
                p.val = rep(unique(p.value),length(Q1)))))
  return(temp_db)
}

ANOVA_sum <- data.frame()
for(i in colnames(MOMcomb_sub)[c(7,8,13,15,17,18,80,81,101,102)]){
  a <- ANOVA_by_group(i)
  ANOVA_sum <- rbind(ANOVA_sum,a)
}

writexl::write_xlsx(ANOVA_sum,"Emul_by_cha.xlsx")

