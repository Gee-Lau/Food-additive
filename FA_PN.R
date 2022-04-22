library(dplyr)

PNdb <- MOMcomb %>% dplyr::filter(Sex == "F")
PNdb <- PNdb %>% dplyr::filter(!Sub_ID %in% laterecruite)
PNoc <- readRDS("Data/PNdb_07092021.rds")
PNdb$ID <- substr(PNdb$Sub_ID,7,9)
PNoc$ID <- substr(PNoc$Sub_ID_Child,7,9)
# PNdb$Weight_child[PNdb$Weight_child<2000] <- NA
PNdb <- merge(PNdb,PNoc,"ID")

PNdb$Weight_child[grepl("lb",PNdb$Weight_unit_child)] <- PNdb$Weight_child[grepl("lb",PNdb$Weight_unit_child)]*453.6
PNdb$Preterm <- ifelse(PNdb$Preterm=="2Yes",1,0)
PNdb$Medication <- ifelse(PNdb$Medication=="2Yes",1,0)
PNdb$`SICU/NICU` <- ifelse(PNdb$`SICU/NICU`=="有",1,0)
PNdb$BMI_Group <- factor(as.factor(PNdb$BMI_Group),
                           levels(as.factor(PNdb$BMI_Group))[c(2,1,3,4)])
PNdb$GComp <- ifelse(PNdb$GComp=="有",1,0)
PNdb$LBW <- ifelse(PNdb$Weight_child<2500,1,0)
PNdb$HBW <- ifelse(PNdb$Weight_child>4000,1,0)
PNdb$NComp <- ifelse(PNdb$NComp=="有",1,0)
PNdb$GBS[PNdb$GBS=="不知道"] <- NA 
PNdb$GBS <- ifelse(PNdb$GBS=="是",1,0)
# PNdb$Weight_child[PNdb$Weight_child<2000|PNdb$Weight_child>5000] <- NA
PNdb[,c(65:78,80:83,104:108)] <- apply(PNdb[,c(65:78,80:83,104:108)], 2, 
                                       function(x) ifelse(grepl("No",x),0,1))

# PNdb <- left_join(PNdb,MOMcomb_sub[,c(1,139:164)],"Sub_ID")

# 8:10,16:19,13,14,72:83,87:89,91:103,106:109
uni_c <- data.frame()
for (i in colnames(PNdb)[c(89:91,98)]){
  if (all(PNdb[,i]==0)|all(PNdb[,i][complete.cases(PNdb[,i])]=="1No")){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("Weight_child~",i)),
                  PNdb[PNdb[,i]!="Low",],family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                    confint(temp_a,level = 0.95),
                    summary(temp_a)$coefficients[,c(4)])
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c[!grepl("Intercep",rownames(uni_c)),]
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
colnames(uni_c) <- c("Est","Lower","Upper","p.value")
uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
# uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))

uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
uni_c$Var <- rownames(uni_c)

for (i in colnames(PNdb)[c(89:91,98)]){
  levels(PNdb[,i]) <- levels(as.factor(PNdb[,i]))[c(3,2,1)]
}

uni_c <- data.frame()
for (i in colnames(PNdb)[c(89:91)]){
  if (all(PNdb[,i]==0)|all(PNdb[,i][complete.cases(PNdb[,i])]=="1No")){
    print(i)
  }
  else {
    print(i)
    for (j in colnames(PNdb)[c(128,131,139,144,151)]){
      temp_a <- glm(as.formula(paste0(j,"~",i)),
                    PNdb[PNdb[,i]!="Low",],family ="binomial")
      temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                      exp(cbind(OR = coef(temp_a), 
                                confint.default(temp_a,level = 0.95))),
                      Outcome = j)
      uni_c <- rbind(uni_c,temp_b)
    }
  }
}
uni_c <- uni_c[!grepl("Intercep",rownames(uni_c)),]
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
uni_c[,c(1:5)] <- apply(uni_c[,c(1:5)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)
# uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.06)
uni_c$Var <- rownames(uni_c)
uni_c$Var <- substr(uni_c$Var,0,3)
uni_c$Var[11:15] <- "CRN"
uni_c <- reshape(uni_c, idvar = "Outcome", timevar = "Var", direction = "wide")
writexl::write_xlsx(uni_c,"PN_outcomes.xlsx")

mglm <- glm(`SICU/NICU`~Age+BMI+Any_disease, PNdb,family = "binomial")
mglmtb <-  cbind(summary(mglm)$coefficients[,c(1,4)],
                 exp(cbind(OR = coef(mglm), 
                           confint.default(mglm,level = 0.95))))
mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)[-1,]



MOMcomb_sub <- MOMcomb %>% dplyr::filter(!Sub_ID %in% laterecruite)
set.seed(0)
uni_c <- data.frame()
for (i in colnames(MOMcomb_sub)[c(7,8,9,118,15:18,12,13,22:59,109:117,119:122,90:102,105:108)]){
  if (all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="M"]==0)|all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub[,i])]=="1No")){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("GComp~",i)),
                  MOMcomb_sub[MOMcomb_sub$Sex=="F",],family ="binomial")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                    confint(temp_a,level = 0.95),
                    summary(temp_a)$coefficients[,c(4)])
    mglmtb_sub <- rbind(mglmtb_sub,temp_b)
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c[!grepl("Intercep",rownames(uni_c)),]
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
# colnames(uni_c) <- c("Est","Lower","Upper","p.value")
uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
uni_c$Pr...z.. <- ifelse(uni_c$Pr...z..<0.001, paste("<0.001"),paste(uni_c$Pr...z..))

uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.06)
uni_c$Var <- rownames(uni_c)

mglmtb <- data.frame()
for (i in colnames(PNdb)[c(89:91,98)]){
  print(i)
  temp_a <- glm(as.formula(paste0("Weight_child ~ Age + GA_Week + BMI +", i)),
                PNdb[PNdb[,i]!="Low",], family ="gaussian")
  temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                  confint(temp_a,level = 0.95),
                  summary(temp_a)$coefficients[,c(4)])
  mglmtb <- rbind(mglmtb,temp_b)
}


mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)
colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                             paste(mglmtb$p.value))

mglmtb <- mglmtb[!grepl("Intercep",rownames(mglmtb)),]
mglmtb$Var <- rownames(mglmtb)
writexl::write_xlsx(mglmtb,"Multi_glm_Emu_BW_14112021.xlsx")

mglmtb <- data.frame()
for (i in colnames(MOMcomb_sub)[c(109:117,119:122)]){
  print(i)
  temp_a <- glm(as.formula(paste0("as.numeric(OVW)~ Age + Alco + PE_Group + Dairy_1_5y +
                                  CarboHfood_1_5y + SoftDrink_5_10y + ProcFruit_1_5y +", i)),
                MOMcomb_sub[MOMcomb_sub$Sex=="M",], family ="binomial")
  temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                  exp(cbind(OR = coef(temp_a), 
                            confint.default(temp_a,level = 0.95))))
  mglmtb <- rbind(mglmtb,temp_b)
}

mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)
# colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$Pr...z.. <- ifelse(mglmtb$Pr...z..<0.001, paste("<0.001"),
                          paste(mglmtb$Pr...z..))
# mglmtb <- mglmtb %>% filter(Pr...z..=="<0.001"|Pr...z..<0.06)
mglmtb <- mglmtb %>% dplyr::filter((!grepl("Intercep",rownames(mglmtb)))&grepl(paste0(colnames(MOMcomb_sub)[c(109:117,119:122)],
                                                                                      collapse = "|"),
                                                                               rownames(mglmtb)))


# Draw table
mglmtb <- mglmtb %>% filter(grepl("High",rownames(mglmtb)))
mglmtb <- mglmtb[-4,]
rownames(mglmtb) <- c("Polysorbate-80", "Carboxymethylcellulose",
                                        "Carrageenan")
colnames(mglmtb) <- c("Estimate","95% CI (lower)",
                      "95%CI (upper)","p-value")

mglmtab <- ggtexttable(mglmtb, rows = rownames(mglmtb),
            theme = ttheme("blank", base_size = 18,
                           padding = unit(c(1,1), "cm"),
                           rownames.style = rownames_style(face = "bold",
                                                           hjust=0, x=0.1,
                                                           size = 18))) %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 5, linetype = 1) %>%
  tab_add_hline(at.row = c(4), row.side = "bottom", linewidth = 5, linetype = 1) %>% 
  tab_add_footnote(text = "*Linear regression, adjsuted for gestational week, materanal age, and pBMI. 
  Norlmal emulsifiers intake group as reference group compare to high intake group.",
                face = "plain", size = 14, hjust = 1)
  

# mglmtb$Var <- str_replace_all(mglmtb$Var,"log2","")
# mglmtb$p.bin <- ifelse(mglmtb$p.value>0.05,"0","1")
# mglmtb$ID <- length(mglmtb$Est):1
# mglmtb$Var <- reorder(mglmtb$Var,mglmtb$ID)

# Birth weight vs FA groups
ggplot(PNdb[complete.cases(PNdb$Sweet_cut),], 
       aes(x=Sweet_cut, y=Weight_child, color = Sweet_cut)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#2c7bb6", "#fdae61","#d7191c")) +
  scale_x_discrete(label = c("Low \n(N = 108)",
                             "Normal \n(N = 228)",
                             "High \n(N = 116)")) +
  theme_classic() +
  stat_compare_means(method = "kruskal",label.x = 2.5,label.y = 63) +
  stat_compare_means(comparisons = statecomp,label = "p.signif",
                     method = "wilcox") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 16),
        axis.title = element_text(size = 20)) +
  xlab("Mother annual P80 intake groups") +
  ylab("Birth weight (g)")


colnames(PNdb)[c(89:91,98)] <- c("polysorbate-80", "carboxymethylcellulose",
                                        "carrageenan", "Overall emulsifier")
colnames(PNdb)[c(89:91,98)] <- c("p80", "carboxymethylcellulose",
                                 "carrageenan", "Overall emulsifier")

# PNdb[,c(89:91,98)] <- apply(PNdb[,c(89:91,98)],2,
#                                    function(x) ifelse(x=="High","2High","1Normal"))
set.seed(100)
plot_list = list()
for (i in c(89:91)) {
  if(i==89){
  my.data <- data.frame(x = PNdb[,i][complete.cases(PNdb[,i])&PNdb[,i]!="Low"],
                        y = PNdb$Weight_child[complete.cases(PNdb[,i])&PNdb[,i]!="Low"])
  j <- colnames(PNdb)[i]
  p <-  ggplot(my.data,  aes(x = x, y = y, color = x)) +
    # geom_violin() +
    geom_boxplot(width = 0.6, size = 1.2,
                 position = position_dodge(0.9)) +
    scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="Normal"&complete.cases(my.data$y))),
                               paste0("High \nN = ",sum(my.data$x=="High"&complete.cases(my.data$y))))) +
    geom_point(shape=16, position=position_jitter(0.3),
               size = 2, alpha = 0.3) +
    scale_color_manual(values=c( "#FFA319FF","#800000FF","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "wilcox", size = 5,
                       label.x = 1,label.y = 1950) +
    # stat_compare_means(comparisons = statecomp,label = "p.signif",
    #                    method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 15)) +
    xlab(paste0("Annual ",j," intake group")) +
    ylab("Children birth weight (g)") 
  plot_list <- c(plot_list, list(p))}
  else{
    my.data <- data.frame(x = PNdb[,i][complete.cases(PNdb[,i])&PNdb[,i]!="Low"],
                          y = PNdb$Weight_child[complete.cases(PNdb[,i])&PNdb[,i]!="Low"])
    j <- colnames(PNdb)[i]
    p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
      # geom_violin() +
      geom_boxplot(width = 0.6, size = 1.2,
                   position = position_dodge(0.9)) +
      scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="Normal"&complete.cases(my.data$y))),
                                 paste0("High \nN = ",sum(my.data$x=="High"&complete.cases(my.data$y))))) +
      geom_point(shape=16, position=position_jitter(0.3),
                 size = 2, alpha = 0.3) +
      scale_color_manual(values=c( "#FFA319FF","#800000FF","#e66101")) +
      theme_classic()+
      stat_compare_means(method = "wilcox", size = 5,
                         label.x = 1,label.y = 1950) +
      # stat_compare_means(comparisons = statecomp,label = "p.signif",
      #                    method = "wilcox",hide.ns = T) +
      theme(legend.position = "none",
            axis.text = element_text(color = "black",size = 12),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title = element_text(size = 16)) +
      xlab(paste0("Annual ",j," intake group")) +
      ylab("Children birth weight (g)") 
    plot_list <- c(plot_list, list(p))}
  }





PNdb$Fa_Sac_cut <- MOMcomb_sub$SAC_cut3[MOMcomb_sub$Sex=="M"][match(substr(PNdb$Sub_ID,6,9),
                                                                   substr(MOMcomb_sub$Sub_ID,6,9))]
PNdb$Fa_Sac_cut <- ifelse(PNdb$Fa_Sac_cut=="High","High","Not high")


ggplot(PNdb[complete.cases(PNdb$Fa_Sac_cut),],  aes(x = Fa_Sac_cut, y = Weight_child, color = Fa_Sac_cut)) +
  geom_violin() +
  geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
  scale_color_manual(values=c("#5e3c99", "#fdb863","#e66101")) +
  geom_hline(aes(yintercept = 2500), 
             color="black", size = 0.5, linetype = "dashed") +
  theme_classic()+
  stat_compare_means(method = "anova",label.x = 1,label.y = 1800) +
  # stat_compare_means(comparisons = statecomp,label = "p.signif",
  #                    method = "t.test",hide.ns = T) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 10),
        axis.title = element_text(size = 12)) +
  xlab("Fathers' saccharin intake") +
  ylab("Children birth weight (g)") 







