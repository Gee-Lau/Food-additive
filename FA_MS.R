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
library("patchwork")

# Correlation between MSS results and FA ques----

Urine_Sweet <- data.frame(read_xlsx("Data/MOMmy_Sweetener_data.xlsx",sheet = "CB"))
Urine_Sweet$Calculated.Amt..ng.ml. <- as.numeric(Urine_Sweet$Calculated.Amt..ng.ml.)
Urine_Sweet$FA_ques <- MomFA$SAC[match(Urine_Sweet$ID,MomFA$Sub_ID)]
Urine_Sweet$FA_ques[Urine_Sweet$Compound=="Aspartame"] <- MomFA$ASP[match(Urine_Sweet$ID[Urine_Sweet$Compound=="Aspartame"],
                                                                          MomFA$Sub_ID)]
Urine_Sweet$FA_ques[Urine_Sweet$Compound=="Sucralose"] <- MomFA$SUC[match(Urine_Sweet$ID[Urine_Sweet$Compound=="Sucralose"],
                                                                          MomFA$Sub_ID)]

plot_list = list()
for (i in unique(Urine_Sweet$Compound)) {
  my.data <- Urine_Sweet %>% filter(Compound==i)
  colnames(my.data)[c(4,5)] <- c("x","y")
  if(all(is.na(my.data$x)))
    next
  a <- cor.test(my.data$x,
                my.data$y,method = "spearman")
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_point(size = 2, color = "#542788") +
      theme_classic() +
      xlab(paste(i," in urine")) +
      ylab(paste(i," from questionnaire")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
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
      xlab(paste(i," in cord blood")) +
      ylab(paste(i," from questionnaire")) +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25, color = "#e08214") +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black")) + 
      annotate("text", x = max(my.data$x,na.rm = T)*3/4,
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3)))
    
    plot_list <- c(plot_list, list(p))
    
  }
}

pdf(file="CB_Ques_sweetener.pdf", width = 12, height = 6)
wrap_plots(plot_list,ncol = 3)
invisible(dev.off())


# Check the sweetener level in different emulsifier intake groups
Samp_list <- data.frame(read_xlsx("Data/FA_Samptest_pilot_list.xlsx"))
Urine_Sweet$Emu_group <- Samp_list$Group[match(Urine_Sweet$ID,Samp_list$STUDY_NO)]

plot_list = list()
for (i in unique(Urine_Sweet$Compound)) {
  my.data <- Urine_Sweet %>% filter(Compound==i)
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
Urine_Sweet <- left_join(Urine_Sweet, MomFA[,c(1,8:10,72:74)],"Sub_ID")
Urine_Sweet <- left_join(Urine_Sweet, MOMcomb_sub[,c(1,13,16,22:62,101,102)],"Sub_ID")
Urine_Sweet <- left_join(Urine_Sweet, PNdb[,c(2,65:88,103:109,118:126,128:134,139,144,150:153)],"Sub_ID")

Urine_Sweet[,c(97:102)] <- apply(Urine_Sweet[,c(97:102)], 2,
                                 function(x) ifelse(is.na(x), NA,
                                                    ifelse(grepl("No", x), 0, 1)))
Urine_Sweet[,c(95,103)] <- apply(Urine_Sweet[,c(95,103)], 2,
                                 function(x) ifelse(is.na(x), NA,
                                                    ifelse(x=="否", 0, 1)))
Urine_Sweet$Del_mode.y <- ifelse(is.na(Urine_Sweet$Del_mode.y), NA,
                                 ifelse(Urine_Sweet$Del_mode.y=="剖腹產",1,0))
colnames(Urine_Sweet)[104] <- "NICU"
Urine_Sweet$GA <- PNdb$GA_Week[match(Urine_Sweet$Sub_ID,PNdb$Sub_ID)]
# Urine_Sweet$Asp_bi <- ifelse(is.na(Urine_Sweet$Aspartame),"0","1")
# con outcomes
uni_c <- data.frame()
for (i in colnames(Urine_Sweet)[c(2:4)]){
  for (j in colnames(Urine_Sweet[c(80:87,95:104,106:108)])) {
    if(all(Urine_Sweet[,j][complete.cases(Urine_Sweet[,j])]==0)){
      print(j)
    }else{
      temp_a <- glm(as.formula(paste0(j,"~",i)),
                    Urine_Sweet,family ="binomial")
      if(temp_a$df.residual==0){
        print(j)
      }else{
        temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                            exp(cbind(OR = coef(temp_a), 
                                      confint.default(temp_a,level = 0.95))),
                            Outcome = j)
      uni_c <- rbind(uni_c,temp_b)}
      
    }
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c[,c(1:5)] <- apply(uni_c[,c(1:5)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)

# uni_c <- uni_c %>% filter(Pr...t..=="<0.001"|Pr...t..<0.05)
# uni_c <- uni_c %>% filter(Pr...z..=="<0.001"|Pr...z..<0.05)

uni_c$Var <- substr(rownames(uni_c),0,3)
uni_c <- reshape(uni_c, idvar = "Outcome", timevar = "Var", direction = "wide")

writexl::write_xlsx(uni_c, "MB_sweetener_outcomes.xlsx")

mglmtb_sub <- data.frame()
for (i in colnames(Urine_Sweet)[c(2:4)]){
  print(i)
  temp_a <- glm(as.formula(paste0("LBW ~ Pre_BMI + Age + SexOB + Preterm + Any_disease +", i)),
                Urine_Sweet, family ="gaussian")
  temp_b <-  cbind(summary(temp_a)$coefficients[,c(1,4)],
                   exp(cbind(OR = coef(temp_a), 
                             confint.default(temp_a,level = 0.95))))
  mglmtb_sub <- rbind(mglmtb_sub,temp_b)
}
mglmtb_sub <- apply(mglmtb_sub, 2, function(x) round(x,4))
mglmtb_sub <- data.frame(mglmtb_sub)
writexl::write_xlsx(mglmtb_sub,"Urine_sweetener_LBW_Mul.xlsx")


uni_c <- data.frame()
for (i in colnames(Urine_Sweet)[c(2:4)]){
  for (j in colnames(Urine_Sweet[c(105,78,79,93,94)])) {
    if(all(Urine_Sweet[,j][complete.cases(Urine_Sweet[,j])]==0)){
      print(j)
    }else{
      temp_a <- glm(as.formula(paste0(j ," ~ ", i)),
                    Urine_Sweet, family ="gaussian")
      temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                      confint(temp_a,level = 0.95),
                      summary(temp_a)$coefficients[,c(4)],
                      Outcome = j)
      uni_c <- rbind(uni_c,temp_b)
    }
  }
}

uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c[,c(1:4)] <- apply(uni_c[,c(1:4)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)
uni_c$Var <- substr(rownames(uni_c),0,3)
uni_c <- reshape(uni_c, idvar = "Outcome", timevar = "Var", direction = "wide")
writexl::write_xlsx(uni_c, "Urine_sweetener_Conoutcomes.xlsx")



# bin outcomes (GComp, GBS, `SICU/NICU`, GDM)
# Boundary: GDM-urine ASP
uni_c <- data.frame()
for (i in colnames(Urine_Sweet)[c(2:4)]){
  if (all(Urine_Sweet[,i]==0)){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("GDM~",i)),
                  Urine_Sweet,family ="binomial")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                    confint(temp_a,level = 0.95),
                    summary(temp_a)$coefficients[,c(4)])
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


mglmtb_sub <- data.frame()
for (i in colnames(Urine_Sweet)[c(2:4)]){
  print(i)
  temp_a <- glm(as.formula(paste0("GDM ~ Pre_BMI + Age +", i)),
                Urine_Sweet, family ="gaussian")
  temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                  confint(temp_a,level = 0.95),
                  summary(temp_a)$coefficients[,c(4)])
  mglmtb_sub <- rbind(mglmtb_sub,temp_b)
}

mglmtb_sub <- apply(mglmtb_sub, 2, function(x) round(x,4))
mglmtb_sub <- data.frame(mglmtb_sub)
colnames(mglmtb_sub) <- c("Est","Lower","Upper","p.value")
mglmtb_sub$p.value <- ifelse(mglmtb_sub$p.value<0.001, paste("<0.001"),
                             paste(mglmtb_sub$p.value))



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

