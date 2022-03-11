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
# Do the similarity plot of FA between sex, dropped----
gender.ds <- vegdist(log2(MomFA[,c(13:22)]),na.rm = T)

gender.md <- monoMDS(gender.ds)
gender.md <- data.frame(cbind(Sex = MomFA$Sex,
                              Fam = substr(MomFA$Sub_ID,7,9),
                              gender.md$points))
gender.md$MDS1 <- as.numeric(gender.md$MDS1)
gender.md$MDS2 <- as.numeric(gender.md$MDS2)
gender.md <- gender.md %>% filter(abs(MDS1)<2) %>% 
  filter(abs(MDS2)<2)

ggsave("FA_MDS_Sex_05082021.eps",width = 6, height = 4.5,units = "in")
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


a <- MomFA %>% filter(grepl("00038",Sub_ID))
                      
set.seed(5)
gender.anosim <- anosim(MomFA[,c(2:10)], 
                        MomFA$Sex, permutations = 9999)
autoplot(gender.anosim, notch = FALSE) +
  ylab("Dissimilarity rank")  +
  ggtitle("Bray-Curtis ANOSIM") +
  theme_classic()

kms = kmeans(MomFA[,3], centers=3,nstart = 1029)
  Position=factor(kms$cluster)
  
# Try HLM----
library(lme4)
library(nlme)
library(merTools)
library(lmerTest)
  
HLMdb <- rbind(combdb[,c(1, 6:8, 14, 19, 38:45, 49, 54)],
               dadcombdb[,c(1, 4, 6 ,7 ,11, 12, 28:33, 35, 36, 34, 8)])
# HLMdb <- HLMdb %>% filter(Group!="Normal")
HLMdb$Diet_type <- as.character(apply(HLMdb[,c(8:12)],1,function(x) paste(colnames(HLMdb)[c(8:12)][x=="Yes"],collapse = ",")))
HLMdb$VM_group <- as.character(lapply(as.list(HLMdb$Diet_type),function(x) strsplit(x,",")[[1]][1]))
HLMdb$VM_group[HLMdb$VM_group=="SpeDie_group"] <- NA

HLMdb$Fam_ID <- substr(HLMdb$Sub_ID,7,9)
HLMdb$FoodA <- ifelse(HLMdb$FoodA=="有","Yes",
                      ifelse(is.na(HLMdb$FoodA),NA,"No"))
HLMdb$Fam_ID <- factor(HLMdb$Fam_ID)
HLMdb <- HLMdb %>% filter(Sub_ID!="PWH100631")

md1 <- gls(Weight~1, data = HLMdb[,-c(1,3)], 
    method = "ML", na.action = "na.omit")
summary(md1)

md2 <- lme(Weight~1, data = HLMdb[,-c(1,3)], 
    method = "ML", na.action = "na.omit", random = ~1|Fam_ID)

summary(md2)
anova(md1,md2)

md3 <- lme(Weight ~ Alco, data = HLMdb[-c(1,3)], method = "ML", 
           na.action = "na.omit", random = ~1|Fam_ID)
summary(md3)
anova(md2,md3)


md4 <-lmer(Weight ~ Group + (1|Fam_ID),
           data = HLMdb[,-c(1,3)], na.action = na.exclude)
summary(md4)
VarCorr(md4)

ICC(outcome = "Weight", group = "Group", data = HLMdb)

HLMdb$Group_bi <- ifelse(HLMdb$Group=="High",1,0)
summary(glm(Group_bi ~ Weight+Age+Alco+Avg_meal+Avg_takeout+VM_group, 
            data = HLMdb[substr(HLMdb$Sub_ID,4,4)==1,], family = "binomial"))
wilcox.test(combdb$Weight[combdb$Group != "Normal"]~combdb$Group[combdb$Group != "Normal"])

test <- HLMdb %>% filter(Alco=="Yes"&substr(HLMdb$Sub_ID,4,4)==1)
summary(glm(Group_bi ~ Weight, 
            data = test, family = "binomial"))

# 1. Early dietary habit and FA----
# 1.1 Distribution of FA early approach 
PFtype <- colnames(PDHdb)[c(9:40)]
PFtype <- as.character(lapply(as.list(PFtype),function(x) strsplit(x,"_")[[1]][1]))
PFtype <- PFtype[!duplicated(PFtype)]


for (i in PFtype){
  PDHdb[,i] <- as.character(apply(PDHdb[,colnames(PDHdb)[grepl(i,colnames(PDHdb))]],
                                  1, function(x)  paste(colnames(PDHdb)[grepl(i,colnames(PDHdb))][x==1][1],collapse = ",")))
  PDHdb[,i] <- as.character(lapply(as.list(PDHdb[,i]),function(x) strsplit(x,"_")[[1]][2]))
  PDHdb[,i][PDHdb[,i]=="NA"] <- NA
}

PDHsum <- gather(PDHdb[,c(44:51)])

PDHsum <- PDHsum%>% group_by(key,value) %>% 
  count()
PDHsum$value[is.na(PDHsum$value)] <- "Not applicable"
colnames(PDHsum) <- c("Category", "Timepoint", "N")
PDHsum$Timepoint <- factor(as.factor(PDHsum$Timepoint),
                           levels(as.factor(PDHsum$Timepoint))[c(5,3,1,4,2)])
PDHsum$Perc <- round(PDHsum$N/426*100,1)

ggsave("PDH_Time_28062021.eps",width = 8, height = 8,units = "in")
ggplot(data = PDHsum, aes( x = Category, y = Perc, 
                           fill = Timepoint)) +
  geom_bar(color = "grey35", stat="identity") +
  theme_classic() +
  scale_x_discrete(name="test", labels=c("CarboHfood"="Carbohydrate \nfoods",
                                         "Dairy"="Processed \ndairy",
                                         "Fastfood"="Fast food",
                                         "ProcFruit"="Processed \nfruit",
                                         "ProcVeg"="Processed \nvegetable",
                                         "Seafood"="Processed \nseafood", 
                                         "Snacks", "SoftDrink"="Soft drink")) +
  scale_fill_manual(name = "Timepoint of first \nfrequently use", 
                    values=c("#f2f0f7",  "#cbc9e2", "#9e9ac8",
                             "#756bb1", "#54278f")) +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 35,hjust = 0.9, vjust = 1),
        legend.position = "bottom",
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 14))+
  ylab("%") 
invisible(dev.off())


# 1.2 MND plot to see the similarity
for (i in c(2,4:7,9:40)){
  PDHdb[,i] <- ifelse(PDHdb[,i] == "Yes", 1,0)
}

PDHdb$Infancy_food_source <- ifelse(PDHdb$Infancy_food_source=="From shops",1,0)
PDHdb$Adolo_food_source <- ifelse(PDHdb$Adolo_food_source=="Preserved",1,0)
dadPDHdb$Infancy_food_source <- ifelse(dadPDHdb$Infancy_food_source==1,1,0)
dadPDHdb$Adolo_food_source <- ifelse(dadPDHdb$Adolo_food_source==2,1,0)

PDHdb <- PDHdb %>% filter(Sub_ID!="PWH100631")
PDHdb$Group <- combdb$Group[match(PDHdb$Sub_ID,combdb$Sub_ID)]
PDHdb$BMI_group <- combdb$BMI_Group[match(PDHdb$Sub_ID,combdb$Sub_ID)]
PDHdb <- PDHdb %>% filter(Complete=="Complete")

PDHdb$Group[is.na(PDHdb$Group)] <- "Normal"

PDH.ds <- rbind(PDHdb[,c(1:40)],dadPDHdb[,c(1:40)])
PDH.ds$Sex <- MomFA$Sex[match(PDH.ds$Sub_ID,MomFA$Sub_ID)]
PDH.ds$Group <- MomFA_W$Group[match(PDH.ds$Sub_ID,MomFA_W$Sub_ID)]
PDH.ds <- PDH.ds %>% filter(complete.cases(Group))
PDH.ds1 <- vegdist(PDH.ds[,c(3:7,9:40)], na.rm = T)

PDH.md <- monoMDS(PDH.ds1)
PDH.md <- data.frame(cbind(Group = PDH.ds$Group,
                           Sex = PDH.ds$Sex,
                           BMI_group =  PDH.ds$BMI_group,
                              PDH.md$points))

PDH.md$MDS1 <- as.numeric(PDH.md$MDS1)
PDH.md$MDS2 <- as.numeric(PDH.md$MDS2)
PDH.md <- PDH.md %>% filter(abs(MDS1)<5&abs(MDS2<5))
# PDH.md$BMI_group[PDH.md$BMI_group=="NA"] <- NA

PDH.md_labels <- data.frame(Group = NA,
                            Sex = c("Female", "Male"), 
                            label = c("Bray-Curtis ANOSIM: \nR2=0.029, p = 0.063",
                                      "Bray-Curtis ANOSIM: \nR2=0.057, p = 0.031"))
PDH.md$Sex <- ifelse(PDH.md$Sex=="F","Female","Male")

ggsave("PDH_MDS_15072021.eps",width = 8, height = 6,units = "in")
ggplot(PDH.md, aes(MDS1, MDS2, col = Group)) +
  geom_point() + 
  scale_color_manual(name = "Food additive group",
                    values=c("#d7191c",  "#2c7bb6", "#fdae61")) +
  
  stat_ellipse(aes(color=Group, group=Group),linetype = 2, size = 1.5) +
  facet_grid(.~Sex) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(color = "black", size = 16)) +
  geom_text(x = -0.015, y = -0.02, color = "black", 
            aes(label = label),
            show.legend = NA,
            data = PDH.md_labels)
  # ggtitle("MDS of annual estimated food additive") +
  # annotate("text", x = -0.052, y = 0.0047, size = 4, #family = "Arial",
  #          label = paste0("Bray-Curtis ANOSIM: \nR2=0.052, p = 0.0072"))
invisible(dev.off())  




set.seed(5)

anosim1 <- function (x, grouping, permutations = 999, distance = "bray", 
                     strata = NULL, parallel = getOption("mc.cores")) 
{
  EPS <- sqrt(.Machine$double.eps)
  if (!inherits(x, "dist")) {
    if ((is.matrix(x) || is.data.frame(x)) && isSymmetric(unname(as.matrix(x)))) {
      x <- as.dist(x)
      attr(x, "method") <- "user supplied square matrix"
    }
    else x <- vegdist(x, method = distance)
  }
  if (any(x[complete.cases(x)] < -sqrt(.Machine$double.eps))) 
    warning("some dissimilarities are negative - is this intentional?")
  sol <- c(call = match.call())
  grouping <- as.factor(grouping)
  if (length(grouping) != attr(x, "Size")) 
    stop(gettextf("dissimilarities have %d observations, but grouping has %d", 
                  attr(x, "Size"), length(grouping)))
  if (length(levels(grouping)) < 2) 
    stop("there should be more than one class level")
  matched <- function(irow, icol, grouping) {
    grouping[irow] == grouping[icol]
  }
  x.rank <- rank(x)
  N <- attr(x, "Size")
  div <- length(x)/2
  irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))
  icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))
  within <- matched(irow, icol, grouping)
  if (!any(within)) 
    stop("there should be replicates within groups")
  aver <- tapply(x.rank, within, mean)
  statistic <- -diff(aver)/div
  cl.vec <- rep("Between", length(x))
  take <- as.numeric(irow[within])
  cl.vec[within] <- levels(grouping)[grouping[take]]
  cl.vec <- factor(cl.vec, levels = c("Between", levels(grouping)))
  ptest <- function(take, ...) {
    cl.perm <- grouping[take]
    tmp.within <- matched(irow, icol, cl.perm)
    tmp.ave <- tapply(x.rank, tmp.within, mean)
    -diff(tmp.ave)/div
  }
  permat <- getPermuteMatrix(permutations, N, strata = strata)
  if (ncol(permat) != N) 
    stop(gettextf("'permutations' have %d columns, but data have %d rows", 
                  ncol(permat), N))
  permutations <- nrow(permat)
  if (permutations) {
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
      if (.Platform$OS.type == "unix" && !hasClus) {
        perm <- unlist(mclapply(1:permutations, function(i, 
                                                         ...) ptest(permat[i, ]), mc.cores = parallel))
      }
      else {
        if (!hasClus) {
          parallel <- makeCluster(parallel)
        }
        perm <- parRapply(parallel, permat, ptest)
        if (!hasClus) 
          stopCluster(parallel)
      }
    }
    else {
      perm <- sapply(1:permutations, function(i) ptest(permat[i, 
      ]))
    }
    p.val <- (1 + sum(perm >= statistic - EPS))/(1 + permutations)
  }
  else {
    p.val <- perm <- NA
  }
  sol$signif <- p.val
  sol$perm <- perm
  sol$permutations <- permutations
  sol$statistic <- as.numeric(statistic)
  sol$class.vec <- cl.vec
  sol$dis.rank <- x.rank
  sol$dissimilarity <- attr(x, "method")
  sol$control <- attr(permat, "control")
  class(sol) <- "anosim"
  sol
}

`getPermuteMatrix` <-
  function(perm, N,  strata = NULL)
  {
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
      perm <- how(nperm = perm)
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
      if (inherits(perm, "how") && is.null(getBlocks(perm)))
        setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
      perm <- shuffleSet(N, control = perm)
    else { # matrix: check that it *strictly* integer
      if(!is.integer(perm) && !all(perm == round(perm)))
        stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
      attr(perm, "control") <-
        structure(list(within=list(type="supplied matrix"),
                       nperm = nrow(perm)), class = "how")
    perm
  }


PDH.anosim <- anosim1(PDH.ds[PDH.ds$Sex=="M",c(3:7,9:40)], 
                     PDH.ds$Group[PDH.ds$Sex=="M"], permutations = 9999)
autoplot(PDH.anosim, notch = FALSE) +
    ylab("Dissimilarity rank")  +
    ggtitle("Bray-Curtis ANOSIM") +
    theme_classic()

# 1.3 Try 0-5 approach (Any sorts of) and food additive
PDHdb$Early_FA <- as.character(apply(PDHdb[,c(9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38)], 1,
                        function(x) if (any(x==1)) {1} else {0}))

PDHdb$All_RP <- combdb$All_RP[match(PDHdb$Sub_ID,combdb$Sub_ID)]
kruskal.test(PDHdb$Early_FA~PDHdb$All_RP)

# Try each FA vs. BMI, weight, height ----
library(corrplot)
library(ggpubr)

HLMdb <- left_join(HLMdb,MomFA_W[,c(1:5,7:15,17:27,29:31)], by = "Sub_ID")
HLMdb$Height <- combdb$Height[match(HLMdb$Sub_ID,combdb$Sub_ID)]
HLMdb$Height[is.na(HLMdb$Height)] <- dadcombdb$Height[match(HLMdb$Sub_ID[is.na(HLMdb$Height)],
                                                            dadcombdb$Sub_ID)]
HLMdb$Income_fam <- combdb$Income_fam_new[match(HLMdb$Sub_ID,combdb$Sub_ID)]
HLMdb$Income_fam[is.na(HLMdb$Income_fam)] <- combdb$Income_fam_new[match(substr(HLMdb$Sub_ID[HLMdb$Sex=="M"],7,9),
                                                                   substr(HLMdb$Sub_ID,7,9))]
HLMdb$Probiotic <- CDHdb$Any_Syn[match(HLMdb$Sub_ID,CDHdb$Sub_ID)]


momall <- HLMdb %>% filter(Sex=="F")
momall$Early_FA <- PDHdb$Early_FA[match(momall$Sub_ID,PDHdb$Sub_ID)]
momall$BF <- PDHdb$BF[match(momall$Sub_ID,PDHdb$Sub_ID)]


for (i in c(13:16,18,36:38)){
  print(colnames(momall)[i])
  print(levels(as.factor(momall[,i])))
}

momall$Avg_takeout <- ifelse(grepl("3",momall$Avg_takeout),">=3", paste(momall$Avg_takeout))
momall$Avg_takeout[momall$Avg_takeout=="NA"] <- NA
momall$Avg_meal <- factor(as.factor(momall$Avg_meal),
                          levels(as.factor(momall$Avg_meal))[c(2,3,1)])
momall$Avg_takeout <- factor(as.factor(momall$Avg_takeout),
                          levels(as.factor(momall$Avg_takeout))[c(3,2,1)])
momall$Stool_type_Group <- factor(as.factor(momall$Stool_type_Group),
                          levels(as.factor(momall$Stool_type_Group))[c(2,3,1)])
momall$Group <- factor(as.factor(momall$Group),
                                  levels(as.factor(momall$Group))[c(2,3,1)])
momall$Group <- as.character(momall$Group)
momall$VM_group <- factor(as.factor(momall$VM_group),
                                  levels(as.factor(momall$VM_group))[c(3,2,1)])
momall$SAC_Group <- factor(as.factor(momall$SAC_Group),
                          levels(as.factor(momall$SAC_Group))[c(2,3,1)])
momall$Early_FA <- as.character(momall$Early_FA)
momall$BF <- as.character(momall$BF)
# momall$FA_log <- log2(momall$All_RP)
momall$Sweet_log <- log2(momall$Sweet_abs)
momall$Emu_log <- log2(momall$Emu_abs)
momall$PresDisp_log <- log2(momall$PresDisp_abs)
# momall <- momall %>% filter(SAC_RP<1000.1)
momall$Sweet_log[momall$Sweet_log=="-Inf"] <- 0
momall$PresDisp_log[momall$PresDisp_log=="-Inf"] <- 0
momall$Alco_con[is.na(momall$Alco_con)] <- 0

# BMI, weight, height, uni-variable

uni_c <- data.frame()
for (i in colnames(momall)[c(4:7,12:14,18,20,22:40,45:47,49,50:55)]){
  if (all(momall[,i]==0)){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~",i)),
                  momall,family ="gaussian")
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
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))

uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
uni_c$Var <- rownames(uni_c)


BMI_Uni <- uni_c

# Weight_uni <- uni_c
# Weight_uni$Upper[Weight_uni$Upper>100] <- NA
# Height_uni <- uni_c

BMI_Uni_plt <- ggplot(BMI_Uni, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+
  scale_y_discrete(name="", labels=c("Avg_takeout1-2"="Average takeout \n1-2/day",
                                         "VM_groupMoverV"="Meat > vegetable",
                                         "SAC"="Saccharin",
                                     "Emu_log"="log2 Emulsifier")) +
  ylab("")+
  xlab("Estimate (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("BMI")

Weight_uni_plt <- ggplot(Weight_uni, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10)) + 
  scale_y_discrete(name="", labels=c("Avg_takeout1-2"="Average takeout \n1-2/day",
                                     "VM_groupMoverV"="Meat > vegetable",
                                     "SAC"="Saccharin",
                                     "Sweet_log" = "log2 Sweetener",
                                     "Emu_log"="log2 Emulsifier",
                                     "TiO2_perKG"="Titanium \nDioxide/year/kg",
                                     "AlcoYes" = "Alcohole drinking")) +
  ylab("")+
  xlab("Estimate (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Weight")


Height_uni_plt <-ggplot(Height_uni, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
    scale_y_discrete(name="", labels=c("Avg_takeout1-2"="Average takeout \n1-2/day",
                                     "Emu_log"="log2  Emusifier",
                                     "BF1" = "Breast feeding")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Height") 
  

ggsave("Outcome_Uni_14072021.eps",width = 13, height = 5,units = "in")
BMI_Uni_plt|Weight_uni_plt|Height_uni_plt
invisible(dev.off())

# 2.2 BMI, weight, height, multi-variate
mglm <- glm(BMI~Age+VM_group+Any_dis+Income_fam + Avg_takeout +
              SAC+Emu_log
            ,momall,family = "gaussian")
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


BMI_mul <- mglmtb
# Weight_mul <- mglmtb
# Weight_mul$Upper[Weight_mul$Upper>100] <- NA
# Height_mul <- mglmtb

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

Weight_mul_plt <- ggplot(Weight_mul, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#54278f") +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+ 
  scale_y_discrete(name="", labels=c("Avg_takeout1-2"="Average takeout \n1-2/day",
                                     "VM_groupMoverV"="Meat > vegetable",
                                     "AlcoYes" = "Alcohole drinking")) +
  ylab("")+
  xlab("Estimate (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Weight")


Height_mul_plt <- ggplot(Height_mul, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#54278f") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  # scale_x_continuous(limits = c(0,2)) +
  scale_y_discrete(name="", labels=c("Emu_log"="log2  Emusifier",
                                     "BF1" = "Breast feeding")) +
  ylab("")+
  xlab("Estimate (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Height") 

ggsave("Outcome_Mul_14072021.eps",width = 12, height = 5,units = "in")
BMI_mul_plt|Weight_mul_plt|Height_mul_plt
invisible(dev.off())


# 3 Maternal
momall$GDM <- as.numeric(combdb$GDM_DM[match(momall$Sub_ID,combdb$Sub_ID)])
momall$Weight_T1BL <- combdb$Weight_T1BL[match(momall$Sub_ID,combdb$Sub_ID)]
momall$Weight_T3BL <- combdb$Weight_T3BL[match(momall$Sub_ID,combdb$Sub_ID)]


uni_c <- data.frame()
for (i in colnames(momall)[c(4:7,12:14,16,18,22:30,39:42,36:38)]){
  
  if (all(momall[,i]==0)){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("Weight_T3BL~",i)),
                  momall,family ="gaussian")
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

GDM_uni <- uni_c
# Weight_T1BL <- uni_c
# Weight_T3BL <- uni_c

GDM_uni_plt <-  ggplot(GDM_uni, aes(x = OR, y = Var))+
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  scale_y_discrete(name="", labels=c( "AlcoYes" = "Alcohole drinking",
                                      "Any_disYes" = "Has co-morbidity",
                                      "FoodAYes" = "Food allergy")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("GDM") 

Weight_T1BL_plt <- ggplot(Weight_T1BL, aes(x = OR, y = Var))+
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  scale_y_discrete(name="", labels=c( "PresDisp_log"="log2  Preservative \nor dispersant")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Weight change (T1 to baseline)") 

Weight_T3BL_plt <- ggplot(Weight_T3BL, aes(x = OR, y = Var))+
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#756bb1") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  # scale_y_discrete(name="")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Weight change (T3 to baseline)") 

ggsave("Secoutcome_Uni_28062021.eps",width = 12, height = 5,units = "in")
GDM_uni_plt|Weight_T1BL_plt|Weight_T3BL_plt
invisible(dev.off())

# 3.2 Maternal, multi-variate
mglm <- glm(GDM~Age+VM_group+Any_dis+Income_fam+
              Avg_takeout+Alco+FoodA
            ,momall,family = "binomial")
mglmtb <- cbind(summary(mglm)$coefficients[,c(1,4)],exp(cbind(OR = coef(mglm), confint(mglm))))
mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)
mglmtb$Pr...t.. <- ifelse(mglmtb$Pr...t..<0.001, paste("<0.001"),
                          paste(mglmtb$Pr...t..))
mglmtb$Pr...z.. <- ifelse(mglmtb$Pr...z..<0.001, paste("<0.001"),
                          paste(mglmtb$Pr...z..))
mglmtb <- mglmtb %>% filter(Pr...t..=="<0.001"|Pr...t..<0.05)
mglmtb <- mglmtb %>% filter(Pr...z..=="<0.001"|Pr...z..<0.05)

mglmtb <- mglmtb %>% filter(!grepl("Intercep",rownames(mglmtb)))

mglmtb$Var <- rownames(mglmtb)
colnames(mglmtb) <- c("Est","p.value","OR","Lower","Upper","Var")
# Weight_T3BL_mul <- mglmtb
GDM_mul <- mglmtb

GDM_mul_plt <-  ggplot(GDM_mul, aes(x = OR, y = Var))+
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#54278f") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  scale_y_discrete(name="", labels=c("Any_disYes" = "Has co-morbidity")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("GDM") 


Weight_T3BL_mul_plt <-  ggplot(Weight_T3BL_mul, aes(x = OR, y = Var))+
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, color = "#54278f") +
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10))+  
  # scale_y_discrete(name="")) +
  ylab("")+
  xlab("Odds ratio (95%-CI)") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="", size = 3.5, hjust = 0) + 
  ggtitle("Weight change (T3 to baseline)") 



ggsave("Secoutcome_Mul_28062021.eps",width = 8, height = 4,units = "in")
GDM_mul_plt|Weight_T3BL_plt
invisible(dev.off())

# Correlation drop
plot_list = list()
for (i in c(22:34)) {
  my.data <- data.frame(x = HLMdb$BMI[HLMdb$Sex=="F"], 
                        y = log2(HLMdb[,i][HLMdb$Sex=="F"]))
  j <- colnames(HLMdb)[i]
  a <- cor.test(my.data$x,
                my.data$y,method = "spearman")
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      scale_color_manual(values=c("lightslateblue")) +
      geom_point(size = 2) +
      theme_classic()+
      ylab(paste("Annual estimated ",j," level")) +
      xlab("Weight") +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25) +
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
      scale_color_manual(values=c("lightslateblue")) +
      geom_point(size = 2) +
      theme_classic() +
      ylab(paste("Annual estimated ",j," level")) +
      xlab("Weight") +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25) +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black")) + 
      annotate("text", x = max(my.data$x,na.rm = T)*3/4,
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3)))
    
    plot_list <- c(plot_list, list(p))
    
  }
}

wrap_plots(plot_list,ncol = 2)

# Try to compare with the mothers with eczema off-springs 

ECZlist <- data.frame(read_excel("Data/Eczema subject's mother ID--16.xlsx"))
momall$EczGroup <- ECZlist$ID[match(momall$Sub_ID,ECZlist$ID)]
momall$EczGroup <- ifelse(is.na(momall$EczGroup),0,1)


plot_list = list()
for (i in c(22:34)) {
  my.data <- data.frame(x = as.character(momall$EczGroup), 
                        y = log2(momall[,i]))
  j <- colnames(momall)[i]
  a <- wilcox.test(my.data$y~my.data$x)
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
      geom_boxplot() +
      scale_color_manual(values=c( "#2b8cbe","#bd0026")) +
      theme_classic()+
      ylab(paste("Annual estimated ",j," level")) +
      stat_compare_means(label.x = 1.5,
                         method = "wilcox", size = 6 ) +
      theme(legend.position = "none",
            axis.text = element_text(color = "black", size = 14),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_blank(),
            plot.background = element_rect(fill = "white", colour = "red3", size = 1)) +
      scale_x_discrete(name="", labels=c( "1" = "Off-spring eczema \nYes (N = 14)",
                                          "0" = "Off-spring eczema \nNo (N = 410)")) 
     
    plot_list <- c(plot_list, list(p))
    
  }else{
    p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
      geom_boxplot() +
      scale_color_manual(values=c( "#2b8cbe","#bd0026")) +
      theme_classic()+
      ylab(paste("Annual estimated ",j," level")) +
      stat_compare_means(label.x = 1.5,
                         method = "wilcox", size = 6 ) +
      theme(legend.position = "none",
            axis.text = element_text(color = "black", size = 14),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_blank()) +
      scale_x_discrete(name="", labels=c( "1" = "Off-spring eczema \nYes (N = 14)",
                                          "0" = "Off-spring eczema \nNo +(N = 410)")) 
    
    
    plot_list <- c(plot_list, list(p))
    
  }
}

ggsave("FA_Ecz_06072021.eps",width = 16, height = 29,units = "in")
wrap_plots(plot_list,ncol = 3)
invisible(dev.off())

# Draw heatmap for proportion in early life dietary habit----
# # PDHsum <- 
# 
# # PDHsum <- reshape(PDHsum[,c(4:6)], idvar = "Item", timevar = "Sex", direction = "wide")
# colnames(PDHsum)[c(2:3)] <- c("M","F")
# rownames(PDHsum) <- PDHsum$Item

ggsave("PDH_gender.eps",width = 20,height = 30,units = "cm")
ggplot(PDHsum, aes(x=Sex, y=Item_od, fill=Perc)) +
  geom_tile(color = "Gray20", aes(width = 1)) +
  scale_fill_gradient2(high = "#a50f15",
                       low = "#ffffb2", 
                       limits = c(1,100)) +
  theme_classic() +
  scale_x_discrete(labels = c("M" = "Male",
                              "F" = "Female")) +
  labs(fill = "Proportion") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black",size = 12),
        axis.text.y = element_text(margin = margin(r = -20)),
        axis.title = element_blank()) 
invisible(dev.off())

# Draw line bar for food additive
MomFA$Group <- ifelse(MomFA$Sex=="F",1,2)
MomFA$Weight <- demodb$Weight[match(MomFA$Sub_ID,demodb$Sub_ID)]
MomFAW <- MomFA
MomFAW[,c(2:10)] <- apply(MomFAW[,c(2:10)],2,function(x) round(x/MomFAW$Weight,1))
MomFAW <- MomFAW %>% filter(complete.cases(Weight))
colnames(MomFAW)[c(3:12)] <- c("P80", "CMC", "CRN", "AlSiO", "X.SO3.2", "TiO2", "ASP", "SUC", "SAC")
colnames(MomFAW)[12] <- "Other_FA"
MomFAsum <- melt(MomFAW[,-c(1,2,12:13)],
     measure.vars = c("P80", "CMC", "CRN", "AlSiO", "X.SO3.2", "TiO2", "ASP", "SUC", "SAC"),
     value.name = "Value")


ggsave("FA_gender.eps",width = 18,height = 12,units = "cm")
ggplot(data = MomFAsum, aes(x = variable, y = log2(Value), color = Sex)) +
  geom_boxplot() +
  scale_color_manual(values=c( "#018571","#0571b0"),
                     labels = c("Female", "Male")) +
  theme_classic() +
  ylab("log 2 food additive value/bnkg") +
  stat_compare_means(method = "t.test") +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.position = c(0.85,0.9),
        legend.direction = "horizontal")
invisible(dev.off())

MomFAsum1 <- melt(MomFA[,c(23:32)],
                 measure.vars = c("log2P80", "log2CMC", 
                                  "log2CRN", "log2AlSiO", "log2X.SO3.2",
                                  "log2TiO2", "log2ASP", "log2SUC", "log2SAC"),
                 value.name = "Value")


ggsave("FA_gender_unadj.eps",width = 18,height = 12,units = "cm")
ggplot(data = MomFAsum1, aes(x = variable, y = log2(Value), color = Sex)) +
  geom_boxplot() +
  scale_color_manual(values=c( "#018571","#0571b0"),
                     labels = c("Female", "Male")) +
  theme_classic() +
  ylab("log 2 food additive value/year") +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.position = c(0.85,0.9),
        legend.direction = "horizontal")
invisible(dev.off())


# Baseline variables and BMI (05082021)----
MOMcomb$Income_fam_new <- demodb$Income_fam_new[match(MOMcomb$Sub_ID,demodb$Sub_ID)]
MOMcomb$Income_fam_new[MOMcomb$Sex=="M"] <- demodb$Income_fam_new[demodb$Sex=="F"][match(substr(MOMcomb$Sub_ID,7,9),
                                                                                          substr(demodb$Sub_ID,7,9))]
MOMcomb$Income_fam_new <- factor(as.factor(MOMcomb$Income_fam_new),
                                levels(as.factor(MOMcomb$Income_fam_new))[c(1,3,4,2)])
colnames(MOMcomb)[c(22:59)] <- str_replace_all(colnames(MOMcomb)[c(22:59)],"-","_")
MOMcomb$OVW <- ifelse(is.na(MOMcomb$OVW),NA
                      ,ifelse(MOMcomb$BMI<23,0,1))

MOMcomb_sub <- MOMcomb %>% dplyr::filter(!Sub_ID %in% laterecruite)
set.seed(10)
uni_c <- data.frame()
for (i in colnames(MOMcomb_sub)[c(139:164)]){
  if (all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"]==0)|all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])]=="1No")){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("Pre_BMI~",i)),
                  MOMcomb_sub[MOMcomb_sub$Sex=="F",],family ="gaussian")
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

# uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
uni_c$Var <- rownames(uni_c)


writexl::write_xlsx(uni_c,"M_demo_BMI_07092021.xlsx")

# Plot log 2 FAs and bmi association
plot_list = list()
for (i in c(61:69,84,85,83,70)) {
  my.data <- data.frame(x = MOMcomb$BMI, y = log2(MOMcomb[,i]), 
                        Sex = MOMcomb$Sex)
  j <- colnames(MOMcomb)[i]
    p <- ggplot(my.data,  aes(x = x, y = y, color = Sex)) +
      scale_color_manual(values=c("#018571","#0571b0")) +
      geom_point(size = 2) +
      theme_classic()+
      ylab(paste("log 2",j)) +
      xlab("BMI") +
      stat_cor(method="spearman",position = "identity",
               r.digits = 2,p.digits = 3, label.x = 30) +
      stat_smooth(method = "lm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25) +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black"))
    plot_list <- c(plot_list, list(p))

}

corlgd <- ggplot(my.data,  aes(x = x, y = y, color = Sex)) +
  scale_color_manual(values=c("#018571","#0571b0"),
                     label = c("Female","Male")) +
  geom_point(size = 2) 

plot_list <- c(plot_list, list(as_ggplot(get_legend(corlgd))))

# library(patchwork)

ggsave("FA_BMI_05082021.eps",width = 16, height = 12,units = "in")
wrap_plots(plot_list,ncol = 3)
invisible(dev.off())


# Uni for bmi groups
# uni_c <- data.frame()
# for (i in colnames(MOMcomb)[c(7,8,9,12,13,15:18,73:82,86:88)]){
#   if (all(MOMcomb[,i]==0)){
#     print(i)
#   }
#   else {
#     print(i)
#     temp_a <- glm(as.formula(paste0("BMI_Congroup~",i)),
#                   MOMcomb[MOMcomb$Sex=="M",],family ="poisson")
#     temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
#                     confint(temp_a,level = 0.95),
#                     summary(temp_a)$coefficients[,c(4)])
#     uni_c <- rbind(uni_c,temp_b)
#   }
# }
# 
# uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# # rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
# colnames(uni_c) <- c("Est","Lower","Upper","p.value")
# uni_c <- apply(uni_c, 2, function(x) round(x,4))
# uni_c <- data.frame(uni_c)
# uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))
# 
# # uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
# uni_c$Var <- rownames(uni_c)
# 
# writexl::write_xlsx(uni_c,"M_demo_BMIgroup_11082021.xlsx")


# Multivariable
# BMI in continuous
# Mother model
mglmtb <- data.frame()
for (i in colnames(MOMcomb_sub)[c(139:164)]){
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~ Age + Edu + SpeDie_group + VM_group +", i)),
                  MOMcomb_sub[MOMcomb_sub$Sex=="F",], family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                    confint(temp_a,level = 0.95),
                    summary(temp_a)$coefficients[,c(4)])
    mglmtb <- rbind(mglmtb,temp_b)
}


# 
# mglm <- glm(BMI ~ Age + Edu + SpeDie_group + VM_group +
#               log2P80_W + log2CMC_W + log2CRN_W + log2AlSiO_W + 
#               log2X.SO3.2_W + log2TiO2_W + log2ASP_W + log2SUC_W +
#               log2SAC_W,
#             MOMcomb[MOMcomb$Sex=="F",],family = "gaussian")
# 
# 
# mglmtb <- cbind(summary(mglm)$coefficients[,c(1)],confint(mglm),
#                 summary(mglm)$coefficients[,c(4)])
mglmtb <- apply(mglmtb, 2, function(x) round(x,4))
mglmtb <- data.frame(mglmtb)
colnames(mglmtb) <- c("Est","Lower","Upper","p.value")
mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))
# mglmtb <- mglmtb %>% filter(p.value=="<0.001"|p.value<0.06)
mglmtb <- mglmtb %>% dplyr::filter((!grepl("Intercep",rownames(mglmtb)))&grepl(paste0(colnames(MOMcomb_sub)[c(139:164)],
                                                                                      collapse = "|"),
                                                                               rownames(mglmtb)))

mglmtb$Var <- str_replace_all(rownames(mglmtb),"_W","")
mglmtb$Var <- str_replace_all(mglmtb$Var,"log2","")
mglmtb$p.bin <- ifelse(mglmtb$p.value>0.05,"0","1")
mglmtb$ID <- length(mglmtb$Est):1
mglmtb$Var <- reorder(mglmtb$Var,mglmtb$ID)

p1 <- ggplot(mglmtb, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, aes(color = p.bin)) + 
  theme_classic() +
  xlab("Estimate") +
  ylab("") +
  scale_color_manual(values = c("gray20", "#756bb1"),
                     label = c("Insignificant","<0.05"),
                     name = "p.value") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10)) 


a <- MOMcomb_sub %>% dplyr::filter(Sex == "F"& VM_group == "MoverV")
mglmtb_sub <- data.frame()
for (i in colnames(a)[c(73:82,86:88)]){
  print(i)
  temp_a <- glm(as.formula(paste0("BMI~ Age + Edu + SpeDie_group +", i)),
                a, family ="gaussian")
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

mglmtb_sub <- mglmtb_sub %>% dplyr::filter((!grepl("Intercep",rownames(mglmtb_sub)))&grepl(paste0(colnames(MOMcomb_sub)[c(73:82,86:88)],
                                                                                                  collapse = "|"),
                                                                                           rownames(mglmtb_sub)))

mglmtb_sub$Var <- str_replace_all(rownames(mglmtb_sub),"_W","")
mglmtb_sub$Var <- str_replace_all(mglmtb_sub$Var,"log2","")
mglmtb_sub$p.bin <- ifelse(mglmtb_sub$p.value>0.05,"0","1")
mglmtb_sub$ID <- length(mglmtb_sub$Est):1
mglmtb_sub$Var <- reorder(mglmtb_sub$Var,mglmtb_sub$ID)

p2 <- ggplot(mglmtb_sub, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, aes(color = p.bin)) + 
  theme_classic() +
  xlab("Estimate") +
  ylab("") +
  scale_color_manual(values = c("gray20", "#756bb1"),
                     label = c("Insignificant","<0.05"),
                     name = "p-value") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10)) 


p1+p2+plot_layout(widths = c(1, 1.08))


# Father model
mglmtb <- data.frame()
for (i in colnames(MOMcomb_sub)[c(139:164)]){
  print(i)
  temp_a <- glm(as.formula(paste0("BMI~ Age +  VM_group + HG + HT +", i)),
                MOMcomb_sub[MOMcomb_sub$Sex=="M",], family ="gaussian")
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
# mglmtb <- mglmtb %>% filter(p.value=="<0.001"|p.value<0.06)
mglmtb <- mglmtb %>% dplyr::filter((!grepl("Intercep",rownames(mglmtb)))&grepl(paste0(colnames(MOMcomb_sub)[c(139:164)],
                                                                                      collapse = "|"),
                                                                               rownames(mglmtb)))

mglmtb$Var <- str_replace_all(rownames(mglmtb),"_W","")
mglmtb$Var <- str_replace_all(mglmtb$Var,"log2","")
mglmtb$p.bin <- ifelse(mglmtb$p.value>0.05,"0","1")
mglmtb$ID <- length(mglmtb$Est):1
mglmtb$Var <- reorder(mglmtb$Var,mglmtb$ID)

p3 <- ggplot(mglmtb, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, aes(color = p.bin)) + 
  theme_classic() +
  xlab("Estimate") +
  ylab("") +
  scale_color_manual(values = c("gray20", "#756bb1"),
                     label = c("Insignificant","<0.05"),
                     name = "p.value") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10)) 


a <- MOMcomb_sub %>% dplyr::filter(Sex == "M"& VM_group == "MoverV")
mglmtb_sub <- data.frame()
for (i in colnames(a)[c(73:82,86:88)]){
  if(i == "log2AlSiO_W"|i == "log2TiO2_W"){
    print(i)
  }else{
  print(i)
  temp_a <- glm(as.formula(paste0("BMI~ Age + Dairy_5_10y + CarboHfood_5_10y + HG + HT +", i)),
                a, family ="gaussian")
  temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                  confint(temp_a,level = 0.95),
                  summary(temp_a)$coefficients[,c(4)])
  mglmtb_sub <- rbind(mglmtb_sub,temp_b)}
}

mglmtb_sub <- apply(mglmtb_sub, 2, function(x) round(x,4))
mglmtb_sub <- data.frame(mglmtb_sub)
colnames(mglmtb_sub) <- c("Est","Lower","Upper","p.value")
mglmtb_sub$p.value <- ifelse(mglmtb_sub$p.value<0.001, paste("<0.001"),
                             paste(mglmtb_sub$p.value))

mglmtb_sub <- mglmtb_sub %>% dplyr::filter((!grepl("Intercep",rownames(mglmtb_sub)))&grepl(paste0(colnames(MOMcomb_sub)[c(73:82,86:88)],
                                                                                                  collapse = "|"),
                                                                                           rownames(mglmtb_sub)))
mglmtb_sub <- rbind(mglmtb_sub[c(1:3),],
                    "AlSiO" = rep(NA,4),
                    mglmtb_sub[4,],
                    "TiO2" = rep(NA,4),
                    mglmtb_sub[c(5:11),])
mglmtb_sub$Var <- str_replace_all(rownames(mglmtb_sub),"_W","")
mglmtb_sub$Var <- str_replace_all(mglmtb_sub$Var,"log2","")
mglmtb_sub$p.bin <- ifelse(mglmtb_sub$p.value>0.05,"0","1")

mglmtb_sub$ID <- length(mglmtb_sub$Est):1
mglmtb_sub$Var <- reorder(mglmtb_sub$Var,mglmtb_sub$ID)

p4 <- ggplot(mglmtb_sub, aes(x = Est, y = Var))+
  geom_vline(aes(xintercept = 0), size = .5, linetype = "dashed" ) +
  geom_errorbarh(aes(xmax = Upper , xmin = Lower), size = 1, height = 
                   .2, color = "gray10") +
  geom_point(size = 3.5, aes(color = p.bin)) + 
  theme_classic() +
  xlab("Estimate") +
  ylab("") +
  scale_color_manual(values = c("gray20", "#756bb1"),
                     label = c("Insignificant","<0.05",""),
                     name = "p-value") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 10)) 
p3+p4+plot_layout(widths = c(1, 1.08))

writexl::write_xlsx(mglmtb,"F_BMI_multi_11082021.xlsx")

# Try in bmi group
library(nnet)
library("AER")
catlm <- multinom(BMI_Group~Age+HL, MOMcomb)
coeftest(catlm)




for (i in colnames(MOMcomb)[c(7,8,9,15:18,12,13,22:59,73:82,86:88,90:102,105:108)]){
  if (all(MOMcomb[,i]==0)){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("BMI~",i)),
                  MOMcomb[MOMcomb$Sex=="F",],family ="gaussian")
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
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))

# uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
uni_c$Var <- rownames(uni_c)

MOMcomb[,c(139:147)] <- log2(MOMcomb[,c(139:147)])
MOMcomb[,c(139:147)][MOMcomb[,c(139:147)]=="Inf"] <- NA
MOMcomb[,c(139:147)][MOMcomb[,c(139:147)]=="Inf"] <- NA

MOMcomb_sub <- MOMcomb %>% dplyr::filter(!Sub_ID %in% laterecruite)

set.seed(10)
uni_c <- data.frame()
for (i in colnames(MOMcomb_sub)[c(22:58,64:76,79:84)]){
  if (all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"]==0)|
      all(MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])]=="1No")|
      i=="Immune"|i=="Cancer"){
    print(i)
  }
  else {
    print(i)
    temp_a <- glm(as.formula(paste0("T2_WC~",i)),
                  MOMcomb_sub[MOMcomb_sub$Sex=="F",],family ="gaussian")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                    exp(cbind(OR = coef(temp_a), 
                              confint.default(temp_a,level = 0.95))))
    uni_c <- rbind(uni_c,temp_b)
  }
}

uni_c <- uni_c[!grepl("Intercep",rownames(uni_c)),]
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]
# colnames(uni_c) <- c("Est","Lower","Upper","p.value")
uni_c <- apply(uni_c, 2, function(x) round(x,4))
uni_c <- data.frame(uni_c)
# uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))

# uni_c <- uni_c %>% filter(p.value=="<0.001"|p.value<0.06)
uni_c$Var <- rownames(uni_c)
uni_c <- uni_c[uni_c$Pr...t..<0.05,]

# BMI_FA_by boxplot tertile group 4 sl present----
# mom TiO2

statecomp <- list(c("Low","Normal"),
                  c("Normal","High"),
                  c("Low","High"))

ggplot(MOMcomb_sub[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub$TiO2_cut3),], 
       aes(x=TiO2_cut3, y=BMI, color = TiO2_cut3)) + 
  geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(fill = TiO2_cut3)) +
  scale_color_manual(values=c("#2c7bb6", "#fdae61","#d7191c")) +
  scale_x_discrete(label = c("Low \n(N = 231)",
                             "Normal \n(N = 260)",
                             "High \n(N = 233)")) +
  theme_classic() +
  stat_compare_means(method = "kruskal",label.x = 2.5,label.y = 53) +
  stat_compare_means(comparisons = statecomp,label = "p.signif",
                     method = "wilcox") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 16),
        axis.title = element_text(size = 20)) +
  xlab("Annual TiO2 intake groups")

ggplot(MOMcomb_sub[MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub$SAC_cut3),], 
       aes(x=SAC_cut3, y=BMI, color = SAC_cut3)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#2c7bb6", "#fdae61","#d7191c")) +
  scale_x_discrete(label = c("Low \n(N = 251)",
                             "Normal \n(N = 251)",
                             "High \n(N = 250)")) +
  theme_classic() +
  stat_compare_means(method = "kruskal",label.x = 2.5,label.y = 63) +
  stat_compare_means(comparisons = statecomp,label = "p.signif",
                     method = "wilcox") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 16),
        axis.title = element_text(size = 20)) +
  xlab("Annual SAC intake groups")

ggplot(MOMcomb_sub[MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub$CMC_cut3),], 
       aes(x=CMC_cut3, y=BMI, color = CMC_cut3)) + 
  geom_boxplot() +
  scale_color_manual(values=c("#2c7bb6", "#fdae61","#d7191c")) +
  scale_x_discrete(label = c("Low \n(N = 251)",
                             "Normal \n(N = 251)",
                             "High \n(N = 250)")) +
  theme_classic() +
  stat_compare_means(method = "kruskal",label.x = 2.5,label.y = 63) +
  stat_compare_means(comparisons = statecomp,label = "p.signif",
                     method = "wilcox") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black",size = 16),
        axis.title = element_text(size = 20)) +
  xlab("Annual SAC intake groups")

plot_list = list()
for (i in c(88:100)) {
  my.data <- data.frame(x = MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])],
                        y = MOMcomb_sub$Pre_BMI[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])])
  j <- colnames(MOMcomb_sub)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    scale_color_manual(values=c("#5e3c99", "#fdb863","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "kruskal",label.x = 1,label.y = 40) +
    stat_compare_means(comparisons = statecomp,label = "p.signif",
                       method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 10),
          axis.title = element_text(size = 12)) +
    xlab(paste0(j)) +
    ylab("BMI") 
    
    
  plot_list <- c(plot_list, list(p))
  
}

ggsave("Mother_FA_Group_24092021.eps",width = 20, height = 10,units = "in")
wrap_plots(plot_list,ncol = 7)
invisible(dev.off())

plot_list = list()
for (i in c(139:151)) {
  my.data <- data.frame(x = MOMcomb_sub[,i][MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub[,i])],
                        y = MOMcomb_sub$BMI[MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub[,i])])
  j <- colnames(MOMcomb_sub)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    scale_color_manual(values=c("#2c7bb6", "#fdae61","#d7191c")) +
    theme_classic()+
    stat_compare_means(method = "kruskal",label.x = 1,label.y = 40) +
    stat_compare_means(comparisons = statecomp,label = "p.signif",
                       method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 10),
          axis.title = element_text(size = 12)) +
    xlab(paste0(j)) +
    ylab("BMI") 
  
  
  plot_list <- c(plot_list, list(p))
  
}

ggsave("Father_FA_Group_24092021.eps",width = 20, height = 10,units = "in")
wrap_plots(plot_list,ncol = 7)
invisible(dev.off())


MOMcomb_sub[,c(94:96,98)] <- apply(MOMcomb_sub[,c(94:96,98)],2,
                                   function(x) ifelse(x=="High","High","Not high"))
plot_list = list()
for (i in c(94:96,98)) {
  my.data <- data.frame(x = MOMcomb_sub[,i][MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub[,i])],
                        y = MOMcomb_sub$BMI[MOMcomb_sub$Sex=="M"&complete.cases(MOMcomb_sub[,i])])
  j <- colnames(MOMcomb_sub)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
    scale_color_manual(values=c("#5e3c99", "#fdb863","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "kruskal",label.x = 1,label.y = 1.2) +
    stat_compare_means(comparisons = statecomp,label = "p.signif",
                       method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 10),
          axis.title = element_text(size = 12)) +
    xlab(paste0(j)) +
    ylab("Male BMI") 
  
  
  plot_list <- c(plot_list, list(p))
  
}

ggsave("Birth_Weight_Emu_14112021.eps",width = 15, height = 7,units = "in")
wrap_plots(plot_list,ncol = 4)
invisible(dev.off())
  


MOMcomb_sub[,c(88:90,97)] <- apply(MOMcomb_sub[,c(88:90,97)],2,
                                   function(x) ifelse(x=="High","2High","1Normal"))

colnames(MOMcomb_sub)[c(88:90,97)] <- c("P80", "carboxymethylcellulose",
                                        "carrageenan", "Overall emulsifier")


plot_list = list()
for (i in c(88,90)) {#88:90,97
  my.data <- data.frame(x = MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])&MOMcomb_sub[,i]!="Low"],
                        y = MOMcomb_sub$T2_WC[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])&MOMcomb_sub[,i]!="Low"])
  j <- colnames(MOMcomb_sub)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    # geom_violin() +
    geom_boxplot(width = 0.6, size = 1.2,
                 position = position_dodge(0.9)) +
    scale_x_discrete(label = c("Normal","High")) +
    geom_point(shape=16, position=position_jitter(0.3),
               size = 2, alpha = 0.3) +
    scale_color_manual(values=c( "#fdb863","#5e3c99","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "wilcox", size = 5,
                       label.x = 0.5,label.y = 0.9) +
    # stat_compare_means(comparisons = statecomp,label = "p.signif",
    #                    method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 14)) +
    xlab(paste0("Annual ",j," intake group")) +
    ylab("Trimester 2 weight gain\n(kg/gestational week)") 
  plot_list <- c(plot_list, list(p))
  
}

# Plot P80 and CRN
set.seed(100)
plot_list = list()
  my.data <- data.frame(x = MOMcomb_sub[,88][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,88])&MOMcomb_sub[,88]!="Low"],
                        y = MOMcomb_sub$T2_WC[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,88])&MOMcomb_sub[,88]!="Low"])
  j <- colnames(MOMcomb_sub)[88]
  p1 <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    # geom_violin() +
    geom_boxplot(width = 0.6, size = 1.2,
                 position = position_dodge(0.9)) +
    scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="Normal"&complete.cases(my.data$y))),
                               paste0("High \nN = ",sum(my.data$x=="High"&complete.cases(my.data$y))))) +
    geom_point(shape=16, position=position_jitter(0.3),
               size = 2, alpha = 0.3) +
    scale_color_manual(values=c( "#fdb863","#5e3c99","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "wilcox", size = 5,
                       label.x = 0.7,label.y = 0.9) +
    # stat_compare_means(comparisons = statecomp,label = "p.signif",
    #                    method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 15)) +
    xlab(paste0("Annual ",j," intake group")) +
    ylab("Trimester 2 weight gain\n(kg/gestational week)") 
  
plot_list <- c(plot_list, list(p1))
set.seed(100)
my.data <- data.frame(x = MOMcomb_sub[,90][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,90])&MOMcomb_sub[,90]!="Low"],
                      y = MOMcomb_sub$T2_WC[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,90])&MOMcomb_sub[,90]!="Low"])
j <- colnames(MOMcomb_sub)[90]
p2 <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
  # geom_violin() +
  geom_boxplot(width = 0.6, size = 1.2,
               position = position_dodge(0.9)) +
  scale_x_discrete(label = c(paste0("Normal \nN = ",sum(my.data$x=="Normal"&complete.cases(my.data$y))),
                             paste0("High \nN = ",sum(my.data$x=="High"&complete.cases(my.data$y))))) +
  geom_point(shape=16, position=position_jitter(0.3),
             size = 2, alpha = 0.3) +
  scale_color_manual(values=c( "#fdb863","#5e3c99","#e66101")) +
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
  xlab(paste0("Annual ",j," intake group")) +
  ylab("Trimester 2 weight gain\n(kg/gestational week)") 
plot_list <- c(plot_list, list(p2))


pdf(file="T2_WC_Emu_24112021..pdf", width = 12, height = 10)
(p1+p2)/
  mglmtab + plot_layout(heights = unit(c(6,2), c('inch', 'null')))
invisible(dev.off())



mglmtb <- data.frame()
for (i in colnames(MOMcomb_sub)[c(88:90)]){
  print(i)
  temp_a <- glm(as.formula(paste0("T2_WC ~ Age + Pre_BMI +", `i`)),
                MOMcomb_sub[MOMcomb_sub[,i]!="Low",], family ="gaussian")
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

# mglmtb <- mglmtb[!grepl("Intercep",rownames(mglmtb)),]3


mglmtb <- mglmtb %>% filter(grepl("High",rownames(mglmtb)))
mglmtb <- mglmtb[-2,]
rownames(mglmtb) <- c("Polysorbate-80", 
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
  tab_add_hline(at.row = c(3), row.side = "bottom", linewidth = 5, linetype = 1) %>% 
  tab_add_footnote(text = "*Linear regression, adjsuted for materanal age, and pBMI. 
  Norlmal emulsifiers intake group as reference group compare to high intake group.",
                   face = "plain", size = 14, hjust = 1) 


plot_list = list()
for (i in c(88:100)) {
  my.data <- data.frame(x = MOMcomb_sub[,i][MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])],
                        y = MOMcomb_sub[MOMcomb_sub$Sex=="F"&complete.cases(MOMcomb_sub[,i])])
  j <- colnames(MOMcomb_sub)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
    scale_color_manual(values=c("#5e3c99", "#fdb863","#e66101")) +
    theme_classic()+
    stat_compare_means(method = "kruskal",label.x = 1,label.y = 1.2) +
    stat_compare_means(comparisons = statecomp,label = "p.signif",
                       method = "wilcox",hide.ns = T) +
    theme(legend.position = "none",
          axis.text = element_text(color = "black",size = 10),
          axis.title = element_text(size = 12)) +
    xlab(paste0(j)) +
    ylab("T2 weight gain (kg/ges week)") 
  
  
  plot_list <- c(plot_list, list(p))
  
}

ggsave("T3_GWG_FA_Group_29092021.eps",width = 20, height = 10,units = "in")
wrap_plots(plot_list,ncol = 7)
invisible(dev.off())

# FA and weight gain stantards

for (i in c(88:100)) {
  print(colnames(MOMcomb_sub)[i])
  print(chisq.test(table(MOMcomb_sub[,i],MOMcomb_sub$IOM_GWG_Over)))
}

MOMcomb_sub$IOM_GWG_Over <- ifelse(MOMcomb_sub$IOM_GWG=="Over",1,0)




# Try pair father and mothers (both high vs others)-----
# paired subgroup
MOMcomb_sub$Sub_ID <- substr(MOMcomb_sub$Sub_ID,0,9)
MOMcomb_pair <- MOMcomb_sub %>% filter(Sex == "M")
MOMcomb_pair <- rbind(MOMcomb_sub[MOMcomb_sub$Sub_ID %in% paste0("PWH1",substr(MOMcomb_pair$Sub_ID,5,9)),],
                      MOMcomb_pair)
ID <- substr(MOMcomb_pair$Sub_ID,5,9)
ID <- ID[duplicated(ID)]
MOMcomb_pair <- MOMcomb_pair %>% filter(Sub_ID %in% paste0("PWH1",ID)|Sub_ID %in% paste0("PWH3",ID))
# regroup each FAs
for(i in 88:100){
  MOMcomb_pair[,i+13] <- NA
  colnames(MOMcomb_pair)[i+13] <- paste0(colnames(MOMcomb_pair)[i],"_Fam")
  for(j in unique(substr(MOMcomb_pair$Sub_ID,5,9))){
    MOMcomb_pair[,i+13][MOMcomb_pair$Sub_ID==paste0("PWH1",j)] <- ifelse(all(c(MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH1",j)],
                                                      MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH3",j)])==3),
                                                  "High",
                                                  ifelse(all(c(MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH1",j)],
                                                               MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH3",j)])==1),"Low","Normal"))
    MOMcomb_pair[,i+13][MOMcomb_pair$Sub_ID==paste0("PWH3",j)] <- ifelse(all(c(MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH1",j)],
                                                                               MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH3",j)])==3),
                                                                         "High",
                                                                         ifelse(all(c(MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH1",j)],
                                                                                      MOMcomb_pair[,i][MOMcomb_pair$Sub_ID==paste0("PWH3",j)])==1),"Low","Normal"))
  }
}



plot_list = list()
for (i in c(101:113)) {
  my.data <- data.frame(x = MOMcomb_pair[,i][MOMcomb_pair$Sex=="F"&complete.cases(MOMcomb_pair[,i])],
                        y = MOMcomb_pair$BMI[MOMcomb_pair$Sex=="F"&complete.cases(MOMcomb_pair[,i])])
  my.data <- my.data %>% filter(x!="Normal")
  j <- colnames(MOMcomb_pair)[i]
  p <- ggplot(my.data,  aes(x = x, y = y, color = x)) +
    geom_violin() +
    geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
    scale_color_manual(values=c("#5e3c99", "#fdb863")) +
    theme_classic()+
    stat_compare_means(method = "wilcox",label.x = 1,label.y = 1.2) +
    xlab(paste0(j)) +
    ylab("BMI") 
  
  
  plot_list <- c(plot_list, list(p))
  
}

wrap_plots(plot_list,ncol = 4)


ggplot(MOMcomb_pair[MOMcomb_pair$Sex=="M"&complete.cases(MOMcomb_pair$SAC_cut3_Fam)&MOMcomb_pair$SAC_cut3_Fam!="Normal",],  
       aes(x = SAC_cut3_Fam, 
           y = BMI,
           color = SAC_cut3_Fam)) +
  geom_violin() +
  geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
  scale_color_manual(values=c("#5e3c99", "#fdb863")) +
  theme_classic()+
  stat_compare_means(method = "wilcox",label.x = 1.3,label.y = 1.2) +
  xlab("Family saccharine intake group") +
  ylab("BMI") +
  ggtitle("Father saccharine and BMI") +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_x_discrete(label = c("High \n(N = 83)","Low \n(N = 72)")) 


ggplot(MOMcomb_pair[MOMcomb_pair$Sex=="F"&complete.cases(MOMcomb_pair$Emu_cut3_Fam)&MOMcomb_pair$Emu_cut3_Fam!="Normal",],  
       aes(x = Emu_cut3_Fam, 
           y = BMI,
           color = Emu_cut3_Fam)) +
  geom_violin() +
  geom_boxplot(width = 0.15, position = position_dodge(0.9)) +
  scale_color_manual(values=c("#5e3c99", "#fdb863")) +
  theme_classic()+
  stat_compare_means(method = "wilcox",label.x = 1.3,label.y = 1.2) +
  xlab("Family emulsifier intake group") +
  ylab("BMI") +
  ggtitle("Mother emulsifier and BMI") +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_x_discrete(label = c("High \n(N = 74)","Low \n(N = 68)")) 


# Try nutrients and clinical outcomes----
NTdb <- data.frame(t(read_excel("Data/Food Additive Calculation_2021.10.08.xlsx",
                                sheet = "Annual Food in kg")))
colnames(NTdb) <- NTdb[1,]
NTdb <- NTdb[-1,]
NTdb <- data.frame(cbind(Sub_ID = rownames(NTdb),
                         NTdb))
MOMcomb <- left_join(MOMcomb,NTdb, "Sub_ID")
MOMcomb_sub[,c(101:125)] <- apply(MOMcomb_sub[,c(101:125)],2,function(x) as.numeric(x))

plot_list = list()
for (i in c(101:125)) {
  my.data <- data.frame(x = MOMcomb_sub$BMI[MOMcomb_sub$Sex=="M"], 
                        y = MOMcomb_sub[,i][MOMcomb_sub$Sex=="M"])
  j <- colnames(MOMcomb_sub)[i]
  a <- cor.test(my.data$x,
                my.data$y,method = "kendall")
  if(a$p.value<0.05){
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      scale_color_manual(values=c("lightslateblue")) +
      geom_point(size = 2) +
      theme_classic()+
      ylab(paste("Annual estimated ",j," level")) +
      xlab("Weight") +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25) +
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
      scale_color_manual(values=c("lightslateblue")) +
      geom_point(size = 2) +
      theme_classic() +
      ylab(paste("Annual estimated ",j," level")) +
      xlab("Weight") +
      stat_smooth(method = "glm",formula = y~x,  
                  se= F, size = 1.3, alpha = 0.25) +
      theme(legend.position = "none",
            axis.text = element_text(colour = "black")) + 
      annotate("text", x = max(my.data$x,na.rm = T)*3/4,
               y = max(my.data$y,na.rm = T), size = 4, 
               label = paste0("R =", round(a$estimate,2),
                              "\np-value = ", round(a$p.value,3)))
    
    plot_list <- c(plot_list, list(p))
    
  }
}

wrap_plots(plot_list,ncol = 5)


ggsave("Nutrients_BMI_Male.eps",width = 30, height = 30,units = "in")
wrap_plots(plot_list,ncol = 5)
invisible(dev.off())

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
