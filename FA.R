library("reshape2")
library("Hmisc")
library("car")
library("dplyr")
library("survival")
library("survminer")
library("tableone")
library("readxl")
library("tidyverse")
library("lubridate")
library("jtools")
library("pROC")
library("psych")
library("ggplot2")
library("gridExtra")
library("corrplot")
library("patchwork")
library("export")
library("rgl")
library("pander")
library("stringi")
library("MicrobiotaProcess")

# Define High/Low FAC group (Using both fa)----
MomFA <- data.frame(read_excel("MOMmy-DietFoodDAT2_DATA_2021.03.24 WL.xlsx", 
                    sheet = "Auunal Additive in mg"))[,-12]
colnames(MomFA)[1] <- "Sub_ID"
MomFA <- cbind(Sub_ID = MomFA[,1],
               Sex = substr(MomFA$Sub_ID,4,4),
               MomFA[,-1])

MomFA$Sex <- ifelse(MomFA$Sex==1,'F',"M")

for (i in 3:12){
  MomFA[,i+10] <- round(MomFA[,i]/median(MomFA[,i])*100,1)
  colnames(MomFA)[i+10] <- paste0(colnames(MomFA)[i],"_RP")
}

for (i in c(2:10)){
  MomFA[,i+38] <- round(MomFA[,i]/median(MomFA[,i])*100,1)
  colnames(MomFA)[i+38] <- paste0(colnames(MomFA)[i],"_RP")
}

for (i in c(37:39)){
  MomFA[,i+12] <- round(MomFA[,i]/median(MomFA[,i])*100,1)
  colnames(MomFA)[i+12] <- paste0(colnames(MomFA)[i],"_RP")
}
MomFA$All_FA_RP <-  round(MomFA$All_FA/median(MomFA$All_FA)*100,1)

for (i in c(40:52)){
  MomFA[,i+13] <- log2(log2(MomFA[,i]))
  colnames(MomFA)[i+13] <- paste0("log2",colnames(MomFA)[i])
}
MomFA$log2All_FA_RP <-  log2(MomFA$All_FA_RP)

# Set cutoffs----

MomFA_F <- MomFA %>% filter(Sex=="F")
for (i in c(2:10)){
  MomFA_F[,i+64] <- ifelse(ntile(MomFA_F[,i],4)==1,"Low",
                           ifelse(ntile(MomFA_F[,i],4)==4,
                                  "High","Normal"))
  colnames(MomFA_F)[i+64] <- paste0(colnames(MomFA_F)[i],"_cut")
}
for (i in c(37:39)){
  MomFA_F[,i+38] <-   ifelse(ntile(MomFA_F[,i],4)==1,"Low",
                             ifelse(ntile(MomFA_F[,i],4)==4,
                                    "High","Normal"))
  colnames(MomFA_F)[i+38] <- paste0(colnames(MomFA_F)[i],"_cut")
}
MomFA_F$All_FA_cut <-  ifelse(ntile(MomFA_F$All_FA,4)==1,"Low",
                              ifelse(ntile(MomFA_F$All_FA,4)==4,
                                     "High","Normal"))
for (i in c(2:10)){
  MomFA_F[,i+77] <- ifelse(ntile(MomFA_F[,i],3)==1,"Low",
                           ifelse(ntile(MomFA_F[,i],3)==3,
                                  "High","Normal"))
  colnames(MomFA_F)[i+77] <- paste0(colnames(MomFA_F)[i],"_cut3")
}
for (i in c(37:39)){
  MomFA_F[,i+51] <-   ifelse(ntile(MomFA_F[,i],3)==1,"Low",
                             ifelse(ntile(MomFA_F[,i],3)==3,
                                    "High","Normal"))
  colnames(MomFA_F)[i+51] <- paste0(colnames(MomFA_F)[i],"_cut3")
}
MomFA_F$All_FA_cut3 <-  ifelse(ntile(MomFA_F$All_FA,3)==1,"Low",
                               ifelse(ntile(MomFA_F$All_FA,3)==3,
                                      "High","Normal"))
# Father
MomFA_M <- MomFA %>% filter(Sex=="M")
for (i in c(2:10)){
  MomFA_M[,i+64] <- ifelse(ntile(MomFA_M[,i],4)==1,"Low",
                           ifelse(ntile(MomFA_M[,i],4)==4,
                                  "High","Normal"))
  colnames(MomFA_M)[i+64] <- paste0(colnames(MomFA_M)[i],"_cut")
}
for (i in c(37:39)){
  MomFA_M[,i+38] <-   ifelse(ntile(MomFA_M[,i],4)==1,"Low",
                             ifelse(ntile(MomFA_M[,i],4)==4,
                                    "High","Normal"))
  colnames(MomFA_M)[i+38] <- paste0(colnames(MomFA_M)[i],"_cut")
}
MomFA_M$All_FA_cut <-  ifelse(ntile(MomFA_M$All_FA,4)==1,"Low",
                              ifelse(ntile(MomFA_M$All_FA,4)==4,
                                     "High","Normal"))
for (i in c(2:10)){
  MomFA_M[,i+77] <- ifelse(ntile(MomFA_M[,i],3)==1,"Low",
                           ifelse(ntile(MomFA_M[,i],3)==3,
                                  "High","Normal"))
  colnames(MomFA_M)[i+77] <- paste0(colnames(MomFA_M)[i],"_cut3")
}
for (i in c(37:39)){
  MomFA_M[,i+51] <-   ifelse(ntile(MomFA_M[,i],3)==1,"Low",
                             ifelse(ntile(MomFA_M[,i],3)==3,
                                    "High","Normal"))
  colnames(MomFA_M)[i+51] <- paste0(colnames(MomFA_M)[i],"_cut3")
}
MomFA_M$All_FA_cut3 <-  ifelse(ntile(MomFA_M$All_FA,3)==1,"Low",
                               ifelse(ntile(MomFA_M$All_FA,3)==3,
                                      "High","Normal"))
MomFA <- rbind(MomFA_F,MomFA_M)

for (i in 66:91){
  MomFA[,i] <- factor(as.factor(MomFA[,i]),
                      levels(as.factor(MomFA[,i]))[c(2,3,1)])
}

MomFA <- MomFA[MomFA$Sub_ID!="PWH100107",]
MomFA <- MomFA[MomFA$Sub_ID!="PWH100046",]



# FA weight adjusted----
MomFA$Weight <- demodb$Weight[match(MomFA$Sub_ID,demodb$Sub_ID)]
MomFA$Weight[is.na(MomFA$Weight)] <- daddemodb$Weight[match(MomFA$Sub_ID[is.na(MomFA$Weight)],
                                                            daddemodb$Sub_ID)]
MomFA_W <- MomFA %>% filter(complete.cases(Weight))

for (i in 3:12){
  MomFA_W[,i+10] <- round(MomFA_W[,i]/MomFA_W$Weight,1)
  colnames(MomFA_W)[i+10] <- paste0(colnames(MomFA_W)[i],"_perKG")
}

MomFA_W$All_RP <- as.numeric(apply(MomFA_W[,c(13:22)],1,function(x) sum(as.numeric(x),na.rm = T)))
colnames(MomFA_W)[23] <- "All_perKG"
MomFA_W <- MomFA_W %>% mutate(Sweet_perKG =  as.numeric(apply(MomFA_W[,c(20:22)],
                                                              1,function(x) sum(as.numeric(x),na.rm = T))),
                          Emu_perKG = as.numeric(apply(MomFA_W[,c(13:15)],
                                                       1,function(x) sum(as.numeric(x), na.rm = T))),
                          PresDisp_perKG = as.numeric(apply(MomFA_W[,c(16:19)],
                                                            1,function(x) sum(as.numeric(x),na.rm = T))),
                          Sweet_abs =  as.numeric(apply(MomFA_W[,c(10:12)],
                                                        1,function(x) sum(as.numeric(x),na.rm = T))),
                          Emu_abs = as.numeric(apply(MomFA_W[,c(3:5)],
                                                     1,function(x) sum(as.numeric(x), na.rm = T))),
                          PresDisp_abs = as.numeric(apply(MomFA_W[,c(6:9)],
                                                          1,function(x) sum(as.numeric(x),na.rm = T))))

MomFA_W$Group <- ifelse(MomFA_W$All_perKG<summary(MomFA_W$All_perKG[MomFA_W$Sex=="F"])[[2]],"Low",
                        ifelse(MomFA_W$All_perKG<summary(MomFA_W$All_perKG[MomFA_W$Sex=="F"])[[5]],
                               "High","Normal"))
MomFA_W$Group[MomFA_W$Sex=="M"] <- ifelse(MomFA_W$All_perKG[MomFA_W$Sex=="M"]<summary(MomFA_W$All_perKG[MomFA_W$Sex=="M"])[[2]],"Low",
                        ifelse(MomFA_W$All_perKG[MomFA_W$Sex=="M"] <summary(MomFA_W$All_perKG[MomFA_W$Sex=="M"])[[5]],
                               "High","Normal"))


pdf(file="Cutoffs_wholedb.pdf", width = 8, height = 6)
for (i in c(3:12,23:26)) {
  my.data <- data.frame(y = MomFA[,i])
  my.data <- data.frame(y = my.data[-tail(order(my.data$y),20),])
  colnames(my.data) <- "y"
  j <- colnames(MomFA)[i]
  p <- gghistogram(data = my.data,x = "y", fill = "#9ecae1") +
    ylab("Frequency") +
    xlab(paste0(j," annual/mg")) +
    geom_vline(aes(xintercept = summary(my.data$y)[2]),
               color="#31a354", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = summary(my.data$y)[5]),
               color="#f03b20", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = tail(my.data$y[order(my.data$y,decreasing = T)],100)[100]), 
               color="#31a354", size = 1) +
    geom_vline(aes(xintercept = head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
               color="#f03b20",  size = 1) +
    scale_x_continuous(breaks = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
                       labels = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]))
  print(p)
  Sys.sleep(2)
print(i)
}
dev.off()

# Mom
momdb <- data.frame(read_excel("MOMmy-DietFoodDAT2_DATA_2021.03.24 WL.xlsx", 
                               sheet = "Auunal Additive in mg"))[,-12]
colnames(momdb)[1] <- "Sub_ID"
momdb <- cbind(Sub_ID = momdb[,1],
               Sex = substr(momdb$Sub_ID,4,4),
               momdb[,-1])

momdb$Sex <- ifelse(momdb$Sex==1,'F',"M")

momdb <- momdb %>% filter(Sex == "F")

for (i in 3:12){
  momdb[,i+10] <- round(momdb[,i]/median(momdb[,i])*100,1)
  colnames(momdb)[i+10] <- paste0(colnames(momdb)[i],"_RP")
}

momdb <- momdb %>% mutate(All_RP = as.numeric(apply(momdb[,c(13:15,17:22)],1,sum)),
                          Sweet_abs =  as.numeric(apply(momdb[,c(8:10)],1,sum)),
                          Emu_abs = as.numeric(apply(momdb[,c(1:3)],1,function(x) sum(as.numeric(x), na.rm = T))),
                          PresDisp_abs = as.numeric(apply(momdb[,c(4:7)],1,sum))) %>% 
  filter(All_RP!=0)

momdb <- momdb[momdb$Sub_ID!="PWH100107",]
momdb <- momdb[momdb$Sub_ID!="PWH100046",]

pdf(file="Cutoffs_mom.pdf", width = 8, height = 6)
for (i in c(3:12,23:26)) {
  my.data <- data.frame(y = momdb[,i])
  my.data <- data.frame(y = my.data[-tail(order(my.data$y),20),])
  colnames(my.data) <- "y"
  j <- colnames(MomFA)[i]
  p <- gghistogram(data = my.data,x = "y", fill = "#9ecae1") +
    ylab("Frequency") +
    xlab(paste0(j," annual/mg")) +
    geom_vline(aes(xintercept = summary(my.data$y)[2]),
               color="#31a354", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = summary(my.data$y)[5]),
               color="#f03b20", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = tail(my.data$y[order(my.data$y,decreasing = T)],100)[100]), 
               color="#31a354", size = 1) +
    geom_vline(aes(xintercept = head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
               color="#f03b20",  size = 1) +
    scale_x_continuous(breaks = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
                       labels = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]))
  print(p)
  Sys.sleep(2)
  print(i)
}
dev.off()

# Dad

daddb <- data.frame(read_excel("MOMmy-DietFoodDAT2_DATA_2021.03.24 WL.xlsx", 
                               sheet = "Auunal Additive in mg"))[,-12]
colnames(daddb)[1] <- "Sub_ID"
daddb <- cbind(Sub_ID = daddb[,1],
               Sex = substr(daddb$Sub_ID,4,4),
               daddb[,-1])

daddb$Sex <- ifelse(daddb$Sex==1,'F',"M")

daddb <- daddb %>% filter(Sex == "M")

for (i in 3:12){
  daddb[,i+10] <- round(daddb[,i]/median(daddb[,i])*100,1)
  colnames(daddb)[i+10] <- paste0(colnames(daddb)[i],"_RP")
}

daddb <- daddb %>% mutate(All_RP = as.numeric(apply(daddb[,c(13:15,17:22)],1,sum)),
                          Sweet_abs =  as.numeric(apply(daddb[,c(8:10)],1,sum)),
                          Emu_abs = as.numeric(apply(daddb[,c(1:3)],1,function(x) sum(as.numeric(x), na.rm = T))),
                          PresDisp_abs = as.numeric(apply(daddb[,c(4:7)],1,sum))) %>% 
  filter(All_RP!=0)

pdf(file="Cutoffs_dad.pdf", width = 8, height = 6)
for (i in c(3:12,23:26)) {
  my.data <- data.frame(y = daddb[,i])
  my.data <- data.frame(y = my.data[-tail(order(my.data$y),20),])
  colnames(my.data) <- "y"
  j <- colnames(MomFA)[i]
  p <- gghistogram(data = my.data,x = "y", fill = "#9ecae1") +
    ylab("Frequency") +
    xlab(paste0(j," annual/mg")) +
    geom_vline(aes(xintercept = summary(my.data$y)[2]),
               color="#31a354", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = summary(my.data$y)[5]),
               color="#f03b20", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = tail(my.data$y[order(my.data$y,decreasing = T)],100)[100]), 
               color="#31a354", size = 1) +
    geom_vline(aes(xintercept = head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
               color="#f03b20",  size = 1) +
    scale_x_continuous(breaks = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]),
                       labels = c(tail(my.data$y[order(my.data$y,decreasing = T)],100)[100],
                                  summary(my.data$y)[2],
                                  summary(my.data$y)[5],
                                  head(my.data$y[order(my.data$y,decreasing = T)],100)[100]))
  print(p)
  Sys.sleep(2)
  print(i)
}
dev.off()

# Print overall mom and dad FA IQR
i = 23
my.data <- data.frame(y = momdb[,i])
my.data <- data.frame(y = my.data[-tail(order(my.data$y),20),])
colnames(my.data) <- "y"
j <- colnames(MomFA)[i]
p <- gghistogram(data = my.data,x = "y", fill = "#9ecae1") +
  ylab("Frequency") +
  xlab(paste0("Mother RAFAC")) +
  geom_vline(aes(xintercept = summary(my.data$y)[2]),
             color="#31a354", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = summary(my.data$y)[5]),
             color="#f03b20", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(summary(my.data$y)[2],
                                summary(my.data$y)[5]),
                     labels = c(summary(my.data$y)[2],
                                summary(my.data$y)[5]))

p
graph2eps(file="Mom_RAFAC_08042021.eps", width = 8, height = 6)

i = 23
my.data <- data.frame(y = daddb[,i])
my.data <- data.frame(y = my.data[-tail(order(my.data$y),20),])
colnames(my.data) <- "y"
j <- colnames(MomFA)[i]
p <- gghistogram(data = my.data,x = "y", fill = "#9ecae1") +
  ylab("Frequency") +
  xlab(paste0("Father RAFAC")) +
  geom_vline(aes(xintercept = summary(my.data$y)[2]),
             color="#31a354", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = summary(my.data$y)[5]),
             color="#f03b20", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(summary(my.data$y)[2],
                                summary(my.data$y)[5]),
                     labels = c(summary(my.data$y)[2],
                                summary(my.data$y)[5]))

p
graph2eps(file="Dad_RAFAC_08042021.eps", width = 8, height = 6)



mhighlist <- substr(momdb$Sub_ID[momdb$All_RP>2625],7,9)
dhighlist <- substr(daddb$Sub_ID[daddb$All_RP>2223.4],7,9)
bothhighlist <- intersect(mhighlist,dhighlist)

bothhigh <- MomFA %>% mutate(Fam_ID = substr(Sub_ID,7,9)) %>% 
  filter(Fam_ID %in% bothhighlist)

bothhigh$FAGroup <- "High"

mlowlist <- substr(momdb$Sub_ID[momdb$All_RP<694],7,9)
dlowlist <- substr(daddb$Sub_ID[daddb$All_RP<625.9],7,9)
bothlowlist <- intersect(mlowlist,dlowlist)

bothlow <- MomFA %>% mutate(Fam_ID = substr(Sub_ID,7,9)) %>% 
  filter(Fam_ID %in% bothlowlist)

bothlow$FAGroup <- "Low"
Seldb <- rbind(bothhigh,bothlow)

Seldb$Sweet_log <- log2(Seldb$Sweet_abs)

Sweetbp <- ggstripchart(Seldb, x = "FAGroup", y = "Sweet_log", size = 3,
             color = "FAGroup", add = c("boxplot"),
             add.params = list(size = 1.5,
                               linetype = "solid", shape="0.5"))+ 
  stat_compare_means(label.x = 0.5,
                     method = "wilcox", size = 6, ) +
  # stat_compare_means(comparisons = compgroup,method = "wilcox.test") +
  scale_color_manual(values=c("#f03b20","#31a354"))  +
  theme_classic() +
  theme(axis.text = element_text(color = "Black", size = 22),
    axis.line = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    axis.ticks.length = unit(.15, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28),
    legend.position = "none") +
  # ylim(0,35) +
  # scale_y_continuous(limits = c(0,400),breaks = c(0,100,200,300,400)) +
  # ggtitle("Samples collected within 14 days of symptom onset")+
  ylab("Log2 Sweetener annual consumption (mg)") +
  scale_x_discrete(labels=c("High (N = 47)",
                            "Low (N = 36)")) 

Seldb$Emu_log <- log2(Seldb$Emu_abs)
Emubp <- ggstripchart(Seldb, x = "FAGroup", y = "Emu_log", size = 3,
             color = "FAGroup", add = c("boxplot"),
             add.params = list(size = 1.5,
                               linetype = "solid", shape="0.5"))+ 
  stat_compare_means(label.x = 0.5,
                     method = "wilcox", size = 6) +
  # stat_compare_means(comparisons = compgroup,method = "wilcox.test") +
  scale_color_manual(values=c("#f03b20","#31a354"))  +
  theme_classic() +
  theme(axis.text = element_text(color = "Black", size = 22),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(.15, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28),
        legend.position = "none") +
  # ylim(0,35) +
  # scale_y_continuous(limits = c(0,400),breaks = c(0,100,200,300,400)) +
  # ggtitle("Samples collected within 14 days of symptom onset")+
  ylab("Log2 Emulsifier annual consumption (mg)") +
  scale_x_discrete(labels=c("High (N = 47)",
                            "Low (N = 36)"))

Seldb$PresDisp_log <- log2(Seldb$PresDisp_abs)
PoDbp <- ggstripchart(Seldb, x = "FAGroup", y = "PresDisp_log", size = 3,
             color = "FAGroup", add = c("boxplot"),
             add.params = list(size = 1.5,
                               linetype = "solid", shape="0.5"))+ 
  stat_compare_means(label.x = 0.5,
                     method = "wilcox", size = 6) +
  # stat_compare_means(comparisons = compgroup,method = "wilcox.test") +
  scale_color_manual(values=c("#f03b20","#31a354"))  +
  theme_classic() +
  theme(axis.text = element_text(color = "Black", size = 22),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(.15, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28),
        legend.position = "none") +
  # ylim(0,35) +
  # scale_y_continuous(limits = c(0,400),breaks = c(0,100,200,300,400)) +
  # ggtitle("Samples collected within 14 days of symptom onset")+
  ylab("Log2 Preservative/Dispersant annual consumption (mg)") +
  scale_x_discrete(labels=c("High (N = 47)",
                            "Low (N = 36)")) 

setEPS(width = 32, height = 12)
postscript("FA_bp_08042021.eps",onefile = FALSE)
Sweetbp|Emubp|PoDbp
invisible(dev.off())

# Define High/Low FAC group (Using bootstrap)----
# Build up matrix
R = 1000                    # number of bootstrap samples
n = nrow(momdb)        # sample size
k = 2                       # number of rows

# set up a empty Rxn matrix B
IQRtest  <-  matrix(nrow = R, ncol = k,
                   dimnames = list(paste("Sample",1:R), 
                                   c("LIQR","UIQR")))

# loop R times

# install.packages("sprof", repos="http://R-Forge.R-project.org")  
 #for stri_list2matrix

set.seed(5)
for(i in 1:R){
  # sample credit data with replacement
  boot.data = MomFA[sample(x = 1:n, size = 200, replace = TRUE),]
  # summary "TP","FP","TN","FN" for sampling
  sumIQR <- paste( summary(boot.data$All_RP)[[2]],
                   summary(boot.data$All_RP)[[5]],sep = ",")
  IQRtest[i,] <- stri_list2matrix(strsplit(sumIQR,","),byrow = T)
}

IQRtest <- data.frame(IQRtest)
IQRtest[,c(1:2)] <- apply(IQRtest[,c(1:2)],2,function(x) as.numeric(x))
hist(IQRtest$LIQR)
mean(IQRtest$LIQR)
mean(IQRtest$UIQR)


# Father cut off
# Build up matrix
R = 1000                    # number of bootstrap samples
n = nrow(momdb)        # sample size
k = 2                       # number of rows

# set up a empty Rxn matrix B
IQRtest  <-  matrix(nrow = R, ncol = k,
                    dimnames = list(paste("Sample",1:R), 
                                    c("LIQR","UIQR")))

# loop R times

# install.packages("sprof", repos="http://R-Forge.R-project.org")  
#for stri_list2matrix

set.seed(5)
for(i in 1:R){
  # sample credit data with replacement
  boot.data = daddb[sample(x = 1:n, size = 200, replace = TRUE),]
  # summary "TP","FP","TN","FN" for sampling
  sumIQR <- paste( summary(boot.data$All_RP)[[2]],
                   summary(boot.data$All_RP)[[5]],sep = ",")
  IQRtest[i,] <- stri_list2matrix(strsplit(sumIQR,","),byrow = T)
}

IQRtest <- data.frame(IQRtest)
IQRtest[,c(1:2)] <- apply(IQRtest[,c(1:2)],2,function(x) as.numeric(x))
hist(IQRtest$LIQR)
mean(IQRtest$LIQR)
mean(IQRtest$UIQR)
dhighlist <- substr(daddb$Sub_ID[daddb$All_RP>mean(IQRtest$UIQR)],7,9)
dlowlist <- substr(daddb$Sub_ID[daddb$All_RP<mean(IQRtest$LIQR)],7,9)

# wilcox test of different groups at different timepoints (Drop)----
PDHdb <- left_join(PDHdb,MomFA[,c(1,22:26)], "Sub_ID")

wiltest <- lapply(PDHdb[,c(4:7,9:40)], function(x) wilcox.test(PDHdb$All_RP ~ x))
wiltest[[i]][["p.value"]]

PDHbygroup <- data.frame()
for (i in c(4:7,9:40)) {
  temp.db <- data.frame(cbind(Var = colnames(PDHdb[i]), Yes = median(PDHdb$All_RP[PDHdb[,i]=="Yes"], na.rm = T), 
                              No = median(PDHdb$All_RP[PDHdb[,i]=="No"], na.rm = T)))
  PDHbygroup <- rbind(PDHbygroup,temp.db)
}

for (i in unique(PDHbygroup$Var)) {
PDHbygroup$p.val[PDHbygroup$Var==i] <- wiltest[[i]][["p.value"]]
}

PDHbygroup <- melt(PDHbygroup,id.vars = c("Var","p.val"))
PDHbygroup$Item <- as.character(lapply(as.list(PDHbygroup$Var),function(x) strsplit(x,"_")[[1]][1]))
PDHbygroup$Timepoint <- as.character(lapply(as.list(PDHbygroup$Var),function(x) strsplit(x,"_")[[1]][2]))
PDHbygroup$Timepoint <- factor(as.factor(PDHbygroup$Timepoint),
                           levels(as.factor(PDHbygroup$Timepoint))[c(3,1,4,2)])
colnames(PDHbygroup)[3:4] <- c("Yes_no","All_RP_median")
PDHbygroup$All_RP_median <- as.numeric(PDHbygroup$All_RP_median)

plot_list = list()
for (i in  unique(PDHbygroup$Item)) {
  temp.db <- PDHbygroup %>% filter(Item == i)
p <- ggplot(data = temp.db, aes(x = Timepoint,
                           y = All_RP_median, color = Yes_no)) +
  scale_color_manual(values=c("#993404","#7fcdbb")) +
  geom_point(aes(color = Yes_no, size = p.val)) +
  geom_path(aes(group = Yes_no)) +
  ggtitle(paste(i)) +
  guides(color = FALSE) +
  geom_hline(aes(yintercept = 1096.65), 
             color="black", size = 0.5, linetype = "dashed") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "bottom")+
  ylab("Estimated FA (median)")

plot_list <- c(plot_list, list(p))
}

pleg <- ggplot(data = temp.db, aes(x = Timepoint,
                                   y = All_RP_median, color = Yes_no)) +
  scale_color_manual(values=c("#993404","#7fcdbb")) +
  geom_point(aes(color = Yes_no)) +
  geom_path(aes(group = Yes_no)) +
  theme(legend.title = element_blank())

plot_list <- c(plot_list, list(as_ggplot(get_legend(pleg))))

setEPS(width = 18, height = 10)
postscript("FA_in groups_30042021.eps",onefile = FALSE)
wrap_plots(plot_list,ncol = 4)
invisible(dev.off())


plot_list = list()
for (i in  unique(PDHbygroup$Item)) {
  temp.db <- cbind(All_RP_median=PDHdb$All_RP,
                   PDHdb[,grepl(i,colnames(PDHdb))])
  temp.db <- melt(temp.db,id.vars = "All_RP_median")
  temp.db$Item <- as.character(lapply(as.list(as.character(temp.db$variable)),
                                      function(x) strsplit(x,"_")[[1]][1]))
  temp.db$Timepoint <- as.character(lapply(as.list(as.character(temp.db$variable)),
                                      function(x) strsplit(x,"_")[[1]][2]))
  temp.db$Timepoint <- factor(as.factor(temp.db$Timepoint),
                                 levels(as.factor(temp.db$Timepoint))[c(3,1,4,2)])
  colnames(temp.db)[3] <- "Yes_no"
p <- ggplot(data = temp.db, aes(x = Timepoint,
                                  y = log2(All_RP_median), color = Yes_no)) +
    scale_color_manual(values=c("#993404","#7fcdbb")) +
    geom_boxplot(aes(color = Yes_no)) +
    ggtitle(paste(i)) +
    stat_compare_means(method = "wilcox", comparisons = temp.db$Yes_no) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          legend.position = "none")+
    ylab("Log2 estimated FA (median)")
  
  plot_list <- c(plot_list, list(p))
}

pleg <- ggplot(data = temp.db, aes(x = Timepoint,
                                   y = All_RP_median, color = Yes_no)) +
  scale_color_manual(values=c("#993404","#7fcdbb")) +
  geom_point(aes(color = Yes_no)) +
  geom_path(aes(group = Yes_no)) +
  theme(legend.title = element_blank())

plot_list <- c(plot_list, list(as_ggplot(get_legend(pleg))))


setEPS(width = 18, height = 10)
postscript("FA_in groups_bp_30042021.eps",onefile = FALSE)
wrap_plots(plot_list,ncol = 4)
invisible(dev.off())

















# Correlation of father and mother FA
FA_Cor <- data.frame(Fam_ID = substr(PDHdb$Sub_ID,7,9),
                     FA_Mom = MomFA$All_RP[substr(MomFA$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)&MomFA$Sex=="F"],
                     FA_Dad = MomFA$All_RP[substr(MomFA$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)&MomFA$Sex=="M"],
                     BMI_Mom = demodb$BMI[substr(demodb$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)],
                     BMI_Dad = daddemodb$BMI[substr(daddemodb$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)],
                     Weight_T1BL = combdb$Weight_T1BL[substr(combdb$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)],
                     Weight_T3BL = combdb$Weight_T3BL[substr(combdb$Sub_ID,7,9) %in% substr(PDHdb$Sub_ID,7,9)])


ggplot(FA_Cor,  aes(x = log2(FA_Cor$FA_Mom), y = log2(FA_Cor$FA_Dad))) +
  geom_point(size = 2) +
  theme_classic()+
  # ylab(paste("Serum ",j," level")) +
  # xlab("Time from onset") +
  stat_cor(method="spearman",position = "identity",
           r.digits = 2,p.digits = 3) +
  stat_smooth(method = "glm",formula = y~x,  
              se= F, size = 1.3, alpha = 0.25) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black")) #,plot.background = element_rect(fill = "white", colour = "red3", size = 1)

ggplot(FA_Cor,  aes(x = log2(FA_Mom), y = Weight_T1BL)) +
  geom_point(size = 2) +
  theme_classic()+
  # ylab(paste("Serum ",j," level")) +
  # xlab("Time from onset") +
  stat_cor(method="pearson",position = "identity",
           r.digits = 2,p.digits = 3) +
  stat_smooth(method = "glm",formula = y~x,  
              se= F, size = 1.3, alpha = 0.25) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"))

# Try Bar distribution of FA 


momdb_t <- data.frame(t(momdb[,c(3:12)]))
colnames(momdb_t) <- momdb$Sub_ID
momdb_t <- apply(momdb_t, 2, function(x) round(x/sum(x)*100,1))
momdb_t <- melt(momdb_t)
momdb_t$Var2 <- substr(momdb_t$Var2,7,9)
momdb_t <- momdb_t %>% filter(Var2 %in% mhighlist)

tiff("Mon_FAtype.tiff", height = 8, width = 20, units = 'in', res=200)
ggplot(data=momdb_t, aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", angle = 45, hjust = 0.5))
invisible(dev.off())



daddb_t <- data.frame(t(daddb[,c(3:12)]))
colnames(daddb_t) <- daddb$Sub_ID
daddb_t <- apply(daddb_t, 2, function(x) round(x/sum(x)*100,1))
daddb_t <- melt(daddb_t)
daddb_t$Var2 <- substr(daddb_t$Var2,7,9)
daddb_t <- daddb_t %>% filter(Var2 %in% dhighlist)


tiff("Dad_FAtype.tiff", height = 8, width = 20, units = 'in', res=200)
ggplot(data=daddb_t, aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", angle = 45, hjust = 0.5))
invisible(dev.off())

ggplot(data=combdb[combdb$Group!="Normal",], aes(x=Group, y=Weight_T3BL)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(color = "black", angle = 45, hjust = 1)) +
  stat_compare_means(method = "t.test")


WeightCh <- aggregate(combdb[,c(6,11,12)], combdb[c("Group")], FUN = function(x) mean(x, na.rm = T))
WeightCh <- melt(combdb[,c(6,11,12,54)])
WeightCh <- WeightCh %>% filter(value<120)

ggplot(data=WeightCh[WeightCh$Group!="Normal",], aes(x=variable, y=value,color = Group)) +
  # geom_line() +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) +
  stat_compare_means(method = "t.test") 

 cWeightCh$variable <- ifelse(WeightCh$variable=="Weight",0,
                            ifelse(WeightCh$variable=="Weight_T1",1,2))
 
ggplot(data=WeightCh, aes(x=variable, y=value, color = Group)) +
  scale_color_manual(values=c("#756bb1","#e34a33","#fdbb84")) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))+
  # ylab("Serum HMGB1 level (ng/mL)") +
  # xlab("Serum IL-1β level (pg/mL)") +
  guides(fill=guide_legend(title="") )+
  stat_cor(method="kendal",position = "identity",
           r.digits = 2,p.digits = 3) +
  stat_smooth(method = "glm",formula = y~x,  
              se= F, alpha = 0.25) 

MOM_GDM <- data.frame(read_xlsx("/media/gee/Data/REDCap/MOMmy-Mother T2 T3 Q_20210716_Wing_20210929_Shilan_20211006.xlsx"))[,c(2,6,30:32)]
colnames(MOM_GDM) <- c("Sub_ID","Tri","Disease","Disease_type","Disease_other")
MOM_GDM$Tri <- ifelse(MOM_GDM$Tri=="第三孕期","T3","T2")
MOM_GDM <- MOM_GDM %>% filter(complete.cases(Tri))
MOM_GDM$DM <- ifelse(MOM_GDM$Disease_type=="糖尿病","Yes","No")
MOM_GDM$DM[grepl("糖尿",MOM_GDM$Disease_other)] <- "Yes"
MOM_GDM$DM[is.na(MOM_GDM$DM)] <- "No"
MOM_GDM <- reshape(MOM_GDM[,-(3:5)],timevar = "Tri",idvar = "Sub_ID",direction = "wide")
colnames(MOM_GDM)[2:3] <- c("T2_DM","T3_DM")
MOMcomb_sub <- merge(MOMcomb_sub,MOM_GDM,"Sub_ID",all.x = TRUE)


