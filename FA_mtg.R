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

# Creating phyloseq db----
FAMT3OTU <- readRDS("Data/Year 1/MT3_M12.rds")
FAMT3OTU <-  subset_samples(FAMT3OTU, grepl("Mother T", Time_point))
FAMT3Taxa <- data.frame(FAMT3OTU@tax_table)
FAMT3Meta <- data.frame(FAMT3OTU@sam_data)
FAMT3Meta <- FAMT3Meta[,c(4,2,8:35,38,39,41,42,44:46,48:50,113:119,121)]
FAMT3Meta <- left_join(FAMT3Meta,MomFA[,c(1:10,37:39,12,66:78)],"Sub_ID")
FAMT3Meta <- FAMT3Meta %>% filter(complete.cases(All_FA))
FAMT3Meta <- left_join(FAMT3Meta,MOMcomb_sub[,c(86,87,1)],"Sub_ID")
FAMT3Meta <- left_join(FAMT3Meta,PNdb[,c(2,106:108,128)],"Sub_ID")
FAMT3Meta$PT <- ifelse(is.na(PT_db$Sub_ID_Child[match(
  substr(FAMT3Meta$Sub_ID,6,9),
  substr(PT_db$Sub_ID_Child,6,9)
)]),"No","Yes")
FAMT3Meta <- FAMT3Meta[!duplicated(FAMT3Meta$Sub_ID),]
FAMT3Meta$T3_WC_Group <- T3$IOM_GWG[match(FAMT3Meta$Sub_ID,
                                          T3$Sub_ID)]
FAMT3Meta$T3_WC_Over <- ifelse(FAMT3Meta$T3_WC_Group=="Over","Over","N_or_Un")
FAMT3Meta$T3_WC_Under <- ifelse(FAMT3Meta$T3_WC_Group=="Under","Under","N_or_Ov")
rownames(FAMT3Meta) <- colnames(FAMT3OTU_sel)[match(
  substr(FAMT3Meta$Sub_ID,7,9),
  substr(colnames(FAMT3OTU_sel),1,3)
)]


FAMT3OTU <- data.frame(FAMT3OTU@otu_table)
colnames(FAMT3OTU) <- str_replace(colnames(FAMT3OTU),"X","")
colnames(FAMT3OTU) <- str_replace(colnames(FAMT3OTU),"\\.","-")
FAMT3OTU <- FAMT3OTU %>% select(colnames(FAMT3OTU)[colnames(FAMT3OTU) %in% FAMT3Meta$sample_name])

FAMT3OTU_sel <- FAMT3OTU %>% filter(grepl("k__Bacteria",rownames(FAMT3OTU)))
FAMT3OTU_sel <- apply(FAMT3OTU_sel, 2, function(x)
  round(x/x[1]*100,2))

FAMT3Meta[,c(62:74)] <- apply(FAMT3Meta[,c(62:74)], 2, function(x)
  ifelse(x=="High","High","Normal"))
rownames(FAMT3Meta) <- FAMT3Meta$sample_name
rownames(FAMT3OTU) <- FAMT3Taxa$Species
rownames(FAMT3Taxa) <- FAMT3Taxa$Species

FAmtgraw <- phyloseq(otu_table(as.matrix(FAMT3OTU_sel), taxa_are_rows = TRUE), 
                     sample_data(FAMT3Meta), 
                     phyloseq::tax_table(as.matrix(FAMT3Taxa)))
FAmtgraw <- prune_taxa(taxa_sums(FAmtgraw) > 0, FAmtgraw)
FAMT3OTU <- data.frame(FAmtgraw@otu_table)

# Select bacteria
# Y1bacmtg <- prune_taxa(grepl("Bac",tax_table(Y1mtgraw)[, "k"]), Y1mtgraw)
# 
# # Select species level and filter out the species <60% in the population, 
# # and recalculate relative abundance
# Y1bacmtgsp <- prune_taxa(complete.cases(tax_table(Y1bacmtg)[, "s"]), Y1bacmtg)
# 
# # spe.sum = tapply(taxa_sums(Y1bacmtgsp), 
# #                  tax_table(Y1bacmtgsp)[, "s"],
# #                  sum, na.rm=TRUE) # Sort by abundance 
# # freqspe = names(sort(spe.sum, TRUE))[1:round(length(spe.sum)*0.8)] # Top 80%
# 
# spe.sum = tapply(apply(Y1bacmtgsp@otu_table,1,function(x) sum(x==0)), 
#                                   tax_table(Y1bacmtgsp)[, "s"],
#                                   sum, na.rm=TRUE) 
# 
# freqspe = names(sort(spe.sum[spe.sum<(0.4*length(spe.sum))], TRUE))
# 
# 
# Y1bacmtgsp = prune_taxa((tax_table(Y1bacmtgsp)[, "s"] %in% freqspe),
#                         Y1bacmtgsp)
# 
# Y1bacmtgsp@otu_table <- otu_table(apply(Y1bacmtgsp@otu_table, 2, function(x)
#   round(x/sum(x)*100,2)), taxa_are_rows = TRUE)
# rownames(Y1bacmtgsp@tax_table@.Data) <- tax_table(Y1bacmtgsp)[, "s"]

# Create original OTU----
Y1mtgraw <- readRDS("Data/Year 1/MT3_M12.rds")
FAMT3Meta <- data.frame(Y1mtgraw@sam_data)
FAMT3OTU <- data.frame(read.table("Data/Year 1/MT3_merged_abundance_20220416.txt",
                       skip = 1, header = T, sep = "\t"))
FAMT3Taxa <- data.frame(matrix(unlist(lapply(strsplit(FAMT3OTU$clade_name,"\\|"), function(x)
  if(length(x)<7){
    c(paste0(x),rep(NA,7-length(x)))
  }
  else{paste0(x)})), 
  ncol = 7, byrow = TRUE))
rownames(FAMT3OTU) <- FAMT3OTU$clade_name 
FAMT3OTU <- FAMT3OTU[,-(1:2)]


FAMT3Meta <- FAMT3Meta[,c(4,2,8:35,38,39,41,42,44:46,48:50,113:119,121)]
FAMT3Meta <- FAMT3Meta %>% filter(startsWith(Sub_ID,"PWH1")) %>% 
  filter(complete.cases(T1_Date))
FAMT3Meta <- left_join(FAMT3Meta,MomFA[,c(1:10,37:39,12,66:78)],"Sub_ID")
FAMT3Meta <- FAMT3Meta %>% filter(complete.cases(All_FA))
FAMT3Meta$sample_name <- str_replace_all(FAMT3Meta$sample_name,"\\-","_")
colnames(FAMT3OTU) <- str_replace(colnames(FAMT3OTU),"X","")
colnames(FAMT3OTU) <- str_replace_all(colnames(FAMT3OTU),"\\.","_")
FAMT3OTU <- FAMT3OTU %>% select(colnames(FAMT3OTU)[which(colnames(FAMT3OTU) %in%
                                                           FAMT3Meta$sample_name)])
FAMT3Meta <- FAMT3Meta %>% filter(sample_name %in% colnames(FAMT3OTU))

FAMT3Meta[,c(62:74)] <- apply(FAMT3Meta[,c(62:74)], 2, function(x)
  ifelse(x=="High","High","Normal"))
rownames(FAMT3Taxa) <- rownames(FAMT3OTU)
rownames(FAMT3Meta) <- FAMT3Meta$sample_name
colnames(FAMT3Taxa) <- c("Kingdom", "Phylum",  "Class", "Order",
                         "Family",  "Genus", "Species")

FAmtgraw <- phyloseq(otu_table(FAMT3OTU, taxa_are_rows = TRUE), 
                     sample_data(FAMT3Meta), 
                     tax_table(as.matrix(FAMT3Taxa)))

FAmtgraw <- prune_taxa(taxa_sums(FAmtgraw) > 0, FAmtgraw)
FAmtgraw <- prune_taxa(grepl("Bac",tax_table(FAmtgraw)[, "Kingdom"]), FAmtgraw)

FAmtgraw@otu_table <- otu_table(apply(FAmtgraw@otu_table, 2, function(x)
  round(x/x[1]*100,2)), taxa_are_rows = TRUE)



# α div----

FAalpha <- cbind(FAMT3Meta,
                 Richness_spe = t(estimateR(round(t(FAmtgraw@otu_table))))[,1],
                 Richness_Chao1 = t(estimateR(round(t(FAmtgraw@otu_table))))[,2],
                 Evenness = diversity(t(FAmtgraw@otu_table)) / log(specnumber(t(FAmtgraw@otu_table))),
                 Shannon = diversity(t(FAmtgraw@otu_table), index = "shannon"),
                 InvSimp = diversity(t(FAmtgraw@otu_table), index = "invsimpson")) 
# FAalpha[,c(62:74)] <- apply(FAalpha[,c(62:74)], 2, function(x)
#   ifelse(x=="High","High","Normal"))
# my.comp <- list( c("Low","High"),
#                  c("Low","Normal"),
#                 c("Normal","High"))
plot_list = list()
for (i in c(62:74)) {
  for(j in c(75:79)){
    my.data <- data.frame(x = FAalpha[,i], 
                          y = FAalpha[,j],
                          lab = colnames(FAalpha)[j])
    x.lab <- str_replace(colnames(FAalpha)[i],"_cut3"," Group")
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_violin(aes(color = x), size = 1,  position = "dodge") +
      geom_point(size = 1,  aes(color = x)) +
      geom_boxplot(width = 0.15, fill="white",aes(color = x)) +
      scale_color_manual(values=c("#BC3C29FF","#00A087FF")) +
      scale_fill_manual(values=c("#BC3C29FF","#00A087FF")) +
      theme_bw() +      
      stat_compare_means(method = "t.test") +
      stat_signif(comparisons = my.comp, label = "p.signif", 
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

pdf("Results/FA_mtg/FA_alpha_sum_11052022.pdf", height = 45, width = 16)
wrap_plots(plot_list,ncol = 5)
invisible(dev.off())

# plot emulsifier and overall FA only

plot_list = list()
for (i in c(71,74)) {
  for(j in c(75:79)){
    my.data <- data.frame(x = FAalpha[,i], 
                          y = FAalpha[,j],
                          lab = colnames(FAalpha)[j])
    # x.lab <- str_replace(colnames(FAalpha)[i],"_cut3"," Group")
    p <- ggplot(my.data,  aes(x = x, y = y)) +
      geom_violin(aes(color = x), size = 1,  position = "dodge") +
      geom_point(size = 1,  aes(color = x)) +
      geom_boxplot(width = 0.15, fill="white",aes(color = x)) +
      scale_color_manual(values=c("#BC3C29FF","#00A087FF")) +
      scale_fill_manual(values=c("#BC3C29FF","#00A087FF")) +
      theme_bw() +      
      stat_compare_means(method = "kruskal") +
      # stat_signif(comparisons = my.comp, label = "p.signif", 
      #             size = 0,
      #             map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, " "=2),
      #             margin_top = 0.05,
      #             method = "wilcox") +
      # xlab(paste(x.lab)) +
      theme(axis.title= element_blank(),
            axis.text = element_text(color = "black", size =13),
            legend.position = "none",
            panel.grid.major = element_line(colour = "gray70")) +
      facet_grid(. ~ lab) +
      theme(strip.background = element_rect(fill="gray40",color = "black"),
            strip.text = element_text(size=15, colour="white"))
    
    plot_list <- c(plot_list, list(p))
  }
}

pdf("Results/FA_mtg/FA_alpha_Emu_all.pdf", height = 10, width = 15)
wrap_plots(plot_list,ncol = 5) +          
  plot_annotation(theme = theme(plot.margin = unit(c(0, 0, 0, 1), "cm")))
grid.draw(textGrob("Emulsifier", 
                   gp = gpar(fontsize=16, fontface = "bold"),
                   y = 0.75, x = 0.015, rot = 90)) 
grid.draw(textGrob("Overall food additive", 
                   gp = gpar(fontsize=16, fontface = "bold"),
                   y = 0.25, x = 0.015, rot = 90)) 
invisible(dev.off())



# β div----
set.seed(100)
FAnmdspre <- vegdist(t(FAmtgraw@otu_table),na.rm = T)
FAnmds <- monoMDS(FAnmdspre)

FAnmds <- data.frame(cbind(FAMT3Meta,
                           FAnmds[["points"]]))

FAnmds$MDS1 <- as.numeric(FAnmds$MDS1)
FAnmds$MDS2 <- as.numeric(FAnmds$MDS2)
# gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
#   filter(abs(MDS2)<0.1)

pdf("Results/FA_mtg/NMDS_CRN.pdf", height = 6, width = 6)
ggplot(FAnmds, aes(MDS1, MDS2, color = CRN_cut)) +
  geom_point() +
  scale_color_manual(name = "CRN intake group",
                     values=c("#BC3C29FF","#00A087FF")) +
  stat_ellipse(aes(color = CRN_cut), na.rm = F) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("") +
  annotate("text",label = paste("PERMANOVA R2 = 0.006\np = 0.004"),
           x = 2, y = 1.8, hjust = 1)
invisible(dev.off())

# 
# set.seed(1)
# FAPCoA <- get_pcoa(obj = FAmtgraw,
#                  method = "hellinger",
#                  distmethod = "euclidean",
#                  taxa_are_rows = TRUE)
# 
# 
# ggordpoint(obj = FAPCoA, biplot=F, speciesannot=F,
#                        factorNames=c("Emu_cut3"), ellipse=F) +
#   theme(legend.position = "bottom") +
#   # scale_fill_manual(name = "Ref",
#   #                    label = c("No","Yes"),
#   #                    values=c("#0072B5FF","#20B54EFF")) +
#   # scale_color_manual(name = "Ref",
#   #                    label = c("No","Yes"),
#   #                    values=c("#0072B5FF","#20B54EFF")) +
#   ggplot2::stat_ellipse() 
#   



FA_dis <- phyloseq::distance(FAmtgraw, 
                             method="bray", weighted=T)
set.seed(1)
Uni_per <- data.frame()
for(i in colnames(FAMT3Meta)[c(8,10:28,30,33,62:74)]){
  print(i)
  temp_db <- FAMT3Meta %>% filter(complete.cases(FAMT3Meta[,i]))
  temp_db <-  subset_samples(FAmtgraw, Sub_ID %in% temp_db$Sub_ID)
  temp_dis <- phyloseq::distance(temp_db, 
                                 method="bray", weighted=T)
  temp_db <-  adonis(as.formula(paste0("temp_dis ~ ", i)), 
                     data = data.frame(temp_db@sam_data),
                     permutation=999)
  Uni_per <- rbind(Uni_per,data.frame(temp_db$aov.tab)[1, 5:6])
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

# Differential species----

library(microbiomeMarker)
FAdif_spe <-  run_lefse(
  FAmtgraw1,
  wilcoxon_cutoff = 0.01,
  group = "GDM",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)

# Or
library(coin)
diff_analysis(obj = FAmtgraw1, classgroup = "GDM",
              mlfun = "lda",
              filtermod = "pvalue",
              firstcomfun = "kruskal_test",
              firstalpha = 0.01,
              strictmod = TRUE,
              secondcomfun = "wilcox_test",
              subclmin = 3,
              subclwilc = TRUE,
              secondalpha = 0.05,
              lda = 2)
                      

# Maaslin
set.seed(100)

for(i in colnames(FAMT3Meta)[c(49:74)]){
  if(all(FAMT3Meta[,i]==0)){
    next
  }else{
    Maaslin2(
      input_data = FAMT3OTU_sel,
      input_metadata = FAMT3Meta %>% select("Sub_ID",
                                           rownames(Sig_Uni_per)[c(1:3,8)],
                                           i),
      output = paste0("Results/FA_mtg/",i,"MaAsLin"),
      min_abundance = 0.05,
      min_prevalence = 0.05,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = "LM",
      max_significance = 0.05,
      fixed_effects = i,
      correction = "BH",
      reference = "Normal",
      standardize = FALSE,
      cores = 1)
  }
}




Maas_FA <- data.frame()
for(i in colnames(FAMT3Meta)[c(49:74)]){
  temp_db <- read.table(paste0("Results/FA_mtg/",i,"MaAsLin/significant_results.tsv"),
                        header = T)
  Maas_FA <- rbind(Maas_FA,temp_db)
  }
Maas_FA$feature[grepl("__",Maas_FA$feature)] <- as.character(lapply(strsplit(Maas_FA$feature[grepl("\\__",Maas_FA$feature)],"\\__"),
                                                                    function(x) paste0(x)[length(x)]))


# Maaslin clinical factor

set.seed(100)
for(i in colnames(FAMT3Meta)[c(75:81,83,84)]){
  Maaslin2(
    input_data = FAMT3OTU_sel,
    input_metadata = FAMT3Meta %>% 
      select(colnames(FAMT3Meta)[c(8,7,10,11:13,15)],
             i),
    output = paste0("Results/FA_mtg/Clinic/",i,"MaAsLin"),
    min_abundance = 0.05,
    min_prevalence = 0.05,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.05,
    fixed_effects = c(i,colnames(FAMT3Meta)[c(8,7,12,13,33)]),
    correction = "BH",
    reference = c("Bachelor/degree","full-time/part-time","MoverV"),
    standardize = FALSE,
    cores = 1)
}

Maas_Clinc <- data.frame()
for(i in colnames(FAMT3Meta)[c(75:81,83,84)]){
  temp_db <- read.table(paste0("Results/FA_mtg/Clinic/",i,"MaAsLin/significant_results.tsv"),
                        header = T)
  Maas_Clinc <- rbind(Maas_Clinc,temp_db)
}
# Maas_Clinc <- Maas_Clinc %>% filter(grepl("g__|s__",feature))
Maas_Clinc$feature[grepl("__",Maas_Clinc$feature)] <- as.character(lapply(strsplit(Maas_Clinc$feature[grepl("\\.",Maas_Clinc$feature)],"\\."),
                                                                    function(x) paste0(x)[length(x)]))
Maas_Clinc$Ass_Group <- ifelse(Maas_Clinc$coef<0, "Neg", "Pos")
Maas_Clinc$Var_fea <- paste0(Maas_Clinc$metadata,"-\n",Maas_Clinc$feature)
# Maas_Clinc$Level <- ifelse(grepl("\\_",Maas_Clinc$feature), "Species", "Genus")

Maas_Clinc <- Maas_Clinc[order(Maas_Clinc$coef,decreasing = T),]
Maas_Clinc$Order <- length(Maas_Clinc$feature):1
Maas_Clinc$feature <- reorder(Maas_Clinc$feature,Maas_Clinc$Order)

# Volcano plot for Maaslin results

ggplot(data = Maas_Clinc, 
       aes(x = coef, y = qval, 
           col = Ass_Group, 
           label = Sig_fea)) +
  geom_point() + 
  xlab("Coef.") +
  ylab("q-value") +
  ggtitle("") +
  labs(title = "Maternal outcomes and differential species",
       subtitle = "MaAsLin results",
       caption = paste0("*Adjusted for: ", paste0(colnames(FAMT3Meta)[c(8,7,12,13,33)],
                                        collapse = ", "))) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none",
    plot.caption = element_text(size = 10, face = "italic")
  ) +
  geom_text_repel() +
  scale_color_manual(values=c("#74add1", "#d73027"),
                     name = "Association",
                     labels = c("Negative", "Positive")) +
  geom_vline(xintercept = c(0), col="black") +
  geom_hline(yintercept = 0.05, col="black")


# Bar plot for Maaslin results

pdf("Results/FA_mtg/MO_MaAsLin_Sum.pdf", height = 9, width = 9)
ggplot(data = Maas_Clinc, 
       aes(x = coef, y = Var_fea, 
           fill = Ass_Group)) +
  geom_col(width = 0.75, position = "dodge") + 
  scale_fill_manual(values=c("#74add1", "#d73027"),
                     name = "Association",
                     labels = c("Negative", "Positive")) +
  xlab("Coef.") +
  ggtitle("") +
  labs(title = "Maternal outcomes and differential species",
       subtitle = "MaAsLin results",
       caption = paste0("*Adjusted for: ", paste0(colnames(FAMT3Meta)[c(8,7,12,13,33)],
                                                  collapse = ", "))) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none",
    plot.caption = element_text(size = 10, face = "italic")
  ) 
invisible(dev.off())
# ANCOM
library(ANCOMBC)
ancombc(phyloseq = FAmtgraw, 
        formula = "P80_cut", 
        p_adj_method = "fdr", 
        zero_cut = 0.90, 
        lib_cut = 1000, 
        group = "P80_cut",
        struc_zero = TRUE, 
        neg_lb = TRUE, 
        tol = 1e-5, 
        max_iter = 100, 
        conserve = TRUE, 
        alpha = 0.05,
        global = FALSE)


# Edger
Edger_FA <- data.frame()
for(i in colnames(FAMT3Meta)[c(62:74)]){
  temp.db <- run_edger(
    FAmtgraw,
    group = i,
    taxa_rank = "Genus",
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
  )
 if(is.null(temp.db@marker_table)){
   next
 }else{
   temp.db <- data.frame(cbind(Var = i,
                               temp.db@marker_table))
   Edger_FA <- data.frame(rbind(Edger_FA,temp.db))
   print(i)
 }
}



# Food additive and bacteria correlation----
# All level OTU
FAMT3OTU <- data.frame(FAmtgraw@otu_table)
FAMT3OTU <- apply(FAMT3OTU)

FAMT3OTU_sel <- FAMT3OTU %>% mutate(sum0 = as.numeric(apply(FAMT3OTU,1,
                                          function(x) sum(x==0)))) %>% 
  filter(sum0<0.95*nrow(FAMT3Meta)) %>% 
  select(-sum0)


# FAMT3OTU_sel <- data.frame(prune_taxa(complete.cases(tax_table(FAmtgraw)[, "Genus"])&
#                              !complete.cases(tax_table(FAmtgraw)[,"Species"]),
#                            FAmtgraw)@otu_table)

write.table(FAMT3OTU_sel,"Data/FA_Mtg/FAMT3_OTU_sel_alllev.txt", sep = "\t")
write.table(data.frame(t(FAMT3Meta[,c(49:57)])),"Data/FA_Mtg/FAMT3_Meta_sel.txt", sep = "\t")

# (Pre-dealing with HALLA)
FA_Bac_Cor <- read.table("Results/FA_mtg/Corr_alllev_pea/all_associations.txt",
                         sep = "\t", header = TRUE)
FA_Bac_Cor$Y_features[grepl("\\|",FA_Bac_Cor$Y_features)] <- as.character(lapply(strsplit(FA_Bac_Cor$Y_features[grepl("\\__",FA_Bac_Cor$Y_features)],"\\|"),
                                                                                 function(x) paste0(x)[length(x)]))


FA_Bac_Cor <- FA_Bac_Cor %>% #filter((association)>0.5) %>% 
  filter(q.values<0.05)


FA_Bac_Cor1 <- rbind(FA_Bac_Cor,FA_Bac_CorMS)
writexl::write_xlsx(FA_Bac_CorMS,"Results/SweMS_Cor.xlsx")


FA_Bac_Cor <- rbind(FA_Bac_Cor,FA_Bac_Cor1)

pdf("test.pdf", height = 15, width = 15)
ggplot(FA_Bac_Cor,aes(x = X_features,
                      y = Y_features,
                      fill = association)) +
  scale_fill_gradient2(name = "R2",
                       midpoint = 0.32,
                       low = "#fee0d2", 
                       high = "#a50f15") 

  

invisible(dev.off())


# Pathway----
library(file2meco)
library(microeco)
library(magsrittr)

FA_PW <- read.table("Data/FA_Mtg/T3_humann_pathabundance_stratified_relab.tsv",
                    sep = "\t", header = T)
colnames(FA_PW) <- str_replace(colnames(FA_PW),"_Abundance","")
colnames(FA_PW) <- str_replace_all(colnames(FA_PW),"\\.","_")
colnames(FA_PW) <- str_replace_all(colnames(FA_PW),"X","")

FA_PW_Sum <- FA_PW
FA_PW_Sum$Pathway <- as.character(lapply(FA_PW_Sum$Pathway,
                                         function(x) strsplit(x,"\\|")[[1]][1]))
FA_PW_Sum <- FA_PW_Sum %>% group_by(Pathway) %>% 
  summarise(across(everything(), list(sum)))
FA_PW_Sum <- data.frame(FA_PW_Sum)
rownames(FA_PW_Sum) <- FA_PW_Sum$Pathway
FA_PW_Sum <- FA_PW_Sum[,-1]
FA_PW_Sum <- FA_PW_Sum %>% mutate(sum0 = as.numeric(apply(FA_PW_Sum,1,
                                                  function(x) sum(x==0, na.rm = T)))) %>% 
  filter(sum0<0.95*(ncol(FA_PW_Sum)-1)) %>% select(-sum0)
FA_PW_Sum <- data.frame(apply(FA_PW_Sum,2,function(x) round(x*100,2)))


# FA_PW <-  FA_PW %>% select(c("Pathway",
#                              colnames(FA_PW)[colnames(FA_PW) %in%
#                                                colnames(FAMT3OTU_sel)]))

rownames(FA_PW) <- FA_PW$Pathway
FA_PW <- FA_PW[,-1]

rownames(FA_PW)[grepl("\t",rownames(FA_PW))] <-  as.character(lapply(rownames(FA_PW)[grepl("\t",rownames(FA_PW))] ,
                                                                           function(x)
                                                                             strsplit(x,"\\\t0")[[1]][1]))

FA_PW <- FA_PW %>% mutate(sum0 = as.numeric(apply(FA_PW,1,
                                                  function(x) sum(x==0, na.rm = T)))) %>% 
  filter(sum0<0.95*(ncol(FA_PW)-1)) %>% select(-sum0)
FA_PW <- data.frame(apply(FA_PW,2,function(x) round(x*100,2)))

write.table(FA_PW,"Data/FA_Mtg/FAMT3_PW_alllev.txt", sep = "\t")
abund_file_path <- "Data/FA_Mtg/FAMT3_PW_alllev.txt"
FA_abutest <- humann2meco(abund_table = abund_file_path, db = "MetaCyc")

FA_PW <- FA_PW %>% filter(rownames(FA_PW) %in% rownames(environment(FA_abutest[["tidy_dataset"]])[["self"]][["tax_table"]]))
a <- data.frame(environment(FA_abutest[["tidy_dataset"]])[["self"]][["tax_table"]])
a <- a %>% filter(rownames(a) %in% rownames(FA_PW))
rownames(FA_PW) <- as.character(apply(a,1, function (x) 
                                        paste0(x[c(2,3,10)],collapse = "_")))


# VT ----

Y1phy <- readRDS("Data/Year 1/MOMmy_phyloseq_20220701.rds")
M1phy <-  subset_samples(Y1phy, grepl("Month 1|Month 2|Mother T3", Time_point)&
                           !grepl("Month 12", Time_point))
# M1phy <- tax_glom(M1phy, taxrank = "Family") # subset
M1OTU <- data.frame(M1phy@otu_table)
M1_Taxa <- data.frame(M1phy@tax_table)
rownames(M1OTU) <- M1_Taxa$Species
M1_Meta <- data.frame(M1phy@sam_data)


# Identify co-occurence species
MB_otu_t <- t(M1OTU)
rownames(M1_Meta) <- paste0("X",rownames(M1_Meta))
MB_otu_t1 <- merge(M1_Meta[,c(3,5)], MB_otu_t, by="row.names")
rownames(MB_otu_t1) <- MB_otu_t1$Row.names
Paired_OTU <- MB_otu_t1[,-1]

Early_pair <- list()
for (i in unique(Paired_OTU$FamilyID)){
  if(all(c("Month 2", "Mother T3") %in% 
         Paired_OTU$Time_point[Paired_OTU$FamilyID==i])|
     all(c("Month 1", "Mother T3") %in% 
         Paired_OTU$Time_point[Paired_OTU$FamilyID==i])
  ){
    Early_pair <- c(Early_pair,i)
  }
}

Early_paired_OTU <- Paired_OTU %>% filter(FamilyID %in% Early_pair) 
Early_paired_OTU <- Early_paired_OTU[order(Early_paired_OTU$FamilyID,
                                           Early_paired_OTU$Time_point),]
Early_paired_OTU_uni <- Early_paired_OTU %>% 
  select(-Time_point) %>% 
  group_by(FamilyID) %>% 
  summarise(across(everything(), list(function(x)
    if(length(x)>2){ # 根据每个familyID内x的数量，进行区分！ 
      ifelse(prod(x[c(1,3)])!=0|prod(x[c(2,3)])!=0,1,0)} # 内部写ifelse进行判断
    else 
    {ifelse(prod(x)!=0,1,0)})))


Early_share_spe <- data.frame(
  apply(Early_paired_OTU_uni[,-1],2,sum)
)
Early_share_spe$Spe <- rownames(Early_share_spe)
colnames(Early_share_spe)[1] <- "Freq"
Early_share_spe <- Early_share_spe %>% filter(Freq != 0)
Early_share_spe <- Early_share_spe[order(Early_share_spe$Freq,decreasing = TRUE),]
Early_share_spe$Order <- 1:length(Early_share_spe$Freq)
Early_paired_OTU_uni <- Early_paired_OTU_uni %>% 
  select(c("FamilyID",colnames(Early_paired_OTU_uni)[colnames(Early_paired_OTU_uni) %in% Early_share_spe$Spe[1:20]]))

Early_paired_OTU_uni <- data.frame(Early_paired_OTU_uni)
Early_paired_OTU_uni$Sub_ID <- Urine_Sweet$Sub_ID[
  match(Early_paired_OTU_uni$FamilyID,substr(Urine_SweetT3_pair$Sub_ID,6,9))
]
Early_paired_OTU_uni <- Early_paired_OTU_uni %>% filter(complete.cases(Sub_ID)) %>% 
  select(-FamilyID)
# strainphlan
FA_VT <- read.table("VT_pair_sum.txt",header = T, sep = "\t")
FA_VT$speceis <- str_remove(FA_VT$speceis,"s__")
FA_VT_sel <- FA_VT %>% filter(Timepoint=="BM1_2") %>% 
  filter(substr(sample,1,3) %in% substr(FAMT3Meta$Sub_ID,7,9)) %>% 
  mutate(Sub_ID = paste0("PWH100",substr(sample,1,3)))
FA_VT_sel <- reshape(FA_VT_sel[,c(6,8,9)],
                     timevar = "speceis", idvar = "Sub_ID",
                     direction = "wide")
# VT_list <- list()
# for(i in colnames(FA_VT_sel)[-1]){
#   if(sum(FA_VT_sel[,i]=="yes", na.rm = T)>10){
#     VT_list <- c(VT_list,i)
#   }
# }
# FA_VT_sel <- FA_VT_sel %>% select(c("Sub_ID",
#                                     colnames(FA_VT_sel)[colnames(FA_VT_sel) %in% VT_list]))
# 
# for(i in colnames(FA_VT_sel)[-1]){
#   print(i)
#   print(table(FA_VT_sel[,i]))
# 
# }

FA_VT_sel$Any_VT <- as.character(apply(FA_VT_sel[,-1],1,
                                       function(x) ifelse(all(is.na(x)),"No","Yes")))
FA_VT_sel <- FA_VT_sel %>% filter(Any_VT=="Yes")


# VIF and emulsifier----
FA_VIF_db <- readRDS("Data/FA_Mtg/MT3_VIF.rds")
VFDB <- data.frame(read_xlsx("~/Dropbox/g/VFDB_Sum/VFDB_Sum.xlsx"))
VFDB$VFDB_lab <- str_remove_all(VFDB$VFDB_lab,"\\>")
VFDB <- data.frame(VFDB[!duplicated(VFDB$VFDB_lab),])
colnames(VFDB) <- "VFDB_lab"
VFDB <- VFDB %>% mutate(VFGID = as.character(lapply(VFDB$VFDB_lab,
                                                    function(x) strsplit(x,"\\(")[[1]][1])),
                        Func = as.character(lapply(VFDB$VFDB_lab,
                                                   function(x) strsplit(x,"\\)")[[1]][3])),
                        Func_sum = as.character(lapply(VFDB$VFDB_lab,
                                                   function(x) strsplit(x,"\\[")[[1]][length(strsplit(x,"\\[")[[1]])-1])),
                        Bac = as.character(lapply(VFDB$VFDB_lab,
                                                  function(x) strsplit(x,"\\[")[[1]][length(strsplit(x,"\\[")[[1]])]))) 
VFDB$Func[!grepl("\\[",VFDB$Func)] <- as.character(lapply(VFDB$VFDB_lab[!grepl("\\[",VFDB$Func)],
                                                          function(x) strsplit(x,"\\)")[[1]][2]))
VFDB$Func[!grepl("\\[",VFDB$Func)] <- as.character(lapply(VFDB$VFDB_lab[!grepl("\\[",VFDB$Func)],
                                                          function(x) str_split_fixed(x,"\\)",3)[3]))
VFDB$Func<- as.character(lapply(VFDB$Func,
                                function(x) strsplit(x,"\\[")[[1]][1]))
# %>% 
#   select(-V1)

VFDB$Bac <- as.character(lapply(VFDB$Bac,function(x) paste0(strsplit(x," ")[[1]][1],"_",
                                                            strsplit(x," ")[[1]][2])))
colnames(VFDB)[2] <- "DB_ID"

VICDB <- data.frame(read_excel("~/Dropbox/g/VFDB_Sum/Victors_db_Sum.xlsx"))
VICDB$ID <- as.character(lapply(VICDB$Vic_lab,
                                function(x) strsplit(x,"\\|")[[1]][4]))

VICDB$ID[!startsWith(VICDB$Vic_lab,"gi")] <- as.character(lapply(VICDB$Vic_lab[!startsWith(VICDB$Vic_lab,"gi")],
                                                                    function(x) strsplit(x,"\\ ")[[1]][1]))
VICDB$Fun <- as.character(lapply(VICDB$Vic_lab,
                                function(x) strsplit(x,"\\|")[[1]][5]))
VICDB$Fun <- as.character(lapply(VICDB$Fun,
                                 function(x) strsplit(x,"\\[")[[1]][1]))
VICDB$Bac <- as.character(lapply(VICDB$Vic_lab,
                                 function(x) strsplit(x,"\\[")[[1]][2]))
VICDB$Bac <- str_remove(VICDB$Bac,"\\]")
VICDB$Fun[!startsWith(VICDB$Vic_lab,"gi")] <- as.character(lapply(VICDB$Vic_lab[!startsWith(VICDB$Vic_lab,"gi")],
                                                                 function(x) strsplit(x,"\\.1")[[1]][2]))
VICDB$Fun[!startsWith(VICDB$Vic_lab,"gi")] <- as.character(lapply(VICDB$Fun[!startsWith(VICDB$Vic_lab,"gi")],
                                                                  function(x) strsplit(x,"\\[")[[1]][1]))
VICDB$ID <- str_replace_all(VICDB$ID,"\\.","_")
colnames(VICDB)[c(2,3)] <- c("DB_ID","Func")
saveRDS(VICDB,"~/Dropbox/g/VFDB_Sum/Victors_db_Sum.rds")

FA_VIF_db <- FA_VIF_db %>% mutate(sum0=as.numeric(apply(FA_VIF_db,1,
                                                        function(x) sum(x!=0)))) %>% 
  filter(sum0<ncol(FA_VIF_db)*0.8) %>% select(-sum0)
FA_VIF <- data.frame(rownames(FA_VIF_db))
colnames(FA_VIF) <- "FA_VIF"

FA_VIF$DB_ID <- NA
FA_VIF$DB_ID[grepl("VFG",FA_VIF$FA_VIF)] <- 
  as.character(lapply(FA_VIF$FA_VIF[grepl("VFG",FA_VIF$FA_VIF)],
                      function(x) strsplit(x,"\\|")[[1]][2]))

FA_VIF$DB_ID<- as.character(lapply(FA_VIF$DB_ID,
                      function(x) strsplit(x,"\\(")[[1]][1]))
FA_VIF$DB_ID[startsWith(FA_VIF$FA_VIF,"Victors")] <- 
  as.character(lapply(FA_VIF$FA_VIF[startsWith(FA_VIF$FA_VIF,"Victors")],
                      function(x) strsplit(x,"\\|")[[1]][5]))
FA_VIF <- left_join(FA_VIF,VFDB[,c(2,3,5)],"DB_ID")
FA_VIF$Func[is.na(FA_VIF$Func)] <- VICDB$Func[match(FA_VIF$DB_ID[is.na(FA_VIF$Func)],
                                                    VICDB$DB_ID)]

FA_VIF$Func[startsWith(FA_VIF$FA_VIF,"virulence")] <- 
  as.character(lapply(FA_VIF$FA_VIF[startsWith(FA_VIF$FA_VIF,"virulence")],
                      function(x) strsplit(x,"\\|")[[1]][length(strsplit(x,"\\|")[[1]])]))





