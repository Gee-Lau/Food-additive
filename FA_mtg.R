library(MicrobiotaProcess)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readxl)
library(phyloseq)
library(vegan)
library(grid)

# Creating phyloseq db----
FAMT3OTU <- readRDS("Data/Year 1/MT3_M12.rds")
FAMT3OTU <-  subset_samples(FAMT3OTU, grepl("Mother T", Time_point))
FAMT3Taxa <- data.frame(FAMT3OTU@tax_table)
FAMT3Meta <- data.frame(FAMT3OTU@sam_data)
FAMT3Meta <- FAMT3Meta[,c(4,2,8:35,38,39,41,42,44:46,48:50,113:119,121)]
FAMT3Meta <- left_join(FAMT3Meta,MomFA[,c(1:10,37:39,12,79:91)],"Sub_ID")
FAMT3Meta <- FAMT3Meta %>% filter(complete.cases(All_FA))

FAMT3OTU <- data.frame(FAMT3OTU@otu_table)
colnames(FAMT3OTU) <- str_replace(colnames(FAMT3OTU),"X","")
colnames(FAMT3OTU) <- str_replace(colnames(FAMT3OTU),"\\.","-")
FAMT3OTU <- FAMT3OTU %>% select(colnames(FAMT3OTU)[colnames(FAMT3OTU) %in% FAMT3Meta$sample_name])
rownames(FAMT3Meta) <- FAMT3Meta$sample_name

FAmtgraw <- phyloseq(otu_table(as.matrix(FAMT3OTU), taxa_are_rows = TRUE), 
                     sample_data(FAMT3Meta), 
                     tax_table(as.matrix(FAMT3Taxa)))

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


# α div----

FAalpha <- cbind(FAMT3Meta,
                 Richness_spe = t(estimateR(round(t(FAmtgraw@otu_table))))[,1],
                 Richness_Chao1 = t(estimateR(round(t(FAmtgraw@otu_table))))[,2],
                 Evenness = diversity(t(FAmtgraw@otu_table)) / log(specnumber(t(FAmtgraw@otu_table))),
                 Shannon = diversity(t(FAmtgraw@otu_table), index = "shannon"),
                 InvSimp = diversity(t(FAmtgraw@otu_table), index = "invsimpson")) 
FAalpha[,c(62:74)] <- apply(FAalpha[,c(62:74)], 2, function(x)
  ifelse(x=="High","High","Normal"))
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

pdf("Results/FA_mtg/FA_alpha_sum.pdf", height = 45, width = 16)
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
FAnmds <- vegdist(t(FAmtgraw@otu_table),na.rm = T)
FAnmds <- monoMDS(FAnmds)

FAnmds <- data.frame(cbind(FAMT3Meta,
                           FAnmds[["points"]]))

FAnmds$MDS1 <- as.numeric(FAnmds$MDS1)
FAnmds$MDS2 <- as.numeric(FAnmds$MDS2)
# gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
#   filter(abs(MDS2)<0.1)

ggplot(FAnmds, aes(MDS1, MDS2, col = Emu_cut3)) +
  geom_point() +
  scale_color_manual(name = "Emulsifier group",
                     values=c("#20B54EFF","#FFDC91FF","#BC3C29FF")) +
  stat_ellipse(linetype = 2, size = 1.5, type = "t") +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("") 

ggplot(Y1nmds, aes(MDS1, MDS2, col = Ref_Hein)) +
  geom_point() +
  scale_color_manual(name = "Ref",
                     label = c("No","Yes"),
                     values=c("#0072B5FF","#20B54EFF")) +
  stat_ellipse(linetype = 2, size = 1.5, aes(color = Disease_M12)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("") 

Y1permanova <- get_dist(Y1bacmtgsp, distmethod ="bray", method="hellinger")

Y1PCA <- get_pcoa(obj = Y1bacmtgsp,
                 method = "hellinger",
                 distmethod = "euclidean",
                 taxa_are_rows = TRUE)


ggordpoint(obj = Y1PCA, biplot=F, speciesannot=F,
                       factorNames=c("Ref_Hein"), ellipse=F) +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Ref",
                     label = c("No","Yes"),
                     values=c("#0072B5FF","#20B54EFF")) +
  scale_color_manual(name = "Ref",
                     label = c("No","Yes"),
                     values=c("#0072B5FF","#20B54EFF")) +
  stat_ellipse() 





