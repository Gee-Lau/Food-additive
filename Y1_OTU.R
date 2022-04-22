library(MicrobiotaProcess)
library(patchwork)
library(ggplot2)
library(dplyr)
library(readxl)
library(phyloseq)
library(vegan)

Y1metadata <- data.frame(read_xlsx("Data/Year 1/Y1_metadata_20042022.xlsx"))
# Y1metadata$Ref_new <- ifelse(Y1metadata$ATB_TF==0&Y1metadata$Weight_child>2500&Y1metadata$Disease_M12==0,
#                              1,0)
Y1metadata[,c(16:28)] <- apply(Y1metadata[,c(16:28)], 2, as.character)
Y1OTU <- data.frame(read.table("Data/Year 1/MOMmy_M12_20220411.txt",
                                  header = T, sep = "\t"))
Y1Taxa <- data.frame(matrix(unlist(lapply(strsplit(Y1OTU$clade_name,"\\|"), function(x)
  if(length(x)<7){
    c(paste0(x),rep(NA,7-length(x)))
  }
  else{paste0(x)})), 
                    ncol = 7, byrow = TRUE))
Y1OTU <- Y1OTU[,-(1:2)]
colnames(Y1OTU) <- Y1metadata$Samp_name
rownames(Y1metadata) <- Y1metadata$Samp_name
colnames(Y1Taxa) <- as.character(lapply(Y1Taxa[9,], function(x)
                                        strsplit(x,"__")[[1]][1]))

Y1mtgraw <- phyloseq(otu_table(Y1OTU, taxa_are_rows = TRUE), 
                     sample_data(Y1metadata), 
                     tax_table(as.matrix(Y1Taxa)))
# Select bacteria
Y1bacmtg <- prune_taxa(grepl("Bac",tax_table(Y1mtgraw)[, "k"]), Y1mtgraw)

# Select species level and filter out the species <60% in the population, 
# and recalculate relative abundance
Y1bacmtgsp <- prune_taxa(complete.cases(tax_table(Y1bacmtg)[, "s"]), Y1bacmtg)

# spe.sum = tapply(taxa_sums(Y1bacmtgsp), 
#                  tax_table(Y1bacmtgsp)[, "s"],
#                  sum, na.rm=TRUE) # Sort by abundance 
# freqspe = names(sort(spe.sum, TRUE))[1:round(length(spe.sum)*0.8)] # Top 80%

spe.sum = tapply(apply(Y1bacmtgsp@otu_table,1,function(x) sum(x==0)), 
                                  tax_table(Y1bacmtgsp)[, "s"],
                                  sum, na.rm=TRUE) 

freqspe = names(sort(spe.sum[spe.sum<(0.4*length(spe.sum))], TRUE))


Y1bacmtgsp = prune_taxa((tax_table(Y1bacmtgsp)[, "s"] %in% freqspe),
                        Y1bacmtgsp)

Y1bacmtgsp@otu_table <- otu_table(apply(Y1bacmtgsp@otu_table, 2, function(x)
  round(x/sum(x)*100,2)), taxa_are_rows = TRUE)
rownames(Y1bacmtgsp@tax_table@.Data) <- tax_table(Y1bacmtgsp)[, "s"]

set.seed(100)
Y1nmds <- vegdist(t(Y1bacmtgsp@otu_table),na.rm = T)
Y1nmds <- monoMDS(Y1nmds)

Y1nmds <- data.frame(cbind(Y1metadata[,c(1,2,29,30,31,23)],
                           Y1nmds[["points"]]))
Y1nmds[,c(3:5)] <- apply(Y1nmds[,c(3:5)], 2, as.character)

Y1nmds$MDS1 <- as.numeric(Y1nmds$MDS1)
Y1nmds$MDS2 <- as.numeric(Y1nmds$MDS2)
# gender.md <- gender.md %>% filter(abs(MDS1)<0.1) %>%
#   filter(abs(MDS2)<0.1)

ggplot(Y1nmds, aes(MDS1, MDS2, col = Disease_Spe_M12)) +
  geom_point() +
  # scale_color_manual(name = "Disease",
  #                    label = c("No","Yes"),
  #                    values=c("#0072B5FF","#BC3C29FF")) +
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





