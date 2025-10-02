library("edgeR")
library("ggplot2")
library("tibble")
library("dplyr")
library("edgeR")
library("ggrepel")
library("readr")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("pheatmap")
library("Cormotif")

##ASHG_RNA#####
###Feature Counts#####
chimp_RNA_pantro5_fc <- read.delim("~/diff_timeline_tes/RNA/featurecounts/c_samples_counts.txt", comment.char="#",row.names=1)
human_hg38_fc <- read.delim("~/diff_timeline_tes/RNA/featurecounts/h_sample_counts.txt", comment.char="#",row.names=1)
human_df <- human_RNA_hg38_fc %>%
  tibble::rownames_to_column("gene")
chimp_df <- chimp_RNA_pantro5_fc %>%
  tibble::rownames_to_column("gene")
RNA_joined_fc <- left_join(human_df, chimp_df, by = "gene")

#To keep only OrthoGenes and sample columns
RNA_fc <- RNA_joined_fc[ , !(names(RNA_joined_fc) %in% c("gene","Chr.x","Start.x","End.x","Strand.x","Length.x","Geneid.x","Chr.y","Start.y","End.y","Strand.y","Length.y","Geneid.y"))]
# 89 columns; 7 x 6 = 42 Human Exp, 7 x 6 = 42 Chimp Exp, 4 Human Replicate
rownames(RNA_fc) <- RNA_joined_fc$gene

#Rename column Names to More Useful Info
col_names<- c("H28126_D0",
              "H28126_D2",
              "H28126_D4",
              "H28126_D5",
              "H28126_D15",
              "H28126_D30",
              "H17_D0",
              "H17_D2",
              "H17_D4",
              "H17_D5",
              "H17_D15",
              "H17_D30",
              "H78_D0",
              "H78_D2",
              "H78_D4",
              "H78_D5",
              "H78_D15",
              "H78_D30",
              "H20682_D0",
              "H20682_D2",
              "H20682_D4",
              "H20682_D5",
              "H20682_D15",
              "H20682_D30",
              "H22422_D0",
              "H22422_D2",
              "H22422_D4",
              "H22422_D5",
              "H22422_D15",
              "H22422_D30",
              "H21792_D0",
              "H21792_D2",
              "H21792_D4",
              "H21792_D5",
              "H21792_D15",
              "H21792_D30",
              "H24280_D0",
              "H24280_D2",
              "H24280_D4",
              "H24280_D5",
              "H24280_D15",
              "H24280_D30",
              "H20682R_D0",
              "H20682R_D2",
              "H20682R_D5",
              "H20682R_D30",
              "C3649_D0",
              "C3649_D2",
              "C3649_D4",
              "C3649_D5",
              "C3649_D15",
              "C3649_D30",
              "C4955_D0",
              "C4955_D2",
              "C4955_D4",
              "C4955_D5",
              "C4955_D15",
              "C4955_D30",
              "C3651_D0",
              "C3651_D2",
              "C3651_D4",
              "C3651_D5",
              "C3651_D15",
              "C3651_D30",
              "C40210_D0",
              "C40210_D2",
              "C40210_D4",
              "C40210_D5",
              "C40210_D15",
              "C40210_D30",
              "C8861_D0",
              "C8861_D2",
              "C8861_D4",
              "C8861_D5",
              "C8861_D15",
              "C8861_D30",
              "C40280_D0",
              "C40280_D2",
              "C40280_D4",
              "C40280_D5",
              "C40280_D15",
              "C40280_D30",
              "C3647_D0",
              "C3647_D2",
              "C3647_D4",
              "C3647_D5",
              "C3647_D15",
              "C3647_D30"
)
colnames(RNA_fc)<-col_names
dim(RNA_fc)

ASHG_RNA_fc <- RNA_fc[ , !(names(RNA_fc) %in% c("H28126_D4",
              "H28126_D5",
              "H17_D5",
              "H78_D5",
              "H20682_D5",
              "H22422_D5",
              "H21792_D5",
              "H24280_D5",
              "H20682R_D5",
              "C3649_D5",
              "C4955_D5",
              "C3651_D5",
              "C40210_D5",
              "C8861_D5",
              "C40280_D5",
              "C3647_D5"))]


#####Unfiltered####
ASHG_RNA_log2cpm <- cpm(ASHG_RNA_fc,log=TRUE)
hist(ASHG_RNA_log2cpm,  main = "Histogram of all counts (unfiltered)",
     xlab =expression("Log"[2]*" counts-per-million"), col =4 )
boxplot(ASHG_RNA_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

timepoint <- c(rep(c("D0","D2","D4","D15","D30"),7),"D0","D2","D30", rep(c("D0","D2","D4","D15","D30"),7))
time_order <- c("D0", "D2","D4","D15","D30")
species <- c(rep("H",35),rep("HR",3),rep("C",35))
smlabel <- (interaction(time,species))

T_ASHG_RNA_log2cpm <- t(ASHG_RNA_log2cpm)
pca_raw <- prcomp(T_ASHG_RNA_log2cpm,scale=TRUE,center=TRUE)
pca_raw_meta <- as.data.frame(pca_raw$x)
pca_raw_meta$sample <- rownames(pca_raw_meta)
pca_raw_meta$species <- species
pca_raw_meta$timepoint <- timepoint
pca_raw_meta$timepoint <- factor(pca_raw_meta$timepoint, levels = time_order)
dim(T_ASHG_RNA_log2cpm)

ggplot(pca_raw_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM ASHG_RNA-Seq Data Unfiltered",
       x = paste0("PC1 (", round(summary(pca_raw)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_raw)$importance[2, 2] * 100, 1), "%)"))

#####RowMu>0####
row_means <- rowMeans(ASHG_RNA_log2cpm)
Filt_RMG0_ASHG_RNA_fc <- ASHG_RNA_fc[row_means >0,]
Filt_RMG0_ASHG_RNA_log2cpm <- cpm(Filt_RMG0_ASHG_RNA_fc,log=TRUE)

hist(Filt_RMG0_ASHG_RNA_log2cpm, main = "Histogram of filtered counts using rowMeans > 0 method",
     xlab =expression("Log"[2]*" counts-per-million"), col =5 )
boxplot(Filt_RMG0_ASHG_RNA_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

######PCA####
T_Filt_RMG0_ASHG_RNA_log2cpm <- t(Filt_RMG0_ASHG_RNA_log2cpm)
pca_RMG0 <- prcomp(T_Filt_RMG0_ASHG_RNA_log2cpm,scale=TRUE,center=TRUE)
pca_RMG0_meta <- as.data.frame(pca_RMG0$x)
pca_RMG0_meta$sample <- rownames(pca_RMG0_meta)
pca_RMG0_meta$species <- species
pca_RMG0_meta$timepoint <- timepoint
pca_RMG0_meta$timepoint <- factor(pca_RMG0_meta$timepoint, levels = time_order)
dim(T_Filt_RMG0_ASHG_RNA_log2cpm)

#PC1 vs PC2
ggplot(pca_RMG0_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM ASHG_RNA-Seq Data RowMean>0",
       x = paste0("PC1 (", round(summary(pca_RMG0)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_RMG0)$importance[2, 2] * 100, 1), "%)"))

######Cor_HeatMap####
cor_Filt_RMG0_ASHG_RNA_log2cpm <- cor(Filt_RMG0_ASHG_RNA_log2cpm, method = "pearson")
time <- factor(timepoint,levels = c("D0","D2","D4","D15","D30"))
metadata <- data.frame(
  sample_cor = colnames(Filt_RMG0_ASHG_RNA_log2cpm),
  species_cor = species,
  timepoint_cor = timepoint
)
ann_colors <- list(
  timepoint_cor = c(
    "D0" = "#883268",   # Purple
    "D2" = "#3E7274",  # blue
    "D4" = "#3E7274",  # blue
    "D15" = "#3E7274",  # blue
    "D30" = "#C03830",  # red
    "D0R" = "#883268",  # pink
    "D2R" = "#3E7274",  # teal
    "D30R" = "#C03830"  # red
  ),
  species_cor = c(
    "H" = "#171717",  # black
    "C" = "#17171717",   # light grey
    "HR" = "#171717B0"#Dark Grey
  )
)
rownames(metadata) <- metadata$sample_cor
pheatmap(cor_Filt_RMG0_ASHG_RNA_log2cpm,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = metadata[, c("species_cor", "timepoint_cor")],
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Sample-Sample Correlation (Pearson-log2CPM)")

###### CorMotif Setup (Pooled by Species + Timepoint) ######


#### Step 1: Prepare Data Set
###Remove specific replicates
samples_to_remove_CorM <- c("H20682R_D0", "H20682R_D2", "H20682R_D30")
Filt_RMG0_ASHG_RNA_fc_NoRep <- Filt_RMG0_ASHG_RNA_fc[, !(colnames(Filt_RMG0_ASHG_RNA_fc) %in% samples_to_remove_CorM)]
dge <- DGEList(counts=Filt_RMG0_ASHG_RNA_fc_NoRep)
dge <- calcNormFactors(dge,method = "TMM")
log2cpm <- cpm(dge,log = TRUE, prior.count=1)

# Set up Variables
groupid <- c(rep(1:3,7),rep(4:6,7))
compid_Cond1 <- c(rep(1,2),rep(4,2))
compid_Cond2 <- c(2,3,5,6)
compid <- cbind(compid_Cond1,compid_Cond2)

#Model Fitting
set.seed(123)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(1:10),max.iter=2000,BIC=TRUE)
motif.fitted$bic
plotIC(motif.fitted)
plotMotif(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.ASHG_RNA<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.ASHG_RNA)
tail(dif.pattern.ASHG_RNA)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)

#Model Fitting
set.seed(91011)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(11),max.iter=2000,BIC=TRUE)
#motif.fitted$bic
#plotIC(motif.fitted)
plotMotif(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.ASHG_RNA<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.ASHG_RNA)
tail(dif.pattern.ASHG_RNA)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)

##ASHG_H3K27ac#####
####samples names####
h_col_names<- c("IgG_H21792_D0",
                "H21792_D0",
                "H21792_D2",
                "H21792_D4",
                "H21792_D5",
                "H21792_D15",
                "H21792_D30",
                "H24280_D0",
                "H24280_D2",
                "H24280_D4",
                "H24280_D5",
                "H24280_D15",
                "H24280_D30",
                "H28126_D0",
                "H28126_D2",
                "H28126_D4",
                "H28126_D5",
                "H28126_D15",
                "H28126_D30",
                "IgG_H21792_D2",
                "H17_3_D0",
                "H17_3_D2",
                "H17_3_D4",
                "H17_3_D5",
                "H17_3_D15",
                "H17_3_D30",
                "H78_1_D0",
                "H78_1_D2",
                "H78_1_D4",
                "H78_1_D5",
                "H78_1_D15",
                "H78_1_D30",
                "H20682_D0",
                "H20682_D2",
                "H20682_D4",
                "H20682_D5",
                "H20682_D15",
                "H20682_D30",
                "H20682R_D0",
                "H20682R_D2",
                "H20682R_D5",
                "H20682R_D30",
                "H22422_D0",
                "H22422_D2",
                "H22422_D4",
                "H22422_D5",
                "H22422_D15",
                "H22422_D30",
                "IgG_H21792_D4",
                "IgG_H21792_D5",
                "IgG_H21792_D15")

c_col_names <- c("IgG_C40280_D0",
                 "C40280_D0",
                 "C40280_D2",
                 "C40280_D4",
                 "C40280_D5",
                 "C40280_D15",
                 "C40280_D30",
                 "C3647_D0",
                 "C3647_D2",
                 "C3647_D4",
                 "C3647_D5",
                 "C3647_D15",
                 "C3647_D30",
                 "C3649_D0",
                 "C3649_D2",
                 "C3649_D4",
                 "C3649_D5",
                 "C3649_D15",
                 "C3649_D30",
                 "IgG_C40280_D2",
                 "C4955_D0",
                 "C4955_D2",
                 "C4955_D4",
                 "C4955_D5",
                 "C4955_D15",
                 "C4955_D30",
                 "C3651_D0",
                 "C3651_D2",
                 "C3651_D4",
                 "C3651_D5",
                 "C3651_D15",
                 "C3651_D30",
                 "C40210_D0",
                 "C40210_D2",
                 "C40210_D4",
                 "C40210_D5",
                 "C40210_D15",
                 "C40210_D30",
                 "C8861_D0",
                 "C8861_D2",
                 "C8861_D4",
                 "C8861_D5",
                 "C8861_D15",
                 "C8861_D30",
                 "IgG_C40280_D4",
                 "IgG_C40280_D5",
                 "IgG_C40280_D15",
                 "IgG_C40280_D30"
)

###Data Upload####
ASHG_c_samples_counts <- read.delim("~/Evol_CM_Diff_TEs/analysis/merged_results/2stp/ASHG_c_samples_counts_2pSTP.txt", comment.char="#")
ASHG_h_samples_counts <- read.delim("~/Evol_CM_Diff_TEs/analysis/merged_results/2stp/ASHG_h_samples_counts_2pSTP.txt", comment.char="#")
ASHG_H3K27ac_joined_fc <- left_join(ASHG_h_samples_counts, ASHG_c_samples_counts, by = "Geneid")
rownames(ASHG_H3K27ac_joined_fc)<-ASHG_H3K27ac_joined_fc$Geneid
head(ASHG_H3K27ac_joined_fc)
dim(ASHG_H3K27ac_joined_fc)

###Union of Peaks####
ASHG_H3K27ac_fc <- ASHG_H3K27ac_joined_fc[ , !(names(ASHG_H3K27ac_joined_fc) %in% c("Geneid","Chr.x","Start.x","End.x","Strand.x","Length.x","Chr.y","Start.y","End.y","Strand.y","Length.y"))]
dim(ASHG_H3K27ac_fc)
col_names<- c("IgG_H21792_D0",
              "H21792_D0",
              "H21792_D2",
              "H21792_D4",
              "H21792_D5",
              "H21792_D15",
              "H21792_D30",
              "H24280_D0",
              "H24280_D2",
              "H24280_D4",
              "H24280_D5",
              "H24280_D15",
              "H24280_D30",
              "H28126_D0",
              "H28126_D2",
              "H28126_D4",
              "H28126_D5",
              "H28126_D15",
              "H28126_D30",
              "IgG_H21792_D2",
              "H17_3_D0",
              "H17_3_D2",
              "H17_3_D4",
              "H17_3_D5",
              "H17_3_D15",
              "H17_3_D30",
              "H78_1_D0",
              "H78_1_D2",
              "H78_1_D4",
              "H78_1_D5",
              "H78_1_D15",
              "H78_1_D30",
              "H20682_D0",
              "H20682_D2",
              "H20682_D4",
              "H20682_D5",
              "H20682_D15",
              "H20682_D30",
              "H20682R_D0",
              "H20682R_D2",
              "H20682R_D5",
              "H20682R_D30",
              "H22422_D0",
              "H22422_D2",
              "H22422_D4",
              "H22422_D5",
              "H22422_D15",
              "H22422_D30",
              "IgG_H21792_D4",
              "IgG_H21792_D5",
              "IgG_H21792_D15",
              "IgG_C40280_D0",
              "C40280_D0",
              "C40280_D2",
              "C40280_D4",
              "C40280_D5",
              "C40280_D15",
              "C40280_D30",
              "C3647_D0",
              "C3647_D2",
              "C3647_D4",
              "C3647_D5",
              "C3647_D15",
              "C3647_D30",
              "C3649_D0",
              "C3649_D2",
              "C3649_D4",
              "C3649_D5",
              "C3649_D15",
              "C3649_D30",
              "IgG_C40280_D2",
              "C4955_D0",
              "C4955_D2",
              "C4955_D4",
              "C4955_D5",
              "C4955_D15",
              "C4955_D30",
              "C3651_D0",
              "C3651_D2",
              "C3651_D4",
              "C3651_D5",
              "C3651_D15",
              "C3651_D30",
              "C40210_D0",
              "C40210_D2",
              "C40210_D4",
              "C40210_D5",
              "C40210_D15",
              "C40210_D30",
              "C8861_D0",
              "C8861_D2",
              "C8861_D4",
              "C8861_D5",
              "C8861_D15",
              "C8861_D30",
              "IgG_C40280_D4",
              "IgG_C40280_D5",
              "IgG_C40280_D15",
              "IgG_C40280_D30"
)
colnames(ASHG_H3K27ac_fc)<-col_names
dim(ASHG_H3K27ac_fc)

ASHG_Subset_H3K27ac_fc <- ASHG_H3K27ac_fc[ , (names(ASHG_H3K27ac_fc) %in% c("IgG_H21792_D0",
                                                                             "H21792_D0",
                                                                             "H21792_D2",
                                                                             "H21792_D4",
                                                                             "H21792_D15",
                                                                             "H21792_D30",
                                                                             "H24280_D0",
                                                                             "H24280_D2",
                                                                             "H24280_D4",
                                                                             "H24280_D15",
                                                                             "H24280_D30",
                                                                             "H28126_D0",
                                                                             "H28126_D2",
                                                                             "H28126_D4",
                                                                             "H28126_D15",
                                                                             "H28126_D30",
                                                                             "IgG_H21792_D2",
                                                                             "H17_3_D0",
                                                                             "H17_3_D2",
                                                                             "H17_3_D4",
                                                                             "H17_3_D15",
                                                                             "H17_3_D30",
                                                                             "H78_1_D0",
                                                                             "H78_1_D2",
                                                                             "H78_1_D4",
                                                                             "H78_1_D15",
                                                                             "H78_1_D30",
                                                                             "H20682_D0",
                                                                             "H20682_D2",
                                                                             "H20682_D4",
                                                                             "H20682_D15",
                                                                             "H20682_D30",
                                                                             "H20682R_D0",
                                                                             "H20682R_D2",
                                                                             "H20682R_D30",
                                                                             "H22422_D0",
                                                                             "H22422_D2",
                                                                             "H22422_D4",
                                                                             "H22422_D15",
                                                                             "H22422_D30",
                                                                             "IgG_H21792_D4",
                                                                             "IgG_H21792_D15",
                                                                             "IgG_C40280_D0",
                                                                             "C40280_D0",
                                                                             "C40280_D2",
                                                                             "C40280_D4",
                                                                             "C40280_D15",
                                                                             "C40280_D30",
                                                                             "C3647_D0",
                                                                             "C3647_D2",
                                                                             "C3647_D4",
                                                                             "C3647_D15",
                                                                             "C3647_D30",
                                                                             "C3649_D0",
                                                                             "C3649_D2",
                                                                             "C3649_D4",
                                                                             "C3649_D15",
                                                                             "C3649_D30",
                                                                             "IgG_C40280_D2",
                                                                             "C4955_D0",
                                                                             "C4955_D2",
                                                                             "C4955_D4",
                                                                             "C4955_D15",
                                                                             "C4955_D30",
                                                                             "C3651_D0",
                                                                             "C3651_D2",
                                                                             "C3651_D4",
                                                                             "C3651_D15",
                                                                             "C3651_D30",
                                                                             "C40210_D0",
                                                                             "C40210_D2",
                                                                             "C40210_D4",
                                                                             "C40210_D15",
                                                                             "C40210_D30",
                                                                             "C8861_D0",
                                                                             "C8861_D2",
                                                                             "C8861_D4",
                                                                             "C8861_D15",
                                                                             "C8861_D30",
                                                                             "IgG_C40280_D4",
                                                                             "IgG_C40280_D15",
                                                                             "IgG_C40280_D30"
))]

dim(ASHG_Subset_H3K27ac_fc)
#####Unfiltered####
ASHG_Subset_H3K27ac_log2cpm <- cpm(ASHG_Subset_H3K27ac_fc,log=TRUE)
hist(ASHG_Subset_H3K27ac_log2cpm,  main = "Histogram of all counts (unfiltered)",
     xlab =expression("Log"[2]*" counts-per-million"), col =4 )
boxplot(ASHG_Subset_H3K27ac_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

timepoint <- c(
  "IgG_D0",
  rep(c("D0","D2","D4","D15","D30"),3),
  "IgG_D2",
  rep(c("D0","D2","D4","D15","D30"),3),
  "D0",
  "D2",
  "D30",
  rep(c("D0","D2","D4","D15","D30"),1),
  "IgG_D4",
  "IgG_D15",
  "IgG_D0",
  rep(c("D0","D2","D4","D15","D30"),3),
  "IgG_D2",
  rep(c("D0","D2","D4","D15","D30"),4),
  "IgG_D4",
  "IgG_D15",
  "IgG_D30")
dim(timepoint)
time_order <- c("D0", "D2","D4","D15","D30")
species <- c(rep("H",32),
             rep("HR",3),
             rep("H",7),
             rep("C",40))

T_ASHG_Subset_H3K27ac_log2cpm <- t(ASHG_Subset_H3K27ac_log2cpm)
pca_raw <- prcomp(T_ASHG_Subset_H3K27ac_log2cpm,scale=TRUE,center=TRUE)
pca_raw_meta <- as.data.frame(pca_raw$x)
pca_raw_meta$sample <- rownames(pca_raw_meta)
pca_raw_meta$species <- species
pca_raw_meta$timepoint <- timepoint
pca_raw_meta$timepoint <- factor(pca_raw_meta$timepoint, levels = time_order)
dim(T_ASHG_Subset_H3K27ac_log2cpm)

ggplot(pca_raw_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM Unfiltered ASHG_Subset_H3K27ac Data",
       x = paste0("PC1 (", round(summary(pca_raw)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_raw)$importance[2, 2] * 100, 1), "%)"))

#####Cor_HeatMap####
cor_ASHG_Subset_H3K27ac_log2cpm <- cor(ASHG_Subset_H3K27ac_log2cpm, method = "pearson")
dim(cor_ASHG_Subset_H3K27ac_log2cpm)
time <- factor(timepoint,levels = c("D0","D2","D4","D15","D30","D0R","D2R","D30R"))
metadata <- data.frame(
  sample_cor = colnames(ASHG_Subset_H3K27ac_log2cpm),
  species_cor = species,
  timepoint_cor = timepoint
)
ann_colors <- list(
  timepoint_cor = c(
    "D0" = "#883268",   # Purple
    "D2" = "#3E7274",  # blue
    "D4" = "#5AAA464D",  # light green
    "D15" = "#C03830",  # red
    "D30" = "#830C05" , # dark red,
    "IgG_D0" = "#EA55B3",  # Pink
    "IgG_D2"= "#0FE8F0",  # teal
    "IgG_D4"= "#57F909",  # green
    "IgG_D15"= "#F51002",  # bright red
    "IgG_D30"= "#171717"  # black
  ),
  species_cor = c(
    "H" = "#171717",  # black
    "C" = "#17171717",   # light grey
    "HR" = "#171717B0"#Dark Grey
  )
)
rownames(metadata) <- metadata$sample_cor
pheatmap(cor_ASHG_Subset_H3K27ac_log2cpm,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = metadata[, c("species_cor", "timepoint_cor")],
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "ASHG_Subset_H3K27ac Unfiltered Sample-Sample Correlation (Pearson-log2CPM)")

###### CorMotif Setup (Pooled by Species + Timepoint) ######


#### Step 1: Prepare Data Set
###Remove specific replicates
samples_to_remove_CorM <- c("H20682R_D0", "H20682R_D2","H20682R_D30")
ASHG_Subset_H3K27ac_fc_NoRep <- ASHG_Subset_H3K27ac_fc[, !(colnames(ASHG_Subset_H3K27ac_fc) %in% samples_to_remove_CorM)]
dge <- DGEList(counts=ASHG_Subset_H3K27ac_fc_NoRep)
dge <- calcNormFactors(dge,method = "TMM")
log2cpm <- cpm(dge,log = TRUE, prior.count=1)

# Set up Variables
groupid <- c(rep(1:3,7),rep(4:6,7))
compid_Cond1 <- c(rep(1,2),rep(4,2))
compid_Cond2 <- c(2,3,5,6)
compid <- cbind(compid_Cond1,compid_Cond2)

#Model Fitting
set.seed(123)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(1:10),max.iter=3000,BIC=TRUE)
motif.fitted$bic
plotIC(motif.fitted)
plotMotif(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.RNA<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.RNA)
tail(dif.pattern.RNA)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)

#####RowMu>0####
row_means <- rowMeans(ASHG_Subset_H3K27ac_log2cpm)
dim(ASHG_Subset_H3K27ac_fc)
Filt_RMG0_ASHG_Subset_H3K27ac_fc <- ASHG_Subset_H3K27ac_fc[row_means >1.8,]
dim(Filt_RMG0_ASHG_Subset_H3K27ac_fc)
Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm <- cpm(Filt_RMG0_ASHG_Subset_H3K27ac_fc,log=TRUE)
dim (Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm)

hist(Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm, main = "Histogram of filtered counts using rowMeans > 1.8 method",
     xlab =expression("Log"[2]*" counts-per-million"), col =5 )
boxplot(Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

######PCA####
T_Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm <- t(Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm)
pca_RMG0 <- prcomp(T_Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm,scale=TRUE,center=TRUE)
pca_RMG0_meta <- as.data.frame(pca_RMG0$x)
pca_RMG0_meta$sample <- rownames(pca_RMG0_meta)
pca_RMG0_meta$species <- species
pca_RMG0_meta$timepoint <- timepoint
pca_RMG0_meta$timepoint <- factor(pca_RMG0_meta$timepoint, levels = time_order)
dim(T_Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm)

#PC1 vs PC2
ggplot(pca_RMG0_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM ASHG_Subset_H3K27ac-Seq Data RowMean>2.5",
       x = paste0("PC1 (", round(summary(pca_RMG0)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_RMG0)$importance[2, 2] * 100, 1), "%)"))

######Cor_HeatMap####
cor_Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm <- cor(Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm, method = "pearson")
time <- factor(timepoint,levels = c("D0","D2","D4","D15","D30"))
metadata <- data.frame(
  sample_cor = colnames(Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm),
  species_cor = species,
  timepoint_cor = timepoint
)
ann_colors <- list(
  timepoint_cor = c(
    "D0" = "#883268",   # Purple
    "D2" = "#3E7274",  # blue
    "D4" = "#5AAA464D",  # light green
    "D15" = "#C03830",  # red
    "D30" = "#830C05" , # dark red,
    "IgG_D0" = "#EA55B3",  # Pink
    "IgG_D2"= "#0FE8F0",  # teal
    "IgG_D4"= "#57F909",  # green
    "IgG_D15"= "#F51002",  # bright red
    "IgG_D30"= "#171717"  # black
  ),
  species_cor = c(
    "H" = "#171717",  # black
    "C" = "#17171717",   # light grey
    "HR" = "#171717B0"#Dark Grey
  )
)
rownames(metadata) <- metadata$sample_cor
pheatmap(cor_Filt_RMG0_ASHG_Subset_H3K27ac_log2cpm,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = metadata[, c("species_cor", "timepoint_cor")],
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Sample-Sample Correlation (Pearson-log2CPM)")

###### CorMotif Setup (Pooled by Species + Timepoint) ######


#### Step 1: Prepare Data Set
###Remove specific replicates
samples_to_remove_CorM <- c("H20682R_D0", "H20682R_D2", "H20682R_D30")
Filt_RMG0_ASHG_Subset_H3K27ac_fc_NoRep <- Filt_RMG0_ASHG_Subset_H3K27ac_fc[, !(colnames(Filt_RMG0_ASHG_Subset_H3K27ac_fc) %in% samples_to_remove_CorM)]
dge <- DGEList(counts=Filt_RMG0_ASHG_Subset_H3K27ac_fc_NoRep)
dge <- calcNormFactors(dge,method = "TMM")
log2cpm <- cpm(dge,log = TRUE, prior.count=1)

# Set up Variables
groupid <- c(rep(1:3,7),rep(4:6,7))
compid_Cond1 <- c(rep(1,2),rep(4,2))
compid_Cond2 <- c(2,3,5,6)
compid <- cbind(compid_Cond1,compid_Cond2)

#Model Fitting
set.seed(123)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(1:10),max.iter=2000,BIC=TRUE)
motif.fitted$bic
plotIC(motif.fitted)
plotMotif(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.ASHG_Subset_H3K27ac<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.ASHG_Subset_H3K27ac)
tail(dif.pattern.ASHG_Subset_H3K27ac)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)

#Model Fitting
set.seed(91011)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(11),max.iter=2000,BIC=TRUE)
#motif.fitted$bic
#plotIC(motif.fitted)
plotMotif(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.ASHG_Subset_H3K27ac<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.ASHG_Subset_H3K27ac)
tail(dif.pattern.ASHG_Subset_H3K27ac)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)

