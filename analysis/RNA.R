#Workflowr setup####
install.packages("workflowr")

library("workflowr")

# Configure Git (only need to do once per computer)
wflow_git_config(user.name = "JohnDHurley", user.email = "jdhurley@utmb.edu")

# Start a new workflowr project
wflow_start("Evol_CM_Diff_TEs")

# Build the site
wflow_build()

#Check status on items
wflow_status()

#To Publish changes (make visible on website)
wflow_publish()

wflow_publish(c("analysis/index.Rmd", "analysis/about.Rmd", "analysis/license.Rmd"),
              "Publish the initial files for myproject")
#Commit message so future you know what you were working on and accomplishing when you published

#To Commit (save them) but not publish
wflow_commit()

#To put your code on github
wflow_use_github("JohnDHurley")

#Load Libraries####
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Cormotif")
}

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

##RNA#####
###Feature Counts#####
chimp_RNA_pantro5_fc <- read.delim("~/diff_timeline_tes/RNA/featurecounts/c_samples_counts.txt", comment.char="#",row.names=1)
human_RNA_hg38_fc <- read.delim("~/diff_timeline_tes/RNA/featurecounts/h_sample_counts.txt", comment.char="#",row.names=1)
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

#####Unfiltered####
RNA_log2cpm <- cpm(RNA_fc,log=TRUE)
hist(RNA_log2cpm,  main = "Histogram of all counts (unfiltered)",
     xlab =expression("Log"[2]*" counts-per-million"), col =4 )
boxplot(RNA_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

timepoint <- c(rep(c("D0","D2","D4","D5","D15","D30"),7),"D0R","D2R","D5R","D30R", rep(c("D0","D2","D4","D5","D15","D30"),7))
time_order <- c("D0", "D2", "D4", "D5", "D15", "D30")
species <- c(rep("H",42),rep("HR",4),rep("C",42))
smlabel <- (interaction(time,species))

T_RNA_log2cpm <- t(RNA_log2cpm)
pca_raw <- prcomp(T_RNA_log2cpm,scale=TRUE,center=TRUE)
pca_raw_meta <- as.data.frame(pca_raw$x)
pca_raw_meta$sample <- rownames(pca_raw_meta)
pca_raw_meta$species <- species
pca_raw_meta$timepoint <- timepoint
pca_raw_meta$timepoint <- factor(pca_raw_meta$timepoint, levels = time_order)
dim(T_RNA_log2cpm)

ggplot(pca_raw_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM RNA-Seq Data Unfiltered",
       x = paste0("PC1 (", round(summary(pca_raw)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_raw)$importance[2, 2] * 100, 1), "%)"))

#####RowMu>0####
row_means <- rowMeans(RNA_log2cpm)
Filt_RMG0_RNA_fc <- RNA_fc[row_means >0,]
Filt_RMG0_RNA_log2cpm <- cpm(Filt_RMG0_RNA_fc,log=TRUE)

hist(Filt_RMG0_RNA_log2cpm, main = "Histogram of filtered counts using rowMeans > 0 method",
     xlab =expression("Log"[2]*" counts-per-million"), col =5 )
boxplot(Filt_RMG0_RNA_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

######PCA####
T_Filt_RMG0_RNA_log2cpm <- t(Filt_RMG0_RNA_log2cpm)
pca_RMG0 <- prcomp(T_Filt_RMG0_RNA_log2cpm,scale=TRUE,center=TRUE)
pca_RMG0_meta <- as.data.frame(pca_RMG0$x)
pca_RMG0_meta$sample <- rownames(pca_RMG0_meta)
pca_RMG0_meta$species <- species
pca_RMG0_meta$timepoint <- timepoint
pca_RMG0_meta$timepoint <- factor(pca_RMG0_meta$timepoint, levels = time_order)
dim(T_Filt_RMG0_RNA_log2cpm)

#PC1 vs PC2
ggplot(pca_RMG0_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM RNA-Seq Data RowMean>0",
       x = paste0("PC1 (", round(summary(pca_RMG0)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_RMG0)$importance[2, 2] * 100, 1), "%)"))

######Cor_HeatMap####
cor_Filt_RMG0_RNA_log2cpm <- cor(Filt_RMG0_RNA_log2cpm, method = "pearson")
time <- factor(timepoint,levels = c("D0","D2","D4","D5","D15","D30","D0R","D2R","D5R","D30R"))
metadata <- data.frame(
  sample_cor = colnames(Filt_RMG0_RNA_log2cpm),
  species_cor = species,
  timepoint_cor = timepoint
)
ann_colors <- list(
  timepoint_cor = c(
    "D0" = "#883268",   # Purple
    "D2" = "#3E7274",  # blue
    "D4" = "#5AAA464D",  # light green
    "D5" = "#94C47D",  # Green
    "D15" = "#C03830",  # red
    "D30" = "#830C05",  # dark red
    "D0R" = "#883268",  # pink
    "D2R" = "#3E7274",  # teal
    "D5R" = "#94C47D",  # bright green
    "D30R" = "#830C05"  # bright red
  ),
  species_cor = c(
    "H" = "#171717",  # black
    "C" = "#17171717",   # light grey
    "HR" = "#171717B0"#Dark Grey
  )
)
rownames(metadata) <- metadata$sample_cor
pheatmap(cor_Filt_RMG0_RNA_log2cpm,
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
samples_to_remove_CorM <- c("H20682R_D0", "H20682R_D2", "H20682_D5", "H20682R_D30")
Filt_RMG0_RNA_fc_NoRep <- Filt_RMG0_RNA_fc[, !(colnames(Filt_RMG0_RNA_fc) %in% samples_to_remove_CorM)]
dge <- DGEList(counts=Filt_RMG0_RNA_fc_NoRep)
dge <- calcNormFactors(dge,method = "TMM")
log2cpm <- cpm(dge,log = TRUE, prior.count=1)

# Set up Variables
groupid <- c(rep(1:6,7),rep(7:12,7))
compid_Cond1 <- c(rep(1,5),rep(7,5))
compid_Cond2 <- c(2,3,4,5,6,8,9,10,11,12)
compid <- cbind(compid_Cond1,compid_Cond2)

#Model Fitting
set.seed(234)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(10),max.iter=2000,BIC=TRUE)
#motif.fitted$bic
plotIC(motif.fitted)
head(motif.fitted$bestmotif$p.post)
dif.pattern.RNA<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.RNA)
tail(dif.pattern.RNA)
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
dif.pattern.RNA<-(motif.fitted$bestmotif$p.post>0.5)
head(dif.pattern.RNA)
tail(dif.pattern.RNA)
topgenelist<-generank(motif.fitted$bestmotif$p.post)
motif_assignments <- apply(motif.fitted$bestmotif$p.post, 1, which.max)
table(motif_assignments)


#####>1.5 >4/S####
species_labels <- c(rep("Human",7*6),
                  rep("HumanR",4),
                  rep("Chimp",7*6))

length(species_labels) == ncol(RNA_log2cpm)
chimp_idx <- which(species_labels == "Chimp")
human_idx <- which(species_labels == "Human")

expr_thresh <- 1.5
min_samples_per_species <- 4

keep_genes_S <- apply(RNA_log2cpm, 1, function(x) {
  sum(x[chimp_idx] > expr_thresh) >= min_samples_per_species ||
    sum(x[human_idx] > expr_thresh) >= min_samples_per_species
})

RNA_log2cpm_S_filtered <- RNA_log2cpm[keep_genes_S, ]
dim(RNA_log2cpm_S_filtered)

filtered_genes_s <- rownames(RNA_log2cpm_S_filtered)
Filt_4pS_RNA_fc <- RNA_fc[filtered_genes_s,]
dim(Filt_4pS_RNA_fc)

Filt_4pS_RNA_log2cpm <- cpm(Filt_4pS_RNA_fc,log=TRUE)
dim (Filt_4pS_RNA_log2cpm)

T_Filt_4pS_RNA_log2cpm <- t(Filt_4pS_RNA_log2cpm)
pca_4pS <- prcomp(T_Filt_4pS_RNA_log2cpm,scale=TRUE,center=TRUE)
pca_4pS_meta <- as.data.frame(pca_4pS$x)
pca_4pS_meta$sample <- rownames(pca_4pS_meta)
pca_4pS_meta$species <- species
pca_4pS_meta$timepoint <- time
pca_4pS_meta$timepoint <- factor(pca_4pS_meta$timepoint, levels = time_order)
dim(T_Filt_4pS_RNA_log2cpm)

ggplot(pca_4pS_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "Log2cpm >1.5 in >4 samples in a species",
       x = paste0("PC1 (", round(summary(pca_4pS)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_4pS)$importance[2, 2] * 100, 1), "%)"))

####>1.5 >4/TPS####
group_labels <- c(rep(c("H_D0","H_D2","H_D4","H_D5","H_D15","H_D30"),7),
                  "HR_D0",
                  "HR_D2",
                  "HR_D5",
                  "HR_D30",
                  rep(c("C_D0","C_D2","C_D4","C_D5","C_D15","C_D30"),7)
)
group_indices <- split (seq_along(group_labels),group_labels)

expr_thresh <- 1.5
min_samples <- 4

keep_genes_TPS <- apply(RNA_log2cpm, 1, function(x) {
  any(sapply(group_indices, function(idx) sum(x[idx] > expr_thresh) >= min_samples))
})

RNA_log2cpm_TPS_filtered <- RNA_log2cpm[keep_genes_TPS, ]
dim(RNA_log2cpm_TPS_filtered)

filtered_genes_TPS <- rownames(RNA_log2cpm_TPS_filtered)
Filt_4pTPS_RNA_Fc <- RNA_fc[filtered_genes_TPS,]
dim(Filt_4pTPS_RNA_Fc)
Filt_4pTPS_RNA_Log2cpm <- cpm(Filt_4pTPS_RNA_Fc,log=TRUE)
dim (Filt_4pTPS_RNA_Log2cpm)

T_Filt_4pTPS_RNA_Log2cpm <- t(Filt_4pTPS_RNA_Log2cpm)
pca_TPS <- prcomp(T_Filt_4pTPS_RNA_Log2cpm,scale=TRUE,center=TRUE)
pca_TPS_meta <- as.data.frame(pca_TPS$x)
pca_TPS_meta$sample <- rownames(pca_TPS_meta)
pca_TPS_meta$species <- species
pca_TPS_meta$timepoint <- time
pca_TPS_meta$timepoint <- factor(pca_TPS_meta$timepoint, levels = time_order)
dim(T_Filt_4pTPS_RNA_Log2cpm)

#PC 1 vs PC2
ggplot(pca_TPS_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "log2cpm is >1.5 in >4/TPS",
       x = paste0("PC1 (", round(summary(pca_TPS)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_TPS)$importance[2, 2] * 100, 1), "%)"))


#####Subset=Pilot####
time_pilot <- c(rep(c("D0","D2","D5","D15"),14))
species_pilot <- c(rep("H",28),rep("C",28))
keep_timepoints <- time %in% c("D0", "D2", "D5", "D15")

RNA_fc_PilotSubset <- RNA_fc[,keep_timepoints]
RNA_log2cpm_PilotSubset <- cpm(RNA_fc_PilotSubset,log=TRUE)
row_means_PilotSubset <- rowMeans(RNA_log2cpm_PilotSubset)
Filt_RMG0_RNA_fc_PilotSubset <- RNA_fc_PilotSubset[row_means_PilotSubset >0,]
dim(Filt_RMG0_RNA_fc_PilotSubset)
Filt_RMG0_RNA_log2cpm_PilotSubset <- cpm(Filt_RMG0_RNA_fc_PilotSubset,log=TRUE)

T_Filt_RMG0_RNA_log2cpm_PilotSubset <- t(Filt_RMG0_RNA_log2cpm_PilotSubset)
pca_pilot <- prcomp(T_Filt_RMG0_RNA_log2cpm_PilotSubset,scale=TRUE,center=TRUE)
pca_pilot_meta <- as.data.frame(pca_pilot$x)
pca_pilot_meta$sample <- rownames(pca_pilot_meta)
pca_pilot_meta$species <- species_pilot
pca_pilot_meta$timepoint <- time_pilot
time_pilot_level <- as.factor(c("D0","D2","D5","D15"))
pca_pilot_meta$timepoint <- factor(pca_pilot_meta$timepoint, levels = time_pilot_level)
dim(T_Filt_RMG0_RNA_log2cpm_PilotSubset)

ggplot(pca_pilot_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM, D0 D2 D5 D15, rowmeans >0",
       x = paste0("PC1 (", round(summary(pca_pilot)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_pilot)$importance[2, 2] * 100, 1), "%)"))

#####Subset=D0-D30####
time_0CM <- c(rep(c("D0","D30"),14),"D0R","D30R")
species_0CM <- c(rep("C",14),rep("H",14),rep("HR",2))
keep_timepoints_0CM <- time %in% c("D0", "D30", "D0R","D30R")

RNA_fc_0CM <- RNA_fc[,keep_timepoints_0CM]
RNA_log2cpm_0CM <- cpm(RNA_fc_0CM,log=TRUE)
dim(RNA_log2cpm_0CM)
row_means_0CM <- rowMeans(RNA_log2cpm_0CM)
Filt_RMG0_RNA_fc_0CM <- RNA_fc_0CM[row_means_0CM >0,]
dim(Filt_RMG0_RNA_fc_0CM)
Filt_RMG0_RNA_log2cpm_0CM <- cpm(Filt_RMG0_RNA_fc_0CM,log=TRUE)

T_Filt_RMG0_RNA_log2cpm_0CM <- t(Filt_RMG0_RNA_log2cpm_0CM)
pca_0CM <- prcomp(T_Filt_RMG0_RNA_log2cpm_0CM,scale=TRUE,center=TRUE)
pca_0CM_meta <- as.data.frame(pca_0CM$x)
pca_0CM_meta$sample <- rownames(pca_0CM_meta)
pca_0CM_meta$species <- species_0CM
pca_0CM_meta$timepoint <- time_0CM
time_0CM_level <- as.factor(c("D0","D30","D0R","D30R"))
pca_0CM_meta$timepoint <- factor(pca_0CM_meta$timepoint, levels = time_0CM_level)
dim(T_Filt_RMG0_RNA_log2cpm_0CM)

ggplot(pca_0CM_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM, D0 D30, rowmeans >0",
       x = paste0("PC1 (", round(summary(pca_0CM)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_0CM)$importance[2, 2] * 100, 1), "%)"))

