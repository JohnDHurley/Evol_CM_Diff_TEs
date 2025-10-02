#CUT&Tag####
##H3K27ac####
###Library Loading####
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
h_samples_counts <- read.delim("~/diff_timeline_tes/H3K27ac/featurecounts/h_samples_counts.txt", comment.char="#")
c_samples_counts <- read.delim("~/diff_timeline_tes/H3K27ac/featurecounts/c_samples_counts.txt", comment.char="#")
H3K27ac_joined_fc <- left_join(h_samples_counts, c_samples_counts, by = "Geneid")
rownames(H3K27ac_joined_fc)<-H3K27ac_joined_fc$Geneid
head(H3K27ac_joined_fc)

###Union of Peaks####
H3K27ac_fc <- H3K27ac_joined_fc[ , !(names(H3K27ac_joined_fc) %in% c("Geneid","Chr.x","Start.x","End.x","Strand.x","Length.x","Chr.y","Start.y","End.y","Strand.y","Length.y"))]
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
colnames(H3K27ac_fc)<-col_names
head(H3K27ac_fc)

H3K27ac_fc_Filt <- H3K27ac_fc[ , !(names(H3K27ac_fc) %in% c(
  "IgG_H21792_D0",
  "IgG_H21792_D2",
  "IgG_H21792_D4",
  "IgG_H21792_D5",
  "IgG_H21792_D15",
  "IgG_C40280_D0",
  "IgG_C40280_D2",
  "IgG_C40280_D4",
  "IgG_C40280_D5",
  "IgG_C40280_D15",
  "IgG_C40280_D30",
  "C40210_D15",
  "C3649_D4"

))]
#####Unfiltered####

H3K27ac_log2cpm <- cpm(H3K27ac_fc_Filt,log=TRUE)
hist(H3K27ac_log2cpm,  main = "Histogram of all counts (unfiltered)",
     xlab =expression("Log"[2]*" counts-per-million"), col =4 )
boxplot(H3K27ac_log2cpm, main = "Boxplots of log cpm per sample",xaxt = "n", xlab= "")
x_axis_labels(labels = label, every_nth = 1, adj=0.7, srt =90, cex =0.4)

timepoint <- c(
  rep(c("D0","D2","D4","D5","D15","D30"),6),
  "D0R",
  "D2R",
  "D5R",
  "D30R",
  rep(c("D0","D2","D4","D5","D15","D30"),1),
  rep(c("D0","D2","D4","D5","D15","D30"),2),
  "D0","D2","D4","D5","D15",
  rep(c("D0","D2","D4","D5","D15","D30"),1),
  rep(c("D0","D2","D4","D5","D15","D30"),1),
  "D0","D2","D4","D5","D30",
  rep(c("D0","D2","D4","D5","D15","D30"),1))
dim(timepoint)
time_order <- c("D0", "D2", "D4", "D5", "D15", "D30")
species <- c(rep("H",36),
             rep("HR",4),
             rep("H",6),
             rep("C",40))

T_H3K27ac_log2cpm <- t(H3K27ac_log2cpm)
pca_raw <- prcomp(T_H3K27ac_log2cpm,scale=TRUE,center=TRUE)
pca_raw_meta <- as.data.frame(pca_raw$x)
pca_raw_meta$sample <- rownames(pca_raw_meta)
pca_raw_meta$species <- species
pca_raw_meta$timepoint <- timepoint
pca_raw_meta$timepoint <- factor(pca_raw_meta$timepoint, levels = time_order)
dim(T_H3K27ac_log2cpm)

ggplot(pca_raw_meta, aes(x = PC1, y = PC2, color = timepoint, shape = species)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "PCA of log2CPM H3K27ac Data Unfiltered",
       x = paste0("PC1 (", round(summary(pca_raw)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_raw)$importance[2, 2] * 100, 1), "%)"))

#####Cor_HeatMap####
cor_H3K27ac_log2cpm <- cor(H3K27ac_log2cpm, method = "pearson")
dim(cor_H3K27ac_log2cpm)
time <- factor(timepoint,levels = c("D0","D2","D4","D5","D15","D30","D0R","D2R","D5R","D30R"))
metadata <- data.frame(
  sample_cor = colnames(H3K27ac_log2cpm),
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
    "D0R" = "#FB019E",  # pink
    "D2R" = "#05EEEE",  # teal
    "D5R" = "#51F701",  # bright green
    "D30R" = "#F71002"  # bright red
  ),
  species_cor = c(
    "H" = "#171717",  # black
    "C" = "#17171717",   # light grey
    "HR" = "#171717B0"#Dark Grey
  )
)
rownames(metadata) <- metadata$sample_cor
pheatmap(cor_H3K27ac_log2cpm,
         fontsize_row = 5,
         fontsize_col = 5,
         annotation_col = metadata[, c("species_cor", "timepoint_cor")],
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "H3K27ac Sample-Sample Correlation (Pearson-log2CPM)")

###### CorMotif Setup (Pooled by Species + Timepoint) ######


#### Step 1: Prepare Data Set
###Remove specific replicates
samples_to_remove_CorM <- c("H20682R_D0", "H20682R_D2", "H20682_D5", "H20682R_D30")
H3K27ac_fc_Filt_NoRep <- H3K27ac_fc_Filt[, !(colnames(H3K27ac_fc_Filt) %in% samples_to_remove_CorM)]
dge <- DGEList(counts=H3K27ac_fc_Filt_NoRep)
dge <- calcNormFactors(dge,method = "TMM")
log2cpm <- cpm(dge,log = TRUE, prior.count=1)

# Set up Variables
groupid <- c(rep(1:6,7),rep(7:12,2),c(7:11),rep(7:12,2),c(8:12),rep(7:12,1))
compid_Cond1 <- c(rep(1,5),rep(7,5))
compid_Cond2 <- c(2,3,4,5,6,8,9,10,11,12)
compid <- cbind(compid_Cond1,compid_Cond2)

#Model Fitting
set.seed(123)
motif.fitted<-cormotiffit(log2cpm,groupid,compid,K=(1:20),max.iter=3000,BIC=TRUE)
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


