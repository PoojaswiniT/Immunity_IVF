#' CELL TYPE PROPORTIONS IN MALE AND FEMALE PBMCs
#' This script uses the pbmc datasets generated from the Oelen et al 2022 paper
#' datasets: https://eqtlgen.org/sc/datasets/1m-scbloodnl.html
#' paper: https://www.nature.com/articles/s41467-022-30893-5
#' The datsets include the singlle cells captured and sequenced using v2 and v3 reagets through the 10x Genomicss controller
#' After downloading, in each folder of v2 and v3 chemistry there must be a barcode file, features file and the matrix file.
#' The cell type info and the pathogem information were also downloaded from the website.
#'
#' The quality control steps follow the standard Seurat pipeline and then the results are anlaysed and visualised to see the differnt cell types found in PBMCs (Peripheral Blood Mononuclear Cells)
#' The cells are assigned sex information from another function that applies cellXY to predict the sex.
#' The analysis is then stratified by sex to see the difference in cell type proportions between the two biological sexes - males and females.
#'
# LOAD LIBRARIES


library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(dittoSeq)
library(speckle)
library(patchwork)
library(ComplexHeatmap)
library(tidyr)
library(colorRamp2)
library(readxl)
library(scCustomize)
library(reshape2)
library(knitr)
library(tibble)
library(cellXY)
library(SingleCellExperiment)
library(caret)
library(randomForest)
library(viridis)
library(ggplotify)
library(lme4)
library(muscatWrapper)
library(writexl)
library(ggpubr)
library(kableExtra)
library(Azimuth)
library(SeuratData)
# SET SEED FOR REPRODUCIBILITY
set.seed(123)

# LOAD INPUT DATA
pbmc_v2_data <- Read10X("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/n_qc_v2_data", gene.column = 1, cell.column = 1)
pbmc_v3_data <-  Read10X("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/n_qc_v3_data", gene.column = 1, cell.column = 1)
#pre_qc <- Read10X("/home/inf-21-2023/thesis/scRNAseq_data_oelen/pre_QC", gene.column = 1, cell.column = 1)

# CREATE THE SEURAT OBJECTS
pbmc_v2 <- CreateSeuratObject(counts = pbmc_v2_data, project = "pbmc_v2")
pbmc_v3 <- CreateSeuratObject(counts = pbmc_v3_data, project = "pbmc_v3")
#pre_qc_so <- CreateSeuratObject(counts = pre_qc, project = "pbmc_pre_qc")
rm(pbmc_v2_data, pbmc_v3_data)
gc()

# READ CELL TYPE INFO
# This version of the TSV file has been updated to include the correct annotations based on the existing cell type annotations in the older file
#cell.info <- fread("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/info/1M_cell_types.tsv")
cell.info <- fread("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/info/1M_cell_types_updated.tsv")

# READ PATHOGEN INFO
pathogen.info <- fread("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/info/1M_assignments_conditions_expid.tsv")

# READ THE META DATA
meta_data <- read_xlsx("/home/inf-21-2023/thesis/scRNAseq_data_oelen/normalised_qc/info/meta_data.xlsx")

# ADD THE INFO TO THE SEURAT OBJECTS
pbmc_v2$cell_type <- cell.info$cell_type_lowerres[match(colnames(pbmc_v2), cell.info$barcode)]
pbmc_v2$cell_subtype <- cell.info$cell_type[match(colnames(pbmc_v2), cell.info$barcode)]
pbmc_v2$subject <- pathogen.info$assignment[match(colnames(pbmc_v2), pathogen.info$barcode)]
pbmc_v2$timepoint <- pathogen.info$timepoint[match(colnames(pbmc_v2), pathogen.info$barcode)]

pbmc_v3$cell_type <- cell.info$cell_type_lowerres[match(colnames(pbmc_v3), cell.info$barcode)]
pbmc_v3$cell_subtype <- cell.info$cell_type[match(colnames(pbmc_v3), cell.info$barcode)]
pbmc_v3$subject <- pathogen.info$assignment[match(colnames(pbmc_v3), pathogen.info$barcode)]
pbmc_v3$timepoint <- pathogen.info$timepoint[match(colnames(pbmc_v3), pathogen.info$barcode)]

# COMBINE BOTH V2 AND V3 CHEMISTRY SEURAT OBJECTS
combined <- merge(
  x = pbmc_v2,
  y = pbmc_v3,
  add.cell.ids = c(
    "pbmc_v2", "pbmc_v3")
)
combined_so <- JoinLayers(combined)
#saveRDS(combined_so, file ="v2_v3_pbmc.rds")
rm(pbmc_v2, pbmc_v3)
gc()

combined_so <- readRDS("/home/inf-21-2023/thesis/scRNAseq_data_oelen/v2_v3_pbmc.rds")

# If needed the cell type annotation can be done using Azimuth
options(future.globals.maxSize = 1000 * 1024^3)

combined_so <- RunAzimuth(combined_so, reference = "pbmcref")

DimPlot(combined_so, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
#saveRDS(combined_so, file = "/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_cell_type_annotates.rds")

combined_so <- readRDS("/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_cell_type_annotates.rds")

# CHECK FOR <NA>

which(is.na(combined_so@meta.data$subject))

# REMOVE THE <NA>
combined_so <- subset(combined_so, cells = setdiff(Cells(combined_so),"pbmc_v3_TTGTTTGCATTCTGTT_190122_lane2"))

# LOAD THE SEX INFO FROM THE META DATA
sex_info <- meta_data$sex

# SORT THE SUBJECTS IN INCREASING ORDER AND ASSIGN THE SEX TO THE META DATA OF THE SEURAT OBJECT
sorted_subjects <- sort(unique(combined_so@meta.data$subject))
if (length(sorted_subjects) == length(sex_info)) {
  combined_so@meta.data$sex <- sex_info[match(combined_so@meta.data$subject, sorted_subjects)]
} else {
  stop("Length mismatch between subjects and sex metadata")
}


# SUBEST THE UNTREATED SAMPLES ONLY
pbmc.UT <- subset(combined_so, subset = timepoint == "UT")
rm(combined_so)
gc()

## QUALITY CONTROL ##
pbmc.UT[["percent.mt"]] <- PercentageFeatureSet(pbmc.UT, pattern = "^MT-")

# Percentage hemoglobin genes
pbmc.UT <- PercentageFeatureSet(pbmc.UT, "^HB[^(P)]", col.name = "percent_hb")

# VISUALIZE QC METRICS BEFORE FILTERING
v.plot <- VlnPlot(pbmc.UT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_hb"), ncol = 4)
v.plot <- v.plot + plot_annotation(title = "Before filtering")
ggsave(filename = "/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot.png", plot = v.plot, width = 10)

# FILTERING AND VISUALISING
pbmc.UT <- subset(pbmc.UT, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 8 & percent_hb < 2)

v.plot2 <- VlnPlot(pbmc.UT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_hb"), ncol = 4)
v.plot2 <- v.plot2 + plot_annotation(title = "After filtering")
ggsave(filename = "/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot2.png", plot = v.plot2, width = 10)

# NORMALIZATION
pbmc.UT <- NormalizeData(pbmc.UT, normalization.method = "LogNormalize", scale.factor = 10000)

# FIND VARIABLE FEATURES
pbmc.UT <- FindVariableFeatures(pbmc.UT, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc.UT)

# SCALING THE DATA
pbmc.UT <- ScaleData(pbmc.UT, features = all.genes)
gc()

pbmc.UT <- readRDS("/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_UT_sex_assigned.rds")

# PCA
pbmc.UT<- RunPCA(pbmc.UT, features = VariableFeatures(object = pbmc.UT))
ElbowPlot(pbmc.UT)

# FIND CLUSTERS
pbmc.UT <- FindNeighbors(pbmc.UT, dims = 1:10)
pbmc.UT<- FindClusters(pbmc.UT, resolution = c(0.5))
#DimPlot(pbmc.UT, group.by = "RNA_snn_res.0.1", label = TRUE)

# NON LINEAR DIMENSIONAL REDUCTION
pbmc.UT <- RunUMAP(pbmc.UT, dims = 1:10)

c29 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gray20", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  "purple", "pink", "cyan", "magenta", "yellow", "gray", "lightblue", "lightgreen"
)

# PLot UMAP for the Azimuth annotation
DimPlot(pbmc.UT, group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, raster = FALSE,  cols = c29) + theme(legend.position = "bottom")

# PLot UMAP for the updated cell type annotation
DimPlot(pbmc.UT, reduction = "umap", group.by = "cell_subtype", label = TRUE, cols = c29) + theme(legend.position = "bottom")
#saveRDS(pbmc.UT, "/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_UT.rds")

# LOAD THE FUNCTION THAT ASSIGNES SEX BASED ON CellXY
source("/home/inf-21-2023/thesis/scRNAseq_data_oelen/sex_assign.R")
pbmc_out <- SexAssign(pbmc.UT)

# you get multiple lists of results, with seurat object, the table with sex assignment and plots, so add the seurat obj back to the original obj name
pbmc.UT <- pbmc_out[[1]]
gc()
#saveRDS(pbmc.UT, "/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_UT_sex_assigned.rds")


pbmc.UT <- readRDS("/home/inf-21-2023/thesis/scRNAseq_data_oelen/pbmc_UT_sex_assigned.rds")

all.genes <- rownames(pbmc.UT)

# NUMBER OF MALES AND FEMALES AFTER SEX ASSIGNMENT
table(unique(pbmc.UT@meta.data[, c("subject", "Predicted_sex")])$Predicted_sex)

# View all the cell types
unique(pbmc.UT@meta.data$cell_type)

# Visualise the number of cells per cell type (also can be seen in umap plots)
cell_counts <- table(pbmc.UT@meta.data$cell_subtype)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("cell_type", "count")

# Plot the histogram
ggplot(cell_counts_df, aes(x = cell_type, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Number of Cells per cell type",
       x = "Cell type",
       y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# there are very few unkown cell types, so remove them
#pbmc.UT <- subset(pbmc.UT, subset = cell_type != "unknown")

# Visualise the  number of cells per cell types for each individual
cell_counts <- table(pbmc.UT@meta.data$subject, pbmc.UT@meta.data$cell_subtype, pbmc.UT@meta.data$Predicted_sex)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("individual", "cell_type","sex", "count")

# Plot the histogram
ggplot(cell_counts_df, aes(x = individual, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  facet_wrap(~cell_type, scales = "free_x") +
  labs(title = "Number of Cells per Individual",
       x = "Individual",
       y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Visualse the cell type distribuion for male and female individuals
dittoBarPlot(pbmc.UT, var = "cell_type", group.by = "subject", main = "Cell type distribution between females and males", xlab = "cell type", ylab = NULL,split.by = "Predicted_sex", split.adjust = list(scales = "free_x"))

# Store the data table
table <- dittoBarPlot(pbmc.UT, var = "cell_subtype", group.by = "subject", main = "Cell type distribution between females and males", xlab = "cell type", ylab = NULL,split.by = "Predicted_sex", split.adjust = list(scales = "free_y"))$data

# Sort the table by the highest cell type count by percent
donors_sorted <- table %>%
  filter(label == "th2 CD4T") %>%
  arrange(percent) %>%
  pull(grouping)

# Visualise the data again to see a pattern
table %>%
  mutate(grouping = factor(grouping, levels = donors_sorted)) %>%
  ggplot(aes(x = grouping, y= percent, fill =label)) + geom_bar(position ="stack", stat = "identity")+
  facet_grid(~ Predicted_sex, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Visualise the differences through box plots and a T test
plot <- table %>%
  mutate(grouping = factor(grouping, levels = donors_sorted)) %>%
  ggplot(aes(x = Predicted_sex, y = percent, fill = Predicted_sex)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7) +
    stat_boxplot(geom = "errorbar", width = 0.25) +     # Add error bars (showing 95% confidence interval)
    geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
    facet_wrap(~ label, scales = "free_y") +
    scale_fill_manual(values = c("M" = "#595959", "F" = "#eb4fb2")) +
    labs(
      title = "Distribution of Cell Types by Sex",
      x = "Sex",
      y = "Percent (%)",
      fill = "Sex"
    ) +
    theme_minimal()+
  stat_compare_means(method = "t.test",label = "p.signif", label.y = 1.1, label.x = 1.2)

ggsave(plot, filename = "/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/distribution_of_cell_types_by_sex.pdf", width = 20, height = 20)

# CALCULATE THE CELL TYPE PROPORTIONS USING SPECKLE
props <- getTransformedProps(clusters = pbmc.UT$cell_subtype,
                             sample = pbmc.UT$subject, transform = "logit")

# CONVERT IT TO A DATAFRAME
freq_df <- props$Proportions %>%
        as.data.frame() %>%
        tidyr::pivot_wider(names_from = "clusters",
                    values_from = "Freq",
                    ) %>%
        tibble::column_to_rownames("sample") %>%
        na.omit()

freq_matrix <- as.matrix(freq_df)


cell_counts <- table(pbmc.UT@meta.data$subject, pbmc.UT@meta.data$cell_subtype, pbmc.UT@meta.data$Predicted_sex)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("individual", "cell_type","sex", "count")

# Compute the total number of cells per individual
total_cell_abundance <- aggregate(count ~ individual, data = cell_counts_df, sum)
total_abundance_vector <- total_cell_abundance$count[match(rownames(freq_matrix), total_cell_abundance$individual)]


subject_sex <- unique(pbmc.UT@meta.data[, c("subject", "Predicted_sex")])
rownames(subject_sex) <- subject_sex$subject

# Reorder the matrix
freq_matrix <- freq_matrix[order(total_abundance_vector, decreasing = TRUE), ]
subject_sex <- subject_sex[order(match(subject_sex$subject, rownames(freq_matrix))), ]

column_ha = rowAnnotation(
  sex = subject_sex$Predicted_sex,
  col = list(sex =  c("M" = "#595959", "F" = "#eb4fb2"))
)
heatmap <- Heatmap(freq_matrix,
    name = "Frequency",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_names_gp = gpar(fontsize=8),
    left_annotation = column_ha,
    col = viridis::viridis(100),
    row_split = subject_sex$Predicted_sex  # Split rows by sex
)

pdf("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/speckle_heatmap_split_by_sex_cell_type.pdf", width = 20, height = 20)
draw(heatmap)
dev.off()


#by both male and female seperate
pbmc.male <- subset(pbmc.UT, subset = Predicted_sex =="M")
pbmc.female <- subset(pbmc.UT, subset = Predicted_sex =="F")

ggplot(cell_counts_df, aes(x = sex, y = count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Abundance of Cell Types by Sex",
       x = "Sex",
       y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



props_m <- getTransformedProps(clusters = pbmc.male$cell_type,
                             sample = pbmc.male$subject, transform = "logit")

# CONVERT IT TO A DATAFRAME
freq_df_m <- props_m$Proportions %>%
        as.data.frame() %>%
        tidyr::pivot_wider(names_from = "clusters",
                    values_from = "Freq",
                    ) %>%
        tibble::column_to_rownames("sample") %>%
        na.omit()
freq_matrix_m <- as.matrix(freq_df_m)

props_f <- getTransformedProps(clusters = pbmc.female$cell_type,
                             sample = pbmc.female$subject, transform = "logit")

# CONVERT IT TO A DATAFRAME
freq_df_f <- props_f$Proportions %>%
        as.data.frame() %>%
        tidyr::pivot_wider(names_from = "clusters",
                    values_from = "Freq",
                    ) %>%
        tibble::column_to_rownames("sample") %>%
        na.omit()

freq_matrix_f <- as.matrix(freq_df_f)
# Rearrange by specific cell type
freq_matrix_m <- freq_matrix_m[order(freq_matrix_m[, "CD8T"], decreasing = TRUE), ]
freq_matrix_f <- freq_matrix_f[order(freq_matrix_f[, "CD4T"], decreasing = TRUE), ]

# cell subtype
freq_matrix_m <- freq_matrix_m[order(freq_matrix_m[, "NKdim"], decreasing = TRUE), ]
freq_matrix_f <- freq_matrix_f[order(freq_matrix_f[, "th2 CD4T"], decreasing = TRUE), ]

freq_matrix_combined <- rbind(freq_matrix_m, freq_matrix_f)

subject_sex <- unique(pbmc.UT@meta.data[, c("subject", "Predicted_sex")])
rownames(subject_sex) <- subject_sex$subject

subject_sex <- subject_sex[order(match(subject_sex$subject, rownames(freq_matrix_combined))), ]
#subject_sex <- subject_sex[order(match(subject_sex$subject, rownames(freq_matrix))), ]

# Rearrange by specific cell type
#freq_matrix <- freq_matrix[order(freq_matrix[, "th2 CD4T"], decreasing = TRUE), ]

# Order by sex
#subject_sex <- subject_sex[order(subject_sex$Predicted_sex), ]
#freq_matrix <- freq_matrix[match(rownames(subject_sex), rownames(freq_matrix)), ]

column_ha = rowAnnotation(
  sex = subject_sex$Predicted_sex,
  col = list(sex =  c("M" = "#595959", "F" = "#eb4fb2"))
)

#  GENERATE THE HEATMAP
heatmap <- Heatmap(freq_matrix_combined,
    name = "Frequency",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_names_gp = gpar(fontsize=8),
    left_annotation = column_ha,
    col = viridis::viridis(100),
    row_split = subject_sex$Predicted_sex
)

pdf("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/speckle_heatmap_by_abundance.pdf", width = 20, height = 20)
draw(heatmap)
dev.off()


###
mean_male <- colMeans(freq_matrix_m, na.rm = TRUE)
mean_female <- colMeans(freq_matrix_f, na.rm = TRUE)

# Calculate log2(male/female) ratio
log2_ratios <- log2(mean_male / mean_female)

# Convert to dataframe for plotting
log2_df <- data.frame(
  cell_type = names(log2_ratios),
  log2_ratio = log2_ratios
)

# Plot
ggplot(log2_df, aes(x = log2_ratio, y = reorder(cell_type, log2_ratio), color = cell_type)) +
  geom_point(size = 5) +
  theme_minimal() +
  labs(title = "Cell Type Frequency Sex Bias",
       x = "Log2(Male/Female)", y = "Cell Type") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Reference line at 0
  scale_color_viridis_d()  # Automatically assigns distinct colors
###


#variance info
cellfreq <- props$Proportions %>%
        as.data.frame() %>%
        tidyr::pivot_wider(names_from = "clusters",
                    values_from = "Freq",
                    ) %>%
        tibble::column_to_rownames("sample")
meta.data <- pbmc.UT@meta.data
sample_info <- meta.data %>%
  distinct(Predicted_sex, subject, orig.ident) %>%
  mutate(project = gsub("-", "_",
                        orig.ident)) %>%
        arrange(subject)

adonis.res <- vegan::adonis2(cellfreq ~ Predicted_sex + orig.ident, data = sample_info, method = "euclidean", permutations = 999)

variance_plot <- data.frame(
  Component = c("Predicted_sex", "orig.ident", "Residual"),
  SumOfSqs = c(adonis.res$SumOfSqs[1], adonis.res$SumOfSqs[2], adonis.res$SumOfSqs[3])
)

# Convert Component to factor for ordering
variance_plot$Component <- factor(variance_plot$Component, levels = c("Residual", "orig.ident", "Predicted_sex"))

# Plot stacked bar plot
ggplot(variance_plot, aes(x = "", y = SumOfSqs, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Sum of Squares", x = "", title = "Variance Partitioning") +
  theme_classic() +
  scale_fill_manual(values = c("Residual" = "#ff4500", "orig.ident" = "#00ff7f", "Predicted_sex" = "#836fff"),
                    labels = c("Residual" = "Residual", "orig.ident" = "Chemistry", "Predicted_sex" = "Sex"))+
                    theme(              axis.text.y =element_blank(),
              axis.ticks.y = element_blank(),
  )








































#########
receptor_genes <- c("^ESR1$", "^ESR2$", "^PGR$","^NR3C3$","^PR$", "^AR$", "^GPER", "^LGR1","^FSHRO","^ODG1", "^FSHR", "^LHR", "^LCGR",
"^LGR2", "^ULG5", "^LHRH","^LHRHR","LSH-R","HHG", "^GNRHR","^GRHR", "^LHRHR", "^GRH$", "^LHCGR" ,"^HHG", "^LCGR", "^LGR2", "LHR", "^ULG5", "^ODG4", "^FSH Receptor", "P23945-1", "P23945", "LH/CG-R","LSH-R","HHG", "ULG5", "LGR2", "LCGR" )

filtered_genes <- grep(paste(receptor_genes,  collapse = "|"), all.genes, value = TRUE)

#cell_types_order <- c("B", "plasma B", "pDC", "mDC", "NKdim", "NKbright", "NK", "th1 CD4T", "th2 CD4T", "naive CD4T", "memory CD4T","reg CD4T", "naive CD8T ","memory CD8T","mono 1", "mono 2" ,"mono 4", "hemapoietic stem", "megakaryocyte", "unkown")
cell_types_order <- c(
  # B Cells
  "B", "plasma B",

  # T Cells (CD4 and CD8 T Cells)
  "th1 CD4T", "th2 CD4T", "naive CD4T", "memory CD4T", "reg CD4T",
  "naive CD4T transitioning to stim", "cyto CD4T",

  # T Cells (CD8 T Cells)
  "naive CD8T", "memory CD8T", "memory CD8T left and naive CD8T right",

  # Double Negative T Cells
  "double negative T",

  # Monocytes
  "mono 1", "mono 2", "mono 3", "mono 4",

  # NK Cells
  "NKdim", "NKbright", "NK",

  # Other
  "hemapoietic stem", "megakaryocyte", "pDC", "mDC"
)
#cell_types_order <- c("B", "plasma B", "DC" , "NK", "CD4T", "CD8T","monocyte", "hemapoietic stem", "megakaryocyte")
receptors_order <- c("ESR1", "ESR2", "GPER1", "AR")


# Annotate the cell types by receotor expression positive or negative
expr_data <- FetchData(pbmc.UT, vars = c("subject", filtered_genes), layer = "data")
expr_data$cell_subtype <- pbmc.UT$cell_subtype
expr_data$subject <- pbmc.UT$subject

# Prepare the data
PrepareBinaryHeatmapData <- function(expr_data, cell_types_order, receptors_order, threshold = 0) {
  # First, aggregate expression by cell subtype and subject
  agg_expr <- expr_data %>%
    group_by(subject, cell_subtype) %>%
    summarise(across(all_of(receptors_order), ~mean(., na.rm = TRUE))) %>%
    ungroup()

  # Then determine positivity
  binary_data <- agg_expr %>%
    mutate(across(all_of(receptors_order), ~as.numeric(. > threshold))) %>%
    group_by(cell_subtype) %>%
    summarise(across(all_of(receptors_order), ~as.numeric(mean(.) > 0))) %>%
    arrange(factor(cell_subtype, levels = cell_types_order)) %>%
    column_to_rownames("cell_subtype")

  return(as.matrix(binary_data))
}

receptor_status <- PrepareBinaryHeatmapData(expr_data, cell_types_order, receptors_order)

# Update receptor_status in the meta.data of pbmc.UT
for (receptor in receptors_order) {
  pbmc.UT@meta.data[[paste0(receptor, "_status")]] <- ifelse(pbmc.UT$cell_subtype %in% rownames(receptor_status[receptor_status[, receptor] == 1, ]), "Positive", "Negative")
}

# Heatmap
CreateBinaryReceptorHeatmap <- function(expr_data, cell_types_order, receptors_order) {

  heatmap_data <- PrepareBinaryHeatmapData(expr_data, cell_types_order, receptors_order)

  custom_breaks <- c(0, 1)
  custom_colors <- c("blue", "red")  # Blue for 0, Red for 1

  # Create heatmap
  pheatmap(
    t(heatmap_data),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = c("blue", "red"),  # Only two colors for binary data
    breaks = c(-0.01, 0.5, 1.01),  # Define breaks tightly around 0 and 1
    legend_breaks = c(0, 1),  # Show only 0 and 1 in the legend
    legend_labels = c("Negative", "Positive"),  # Meaningful labels for 0 and 1
    main = "Heatmap of Receptor Positivity",
    fontsize = 10,
    cellwidth = 20,
    cellheight = 15)
}

# Generate the heatmap
CreateBinaryReceptorHeatmap(expr_data, cell_types_order, receptors_order)

# Convert the matrix to long format
expr_long <- melt(expr_data, id.vars = c("cell_subtype", "subject"),
                   variable.name = "gene", value.name = "expression")

expr_long <- expr_long %>%
  group_by(cell_subtype, subject, gene) %>%
  summarise(percent_expression = sum(expression > 0) / n() * 100) %>%
  ungroup()

expr_long$subject <- factor(expr_long$subject, levels = unique(expr_long$subject))


violin_plot <- ggplot(expr_long, aes(x = cell_subtype, y = percent_expression, fill = subject)) +
  geom_violin() +
  labs(x = "Cell Type", y = "Percent Expression per Subject", title = "Violin Plot of Percent Expression per Subject by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot_percent_expression1.pdf", violin_plot, width = 10, height = 8)


heatmap <- ggplot(expr_long, aes(x = subject, y = gene, fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Rotate for visibility
  labs(x = "Subject", y = "Receptor Genes", fill = "Expression")

ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/receptor_expression_by_cell_type2.pdf", heatmap, width = 20, height = 20)

# Split the data to male and female data sets to see sex differences
pbmc_male <- subset(pbmc.UT, subset = Predicted_sex == "M")
pbmc_female <- subset(pbmc.UT, subset = Predicted_sex == "F")

expr_data_male <- FetchData(pbmc_male, vars = c("subject", filtered_genes), layer = "data")
expr_data_male$cell_subtype <- factor(pbmc_male$cell_subtype, levels = cell_types_order)
expr_data_male$subject <- pbmc_male$subject

expr_long_male <- melt(expr_data_male, id.vars = c("cell_subtype", "subject"),
                   variable.name = "gene", value.name = "expression")

# For heatmap
expr_long_male <- expr_long_male %>%
  group_by(cell_subtype, subject, gene) %>%
  summarise(
    expressing_cells = sum(expression > 0),  # Count cells expressing the gene
    total_cells = n(),
  ) %>%
  ungroup()

heatmap_male <- ggplot(expr_long_male, aes(x = factor(subject, levels = unique(subject)), y = factor(gene, levels = receptors_order), fill = expressing_cells)) +
  geom_tile() +
  scale_fill_viridis_c(option = "cividis", na.value = "white") +
  facet_wrap(~ cell_subtype, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Subject", y = "Receptor Genes", fill = "Expression")


# only for boxplots and violin plot, not for heatmaps
expr_long_male <- expr_long_male %>%
  group_by(cell_subtype, subject, gene) %>%
  summarise(
    expressing_cells = sum(expression > 0),
    total_cells = n(),
    percent_expression = (expressing_cells / total_cells) * 100
  ) %>%
  ungroup()

expr_long_male$subject <- factor(expr_long_male$subject, levels = unique(expr_long_male$subject))

expr_data_female <- FetchData(pbmc_female, vars = c("subject", filtered_genes), layer = "data")
expr_data_female$cell_subtype <- factor(pbmc_female$cell_subtype, levels = cell_types_order)
expr_data_female$subject <- pbmc_female$subject


expr_long_female <- melt(expr_data_female, id.vars = c("cell_subtype", "subject"),
                   variable.name = "gene", value.name = "expression")

expr_long_female <- expr_long_female %>%
  group_by(cell_type, subject, gene) %>%
  summarise(
    expressing_cells = sum(expression > 0),
    total_cells = n()
  ) %>%
  ungroup()


# only for box plots and violin plots, not for heatmaps
expr_long_female <- expr_long_female %>%
  group_by(cell_subtype, subject, gene) %>%
  summarise(
    expressing_cells = sum(expression > 0),
    total_cells = n(),
    percent_expression = (expressing_cells / total_cells) * 100
  ) %>%
  ungroup()

expr_long_female$subject <- factor(expr_long_female$subject, levels = unique(expr_long_female$subject))


boxplot_male <- ggplot(expr_long_male, aes(x = gene, y = expression)) +
  geom_boxplot() +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Gene", y = "Expression", title = "Boxplot of Gene Expression in Males")
ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/box_plot_male.pdf", boxplot_male, width = 25, height = 20)

# Boxplot for female data to check for outliers
boxplot_female <- ggplot(expr_long_female, aes(x = gene, y = expression)) +
  geom_boxplot() +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Gene", y = "Expression", title = "Boxplot of Gene Expression in Females")
ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/box_plot_female.pdf", boxplot_female, width = 25, height = 20)


heatmap_male <- ggplot(expr_long_male, aes(x = factor(subject, levels = unique(subject)), y = factor(gene, levels = receptors_order), fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Subject", y = "Receptor Genes", fill = "Expression")

ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/receptor_expression_by_cell_type_male2.pdf", heatmap_male, width = 25, height = 20)

heatmap_female <- ggplot(expr_long_female, aes(x = factor(subject, levels = unique(subject)), y = factor(gene, levels = receptors_order), fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Subject", y = "Receptor Genes", fill = "Expression")

ggsave("/home/inf-21-2023/thesis/scRNAseq_data_oelen/receptor_expression_by_cell_type_female2.pdf", heatmap_female, width = 25, height = 20)

for (cell_type in unique(expr_long_male$cell_type)) {
  cell_data <- expr_long_male %>% filter(cell_type == !!cell_type) %>% mutate(gene = factor(gene, levels = receptors_order))
  violin_plot <- ggplot(cell_data, aes(x = gene, y = percent_expression, fill = gene)) +
    geom_violin() +
    geom_point(position = position_jitter(width = 0.2), size = 1) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    labs(x = "Receptor gene", y = "Percent Expression", title = paste("Plot of Percent Expression per Cell Type - Male",cell_type, "cell")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot_percent_expression_m", cell_type, ".pdf"), violin_plot, width = 10, height = 8)
}


for (cell_type in unique(expr_long_female$cell_type)) {
  cell_data <- expr_long_female %>% filter(cell_type == !!cell_type) %>% mutate(gene = factor(gene, levels = receptors_order))
  violin_plot <- ggplot(cell_data, aes(x = gene, y = percent_expression, fill = gene)) +
    geom_violin() +
    geom_point(position = position_jitter(width = 0.2), size = 1) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    labs(x = "Receptor gene", y = "Percent Expression", title = paste("Plot of Percent Expression per Cell Type - Female", cell_type, "cell")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot_percent_expression_f", cell_type, ".pdf"), violin_plot, width = 10, height = 8)
}

for (gene in unique(expr_long_male$gene)) {
  gene_data <- expr_long_male %>% filter(gene %in% filtered_genes)
  violin_plot <- ggplot(gene_data, aes(x = percent_expression, y = cell_type, fill = gene)) +
    geom_bar(position ="dodge", stat = "identity") +
    geom_point(aes(group = gene), position = position_dodge(width = 0.75), size = 1) +
    #stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
    labs(x = "Cell Type", y = "Percent Expression", title = paste("Plot of Percent Expression per Cell Type - Male")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot_percent_expression_m", gene, ".pdf"), violin_plot, width = 10, height = 8)
}

for (cell_type in unique(expr_long_female$cell_type)) {
  cell_data <- expr_long_female %>% filter(cell_type == !!cell_type)
  box_plot <- ggplot(cell_data, aes(x = gene, y = percent_expression, fill = gene)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 95, size = 4, fill = "white") +  # Mean
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 1, color = "black") +  # SE bars
    labs(title = paste("Plot of Percent Expression with Mean ± SEper Cell Type - Female", cell_type, "cell") ,
         y = "Percent Expression", x = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/box_plot_percent_expression_f", cell_type, ".pdf"), box_plot, width = 10, height = 8)
}

# Annotate the cell types by receotor expression positive or negative
expr_data <- FetchData(pbmc.UT, vars = c("subject", filtered_genes, "Predicted_sex"), layer = "data")
expr_data$cell_subtype <- pbmc.UT$cell_subtype
expr_data$subject <- pbmc.UT$subject

# Convert the matrix to long format
expr_long <- melt(expr_data, id.vars = c("cell_subtype", "subject", "Predicted_sex"),
                   variable.name = "gene", value.name = "expression")

expr_long <- expr_long %>%
  group_by(cell_subtype, subject, gene, Predicted_sex) %>%
  summarise(percent_expression = sum(expression > 0) / n() * 100) %>%
  ungroup()

# Create a new column for combined labels (Sex_Receptor)
combined_expr_data <- expr_long %>%
  mutate(sex_receptor = factor(paste(Predicted_sex, gene, sep = "_"), levels = c(rbind(paste("F", receptors_order, sep = "_"), paste("M", receptors_order, sep = "_")))))

# Plot violin plot with Sex_Receptor on x-axis

violin_plot_combined <- ggplot(combined_expr_data, aes(x = sex_receptor, y = percent_expression, fill = Predicted_sex)) +
  geom_violin(trim = TRUE) +  # Violin plot for Male and Female
  geom_point() +
 # geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +  # Add median markers
  facet_wrap(~ cell_subtype, scales = "free_y") +  # Facet by cell type
  labs(x = "Receptor (Sex)", y = "Percent Expression", title = "Percent Expression Across Receptors by Sex and Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )
# Save the updated plot
ggsave(
  "/home/inf-21-2023/thesis/scRNAseq_data_oelen/results/violin_plot_percent_expression_combined_sex_receptor.pdf",
  violin_plot_combined,
  width = 12,
  height = 10
)
