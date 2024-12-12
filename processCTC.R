# rmf 6.13.2023, last modified 12.12.2024

# script assumes the following directories exist:
# plots/

#################
### FUNCTIONS ###
#################

AddMetadataToObj <- function(obj, barcodes, mito_ids){
  obj$sample <- barcodes$sample
  obj$dataset <- barcodes$dataset
  obj$cell <- barcodes$cell
  
  print("Adding mitochondrial gene percentage...")
  # get mito ids that are found in the object (they are not all guaranteed to be)
  mito_ids_in_obj <- mito_ids[which(mito_ids %in% row.names(obj))]
  obj$percent_mito <- PercentageFeatureSet(obj, features = mito_ids_in_obj)
  
  return(obj)
}

############
### MAIN ###
############

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(cowplot)

# data files
dgeCountsFile <- "zUMIs_output/expression/CTCsmartseq3.dgecounts.rds"  # output of running zUMIs
geneNameFile <- "zUMIs_output/expression/CTCsmartseq3.gene_names.txt"  # output of running zUMIs
barcodeFile <- "reads_for_zUMIs.samples.txt"  # output of remultiplexing fastqs
umiCountsFile <- "zUMIs_output/stats/CTCsmartseq3.UMIcounts.txt"  # output of running zUMIs
mitoMapFile <- "ensembleIDtoGeneNameMap.txt"  # provided in github repository

# read in data
dge_counts <- readRDS(dgeCountsFile)
gene_names <- read.table(geneNameFile, header = T)
barcodes <- read.table(barcodeFile, header = T)
raw_umi_counts <- read.table(umiCountsFile, header = T)
mito_map <- read.table(mitoMapFile)
colnames(mito_map) <- c("ensembl_id","gene_name")

######################
### pre-processing ###
######################

# add info to the barcodes table
barcodes <- separate(barcodes, col = sample, into = c("cell","dataset"), 
                     sep = "_", extra = "merge", remove = FALSE)

# check that all cells are in the DGE object
check <- row.names(barcodes[which(!(barcodes$BC %in% colnames(dge_counts$umicount$inex$all))),]) 
length(check) == 0  # TRUE if all cells are in the object

# add some info to the umi counts table
names(raw_umi_counts)[names(raw_umi_counts) == 'SampleID'] <- 'BC'
raw_umi_counts <- merge(raw_umi_counts, barcodes, by = "BC")

# subset by type of cancer cell line
raw_umi_melanoma <- subset(raw_umi_counts, subset = dataset == "melanoma_cell_line")
raw_umi_prostate <- subset(raw_umi_counts, subset = dataset == "prostate_cell_line")

# plot the number of reads per cell from zUMIs output
umi_counts_exon <- raw_umi_counts[raw_umi_counts$type == "Exon",]
umi_counts_intron <- raw_umi_counts[raw_umi_counts$type == "Intron",]
umi_counts_inex <- raw_umi_counts[raw_umi_counts$type == "Intron+Exon",]

# exons
ggplot(data = umi_counts_exon, aes(x = BC, y = Count, fill = dataset)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Reads Overlapping Exons")
ggsave("plots/barplot_readsPerCell_exons.png")
ggsave("plots/barplot_readsPerCell_exons.pdf")

# inex
ggplot(data = umi_counts_inex, aes(x = BC, y = Count, fill = dataset)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(trans='log10') +
  ggtitle("Reads Overlapping Introns and Exons")
ggsave("plots/barplot_readsPerCell_inex.png")
ggsave("plots/barplot_readsPerCell_inex.pdf")

# prep for creating seurat object
barcodes_melanoma <- subset(barcodes, subset = dataset == "melanoma_cell_line")
barcodes_prostate <- subset(barcodes, subset = dataset == "prostate_cell_line")

umi_count <- as.matrix(dge_counts$umicount$inex$all)

umi_melanoma <- umi_count[, which(colnames(umi_count) %in% barcodes_melanoma$BC)]
umi_prostate <- umi_count[, which(colnames(umi_count) %in% barcodes_prostate$BC)]

# create objects
obj_melanoma <- CreateSeuratObject(counts = umi_melanoma, project = "melanoma")
obj_prostate <- CreateSeuratObject(counts = umi_prostate, project = "prostate")

# add metadata
obj_melanoma <- AddMetadataToObj(obj_melanoma, barcodes_melanoma, mito_map$ensembl_id)
obj_prostate <- AddMetadataToObj(obj_prostate, barcodes_prostate, mito_map$ensembl_id)

############################
### visualize QC metrics ###
############################

# legend code adapted from: https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd

obj_qc <- merge(obj_melanoma, y = obj_prostate)

features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
plot_titles <- c("Number of\nGenes Detected", "Number of\nMolecules Detected", "Percent Mitochondrial\nReads")

plist <- list()
i <- 1
for (feature in features_to_plot){
  p <- VlnPlot(obj_qc, features = features_to_plot[[i]], group.by = "dataset") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(fill=guide_legend(title="Tissue    ")) +
    scale_fill_discrete(labels=c('Melanoma    ', 'Prostate')) +
    theme(legend.box = "horizontal", legend.direction = "horizontal") +
    ggtitle(plot_titles[[i]]) +
    theme(plot.title = element_text(size = 14))
  plist[[i]] <- p
  i <- i + 1
}

# now go through and hide the legend for each plot
plist_nolegend <- list()
i <- 1
for (p in plist){
  p_nolegend <- p + theme(legend.position="none")
  plist_nolegend[[i]] <- p_nolegend
  i <- i+1
}
myplot <- plot_grid(plotlist = plist_nolegend, ncol = 3)

# extract legend from first plot list (horizontal)
legend <- get_legend(
  plist[[1]] + 
    guides(color = guide_legend(ncol = 2)) +
    theme(legend.box.margin = margin(0, 0, 18, 0)) #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
)

# add the legend underneath the row we made earlier
# give it 10% of the height of one plot (via rel_heights)
plot_grid(myplot, legend, nrow = 2, rel_heights = c(1, .1))
ggsave("plots/violin_QC_cellLines.png")
ggsave("plots/violin_QC_cellLines.pdf")

# plot melanoma only
plist <- list()
i <- 1
for (feature in features_to_plot){
  p <- VlnPlot(obj_melanoma, features = features_to_plot[[i]], pt.size = 2) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(fill=guide_legend(title="Tissue    ")) +
    ggtitle(plot_titles[[i]]) +
    theme(plot.title = element_text(size = 12))
  plist[[i]] <- p
  i <- i + 1
}

# now go through and hide the legend for each plot
plist_nolegend <- list()
i <- 1
for (p in plist){
  p_nolegend <- p + theme(legend.position="none")
  plist_nolegend[[i]] <- p_nolegend
  i <- i+1
}
myplot <- plot_grid(plotlist = plist_nolegend, ncol = 3)
myplot
ggsave("plots/violin_QC_melanoma.png")
ggsave("plots/violin_QC_melanoma.pdf")

#############################
### now plot individually ###
#############################

features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
plot_titles <- c("Number of Genes Detected", "Number of Molecules Detected", "Percent Mitochondrial Reads")

i <- 1
for (feature in features_to_plot){
  VlnPlot(obj_qc, features = features_to_plot[[i]], group.by = "dataset") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(fill=guide_legend(title="Tissue")) +
    scale_fill_discrete(labels=c('Melanoma', 'Prostate')) +
    ggtitle(plot_titles[[i]])
  ggsave(paste("plots/violin_QC_cellLines_",features_to_plot[[i]],".png", sep = ""))
  ggsave(paste("plots/violin_QC_cellLines_",features_to_plot[[i]],".pdf", sep = ""))
  i <- i + 1
}

#####################################
### print min/max feature content ###
#####################################

min(obj_melanoma$nCount_RNA)
min(obj_prostate$nCount_RNA)

min(obj_melanoma$nFeature_RNA)
min(obj_prostate$nFeature_RNA)

max(obj_melanoma$percent_mito)
max(obj_prostate$percent_mito)

##############################
### initial analysis merge ###
##############################

obj_melanoma <- NormalizeData(obj_melanoma)
obj_prostate <- NormalizeData(obj_prostate)

# melanoma only
obj_melanoma <- FindVariableFeatures(obj_melanoma)
obj_melanoma <- ScaleData(obj_melanoma)
obj_melanoma <- RunPCA(obj_melanoma, npcs = 10)
obj_melanoma <- RunUMAP(obj_melanoma, dims = 1:5, n.neighbors=9L)
obj_melanoma <- FindNeighbors(obj_melanoma, dims = 1:5)
obj_melanoma <- FindClusters(obj_melanoma, res = 0.5)
#saveRDS(obj_melanoma, file = "obj_analyzed_melanoma.RDS")

# merged
obj_merged <- merge(obj_melanoma, y = obj_prostate, merge.data = TRUE)
#saveRDS(obj_merged, file = "obj_merged.RDS")

obj_merged <- FindVariableFeatures(obj_merged)
obj_merged <- ScaleData(obj_merged)
obj_merged <- RunPCA(obj_merged, npcs = 30)
ElbowPlot(obj_merged)
obj_merged <- RunUMAP(obj_merged, dims = 1:5)
obj_merged <- FindNeighbors(obj_merged, dims = 1:5)
combined_obj <- FindClusters(obj_merged, res = 0.5)

#saveRDS(combined_obj, file = "obj_analyzed.RDS")
combined_obj <- readRDS("obj_analyzed.RDS")

# visualize PCA
DimPlot(combined_obj, reduction = "pca", group.by = "dataset") +
  scale_color_hue(labels=c('Melanoma', 'Prostate')) +
  ggtitle("")
ggsave("plots/dimplot_PCA_byDataset.png")
ggsave("plots/dimplot_PCA_byDataset.pdf")

# visualize UMAP
DimPlot(combined_obj, reduction = "umap", group.by = "dataset") +
  scale_color_hue(labels=c('Melanoma', 'Prostate')) +
  ggtitle("")
ggsave("plots/dimplot_UMAP_byDataset.png")
ggsave("plots/dimplot_UMAP_byDataset.pdf")

##############################
### visualize marker genes ###
##############################

# to find marker gene IDs:
# google the name, eg NGFR to get ENSEMBL ID without decimal
# use the following line to get the ID with decimal in this dataset:
# row.names(combined_obj@assays$RNA@counts)[grep("ENSG00000064300", row.names(combined_obj@assays$RNA@counts))]
# add to the markers list, add gene name to labels list

# MIA, MCSP (aka CSPG4), MCAM, MART-1 (aka MLANA), and PTPRC
# plus: PRAME, c-kit (not in the dataset), and PD-L1 (aka CD274-- not expressed)
labels <- c("MIA","MCSP","MCAM",
            "MART1","PTPRC","PRAME")
markers <- c("ENSG00000261857.7","ENSG00000173546.7","ENSG00000076706.17",
             "ENSG00000120215.10","ENSG00000081237.20", "ENSG00000185686.18")

# label source code: https://github.com/wilkelab/cowplot/blob/master/vignettes/shared_legends.Rmd
# violin plot, first with legend so we can extract later
plist <- list()
i <- 1
for (m in markers){
  p <- VlnPlot(obj_melanoma, assay = "RNA", features = markers[[i]], pt.size = 2) +
    ggtitle(labels[[i]]) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size = 12)) +
    theme(legend.box = "horizontal", legend.direction = "horizontal") 
  plist[[i]] <- p
  i <- i+1
}

# now go through and hide the legend for each plot
plist_nolegend <- list()
i <- 1
for (p in plist){
  p_nolegend <- p + theme(legend.position="none")
  plist_nolegend[[i]] <- p_nolegend
  i <- i+1
}
myplot <- plot_grid(plotlist = plist_nolegend)
myplot
ggsave("plots/violin_markergenes_melanoma.pdf")
ggsave("plots/violin_markergenes_melanoma.png")
