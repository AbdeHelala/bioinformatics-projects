# Name: Abdelsalam Helala
# ID: 7056985
# Name: Ahmed Lamloum
# ID: 7003029


#Function to make plots

save_and_display_plot <- function(plot, filename, width = 10, height = 5, plot_type = "ggplot") {
  if (plot_type == "ggplot") {
    # Save ggplot as PNG
    ggsave(filename, plot = plot, device = "png", width = width, height = height)
    print(plot)  # Print plot inline
  } else if (plot_type == "base") {
    # Save base R plot as PNG
    png(filename, width = width * 100, height = height * 100)
    plot()
    dev.off()
    print(paste("Plot saved as", filename))  
  }
}

setwd("~/Single-cell/second-project/single-cell/scbi_p2/project2")
set.seed(50)

library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ArchR)
library(GenomicRanges)
library(SummarizedExperiment)
library(Cairo)
library(hexbin)
library(pheatmap)

ArchR::installExtraPackages()

# Setup ArchR environment
addArchRThreads(threads = 1)
addArchRGenome("hg38")       

# Define fragment files and output directory
fragment_files <- c(
  "hft_ctx_w21_dc1r3_r1_atac_fragments.tsv.gz",
  "hft_ctx_w21_dc2r2_r1_atac_fragments.tsv.gz",
  "hft_ctx_w21_dc2r2_r2_atac_fragments.tsv.gz"
)
sample_names <- c("dc1r3_r1", "dc2r2_r1", "dc2r2_r2")

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = fragment_files,
  sampleNames = sample_names,
  minFrags = 500,    
  minTSS = 4,       
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 1000),  
  addGeneScoreMat = TRUE  
)

# Identify doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

# Create ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE
)

# Retrieve doublet scores and plot distribution
doublet_scores <- as.data.frame(getCellColData(proj, select = "DoubletScore"))
doublet_hist_plot <- ggplot(doublet_scores, aes(x = DoubletScore)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  ggtitle("Doublet Score Distribution") +
  xlab("Doublet Score") +
  theme_minimal()
save_and_display_plot(doublet_hist_plot, "Doublet_Score_Distribution.png")

# Filter out doublets
proj <- filterDoublets(proj)


# Inspect cell metadata
num_cells <- nCells(proj)
median_tss <- median(getCellColData(proj, "TSSEnrichment")$TSSEnrichment)
median_fragments <- median(getCellColData(proj, "nFrags")$nFrags)

# Print results
cat("Total number of cells in the project:", num_cells, "\n")
cat("Median TSS Enrichment value:", median_tss, "\n")
cat("Median number of fragments:", median_fragments, "\n")

# Dimension of the peak set
if (!is.null(proj$peakSet)) {
  peak_set_dim <- dim(proj$peakSet)
  cat("Dimension of the peak set:", paste(peak_set_dim, collapse = " x "), "\n")
} else {
  cat("Peak set has not been generated yet.\n")
}


# QC Metrics and Summary
num_cells <- nCells(proj)  # Total number of cells
median_tss <- median(getCellColData(proj, "TSSEnrichment")$TSSEnrichment)
median_fragments <- median(getCellColData(proj, "nFrags")$nFrags)

# Print results
cat("Total number of cells in the project:", num_cells, "\n")
cat("Median TSS Enrichment value:", median_tss, "\n")
cat("Median number of fragments:", median_fragments, "\n")


# 1. Fragment Length Distribution
fragment_length_plot <- plotFragmentSizes(proj)
save_and_display_plot(fragment_length_plot, "Fragment_Length_Distribution.png")

# 2. TSS Enrichment Distribution
tss_enrichment_plot <- plotTSSEnrichment(proj)
save_and_display_plot(tss_enrichment_plot, "TSS_Enrichment_Distribution.png")

# 3. Fragments vs TSS Enrichment Scatter Plot
scatter_plot <- ggplot(getCellColData(proj), aes(x = nFrags, y = TSSEnrichment, color = Sample)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +  # Log scale for better visualization
  ggtitle("Fragments vs. TSS Enrichment") +
  theme_minimal()

save_and_display_plot(scatter_plot, "Fragments_vs_TSS_Enrichment.png")


# QC Metrics and Plotting Functions
plot_QC_metrics <- function(proj) {
  # Median TSS Enrichment and Fragment count
  numCells <- nCells(proj)
  medianTSS <- median(getCellColData(proj, "TSSEnrichment")$TSSEnrichment)
  medianFragments <- median(getCellColData(proj, "nFrags")$nFrags)
  
  # Print summary stats
  print(paste("Total number of cells:", numCells))
  print(paste("Median TSS Enrichment:", medianTSS))
  print(paste("Median number of fragments:", medianFragments))
  
  # Generate QC plots
  p1 <- plotFragmentSizes(proj)
  p2 <- plotTSSEnrichment(proj)
  combined <- p1 + p2
  
  # Save and display combined plot
  save_and_display_plot(combined, "Combined_Fragment_TSS_Plots.png")
}

# Run QC Plotting
plot_QC_metrics(proj)


# Apply stricter filtering criteria
filtered_cells <- rownames(proj@cellColData)[
  proj@cellColData$nFrags > 3200 &                
    proj@cellColData$nFrags < 100000 &              
    proj@cellColData$TSSEnrichment > 10 &           
    proj@cellColData$DoubletScore < 50              
]

# Subset the project based on stricter filters
proj_filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = filtered_cells,
  outputDirectory = "FilteredProject_Strict",
  force = TRUE
)

# Save the filtered project
saveArchRProject(ArchRProj = proj_filtered, overwrite = TRUE)

# Recalculate QC metrics
num_cells_filtered <- nCells(proj_filtered)  # Total number of filtered cells
median_tss_filtered <- median(getCellColData(proj_filtered, "TSSEnrichment")$TSSEnrichment)
median_fragments_filtered <- median(getCellColData(proj_filtered, "nFrags")$nFrags)


cat("After filtering:\n")
cat("Total number of cells:", num_cells_filtered, "\n")
cat("Median TSS Enrichment:", median_tss_filtered, "\n")
cat("Median number of fragments:", median_fragments_filtered, "\n")



#-------------------------------------------------------------------





# Perform Iterative LSI for dimensionality reduction
proj_filtered <- addIterativeLSI(
  ArchRProj = proj_filtered,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = 0.2,       
    sampleCells = 10000,    
    n.start = 10),          
  varFeatures = 25000,      
  dimsToUse = 1:30          
)


# Compute UMAP coordinates
proj_filtered <- addUMAP(
  ArchRProj = proj_filtered,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

# Plot UMAP colored by Sample
umap_sample_plot <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)
save_and_display_plot(umap_sample_plot, "UMAP_Sample_Plot.png")

# Plot UMAP colored by TSS Enrichment
umap_tss_plot <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "TSSEnrichment",
  embedding = "UMAP"
)
save_and_display_plot(umap_tss_plot, "UMAP_TSS_Enrichment.png")

# Plot UMAP colored by Number of Fragments
umap_fragments_plot <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "nFrags",
  embedding = "UMAP"
)
save_and_display_plot(umap_fragments_plot, "UMAP_Fragments.png")




# Add Harmony for batch correction
proj_filtered <- addHarmony(
  ArchRProj = proj_filtered,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

# Recalculate UMAP with batch correction
proj_filtered <- addUMAP(
  ArchRProj = proj_filtered,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

# Plot batch-corrected UMAP
umap_harmony_plot <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)
save_and_display_plot(umap_harmony_plot, "UMAP_Harmony_BatchCorrected.png")


# Apply Louvain clustering
proj_filtered <- addClusters(
  input = proj_filtered,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5
)

# Plot UMAP colored by Clusters
umap_cluster_plot <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)
save_and_display_plot(umap_cluster_plot, "UMAP_Clusters.png")

# Count the number of cells per cluster
cluster_counts <- table(getCellColData(proj_filtered, "Clusters"))
print(cluster_counts)


# Add group coverages
proj_filtered <- addGroupCoverages(
  ArchRProj = proj_filtered,
  groupBy = "Clusters"
)

# Perform peak calling using MACS2
proj_filtered <- addReproduciblePeakSet(
  ArchRProj = proj_filtered,
  groupBy = "Clusters",
  pathToMacs2 = "/home/abkh00004/miniconda3/envs/single-cell/bin/macs2"
)

# Add the peak matrix to the project
proj_filtered <- addPeakMatrix(proj_filtered)

# Identify marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Generate heatmap for marker peaks
markerHeatmap <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

# Save the heatmap
png("Marker_Peak_Heatmap.png", width = 1000, height = 800)
draw(markerHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# Define genes of interest
genes_of_interest <- c("TOP2A", "MKI67", "AURKA", "SATB2", "SLC12A7")

# Generate browser tracks for the selected genes
plotGenes <- plotBrowserTrack(
  ArchRProj = proj_filtered,
  groupBy = "Clusters",
  geneSymbol = genes_of_interest,
  upstream = 5000,
  downstream = 5000
)

# Save each gene's plot individually
for (i in seq_along(plotGenes)) {
  # Extract individual plot
  p <- plotGenes[[i]]
  
  # Define the gene name and file name
  gene_name <- genes_of_interest[i]
  filename <- paste0("Gene_Browser_Track_", gene_name, ".png")
  
  # Save the plot using ggsave
  ggsave(filename, plot = p, width = 8, height = 6)
  cat("Saved:", filename, "\n")
}


#-------------------------------------------------------------

# Compute Gene Activity Scores
proj_filtered <- GeneScoreMatrix(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",  
  matrixName = "GeneActivity"
)


# Identify marker genes
markerGenes <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


# Plot a heatmap of top marker genes
heatmap_genes <- plotMarkerHeatmap(
  seMarker = markerGenes,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

png("Marker_Gene_Heatmap.png", width = 1000, height = 800)
ComplexHeatmap::draw(heatmap_genes, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

#----------------------------------------------------


library(Rmagic)

topGenes <- rowData(markerGenes)$name[1:5]  # Extract first 5 marker genes

proj_filtered <- addImputeWeights(proj_filtered)

for (gene in topGenes) {
  
  # UMAP plot without MAGIC smoothing
  plot_without_magic <- plotEmbedding(
    ArchRProj = proj_filtered,
    colorBy = "GeneScoreMatrix",
    name = gene,
    embedding = "UMAP",
    imputeWeights = NULL  
  )
  
  # UMAP plot with MAGIC smoothing
  plot_with_magic <- plotEmbedding(
    ArchRProj = proj_filtered,
    colorBy = "GeneScoreMatrix",
    name = gene,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_filtered)  
  )
  
  # Save the plots
  save_and_display_plot(plot_without_magic, paste0(gene, "_UMAP_Without_MAGIC.png"), plot_type = "ggplot")
  save_and_display_plot(plot_with_magic, paste0(gene, "_UMAP_With_MAGIC.png"), plot_type = "ggplot")
}

#-------------------------------------------------------


# Add Motif Annotations
proj_filtered <- addMotifAnnotations(
  ArchRProj = proj_filtered,
  motifSet = "cisbp"  # Using the CIS-BP motif database
)

# Compute Motif Enrichment and Activity
proj_filtered <- addDeviationsMatrix(
  ArchRProj = proj_filtered,
  peakAnnotation = "Motif",   # Annotation of motifs in peak regions
  matrixName = "MotifMatrix"  # Name of the deviations matrix
)

# Extract the "name" and "combinedVars" from motifVar$data
motifVar_df <- motifVar$data  # Access the data frame directly


# Extract the top 2 most variable motifs
if (!is.null(motifVar_df$name) & !is.null(motifVar_df$combinedVars)) {
  topTFs <- motifVar_df$name[order(motifVar_df$combinedVars, decreasing = TRUE)[1:2]]
} 

print(topTFs)

# Plot UMAP embeddings for the top TF motifs
for (tf in topTFs) {
  plotTF <- plotEmbedding(
    ArchRProj = proj_filtered,
    colorBy = "MotifMatrix",
    name = tf,
    embedding = "UMAP"
  )
  # Save the plot
  save_and_display_plot(plotTF, paste0(tf, "_UMAP.png"), plot_type = "ggplot")
}

# Plot motif activity distribution across clusters
plotMotif <- plotGroups(
  ArchRProj = proj_filtered,
  groupBy = "Clusters",      # Group by clusters
  colorBy = "MotifMatrix",   # Use MotifMatrix for plotting
  name = topTFs             
)

# Save the motif activity plot
save_and_display_plot(plotMotif, "Motif_Activity_Distribution.png", plot_type = "ggplot")

getFeatures(proj_filtered, useMatrix = "MotifMatrix")

# Get the full feature names for top TFs
features <- getFeatures(proj_filtered, useMatrix = "MotifMatrix")
topTFs_full <- features[grep(paste(topTFs, collapse = "|"), features)]
print(topTFs_full)

# Filter to include only deviations
topTFs_full <- topTFs_full[grep("^deviations:", topTFs_full)]
print(topTFs_full)  # Confirm the filtered names

for (tf in topTFs_full) {
  plotTF <- plotEmbedding(
    ArchRProj = proj_filtered,
    colorBy = "MotifMatrix",
    name = tf,
    embedding = "UMAP"
  )
  # Save the plot
  save_and_display_plot(plotTF, paste0(tf, "_UMAP.png"), plot_type = "ggplot")
}


plotMotif <- plotGroups(
  ArchRProj = proj_filtered,
  groupBy = "Clusters",      # Group by clusters
  colorBy = "MotifMatrix",   # Use MotifMatrix for plotting
  name = topTFs_full         
)

# Save each plot in the list
if (is.list(plotMotif)) {
  for (i in seq_along(plotMotif)) {
    # Generate a filename for each plot
    filename <- paste0("Motif_Activity_Distribution_", i, ".png")
    
    # Save each plot
    save_and_display_plot(plotMotif[[i]], filename, plot_type = "ggplot")
  }
} else {
  save_and_display_plot(plotMotif, "Motif_Activity_Distribution.png", plot_type = "ggplot")
}

