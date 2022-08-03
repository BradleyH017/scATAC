##### Bradley June 2022
##### scATAC analysis using archR https://www.archrproject.com/ https://www.archrproject.com/articles/Articles/tutorial.html
##### Requires fragments as input only
##### OLD conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env r-base r-seurat signac r-matrix bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase
##### OLD conda install -c conda-forge r-devtools
##### conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_macs2_env r-base r-matrix bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase macs2 r-devtools r-harmony bioconductor-qusage bioconductor-fgsea
##### R
##### Sys.setenv(CONDA_BUILD_SYSROOT="/")
##### devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env


# Set up
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/data_archR")
library(ArchR)
library(parallel)
set.seed=1
addArchRThreads(threads = 25)
samples <- readLines("samples_use.txt")
inputFiles <- paste0("filtered_", samples)


# Add genome version
addArchRGenome("hg38")

# Reformat the fragment files, overwrite the variable
reformatFragmentFiles(inputFiles)
inputFiles <- paste0("filtered_", samples, ".gz")

# Load data
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = samples,
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## Find doublets -- NOT RUN as using only the 'filtered' cells from the authors
#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#  LSIMethod = 1
#)

# Now make the project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "../results_archR/",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

# PLotting some QC metrics
p1 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

p3 <- plotGroups(
    ArchRProj = proj,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
  )

p4 <- plotFragmentSizes(ArchRProj = proj)
p4

p5 <- plotTSSEnrichment(ArchRProj = proj)
p5

# Save together
plotPDF(p1,p3,p4,p5, name = "QC-Sample-All-Profile.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)


## Filter doublets - Something wrong with 'Cairo' installation
#proj <- filterDoublets(ArchRProj = proj)

# Dimensionallity reduction and clustering
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# Batch correction using HArmony. Requires installation into conda env
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Add clusters - Area that could be tweaked
proj <- addClusters(input = proj, reducedDims = "Harmony", force=T)

# Visualise
proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", force=T)

# Plot with Cluster
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
#ggAlignPlots(p1, p2, type = "h")

# Save Dimplot for Thesis
png(file="../../Fig_plots/Thesis/Chapter2/scATAC_umap_clusters.png", width=7.5*ppi, heigh=7.5*ppi, res=ppi)
plotEmbedding(ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP",
    size=0.4,
    plotAs="points",
    keepAxis=T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) + ggtitle("")
dev.off()

# Assigning clusters with GeneScores
proj <- addImputeWeights(proj)

# Plot with marker genes - This is the 3 cis-eQTLs, + 12 genes which are correlated with one another in cluster 11 of the scRNAseq specifically
markerGenes  <- c(
    "C11orf53", "COLCA1", "COLCA2",
    "LRMP", "SH2D6", "HTR3E", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PTGS1", "IL17RB", "AZGP1", "GNG13"
  )


pmarks <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

#Rearrange for grid plotting
pmarks2 <- lapply(pmarks, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 2),pmarks2))


# Doing this myself for thesis
pmarks <- vector("list", length = length(markerGenes))
for(g in seq_along(markerGenes)){
  temp <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "GeneScoreMatrix",
      name = markerGenes[g],
      embedding = "UMAP",
      size=0.1,
      imputeWeights = getImputeWeights(proj),
      plotAs="points",
      keepAxis=T) +
      ggtitle(markerGenes[g]) +
      xlab("UMAP_1") +
      ylab("UMAP_2") +
      theme(plot.title=element_text(face="bold", size = 14)) +
      theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10))
  ppi=300
  png(file=paste0("../../Fig_plots/Thesis/Chapter2/", markerGenes[g], "_scATAC_umap.png"), height=4*ppi, width=4*ppi, res=ppi)
  print(temp)
  dev.off()
  print("plotted")
}



trans_genes <- c("POU2F3","ITPRID1","BMX","SH2D6","SH2D7","CHAT","HTR3E","TRPM5","AZGP1","OGDHL","AVIL","KLK13","PSTPIP2","PIK3CG","PLCG2","TAS1R3")
p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = trans_genes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 8),p2))


# Visualising marker genome tracks
all_int <- unique(c(markerGenes, trans_genes))
p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Clusters",
    geneSymbol = all_int,
    upstream = 50000,
    downstream = 50000
)

plotPDF(plotList = p,
    name = "Plot-Tracks-All-Int-Genes.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)

grid::grid.newpage()
grid::grid.draw(p$C11orf53)

# Saving the ArchRProject
setwd("../results_archR")
proj <- saveArchRProject(ArchRProj = proj)


# Plotting browser track for


# Getting markers
markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
c5 <- markerList$C5
for(c in seq_along(markerList)){
  write.csv(as.data.frame(markerList[[c]]), paste0("Markers/cluster_", c, ".csv"))
}

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
c5 <- markerList$C5
for(c in seq_along(markerList)){
  write.csv(as.data.frame(markerList[[c]]), paste0("Markers/relax_cluster_", c, ".csv"))
}


# Enrichment of the markers of each clusters against the Smillie cluster markers.
# Requires installation of qusage and fgsea
library(fgsea)
library(qusage)
# Load in the smillie sets
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")
gene_list_dir <- paste(pathOut, "tables/markers/Authors", sep ="/")
sets <- list.files(gene_list_dir)
sets <- paste(gene_list_dir, sets, sep = "/")
smil_sets <- lapply(sets, read.gmt)
# reformat
for(x in seq(1:length(smil_sets))){
  names(smil_sets)[x] <- names(smil_sets[[x]]);
  smil_sets[[x]] <- smil_sets[[x]][[1]]
}

# Add trans-eQTL gene list and cluster 11 markers
explore_sets <- list()
explore_sets$rs3087967_trans_eQTLs <- c("POU2F3","ITPRID1","BMX","SH2D6","SH2D7","CHAT","HTR3E","TRPM5","AZGP1","OGDHL","ACTG1P22","AVIL","KLK13","IRAG2","PSTPIP2","PIK3CG","PLCG2","TAS1R3")
c11 <- read.csv("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/tables/markers_unmerged_clusters/Cluster11_markers.csv")
explore_sets <- c(explore_sets, smil_sets)

# Reformat our own cluster markers. Make geneList
clusters <- paste0("cluster_", seq(1:5))
c_list <- vector("list", length=length(clusters))
c_gl <- vector("list", length=length(clusters))
for(c in seq_along(clusters)){
  c_list[[c]] <- as.data.frame(markerList[[c]])
  c_gl[[c]] <- c_list[[c]]$Log2FC
  names(c_gl[[c]]) <- c_list[[c]]$name
}


# Now do the enrichment
fgsea_Res <- vector("list", length = length(clusters))
names(fgsea_Res) <- clusters
for(q in seq(1:length(fgsea_Res))){
  print(paste("testing cluster ", names(fgsea_Res)[q], sep = ""));
  fgsea_Res[[q]] <-  fgseaMultilevel(explore_sets, c_gl[[q]], minSize = 0, maxSize = 500, scoreType = "pos", eps = 0)
  fgsea_Res[[q]]$cluster <- rep(names(fgsea_Res)[q], nrow(fgsea_Res[[q]]));
  fgsea_Res[[q]] <- fgsea_Res[[q]][order(fgsea_Res[[q]]$padj),];
  fgsea_Res[[q]]$BH_p.adj <- p.adjust(fgsea_Res[[q]]$pval, method = "BH")
}

# Save this
fgsea_Res_df <- do.call(rbind, fgsea_Res)
write.csv(as.data.frame(fgsea_Res_df[,-8]), "Markers/fgsea_all_markers.csv")

# PLot This
trans_anno <- paste0("NES=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "rs3087967_trans_eQTLs",]$NES, 2), "\n", "p=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "rs3087967_trans_eQTLs",]$pval, 2))
Trans <- plotEnrichment(explore_sets[['rs3087967_trans_eQTLs']], c_gl[[5]]) + theme_bw() + ggtitle("Trans-eQTL targets in Cluster 5") + xlab("Rank") + ylab("Enrichment score") + theme(plot.title=element_text(face="bold", size = 14)) + annotate(geom="text", x=220, y=0.7, label=trans_anno, color="black", size = 6, hjust=1)

c11_anno <- paste0("NES=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "cluster11_scRNA",]$NES, 2), "\n", "p=",signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "cluster11_scRNA",]$pval, 2))
Cluster11 <- plotEnrichment(explore_sets[['cluster11_scRNA']], c_gl[[5]]) + theme_bw() + ggtitle("Cluster11_scRNAseq markers in Cluster 5") + xlab("Rank") + ylab("Enrichment score") + theme(plot.title=element_text(face="bold", size = 14)) + annotate(geom="text", x=220, y=0.6, label=c11_anno, color="black", size = 6, hjust=1)

# Put together and save
library(patchwork)
ppi=300
png(file="Plots/Enrichment_GS_cluster5.png", width=12*ppi, height=5*ppi, res=ppi)
Trans + Cluster11
dev.off()

# Draw heatmap
labels <- intersect(explore_sets[['cluster11_scRNA']], c5$name)
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  transpose = TRUE,
  labelMarkers=labels
)
png(file="Plots/Heatmap_cluster11_in_cluster5.png", width=12*ppi, height=5*ppi, res=ppi)
draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

##### Co-accessibility analysis
proj <- loadArchRProject(path = "./")

# Add group coverage
proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Clusters",
    useLabels = TRUE
)

# Add peak matrix
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy="Clusters",
    peakMethod = "Macs2"
#    pathToMacs2="/home/s1987229/macs2"
)

proj <- addPeakMatrix(
  ArchRProj = proj
)


proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "Harmony"
)

# Obtain
cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = T
)

# Calculate marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Plot markers
markerGenes  <- c(
    "MUC2",  #Goblet
    "CLCA1", #Goblet
    "TRPM5", "POU2F3", "HTR3E", "SH2D6", "SH2D7", "IL17RB", # Tuft
    "C11orf53", "COLCA1", "COLCA2"
  )

# Calc markers
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1", returnGR=T)
markerList

# Calc heatmap
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1",
  transpose = TRUE
)

# Draw
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# PLot example peak
p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Clusters",
    geneSymbol = c("FCGBP"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1", returnGR = TRUE)["C1"],
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$FCGBP)


# Calculate DA peaks between C2 and C5
markerTest <- getMarkerFeatures(ArchRProj=proj,
  useMatrix="PeakMatrix",
  testMethod = "wilcoxon",
  groupBy="Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups="C5",
  bgdGroups="C2")

res <- as.data.frame(getMarkers(markerTest, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR=T))
resup <- res[res$Log2FC > 0,]

pma <- markerPlot(seMarker = markerTest, name = "C5", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma



p <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = "Clusters",
      geneSymbol = markerGenes,
      upstream = 50000,
      downstream = 50000,
      loops = getCoAccessibility(proj)
  )

grid::grid.newpage()
grid::grid.draw(p$CLCA1)

# save
saveArchRProject(ArchRProj = proj)


#### Motif enrichment - non Custom. Using JASPAR as this is where the C11orf53/POU2F3 enrichment came from
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "Motif", force=T)
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1"
  )
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

# plot
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

##### Custom enrichment. https://www.archrproject.com/bookdown/custom-enrichment.html
# ChIP used hg38 genome
# As far as I can tell, the ATACseq also used hg38 (although they describe GRCh38, which I think is the same)

C11orf53_peaks <- c(
  OCAT_NCIH211_1 <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657742/suppl/GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_peaks.narrowPeak.gz",
  OCAT_NCIH211_2 <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657743/suppl/GSM5657743_NCIH211_OCA-T_Ab2_S7_R1_001.fastq.gz_peaks.narrowPeak.gz",
  OCAT_NCIH526_1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657751/suppl/GSM5657751_NCIH526_OCA-T_Ab1_S6_R1_001.fastq.gz_bt2_snstv.sorted_trimmed_rd.bam_peaks.narrowPeak.gz",
  OCAT_NCIH526_2 <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657752/suppl/GSM5657752_NCIH526_OCA-T_Ab2_S5_R1_001.fastq.gz_bt2_snstv.sorted_trimmed_rd.bam_peaks.narrowPeak.gz"
)
proj <- addPeakAnnotations(ArchRProj = proj, regions = C11orf53_peaks, name = "C11orf53_ChIP", force=T)

# may need to recalculate markers if loaded
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Calculate the enrichment
enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "C11orf53_ChIP"
  )

#rename
rownames(enrichRegions)=enrichRegions@NAMES

# Have a peak
enrichRegions


heatmapRegions <- plotEnrichHeatmap(enrichRegions, cutOff=0, transpose = TRUE)
ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", annotation_legend_side = "bot")


# Extract results
df <- data.frame(TF = rownames(enrichRegions),
                mlog10Padj = assay(enrichRegions)[,1],
                mlog10p=assay(enrichRegions)[,2],
                Enrichment=assay(enrichRegions)[,3],
                BackgroundProporition=assay(enrichRegions)[,4],
                nBackground=assay(enrichRegions)[,5]
                )
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))



#### TEsting motif deviations. This calculates the motif variability
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

# Extracting a subset of motifs for downstream analysis
motifs <- c("POU2F3", "OCT4", "OCT11", "OCT2")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

# Plot as a ridge plot
p <- plotGroups(ArchRProj = proj,
  groupBy = "Clusters",
  colorBy = "MotifMatrix",
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj)
)
p

# PLot on UMAP
p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)
p2

p+p2

#### Motif deviations using the custom input - THIS IS GREAT
C11orf53_peaks <- c(
  "OCAT_NCIH211_1" <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657742/suppl/GSM5657742_NCIH211_OCA-T_Ab1_S8_R1_001.fastq.gz_peaks.narrowPeak.gz",
  "OCAT_NCIH211_2" <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657743/suppl/GSM5657743_NCIH211_OCA-T_Ab2_S7_R1_001.fastq.gz_peaks.narrowPeak.gz",
  "OCAT_NCIH526_1" <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657751/suppl/GSM5657751_NCIH526_OCA-T_Ab1_S6_R1_001.fastq.gz_bt2_snstv.sorted_trimmed_rd.bam_peaks.narrowPeak.gz",
  "OCAT_NCIH526_2" <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5657nnn/GSM5657752/suppl/GSM5657752_NCIH526_OCA-T_Ab2_S5_R1_001.fastq.gz_bt2_snstv.sorted_trimmed_rd.bam_peaks.narrowPeak.gz"
)

proj <- addPeakAnnotations(ArchRProj = proj, regions = C11orf53_peaks, name = "ChIP_C11orf53", force=T)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "ChIP_C11orf53",
  matrixName="ChIP_C11orf53_Matrix",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, plot = TRUE, name = "ChIP_C11orf53_Matrix")

# Plot deviation z-scores over UMAP embedding
markerChIP <- getFeatures(proj, useMatrix = "ChIP_C11orf53_Matrix")
markerChIP <- sort(grep("z:", markerChIP, value = TRUE))
markerChIP
p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "ChIP_C11orf53_Matrix",
    name = markerChIP,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)
p

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 2), p2))

ppi=300
png(file="Plots/C11orf53_motif_accessibility_UMAPs_from_ChIP.png", width=8*ppi, height=8*ppi, res=ppi)
do.call(cowplot::plot_grid, c(list(ncol = 2),p2))
dev.off()

# make this slightly nicer for Thesis
ac_plots <- vector("list", length=4)
names(ac_plots) <- c("NCIH211.1", "NCIH211.2", "NCIH526.1", "NCIH526.2")
for(rep in seq_along(ac_plots)){
  temp <- plotEmbedding(
      ArchRProj = proj,
      colorBy = "ChIP_C11orf53_Matrix",
      name = markerChIP[rep],
      embedding = "UMAP",
      keepAxis=T,
      plotAs="points",
      size=0.1,
      imputeWeights = getImputeWeights(proj)
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(plot.title=element_text(face="bold", size = 14), legend.position="none") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) +
  ggtitle(names(ac_plots)[rep])
  ac_plots[[rep]] <- temp
}

ppi=300
png(file="../../Fig_plots/Thesis/Chapter2/C11orf53_ChIP_accessibility_umaps.png", width=16*ppi, height=5*ppi, res=ppi)
ac_plots[[1]] + ac_plots[[2]] + ac_plots[[3]] + ac_plots[[4]] + plot_layout(ncol=4)
dev.off()

### Motif footprinting
motifPositions <- getPositions(proj)
motifs <- c("POU2F3", "OCT4", "OCT11", "OCT2")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

# Calculate
seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "Clusters"
)

# Plot
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)


# save
proj <- saveArchRProject(ArchRProj = proj)

###### Identify Correlated TF Motifs and TF Gene Score/Expression
# Want to do this with the features in the ChIP_C11orf53_Matrix? This seems to be causing an error

# Using published motifs
# Calculate variability
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Clusters")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
# Calculate the correlations
corGSM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
)
# Add max Delta
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
# PLot
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p

# Doing this manually
# Want to look at the correlation between C11orf53 motif accessibility and expression of 11q23.1 targets
# Must be within the object somewhere
GS <- getMatrixFromProject(ArchRProj=proj,
  useMatrix="GeneScoreMatrix")
GS_df <- assays(GS)$GeneScoreMatrix
rownames(GS_df) <- rowData(GS)$name
GS_df <- GS_df[rownames(GS_df) %in% trans_genes,]

# extract ChIP score
ChIP_C11orf53 <- getMatrixFromProject(ArchRProj=proj,
  useMatrix= "ChIP_C11orf53_Matrix")
C53_df <- assays(ChIP_C11orf53)$z

# Test using second scores
# Is there a correlation between accessibility of C11orf53 and 11q23.1 genes?
clusters <- levels(factor(proj@cellColData$Clusters))
cor_clusters <- vector("list", length = length(clusters))
p_clusters <- vector("list", length=length(clusters))
for(c in seq_along(clusters)){
  cells.use <- rownames(proj@cellColData[proj@cellColData$Clusters == clusters[c],])
  temp <- GS_df[,colnames(GS_df) %in% cells.use]
  temp <- t(temp)
  temp_c53 <- C53_df[,colnames(C53_df) %in% cells.use]
  temp_c53 <- t(temp_c53)
  cor_clusters[[c]] <- matrix(nrow=ncol(temp), ncol=ncol(temp_c53))
  rownames(cor_clusters[[c]]) <- colnames(temp)
  colnames(cor_clusters[[c]]) <- colnames(temp_c53)
  p_clusters[[c]] <- matrix(nrow=ncol(temp), ncol=ncol(temp_c53))
  for(r in 1:ncol(temp)){
    for( s in 1:ncol(temp_c53)){
      cor_clusters[[c]][r,s] <- cor(temp[,r], temp_c53[,s])
      p_clusters[[c]][r,s] <- cor.test(temp[,r], temp_c53[,s], method="pearson")$p.value
      if(p_clusters[[c]][r,s] > 0.05){
        cor_clusters[[c]][r,s] <- 0
      }
    }
  write.csv(cor_clusters[[c]], paste0("tables/", "sig_cor_trans_C53_ChIP_clusters_", c, ".csv"))
  }
}

# Plot
ppi=300
png(file="Plots/correlation_plots_C53_ChIP_clusters.png", width=10*ppi, height=6*ppi, res=ppi)
par(mfrow=c(3,2))
for(cluster in seq_along(cor_clusters)){
  if(nrow(cor_clusters[[cluster]])>ncol(cor_clusters[[cluster]])){
    cor_clusters[[cluster]] <- t(cor_clusters[[cluster]])
  }
  rownames(cor_clusters[[cluster]]) <- c("OCAT_NCIH211_1", "OCAT_NCIH211_2", "OCAT_NCIH526_1", "OCAT_NCIH526_2")
  corrplot::corrplot(cor_clusters[[cluster]], cl.pos = 'n', is.corr = T , na.label = "NA", na.label.col = "white", mar=c(0,0,0,0), tl.cex=1.5, tl.col="black", bg=NA, tilte="", col=colorRampPalette(c("blue","white","red"))(200))
}
dev.off()

##################### BELOW IS SCRAP



# Exploratry extraction of matrices
GS <- getMatrixFromProject(ArchRProj=proj,
  useMatrix="GeneScoreMatrix")

# Exploratory loading of pseudo-bulk calls for samples per cell
file="/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/results_archR/GroupCoverages/Clusters/C5._.colon_transverse_SM.CSSDA_rep1_fragments.bed.insertions.coverage.h5"
h5ls(file)
mydata <- h5read(file, "/Coverage/chr11")





##### Motif enrichment analysis in trans-eQTL targets that are enriched in tuft-like cells (using homer to match ChIPSeq)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "homer", name = "Motif")
motifsUp <- peakAnnoEnrichment(
    seMarker = markersGS,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
