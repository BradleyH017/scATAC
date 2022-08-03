###### Bradley July 2022
###### scATAC ssGSEA for thesis and plots
###### Temp using different environmnent with the correct packages installed
###### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scRNA_mouse_env

###### Set up
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/results_archR/")
library(ArchR)
library(parallel)
set.seed=1
addArchRThreads(threads = 25)

# Load project
proj <- loadArchRProject(path = "./")

# Load the markers
markerList <- vector("list", length=5)
for(c in seq_along(markerList)){
  markerList[[c]] <- read.csv(paste0("Markers/relax_cluster_", c, ".csv"), row.names=1)
}

# Load the target signatures we want to test against
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
explore_sets$rs3087967_trans_eQTLs_RNAseq <- c("POU2F3","ITPRID1","BMX","SH2D6","SH2D7","CHAT","HTR3E","TRPM5","AZGP1","OGDHL","ACTG1P22","AVIL","KLK13","IRAG2","PSTPIP2","PIK3CG","PLCG2","TAS1R3")
explore_sets$FDR_sig_Vaughan_Shaw_HT12 <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "EIF4EBP", "GIN1", "SPAG6")
c11 <- read.csv("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/tables/markers_unmerged_clusters/Cluster11_markers.csv")
explore_sets <- c(explore_sets, smil_sets)
explore_sets$cluster11_scRNA <- c11$gene

# Reformat our own cluster markers. Make geneList
clusters <- paste0("cluster_", seq(1:5))
c_list <- vector("list", length=length(clusters))
c_gl <- vector("list", length=length(clusters))
for(c in seq_along(clusters)){
  c_list[[c]] <- as.data.frame(markerList[[c]])
  c_gl[[c]] <- c_list[[c]]$Log2FC
  names(c_gl[[c]]) <- c_list[[c]]$name
}

# Marker Enrichment
fgsea_Res <- vector("list", length = length(clusters))
names(fgsea_Res) <- clusters
for(q in seq(1:length(fgsea_Res))){
  print(paste("testing cluster ", names(fgsea_Res)[q], sep = ""));
  fgsea_Res[[q]] <-  fgseaMultilevel(explore_sets, c_gl[[q]], minSize = 0, maxSize = 500, scoreType = "pos", eps = 0)
  fgsea_Res[[q]]$cluster <- rep(names(fgsea_Res)[q], nrow(fgsea_Res[[q]]));
  fgsea_Res[[q]] <- fgsea_Res[[q]][order(fgsea_Res[[q]]$padj),];
  fgsea_Res[[q]]$BH_p.adj <- p.adjust(fgsea_Res[[q]]$pval, method = "BH")
}

fgsea_Res_df <- do.call(rbind, fgsea_Res)
trans_anno <- paste0("NES=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "FDR_sig_Vaughan_Shaw_HT12",]$NES, 2), "\n", "p=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "FDR_sig_Vaughan_Shaw_HT12",]$pval, 2))
Trans <- plotEnrichment(explore_sets[['FDR_sig_Vaughan_Shaw_HT12']], c_gl[[5]]) + theme_bw() + ggtitle("11q23.1 trans-eQTL targets in Cluster 5") + xlab("Rank") + ylab("Enrichment score") + theme(plot.title=element_text(face="bold", size = 14)) + annotate(geom="text", x=220, y=0.6, label=trans_anno, color="black", size = 6, hjust=1)

c11_anno <- paste0("NES=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "cluster11_scRNA",]$NES, 2), "\n", "p=",signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "cluster11_scRNA",]$pval, 2))
Cluster11 <- plotEnrichment(explore_sets[['cluster11_scRNA']], c_gl[[5]]) + theme_bw() + ggtitle("Enrichment of Cluster11_scRNAseq markers in Cluster 5") + xlab("Rank") + ylab("Enrichment score") + theme(plot.title=element_text(face="bold", size = 14)) + annotate(geom="text", x=220, y=0.6, label=c11_anno, color="black", size = 6, hjust=1)

Tuft_anno <- paste0("NES=", signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "Tuft",]$NES, 2), "\n", "p=",signif(fgsea_Res[[5]][fgsea_Res[[5]]$pathway == "Tuft",]$pval, 2))
Tuft <- plotEnrichment(explore_sets[['Tuft']], c_gl[[5]]) + theme_bw() + ggtitle("Tuft cell markers in Cluster 5") + xlab("Rank") + ylab("Enrichment score") + theme(plot.title=element_text(face="bold", size = 14)) + annotate(geom="text", x=220, y=0.5, label=c11_anno, color="black", size = 6, hjust=1)



## SS-GSEA
# conda install -c bioconda bioconductor-escape
library(escape)
#load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))
GS <- c(GSEABase::GeneSet(explore_sets[['FDR_sig_Vaughan_Shaw_HT12']]), GSEABase::GeneSet(smil_sets[['Tuft']][1:50]), GSEABase::GeneSet(c11$gene))
GS[[1]]@setName <- "FDR_sig_Vaughan_Shaw_HT12"
GS[[2]]@setName <- "Tuft"
GS[[3]]@setName <- "scRNAseq_cluster11"
gs <- getMatrixFromProject(ArchRProj=proj, useMatrix="GeneScoreMatrix")
gs_mat <- assay(gs_mat)
rownames(gs_mat) <- gs@elementMetadata$name
ES <- enrichIt(obj = as.matrix(gs_mat),
                      gene.sets = GS,
                      groups = 5000, cores = 2)
ES <- ES[match(rownames(proj@cellColData), rownames(ES)),]
all(rownames(proj@cellColData) == rownames(ES))
proj@cellColData <- cbind(proj@cellColData, ES)

# Plot this
gsea_pl <- vector("list", length= length(GS))
names(gsea_pl) <- c(GS[[1]]@setName, GS[[2]]@setName, GS[[3]]@setName)
for(set in seq_along(GS)){
  temp <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = names(gsea_pl)[set], embedding = "UMAP", size=0.1,
    imputeWeights = getImputeWeights(proj),
    plotAs="points",
    keepAxis=T) +
    ggtitle("") +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme(plot.title=element_text(face="bold", size = 14), legend.position="none") +
    theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10))
  gsea_pl[[set]] <- temp
}


# Pot all together (along with UMAP)
library(pathwork)
pclust <- plotEmbedding(ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP",
    size=0.4,
    plotAs="points",
    keepAxis=T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size = 10)) + ggtitle("")

ppi=300
png(file="../../Fig_plots/Thesis/Chapter2/umap_cluster_enrichments.png", width=16*ppi, height=8*ppi, res=ppi)
pclust + ((Trans + gsea_pl[[1]]) / (Tuft + gsea_pl[[2]]))
dev.off()

pclust + ((Trans + gsea_pl[[1]]) / (Tuft + gsea_pl[[2]]))
