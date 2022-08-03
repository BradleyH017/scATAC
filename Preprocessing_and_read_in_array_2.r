# Bradley April 2022
# Processing of scATAC files
# cd /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/
# conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/scATAC_env r-base=4.1.2 r-seurat signac r-matrix bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase
# source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/scATAC_env
# See ATAC demultiplex for how these were generated

# Have adjusted this so runs on each sample individually in an array, then can combine later
library(optparse)
option_list = list(make_option(c("-s", "--sample"), action = "store", default = NA, type ="character", help="which sample to use"))
opt = parse_args(OptionParser(option_list=option_list))
sample = opt$s;
print(paste("The option for sample is", sample, sep = " "))


# Load in the sparse matrices for this sample
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC")


# 1. ~~~~~~~~~~~~ Load in the matrices to a list of sparse matrices
library("Matrix")
res <- readMM(paste("data/processed/", sample, "_rep1_fragments.mtx", sep = ""))
peaks <- read.delim(paste("data/processed/", sample, "_rep1_fragments_coo.bin.ygi", sep = ""), sep = "\t", header =F)
peak_tog <- paste(paste(peaks[,1], peaks[,2], sep = ":"), peaks[,3], sep = "_")
colnames(res) <- peak_tog
rownames(res) <- paste(sample, readLines(paste("data/", "barcodes_final_raw_temp_", sample, "_rep1_fragments.bed.xgi", sep = "")), sep = ".")
res <- t(res)
print(paste("done loading", sample))


# 2. ~~~~~~~~~~~~ Initial filtration

# Paper says that cells with < 1000 total calls were excluded
res <- res[,colSums(res) > 1000]



# Am now following seurat workflow - would MUCH rather do this in R and seems relatively straightforward with this
# https://satijalab.org/signac/articles/pbmc_vignette.html#integrating-with-scrna-seq-data-1
# Generate a chromatin assay object / seurat object for each sample - originally had used GRCh38, now trying hg38
library(BiocGenerics)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

# Generate annotations
edb = EnsDb.Hsapiens.v86
annotations <- GetGRangesFromEnsDb(ensdb = edb, standard.chromosomes = TRUE)
# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'

# 3. ~~~~~~~~~~~~~ Proper processing. Generate all of the required matrices for the analsysis
# Now process this sample - Error: Overlapping ranges supplied. Ranges should be non-overlapping
#temp <- res[1:1000,1:1000]
colnames(res) <- unlist(strsplit(colnames(res), "\\."))[c(F,T)]
chrom_test <- CreateChromatinAssay(
  counts=res,
  sep = c(":","_"),
  genome = 'hg38',
  fragments=paste("data/frags_", sample, "_rep1_fragments.bed.gz", sep = ""),
)
print("Done chromatin assay")


# Now create seurat object
seur <- CreateSeuratObject(
  counts = chrom_test,
  assay = "fragments",
  min.cells = 0,
  min.features = 0
)
seur$dataset <- sample
print("Done seurat object")


# annotate
Annotation(seur) <- annotations
print("Annotated")
g=genes(edb)


# ~~~~~~~~~~~~ Now generate the peak matrix (for real this time)
peaks.gr <- CallPeaks(
  object = seur,
  assay="fragments",
  group.by = NULL,
  macs2.path = "/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/snapatac_env2/bin/macs2"
)
print("Peaks called")

# Also count fragments for each sample
frags.list <- CreateFragmentObject(
path = paste("data/frags_", sample,
  "_rep1_fragments.bed.gz", sep = ""),
  cells = colnames(seur@assays$fragments))
print("Created fragment object")

# Load the combined peaks (see the next script 'combine_seur_objects_v2' for how this was done)
combined.peaks <- read.csv("data/seur_r_objects_tog/combined.peaks.all.samps.csv", row.names=1)
combined.peaks.gr <- makeGRangesFromDataFrame(combined.peaks)
print("Loaded combined peaks")

# Can now create a peak x cell matrix for each sample
counts.list <- FeatureMatrix(
  fragments = frags.list,
  features = combined.peaks.gr,
  cells = colnames(seur@assays$fragments))
print("Generated the peak x cell matrix")
counts.list <- counts.list[,colnames(counts.list) %in% colnames(seur@assays$fragments)]

# Add this to the seurat object
atac_assay <- CreateChromatinAssay(counts.list,
  fragments = frags.list,
  genome = 'hg38',
  min.cells = 0,
  min.features = 0)
seur[['ATAC']] <- atac_assay
print("Added the cell x peak matrix")

# Count frags - for QC in next script
frags.res <- CountFragments(paste("data/frags_", sample,
  "_rep1_fragments.bed.gz", sep = ""),
  cells = colnames(seur@assays$fragments))
seur@meta.data <- cbind(seur@meta.data, frags.res)
print("Counted fragments, added to meta")

# Calculate FRiP
seur <- FRiP(
  object = seur,
  assay = 'ATAC',
  total.fragments = 'nCount_fragments'
)
print("Calculated FRiP")

# Calculate fragments in black list region
seur$blacklist_fraction <- FractionCountsInRegion(
  object = seur,
  assay = 'fragments',
  regions = blacklist_hg38_unified)

# save
save(seur, file=paste("objects/pp_in_samps/", sample, "_pp2.Rds", sep = ""))
print("Saved the processed seur object for this sample - all the trimmings")


# Do some TSS and Nucleosome signal plots
# TSS
DefaultAssay(seur) <- "ATAC"
Annotation(seur) <- annotations
seur <- TSSEnrichment(object = seur, fast = FALSE, assay="ATAC")
seur$high.tss <- ifelse(seur$TSS.enrichment > 2, 'High', 'Low')


pathOut=paste("results/QC/seur_processed_QC", sample, sep= "/")
if(file.exists(pathOut) !=T){
  dir.create(pathOut)
}
pdf(file=paste(pathOut, "TSS_ATAC_sep_2.pdf", sep = "/"))
TSSPlot(seur, group.by = 'high.tss') + NoLegend() + ggtitle(sample)
dev.off()
print("Done TSS plot")
# NS
seur <- NucleosomeSignal(object = seur)
seur$nucleosome_group <- ifelse(seur$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

pdf(file=paste(pathOut, "NucleosomeSignal_ATAC_sep_4.pdf", sep = "/"))
FragmentHistogram(object = seur, group.by = 'nucleosome_group') + ggtitle(sample)
dev.off()
print("Done NS plot")
