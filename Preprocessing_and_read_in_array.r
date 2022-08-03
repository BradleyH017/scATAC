# Bradley March 2022
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


# Subset for common peaks - see Processing and read in r script for this
common_peaks <- read.csv("common_peaks.csv")
res <- res[rownames(res) %in% common_peaks$x,]


# 2. ~~~~~~~~~~~~ Filtration

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
  assay = "peaks",
  min.cells = 0,
  min.features = 0
)
seur$dataset <- sample
print("Done seurat object")


# annotate
Annotation(seur) <- annotations
print("Annotated")
g=genes(edb)
# Generate the peak matrix
seur@assays$peak_matrix <- FeatureMatrix(
  fragments = Fragments(seur),
  features = g
)
print(dim(seur))
print("Peak matrix generated")

# Explore alternative way to create peaks with seurat
explore_peaks <- F
if(explore_peaks == T){
  seur <- CallPeaks(
    object = seur,
    assay="peaks",
    group.by = NULL,
    macs2.path = "/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/snapatac_env2/bin/macs2"
  )
}

# Now generate the bin matrix
seur@assays$bin_matrix <- GenomeBinMatrix(
  fragments = Fragments(seur),
  genome = seqlengths(seur),
  binsize = 5000
)
print("Bin matrix generated")

# save
save(seur, file=paste("objects/pp_in_samps/", sample, "_pp.Rds", sep = ""))

# This seems to work within sample, but when combined together was causing problems. Need to properly examine the way in which these objects are combined.
