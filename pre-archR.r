###### Bradley Jun 2022
###### Pre-ArchR script
##### conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env r-base r-seurat signac r-matrix bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase
##### conda install -c conda-forge r-devtools
##### R
##### Sys.setenv(CONDA_BUILD_SYSROOT="/")
##### devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env

##### Filtering fragment files for the cells included in the 'Colon-Epithelial' annotations of the authors
# Move raw frags to archR directory
cd /gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/data
while read f
do
echo "$f"
cp $f ../data_archR/
done <samples.txt

# Remove the GSM.. names for ease
cd ../data_archR/
for x in *; do
mv $x `echo $x | cut -c 12-`
done


# Manipulate meta file in R
R
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/data")
meta <- read.table("catlas/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
meta <- meta[grep("colon", rownames(meta)),]
meta <- meta[grep("Colon|Tuft", meta$cell_type),]
smooth <- grep("Smooth", meta$cell_type)
meta <- meta[-smooth,]
samples.use <- names(table(meta$sample)[table(meta$sample)>2]) # 5 samples
meta.use <- meta[meta$sample %in% samples.use,] #21626 cells
samp.cell <- vector("list", length = length(samples.use))

# Save a list of the epithelial, colon cells in directory
for(samp in seq_along(samples.use)){
  temp <- rownames(meta.use[meta.use$sample == samples.use[samp],])
  samples.use[samp] <- gsub("_1", "", samples.use[samp])
  temp <- unlist(strsplit(temp, "\\+"))[c(F,T)]
  write.table(temp, paste0("../data_archR/", "use_cells_", samples.use[samp], ".txt"), sep = "\t", quote=F, row.names=F, col.names=F)
}

# Save a list of the usefule samples
samples.use <- paste0(samples.use, "_rep1_fragments.bed")
write.table(samples.use, "../data_archR/samples_use.txt", sep = "\t", quote=F, row.names=F, col.names=F)

# Run the file intersection script
