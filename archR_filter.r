##### Bradley June 2022
#### Intersecting the raw bed files based on their intersection with intersecting cells as annotated by the authors
##### conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env r-base r-seurat signac r-matrix bioconductor-genomeinfodb  bioconductor-ensdb.hsapiens.v86 r-ggplot2 r-patchwork  bioconductor-biovizbase
##### conda install -c conda-forge r-devtools
##### R
##### Sys.setenv(CONDA_BUILD_SYSROOT="/")
##### devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
##### source activate /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/scATAC/archR_env

# Set up
setwd("/gpfs/igmmfs01/eddie/CCGG-tumour-WGS/BradTemp/scATAC/data_archR")
library(optparse)
option_list = list(make_option(c("-s", "--sample"), action = "store", default = NA, type ="character", help="which sample to use"))
opt = parse_args(OptionParser(option_list=option_list))
sample = opt$s;
print(paste("The option for sample is", sample, sep = " "))

# Load the fragment file
frags <- read.table(sample, header=F)

# Subset
use_cells <- readLines(paste0("use_cells_", sample))
frags <- frags[frags$V4 %in% use_cells,]

# Save
write.table(frags, paste0("filtered_", sample), sep = "\t", quote=F, row.names=F, col.names=F)
