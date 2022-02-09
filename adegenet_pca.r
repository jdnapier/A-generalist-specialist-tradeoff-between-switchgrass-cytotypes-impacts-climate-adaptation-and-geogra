# Script for PCA analysis using adegenet

# Script uses some custom functions from `custom_sg_analysis_function.r`;
#  the path to this file must be set in `custom_func_file` below

########

# args[1] = data_file = Full path and file name for .rds file containing 
#                         genlight genotype object

# args[2] = samp_file = Full path and file name for list of samples to 
#			  be included in analysis. If `all`, then all samples
#			  in the genlight object are included

# args[3] = snp_file = Full path and file name for list of SNPs to be
#			 included in the analysis. If `all`, then all SNPs
#			 in the genlight object are included

# args[4] = out_dir = The directory in which results will be saved

# args[5] = out_pre = Filename prefix for the results file

# args[6] = maf_cut = the minor allele cutoff for SNPs to be included in the
#			analysis. ex: 0.05

# args[7] = max_miss = the minimum percentage of samples with a genotype for
#			 a SNP to be included - appologies for the confusing
#			 variable name. ex: 0.8

args = commandArgs(trailingOnly = TRUE)

custom_func_file <- '/PATH/TO/custom_sg_analysis_functions.'
source(custom_func_file)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)
library(data.table)

### INPUT DATA ###
data_file <- args[1]
tot_gl_0 <- readRDS(data_file)

samp_file <- args[2]
if(samp_file == 'all'){
  keep_samps <- indNames(tot_gl_0)} else {
    keep_samps <- fread(samp_file, header = F)[, V1]
}

snp_file <- args[3]
if(snp_file == 'all'){
  keep_snps <- locNames(tot_gl_0)} else {
    keep_snps <- fread(snp_file, header = F)[, V1]
}

### SET OUTPUT ###
out_dir <- args[4]
out_dir <- add_slash(out_dir)

#base_in <- basename(data_file)
out_pre <- args[5]
out_file_short <- paste(out_pre, '.PCAresults.rds', sep = '')
out_file <- paste(out_dir, out_file_short, sep = '')

### SET VARIABLES ###
maf_cut <- as.numeric(args[6])

max_miss <- as.numeric(args[7])

############
# Filter by SNPs and Samples
samp_inds <- which(indNames(tot_gl_0) %in% keep_samps)
snp_inds <- which(locNames(tot_gl_0) %in% keep_snps)

tot_gl <- tot_gl_0[samp_inds, snp_inds]

# calculate MAF
tet_inds <- which(ploidy(tot_gl) == 2)
oct_inds <- which(ploidy(tot_gl) == 4)

n_tet_na <- apply(as.matrix(tot_gl[tet_inds,]), 2, function(x) sum(is.na(x)))
n_tet_genos <- length(tet_inds) - n_tet_na
n_oct_na <- apply(as.matrix(tot_gl[oct_inds,]), 2, function(x) sum(is.na(x)))
n_oct_genos <- length(oct_inds) - n_oct_na

# filter by missing data
miss_vec <- (n_tet_na + n_oct_na) / (nInd(tot_gl))

keep_miss_inds <- which(miss_vec <= max_miss)

# calculate allele frequencies
sum_tet <- apply(as.matrix(tot_gl[tet_inds]), 2, function(x) 
  sum(x, na.rm = T)/2)
sum_oct <- apply(as.matrix(tot_gl[oct_inds]), 2, function(x) 
  sum(x, na.rm = T)/4)

tot_af <- (sum_tet + sum_oct) / (n_tet_genos + n_oct_genos)
over_inds <- which(tot_af > 0.5)

# filter by MAF
maf_vec <- tot_af
maf_vec[over_inds] <- 1-tot_af[over_inds]
keep_maf_inds <- which(maf_vec >= maf_cut)

# select SNPs that pass both filters
keep_snp_inds <- intersect(keep_miss_inds, keep_maf_inds)

filt_gl <- tot_gl[, keep_snp_inds]

# fun PCA
n_eig <- nInd(filt_gl) - 1

tot_pca <- glPca(filt_gl, nf = n_eig, loadings = F, alleleAsUnit = F, useC = F)

saveRDS(tot_pca, out_file)

quit(save = 'no')

