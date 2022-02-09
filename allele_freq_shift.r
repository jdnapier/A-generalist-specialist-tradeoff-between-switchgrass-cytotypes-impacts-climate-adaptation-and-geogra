# Script for calculating tha Allele Frequency Shift from in a test set of
#  samples from frequency in one training set towards the frequenciy in
#  another training set

# For example, can be used for looking at the allele frequency shift from 
#  POP_A towards POP_B in a subset of samples from POP_A (this notation is used
#  in the explanation of the arguements below)

# This script uses some custom function from `custom_sg_analysis_functions.r`
#   The path to this file must be set in `custom_func_file` in the script

#########

# args[1] = vcf_file = Full path and file name for the VCF file containing the
#			target SNPs - generally these are SNPs with high Fst
#			between 2 focal populations (ex POP_A and POP_B) and
#			are the same SNPs present in the files in args[2] and
#			args[3]

# args[2] = allele_state_file = the allele states of the SNPs in the analysis,
#				  indicating which alleles are indicative of
#				  which populations at each SNP. This file can 
#				  be generated using the `set_allele_states.r` 
#				  comparing POP_A and POP_B training sets

# args[3] = training_freq_file = file containing the REF allele frequency in
#				   the training sets of both populations of
#				   interest in this analysis. This file can be
#				   generated using the `set_allele_states.r`
#			 	   script comparing POP_A and POP_B training
#				   sets

# args[4] = samp_name_file = Full path and file name of file containing list
#				of samples to be examined, in same format
#				as the names in the VCF. For example, the
#				names of samples-of-interest from POP_A

# args[5] = out_dir = Path to directory where output file will be saved

# args[6] = out_pre = Desired prefix for filenames of output file

# args[7] = max_miss = the maximum percentage of missing data allowed at a 
#			 SNP to be included in the analysis. ex: 0.2

# Explanation of output:
# Table with info about each SNP, including a column that contains the percent
#   shift from POP_A to POP_B in the samples-of-interest. This value depends
#   on the order of the populations in the files used for args[2] and args[3],
#   so it is recommended to compare allele frequencies in the samples-of-
#   interest to each training set to make sure results make sense. If the
#   order of populations is reversed, then need to swap columns 3 and 4 in 
#   args[3]

###########

args <- commandArgs(trailingOnly = T)

### LOAD PACKAGES ###
library(data.table)

custom_func_file <- 'PATH/TO/custom_sg_analysis_functions.r'
source(custom_func_file)

### INPUT DATA ###

vcf_file <- args[1]
vcf_in <- read_vcf(vcf_file)

allele_state_file <- args[2]
allele_states <- fread(allele_state_file)

training_freq_file <- args[3]
train_freq <- fread(training_freq_file)

samp_name_file <- args[4]
samp_names <- fread(samp_name_file, header = F)

### SET OUTPUTS ###
out_dir <- args[5]
out_dir <- add_slash(out_dir)
out_pre <- args[6]
out_file <- paste(out_dir, out_pre, '.Fvals.txt', sep = '')

### SET VARIABLES ###
max_miss <- as.numeric(args[7])

##########

# Calculate the REF allele frequency in the test samples
snp_names <- paste(allele_states$CHR, allele_states$POS, sep = '_')

genos_1 <- vcf_many_SNP_genotypes(vcf = vcf_in, snp_name_vec = snp_names,
  samp_vec = samp_names$V1)

genos_1[genos_1 == './.'] <- NA
genos_1[genos_1 == '0/0'] <- 2
genos_1[genos_1 == '0/1'] <- 1
genos_1[genos_1 == '1/1'] <- 0

# generate matrix in format of rows = Samples, columns = SNPs
genos_1_mat <- as.matrix(genos_1[, 2:ncol(genos_1)])

# calculate allele frequencies in test sample set
n_miss <- apply(genos_1_mat, 1, function(x) sum(is.na(x)))
n_genos <- ncol(genos_1_mat) - n_miss

max_miss_n <- length(samp_names$V1) * max_miss
max_miss_inds <- which(n_miss > max_miss_n)

allele_count <- apply(genos_1_mat, 1, function(x) sum(as.numeric(x),
  na.rm = T))

ref_af <- allele_count / (n_genos * 2)
ref_af[max_miss_inds] <- NA

# Calculate the REF freq difference between pop2 and pop1 and convert
#  to their range
ref_freq_dif <- unlist(train_freq[,4] - train_freq[,3])
ref_freq_range <- abs(ref_freq_dif)

# find SNPs where pop1 has higher REF freq
pop1_ref_hi <- which(ref_freq_dif < 0)

# calc diff in REF freq between subgrp and pop1
sub_v_pop1_ref <- ref_af - unlist(train_freq[,3])

# adjust diff for SNPs where pop1 has higher REF freq than pop2 because
#  want this to represent difference in direction of pop2
sub_v_pop1_ref_2 <- sub_v_pop1_ref
sub_v_pop1_ref_2[pop1_ref_hi] <- sub_v_pop1_ref[pop1_ref_hi]*-1

F_sub_v_pop1 <- sub_v_pop1_ref_2 / ref_freq_range

# Generate final File
# add "cumulative positions" for plotting
train_freq <- add_cumulative_pos(train_freq)

train_freq[, REF_STATE := allele_states[, list(REF_STATE)]]
train_freq[, ALT_STATE := allele_states[, list(ALT_STATE)]]

train_freq[, subgrp_ref_freq := ref_af]

train_freq[, F_subgrp_v_pop1 := F_sub_v_pop1]

fwrite(train_freq, out_file, sep = '\t')

quit(save = 'no')


