# Assign likely ancestry origin to the alleles that can be used to look for 
#  evidence of introgressions

# INPUT EXPLANATIONS
# pop1_hw_file = input [1]; file path for vcfTools hardy-weinberg output for 
#			training population 1 (pop1)
# pop2_hw_file = input [2]; file path for vcfTools hardy-weinberg output for 
#			training population 2 (pop2)
# pop1_name = input [3]; name/abbreviation for pop1, used for file names
# pop2_name = input [4]; name/abbreviation for pop2, used for file names
# out_dir = input [5]; the directory where output files will be saved
# missing_cut = input [6]; the minimum percentage of samples in both 
#			training set that SNP but be present in to be retained
#			ex: 0.8
# n_training = input [7]; the number of samples in a training set; currently
#			assumes training sets use same number, but can
#			update in the future in need be; ex: 40
# min_keep_freq = input [8]; the minimum allele frequency for an allele to be
#			retained as informative; ex: 0.5
# prob_allele_ratio = input [9]; the probability allele ratio of popA/popB
#			for an allele to be considered informative

# Explanation of outputs
# allele_freq_file = table with REF allele freq for both pops at each 
#   retained SNP position
# keep_position_file = table with CHR and POS for retained SNPs; to be used
#   for future vcfTools needs
# allele_state_file = table REF and ALT allele state designation at each
#   retained SNP

### LOAD ENVIRONMENT ###
args <- commandArgs(trailingOnly=T)

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
pop1_hw_file <- args[1]
pop1_hw <- fread(pop1_hw_file)

pop2_hw_file <- args[2]
pop2_hw <- fread(pop2_hw_file)

### SET OUTPUT ###
pop1_name <- args[3]
pop2_name <- args[4]

out_dir <- args[5]
tmp_string <- rev(unlist(strsplit(out_dir, split = '')))
if(tmp_string[1] != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

allele_freq_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_ref_freq.txt', sep = '')

keep_position_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_keep_pos.txt', sep = '')

allele_state_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_allele_states.txt', sep = '')

### SET VARIABLES ###
missing_cut <- as.numeric(args[6])

n_training <- as.numeric(args[7])

min_keep_freq <- as.numeric(args[8])

prob_allele_ratio <- 2.5

#####################

### Filter out sites with high missing data

pop1_geno_count_pre <- lapply(strsplit(
  unlist(pop1_hw[, c('OBS(HOM1/HET/HOM2)')]), 
  split = '/', fixed = T), function(x) as.numeric(x))

pop1_n_genos <- unlist(lapply(pop1_geno_count_pre, function(x) sum(x)))
pop1_lowcount_inds <- which(pop1_n_genos < (n_training * missing_cut))

pop2_geno_count_pre <- lapply(strsplit(
  unlist(pop2_hw[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

pop2_n_genos <- unlist(lapply(pop2_geno_count_pre, function(x) sum(x)))
pop2_lowcount_inds <- which(pop2_n_genos < (n_training * missing_cut))

lowcount_inds <- sort(union(pop1_lowcount_inds, pop2_lowcount_inds))

pop1_hw_filt <- pop1_hw[-lowcount_inds]
pop2_hw_filt <- pop2_hw[-lowcount_inds]

#######
# tally the number of samples with each genotype at each SNP and the
#   REF and ALT frequency at each SNP in each pop
pop1_geno_count <- lapply(strsplit(
  unlist(pop1_hw_filt[, c('OBS(HOM1/HET/HOM2)')]),    
  split = '/', fixed = T), function(x) as.numeric(x))

pop1_ref_freq <- unlist(lapply(pop1_geno_count, function(x) 
  (x[1] * 2 + x[2])/(sum(x)*2)))
pop1_alt_freq <- unlist(lapply(pop1_geno_count, function(x) 
  (x[3] * 2 + x[2])/(sum(x)*2)))

pop2_geno_count <- lapply(strsplit(
  unlist(pop2_hw_filt[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

pop2_ref_freq <- unlist(lapply(pop2_geno_count, function(x) 
  (x[1] * 2 + x[2])/(sum(x)*2)))
pop2_alt_freq <- unlist(lapply(pop2_geno_count, function(x) 
  (x[3] * 2 + x[2])/(sum(x)*2)))

summary(abs(pop2_ref_freq - pop1_ref_freq))

freq_table <- data.table(CHR = pop1_hw_filt$CHR, POS = pop1_hw_filt$POS,
  pop1_ref_freq = pop1_ref_freq,
  pop2_ref_freq = pop2_ref_freq)

pop_freq_names <- paste(c(pop1_name, pop2_name), '_ref_freq', sep = '')
colnames(freq_table)[c(3,4)] <- pop_freq_names

fwrite(freq_table, file = allele_freq_file)

###
# calculate the probabilty of detecting the REF of ALT allele at least
#  once in a sample given the training sets allele frequencies
pop1_prob_ref_allele <- sapply(pop1_ref_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))
pop1_prob_alt_allele <- sapply(pop1_alt_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))

pop2_prob_ref_allele <- sapply(pop2_ref_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))
pop2_prob_alt_allele <- sapply(pop2_alt_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))

pop2_only_REF <- which(pop1_ref_freq == 0)
pop1_only_REF <- which(pop2_ref_freq == 0)
pop1_mainly_REF <- which(
  pop1_ref_freq > min_keep_freq &
  (pop1_prob_ref_allele/pop2_prob_ref_allele) > prob_allele_ratio &
  pop2_ref_freq > 0)
pop2_mainly_REF <- which(
  pop2_ref_freq > min_keep_freq &
  (pop2_prob_ref_allele / pop1_prob_ref_allele) > prob_allele_ratio &
  pop1_ref_freq > 0)
tmp_REF_tot <- sort(c(pop2_only_REF, pop1_only_REF, 
  pop1_mainly_REF, pop2_mainly_REF))
noinfo_REF <- setdiff(seq(length(pop1_ref_freq)), tmp_REF_tot)

length(tmp_REF_tot)
# 35378
sum(duplicated(tmp_REF_tot))
# 0
length(noinfo_REF)
# [1] 23451

pop2_only_ALT <- which(pop1_alt_freq == 0)
pop1_only_ALT <- which(pop2_alt_freq == 0)
pop1_mainly_ALT <- which(
  pop1_alt_freq > min_keep_freq &
  (pop1_prob_alt_allele/pop2_prob_alt_allele) > prob_allele_ratio &
  pop2_alt_freq > 0)
pop2_mainly_ALT <- which(
  pop2_alt_freq > min_keep_freq &
  (pop2_prob_alt_allele / pop1_prob_alt_allele) > prob_allele_ratio &
  pop1_alt_freq > 0)
tmp_ALT_tot <- sort(c(pop2_only_ALT, pop1_only_ALT,
  pop1_mainly_ALT, pop2_mainly_ALT))
noinfo_ALT <- setdiff(seq(length(pop1_alt_freq)), tmp_ALT_tot)

length(tmp_ALT_tot)
# 43724
sum(duplicated(tmp_ALT_tot))
# 0
length(noinfo_ALT)
# 15105 

tmp_BOTH <- c(tmp_REF_tot, tmp_ALT_tot)
length(setdiff(seq(length(pop1_alt_freq)), tmp_BOTH))
# 0 not included in any category 0
# 45 not included in any category if use prob_allele_ratio = 2.67
# 332 not included in any category if use prob_allele_ratio = 3

sum(duplicated(tmp_BOTH))
# 20273 SNPs informative for both alleles

allele_assign_tab <- data.table(CHR = pop1_hw_filt$CHR,
  POS = pop1_hw_filt$POS,
  REF_STATE = as.character(NA),
  ALT_STATE = as.character(NA))

allele_assign_tab[pop2_only_REF, REF_STATE := paste(pop2_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop1_only_REF, REF_STATE := paste(pop1_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop2_mainly_REF, REF_STATE := paste(pop2_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[pop1_mainly_REF, REF_STATE := paste(pop1_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[noinfo_REF, REF_STATE := 'NO_INFO']

allele_assign_tab[pop2_only_ALT, ALT_STATE := paste(pop2_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop1_only_ALT, ALT_STATE := paste(pop1_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop2_mainly_ALT, ALT_STATE := paste(pop2_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[pop1_mainly_ALT, ALT_STATE := paste(pop1_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[noinfo_ALT, ALT_STATE := 'NO_INFO']

fwrite(allele_assign_tab, file = allele_state_file,
  sep = '\t')

fwrite(allele_assign_tab[, c('CHR', 'POS')], 
  file = keep_position_file, sep = '\t')

quit(save = 'no')

