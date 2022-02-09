# Custom functions required by some of the scripts of this analysis

add_slash <- function(dir_string){
  # Add '/' to end of a string, particularly directory names that will be
  #  used with paste() to generate file names
  # INPUTS
  # dir_string = character string; ex: '/FULL/DIR/STRING'
  # OUPUT
  # character string
  #   if dir_string already ends with '/', then returns dir_string
  #   if dir_string does not end with '/', returns dir_string ending with '/'
  #########
  final_string <- dir_string
  last_char <- rev(unlist(strsplit(dir_string, split = '')))[1]
  if(last_char != '/'){
    final_string <- paste(final_string, '/', sep = '')
  }
  return(final_string)
}

gen_gl_preobj <- function(vcf, oct_libs = c(), tet_libs = c(),
  maf_cut = 0.001){
  # Function for filtering vcf by MAF freq and generating necessary
  #  objects to be used for making a "genlight" object
  # INPUTS
  # vcf = inputted vcf with header; needs following columns:
  #        ID with snp names
  #        CHROM with chromosome name
  #        POS with SNP position
  #        REF with reference allele
  #        ALT with alternate allele
  #        Sample library names as column names
  # oct_libs = names of 8X libraries
  # tet_libs = names of 4X libraries
  # maf_cut = minor allele frequency cutoff
  # OUTPUT
  # list with following elements, elements of list are based on vcf with 
  #   SNPs with MAF below the cutoff removed
  # 'geno_mat' = matrix of genotypes showing the number of ALT alleles;
  #                maximum number is 2 for 4X and 4 for 8X libraries
  # 'loc.names' =  vector of the names of the retained SNPs
  # 'chromosome' = vector of the chromosome name(s)
  # 'position' = vector of the position of the SNPs
  # 'ploidy' = vector of the ploidy for each sample
  #              2 = 4X switchgrass (because of disomic inheritance)
  #              4 = 8X switchgrass (because of tetrasomic inheritance)
  # 'loc.all' = vector of the REF and ALT alleles, ex: 'g/a'
  #############
  geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
  oct_alt_dose_vec <- c('0', '1', '2', '3', '4')
  tet_alt_dose_vec <- c('0', '1', '1', '1', '2')
  #
  if(length(oct_libs) > 0){
    oct_df_tmp <- vcf[, oct_libs]
    oct_df <- data.frame(apply(oct_df_tmp, 2, function(y) 
      unlist(lapply(strsplit(y, split = ':'), function(x) x[1]))), 
      stringsAsFactors = F)
    oct_df[oct_df == './.'] <- NA
    for(i in seq(length(geno_vec))){
      oct_df[oct_df == geno_vec[i]] <- oct_alt_dose_vec[i]
    }
    for(j8 in seq(ncol(oct_df))){
      oct_df[, j8] <- as.numeric(oct_df[, j8])
    }
    sum_alt_8 <- apply(oct_df, 1, function(x) sum(unlist(x), na.rm = T))
  } else{sum_alt_8 <- c()}
  #
  if(length(tet_libs) > 0){
    tet_df_tmp <- vcf[, tet_libs]
    tet_df <- data.frame(apply(tet_df_tmp, 2, function(y)
      unlist(lapply(strsplit(y, split = ':'), function(x) x[1]))), 
      stringsAsFactors = F)
    tet_df[tet_df == './.'] <- NA
    for(i in seq(length(geno_vec))){
      tet_df[tet_df == geno_vec[i]] <- tet_alt_dose_vec[i]
    }
    for(j4 in seq(ncol(tet_df))){
      tet_df[, j4] <- as.numeric(tet_df[, j4])
    }
    sum_alt_4 <- apply(tet_df, 1, function(x) sum(unlist(x), na.rm = T))
  } else{sum_alt_4 <- c()}
  #
  if(length(oct_libs) == 0){
    tot_df <- tet_df
    sum_alt_total <- sum_alt_4
  } else if(length(tet_libs) == 0){
    tot_df <- oct_df
    sum_alt_total <- sum_alt_8/2
  } else{
    tot_df <- cbind(oct_df, tet_df)
    sum_alt_total <- (sum_alt_8/2) + sum_alt_4
  }
  #
  n_nas <- apply(tot_df, 1, function(x) sum(is.na(unlist(x))))
  n_gsamps <- ncol(tot_df) - n_nas
  allele_freq <- sum_alt_total / (2*n_gsamps)
  allele_freq_1 <- 1 - allele_freq
  minor_af <- apply(cbind(allele_freq, allele_freq_1), 1, function(x) min(x)[1])
  remove_inds <- which(minor_af < maf_cut)
  #
  keep_inds <- setdiff(seq(nrow(vcf)), remove_inds)
  #
  preobj_ls <- list()
  preobj_ls[['geno_mat']] <- tot_df[keep_inds, ]
  preobj_ls[['loc.names']] <- vcf$ID[keep_inds]
  preobj_ls[['chromosome']] <- vcf$CHROM[keep_inds]
  preobj_ls[['position']] <- vcf$POS[keep_inds]
  preobj_ls[['ploidy']] <- c(rep(4, times = length(oct_libs)),
    rep(2, times = length(tet_libs)))
  preobj_ls[['loc.all']] <- paste(tolower(vcf$REF[keep_inds]),
    tolower(vcf$ALT[keep_inds]), sep = '/')
  #
  return(preobj_ls)
}

gen_gl_object <- function(preobj_list){
  # Function to generate 'genlight' object from the 'preobj_list' generated
  #   using the 'gen_gl_preobj' function
  # INPUTS
  # preobj_list = list generated by 'gen_gl_preobj' that contains elements
  #                for making the 'genlight' object
  # OUTPUT
  # 'genlight' object that includes genotypes, sample ploidy (2 = 4X 
  #    switchgrass, 4 = 8X switchgrass), SNP name, SNP chromosome, 
  #   SNP position, and alleles at each SNP
  ######
  gl <- new('genlight', gen = t(as.matrix(preobj_list[['geno_mat']])),
    ploidy = preobj_list[['ploidy']],
    loc.names = preobj_list[['loc.names']],
    chromosome = preobj_list[['chromosome']],
    position = preobj_list[['position']],
    loc.all = preobj_list[['loc.all']]
  )
  return(gl)
}

read_vcf <- function(vcf_file){
  # Function for inputing a VCF into R for analysis
  #   mainly for inputting VCFs generated for switchgrass analysis; not sure
  #   about the universality of these functions
  #   - does NOT work for .gz files - those need to be uncompressed first
  # INPUTS #
  # vcf_file = full path to the vcf file; must already be uncompressed
  # OUTPUT #
  # data.table with VCF info, each column still in the original formatting
  ###########
  vcf_in <- read.table(vcf_file, header = F, stringsAsFactors = F)
  vcf_head_tmp <- system(paste('grep CHR ', vcf_file, sep = ''),
    intern = T)
  vcf_head_2 <- sub('#CHROM', 'CHROM', vcf_head_tmp)
  vcf_head_3 <- unlist(strsplit(vcf_head_2, split = '\t'))
  colnames(vcf_in) <- vcf_head_3
  vcf_in <- data.table(vcf_in)
  return(vcf_in)
}

vcf_many_SNP_genotypes <- function(vcf, snp_name_vec, samp_vec = c()){
  # Get genotypes for samples at multiple SNPs
  # INPUTS #
  # vcf = vcf loaded into R
  # snp_name_vec = vector of snp_names to be included, ; usually in format 
  #  of 'CHRNAME_CHRPOS'
  # samp_vec = vector of sample names. If blank, then returns genotypes for
  #  all samples
  # OUTPUT #
  # data.table of genotypes, first column = ID (snp_name), other columns
  #  are each sample in samp_vec
  ##########
  # check samp_vec
  if(length(samp_vec) == 0){
    format_col <- which(colnames(vcf) == 'FORMAT')
    samp_vec <- colnames(vcf)[(format_col + 1):ncol(vcf)]
  }
  if(length(samp_vec) > 0){
    miss_names <- setdiff(samp_vec, colnames(vcf))
    if(length(miss_names) > 0){
      stop(paste(paste(miss_names, collapse = ','), 'missing from VCF'))
    }
  }
  miss_snps <- setdiff(snp_name_vec, vcf$ID)
  if(length(miss_snps) > 0){
    stop(paste(paste(miss_snps, collapse = ','), 'missing from VCF'))
  }
  vcf_inds <- which(vcf$ID %in% snp_name_vec)
  geno_df <- apply(vcf[vcf_inds, samp_vec, with = F], 2, function(x)
    unlist(lapply(strsplit(x, split = ':'), function(x) x[1])))
  geno_tab <- data.table(ID = vcf$ID[vcf_inds], geno_df)
  return(geno_tab)
}

add_cumulative_pos <- function(ref_freq_tab, chrom_space = 1e5){
  #######
  # calculate the cumulative position of SNPs across all chromosome in 
  #  'ref_freq_tab' so can plot linearly
  ### INPUTS ###
  # ref_freq_tab = table with CHR, POS, and REF_freq columns for different
  #                   populations
  # chrom_space = the space to put between the last postion of a chromosome
  #                 and  position 1 of the next chromosome
  ### OUTPUT ###
  # same date table as 'ref_freq_tab' but with an additional column, POS_CUM
  #########
  max_chrom_pos <- c()
  for(chrm in unique(ref_freq_tab$CHR)){
    max_pos <- max(ref_freq_tab$POS[ref_freq_tab$CHR == chrm])
    max_chrom_pos[[chrm]] <- max_pos
  }
  max_chrom_pos <- unlist(max_chrom_pos)
  # figure out the cumulative amount to add to each chromosome
  cum_chrom_amount <- rep(0, times = length(max_chrom_pos))
  for(i in c(2:length(cum_chrom_amount))){
    cum_chrom_amount[i] <- sum(max_chrom_pos[1:(i-1)])+(chrom_space*(i-1))
  }
  names(cum_chrom_amount) <- names(max_chrom_pos)
  #
  ref_freq_tab[, POS_CUM := as.numeric(NA)]
  for(chrm in unique(ref_freq_tab$CHR)){
    tmp_inds <- which(ref_freq_tab$CHR == chrm)
    tmp_new_pos <- ref_freq_tab$POS[tmp_inds] + cum_chrom_amount[chrm]
    ref_freq_tab[tmp_inds, POS_CUM := tmp_new_pos]
  }
  return(ref_freq_tab)
}





