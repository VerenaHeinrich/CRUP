

##################################################################
#definition: headline
##################################################################

headline <- "\n*****************************************************************************************

  C                               *** This is CRUP ***                                C
  *             (C)ondition-specific (R)egulatory (U)nits (P)rediction                *	
  R                                                                                   R	
  *                                                                                   *	
  U                           (version 1.1 - 15/12/2020)                              U	
  *                                                                                   * 
  P                                     contact:                                      P 
  *                 		 heinrich@molgen.mpg.de          	              * 

*****************************************************************************************\n\n"

##################################################################
# definition: parameter spectrum:
##################################################################
        
spectrum <- matrix(c( 'norm',           'N', 0, "logical",       "computes normalized count values for .bam files",
                      'prediction',     'P', 0, "logical",       "runs CRUP - EP: (E)nhancer (P)rediction from histone modification",
                      'dynamics',       'D', 0, "logical",       "runs CRUP - ED: assigns (E)nhancer to (D)ynamic conditions",
                      'targets',        'T', 0, "logical",       "runs CRUP - ET: correlates (E)nhancer to (T)arget genes",
		      'Inputfree',	'I', 0, "logical",	 "runs CRUP - N without additional ChIP-seq Input normalization",
                      'cores',          'x', 1, "integer",       "number of cores to use (DEFAULT:1)",
                      'input',          'i', 1, "character",     "summary input for ChIP-seq experiments (in .txt format)",
                      'matrix',         'm', 1, "character",     "normalized data matrix (.rds file format)",
                      'classifier',     'c', 1, "character",     "directory of enhancer classifier (DEFAULT: DATA/CLASSIFIER/)",
                      'cutoff',         'u', 1, "double",        "cutoff for probabilities [0,1] (DEFAULT: 0.5)",
                      'distance',       'd', 1, "integer",       "maximum distance (bp) for peak clustering (DEFAULT: 12500)",
                      'genome',         'g', 1, "character",     "genome used in alignments ('hg19', 'mm10', 'mm9' or 'hg38'). Alternatively: location of a fasta index in .fai format with chromosome lengths.",
                      'sequencing',     's', 1, "character",     "type of sequencing ('paired' or 'single')",
                      'outdir',         'o', 1, "character",     "output directory (DEFAULT: same as 'file' directory)",
                      'probabilities',  'p', 1, "character",     "probabilities in rds format. list: delimiter samples: ':', delimiter conditions: ','",
                      'label',          'l', 1, "character",     "aternative labels for conditions (DEFAULT: cond1,cond2, ..)",
                      'w_0',            'w', 1, "double",        "minimum difference between group means [0,1]. (DEFAULT: 0.5)",
                      'threshold',      't', 1, "double",        "threshold for p-values in [0,1]. (DEFAULT: 0.05)",
                      'threshold_c',    'C', 1, "double",        "threshold for correlation in [0.5,1]. (DEFAULT: 0.9)",
                      'flank',          'f', 1, "integer",       "length of flanking region for summarizing differential regions. (DEFAULT: 1000)",
                      'n.flank',        'F', 1, "integer",       "number of (+/-) flanking region for summarizing differential peaks. (DEFAULT: 2)",
                      'regions',        'r', 1, "character",     "text file with condition-specific regions in txt format",
                      'RNA',            'E', 1, "character",     "RNA-seq experiments in bam format. list: delimiter samples: ':', delimiter conditions: ','",
                      'expression',     'e', 1, "character",     "gene expression counts for all samples and conditions",
                      'TAD',            'b', 1, "character",     ".bed file with TADs (DEFAULT: DATA/mESC_mapq30_KR_all_TADs.bed)",
                      'mapq',           'q', 1, "integer",       "minimum mapping quality (DEFAULT:10)",
		      'nearest',	'n', 0, "logical",	 "logical: if set, takes the nearest gene to build regulatory regions",
                      'help',           'h', 0, "logical",       "this help message"
),ncol = 5,byrow = T)


##################################################################
# definition: valid parameter values:
##################################################################

sequencing_values <- c('paired', 'single')
genome_values <- c('hg19', 'mm10', 'mm9', 'hg38')

##################################################################
#definition: output messages
##################################################################

done <- function() {cat(".. done.\n")}
skip <- function() {cat("\t ..")}
startPart <- function(m) {cat(paste0("\n--- ",m," ---\n\n"))}
endPart <- function() {cat("\n\t>>> All done!\n")}

##################################################################
#definition: standard GRanges/DataFrame header:
##################################################################

GR_header <- c("seqnames", "start","end","width")
GR_header_short <- c("seqnames", "start","end")
DF_header <- c("chr", "start","end")

##################################################################
# function: re-header
##################################################################

reheader_DF <- function(DF, header){
  colnames(DF)[1:length(header)] <- header
  return(DF)
}

##################################################################
# function: check if file exists
##################################################################

check_file <- function(f){
  if (!(file.exists(f))) { 
    cat("File ",f," does not exist.\n");
    q();
  }
}

##################################################################
# function: check if outdir exists
#
#(d = outdir, alt_d = alternative dir)
##################################################################

check_outdir <- function(d, alt_d){
  if (is.null(d)) d <- dirname(alt_d)
  if (!dir.exists(d)) {
    cat(paste0("Output directory '",d,"'is not a valid directory. \n
              Output directory is set to ",dirname(alt_d)));
    d = paste0(dirname(alt_d),"/")
  }
  return(d)
}

##################################################################
# function: stop_script
##################################################################

stop_script <- function(parameter){
  cat(paste(getopt(spectrum[(spectrum[,1] %in% c(parameter, "help")),], usage = T),"\n"));
  q();
}


##################################################################
# function: install packages
##################################################################

pkgInstall <- function(pkg) { 

  if(! isTRUE(pkg %in% .packages(all.available=T))) { 
    
    x <- tryCatch({eval(parse(text = sprintf("install.packages(\"%s\", dependencies = T)", pkg)))})

    if(is.null(x)) {

	
  	if (!'BiocManager' %in% installed.packages()){
		if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")
		BiocManager::install()
	}
	x <- tryCatch({eval(parse(text = sprintf("BiocManager::install(\"%s\")", pkg)))})

	# old: 
	# x <- tryCatch({source("http://bioconductor.org/biocLite.R")
    	# eval(parse(text = sprintf("biocLite(\"%s\")", pkg)))})
    }

    if(is.null(x)) {
	cat(paste0("Unable to install package: ",pkg,".\n"));
  	q();
    }else{
	    eval(parse(text = sprintf("require(\"%s\")", pkg)))
    }
  }
}

##################################################################
# function: check if package is installed
##################################################################

pkgTest <- function(pkg){
  if (!pkg %in% installed.packages()) {
	cat(paste0("package ",pkg, " not found. Install now."))
	return(F)
  }else{
	return(T)	
  }
} 

##################################################################
# function: load package
##################################################################

pkgLoad <- function(pkg) {
  if(pkgTest(pkg) == F){
 	pkgInstall(pkg)
  }else{
	#cat(paste0(skip(), "load package ", pkg))
  	suppressMessages(library(pkg, character.only = TRUE))
  	#done()
  }
}

##################################################################
# function: partition genome into bins
##################################################################

get_binned_genome_old_delete_me <- function(txdb, width){
  
  # get binned genome from txdb object
  binned <- tileGenome(seqinfo(txdb),
                       tilewidth = width,
                       cut.last.tile.in.chrom = TRUE)
  
  # only take autosomes and X chromosome
  seqlevels(binned, pruning.mode = "coarse") <- seqlevels(txdb)[grep("^chr[0-9]{,2}$|chrX$", seqlevels(txdb))]
  
  # delete last region in each chromosome that is smaller than 100
  return(binned[-which(width(binned) != 100)])
}


get_binned_genome <- function(seqinfo, width){
  
  # get binned genome from txdb object
  binned <- tileGenome(seqinfo,
	                         tilewidth = width, 
	                         cut.last.tile.in.chrom = TRUE)
  
  # only take autosomes and X chromosome
  if(grepl('^chr',seqlevels(seqinfo)[1]) == T) seqlevels(binned, pruning.mode = "coarse") <- seqlevels(seqinfo)[grep("^chr[0-9]{,2}$|chrX$", seqlevels(seqinfo))]
  
  # delete last region in each chromosome that is smaller than 100
  return(binned[-which(width(binned) != 100)])
}



##################################################################
# function: get summarized counts in defined bins
##################################################################

get_bamProfile <- function(bampath, gr, mapq, sequencing){
  
  # gives information about path of bam and bai file (Object)
  bf = BamFile(bampath)
  
  # gives length of different chromosomes (SeqInfo class)
  si = seqinfo(bf)
  
  # fix chromosome prefix
  seqlevels(gr) = gsub("chr","",seqlevels(gr))
  if (grepl("chr",seqlevels(si)[1])) {
    seqlevels(gr) = paste0("chr", seqlevels(gr))
  }
  
  # filter chromosomes
  gr = gr[seqnames(gr) %in% seqnames(si)]
  
  # count and summarize
  if (sequencing == "paired") {
    sapply( bamProfile( bampath = bampath, 
                        gr = gr, 
                        binsize = 100, 
                        mapqual = mapq, 
                        ss = FALSE, 
                        paired.end = "midpoint", 
                        filteredFlag = 1024,
                        verbose = FALSE)@signals, 
            function(x) x)
    
  } else if (sequencing == "single") {
    sapply( bamProfile( bampath = bampath, 
                        gr = gr, 
                        binsize = 100, 
                        mapqual = mapq, 
                        ss = FALSE, 
                        shift = 100,
                        filteredFlag = 1024,
                        verbose = FALSE)@signals, 
            function(x) x)
  }
}

##################################################################
# function: quantile normalize with target
##################################################################

get_targetQuantileNorm <- function(ecdf){
  
  ### recreate HM count vector from ecdf of mESC
  x.ref <- knots(ecdf)
  y.ref <- ecdf(x.ref)
  
  # 26337756 = nb of 100 bp bins in mm10
  temp <- c(0, round(y.ref*26337756, digits = 0))
  ret <- unlist(lapply(1:length(x.ref), function(x) rep(x.ref[x], temp[x + 1] - temp[x])))
  
  return(ret)
}


##################################################################
# function: create extended data matrix (plus/minus) bins
##################################################################

extend_dataMatrix <- function(N, df, f){
  
  N_regions <- nrow(df)
  
  # make extended feature vector
  f_ext <- NULL
  for (i in 1:N) {
    f_ext <- c(f_ext, paste0(f, "_left", i ))
    f_ext <- c(f_ext, paste0(f, "_right", i ))
  }
  
  # prepare new matrix
  df_ext <- cbind(df[,c(GR_header_short, f)], matrix(0, nrow = N_regions, ncol = length(f_ext)))
  colnames(df_ext) <- c(DF_header, f, f_ext)
  
  # make extended data matrix
  region_ext <- (N + 1):(N_regions - N)
  
  for (i in 1:N) {
    df_ext[region_ext, f_ext[((2*i - 2)*length(f) + 1):((2*i - 1)*length(f))]] <- df[region_ext - i, f]
    df_ext[region_ext, f_ext[((2*i - 1)*length(f) + 1):(2*i*length(f))]] <-  df[region_ext + i, f]
  }
  
  return(df_ext)
}


#############################################################################
# function: sort peak candidates and remove overlapping peaks
# used in peak_calling function
#############################################################################

sort_peaks <- function(peaks){
  
  # sort peaks according to probability score
  peaks <- peaks[sort(score(peaks), decreasing = T, index.return = T)$ix]
  
  count <- 0
  while (length(peaks) > (count + 1)) {
    
    count <- count + 1
    overlap.to <- findOverlaps(query = peaks[count], subject = peaks)@to
    
    if (length(overlap.to) == 1) next
    
    delete.index <- sort(overlap.to, decreasing = F)[-1]
    peaks <- peaks[-delete.index]
  }
  return(peaks)
}


#############################################################################
# function: call peaks (from probabilities)
#############################################################################

get_enhancerPeaks <- function(gr, cutoff, cores){
  
  # define and merge all regions under background probability 
  candidates <- reduce(gr[which(score(gr) > cutoff)])
  
  # possible 100 bp peaks
  peaks <- gr[findOverlaps(gr, candidates)@from]
  
  # create 1100 bp peak regions
  start(peaks) <- start(peaks) - 500
  width(peaks) <- 1100
  
  # call peaks for each chromosome separately
  out <- mclapply(split(peaks, seqnames(peaks)), 
                  sort_peaks, 
                  mc.cores = cores)
  
  return(do.call("c", unname(out)))
}


##################################################################
# function: call super enhancer candidates from called peaks
##################################################################

get_enhancerCluster <- function(peaks, peak.gap, cores){
  
  # cluster candidates (neighboring enhancers within 12.5 kb)
  red.peaks <- reduce(peaks, min.gapwidth = peak.gap)
  cluster <- red.peaks[which(width(red.peaks) > 1100)]
  
  # sort cluster accroding to maximum probability
  if (length(cluster) > 0) {
    sort.max <- unlist(mclapply(1:length(cluster), 
                                function(x) max(peaks$score[findOverlaps(cluster[x], peaks)@to]),
                                mc.cores = cores))
    cluster <- cluster[sort(sort.max, index.return = T, decreasing = T)$ix]
  }  
  return(cluster)
}


##################################################################
# function: create probabilty matrix for all given files:
##################################################################
get_probabiltyMatrix <- function(files, IDs){
  
  probs <- GRanges()
  for (i in 1:length(files)) {
    for (j in 1:length(files[[i]])) {
      
      # read rds file:
      this <- readRDS(files[[i]][j])
      colnames(elementMetadata(this)) <- IDs[[i]][j]
 
      # combine:
      if(length(probs) > 0 ){
      	probs <- merge(probs, this, all = TRUE)
      }else{
	probs <- this
	}
    }
  }

  return(probs)
}



##################################################################
# function: compute t test statistic
#
# #mean_a -mean_b <= w_0 == differenz is smaller or equal than w_0
##################################################################

get_tStatistic <- function(a, b, w_0){
  
  a <- as.matrix(a)
  b <- as.matrix(b)
  
  len <- dim(a)[1] 
  
  n_a <- dim(a)[2]
  n_b <- dim(b)[2]
  
  mean_a <- rowMeans(a,na.rm = T)
  mean_b <- rowMeans(b,na.rm = T)
  
  if (n_a == 1) {
    var_a <- rep(0,len)
  }else{
    var_a  <- rowVars(as.matrix(a),na.rm = T)
  }
  
  if (n_b == 1) {
    var_b <- rep(0,len)
  }else{
    var_b <- rowVars(as.matrix(b),na.rm = T)
  }
  
  if (n_a == 1 & n_b == 1) {
    S <- rep(1/1000000, len) # a small number
  }else{
    S_2 <- ((n_a - 1)*var_a + (n_b - 1)*var_b)/(n_a + n_b - 2)
    S   <-  sqrt(S_2)*sqrt(1/n_a + 1/n_b)
  }
  
  diff = (mean_a-mean_b)
  S[which(S == 0)] <- 1/1000000 # a small number
  res  <- sqrt((n_a*n_b)/(n_a + n_b))*(abs(diff) - w_0)/S 
  
  #generate complete cases:
  res[is.infinite(res)]  <-  NA
  #res[which(res<0)]  <-  NA

  return(res)
}  


##################################################################
# function: get empirical p values
##################################################################

get_empPvalues <- function(i, comb, w_0, IDs, probs, probs.shuffled){

	# this pairwise samples:
	IDs.1 = IDs[[comb[,i][1]]]
  	IDs.2 = IDs[[comb[,i][2]]]

	sign <- c()
	null_statistic <-c()
	statistic <- c()

	if(IDs.1[1] != 'background' & IDs.2[1] != 'background'){

		# true statistic:
		statistic <- get_tStatistic(  	probs[,IDs.1], 	
	                            		probs[,IDs.2], 
						w_0
		)

		# random statistic:
		null_statistic <- get_tStatistic( probs.shuffled[,IDs.1], 
						  probs.shuffled[,IDs.2], 
			                    	  w_0=0
		)

		# get sign (+1 or -1 or NA):
		sign <- sign(rowMeans(as.matrix(probs[,IDs.1]),na.rm = T) - rowMeans(as.matrix(probs[,IDs.2]),na.rm = T))
		sign[which(sign==0)] = NA

	}else{
		if(IDs.1[1] == 'background'){

			# true statistic:
			statistic <- get_tStatistic(  	rep(0,nrow(probs)), 	
		                            		probs[,IDs.2], 
							w_0=0.5
			)

			# random statistic:
			null_statistic <- get_tStatistic( rep(0,nrow(probs)), 
							  probs.shuffled[,IDs.2], 
				                    	  w_0=0
			)
		
			# get sign (+1 or -1 or NA):
			sign <- sign(rowMeans(as.matrix(rep(0,nrow(probs))),na.rm = T) - rowMeans(as.matrix(probs[,IDs.2]),na.rm = T))
			sign[which(sign==0)] = NA
		}

		if(IDs.2[1] == 'background'){

			# true statistic:
			statistic <- get_tStatistic(	probs[,IDs.1], 
							rep(0,nrow(probs)),
							w_0=0.5
			)

			null_statistic <- get_tStatistic( probs.shuffled[,IDs.1], 
							  rep(0,nrow(probs)), 
				                    	  w_0=0
			)

			# get sign (+1 or -1 or NA):
			sign <- sign(rowMeans(as.matrix(probs[,IDs.1]),na.rm = T) - rowMeans(as.matrix(rep(0,nrow(probs))),na.rm = T))
			sign[which(sign==0)] = NA
		}
	}


	# get two-sided pvalues:
  	ECDF <- ecdf(null_statistic)
  	p.values <- 1-ECDF(statistic)

	#p.values[which(statistic<0)] = 1
	p.values[is.na(sign)] = 1
	p.values[is.na(p.values)] = 1

	return(list(p.values=p.values, sign=sign))
}



##################################################################
# function: get pairwise p values
##################################################################

get_pairwisePvalues <- function(probs, IDs, w_0, cores){

  # define all pairwise combinations:
  comb <- combn(seq(length(IDs)),2)

  # generate null distribution:
  probs.shuffled <- sample(mcols(probs))
 
  # get empirical p.values:
  p.values <- mclapply( as.list(seq(dim(comb)[2])), mc.cores = cores,
			function(i) get_empPvalues( 	i,
						 	comb,
							w_0, 
							IDs,
							mcols(probs),
							probs.shuffled)
                        )
  names(p.values) = apply(comb,2, function(x) paste(x, collapse = ","))

  return(p.values)
}

##################################################################
# function: function to translate significant peak pattern
#
# output: vector \in [0,1]
##################################################################

get_positionPattern <- function(p, i, threshold, pattern_empty, comb, n.flank){
  
  # define region with flanking positions
  pos <- c(rev(i - 1:n.flank), i, (i + 1:n.flank))
  
  # all pairwise p-values at position i and +/- 2 positions:
  this.pvalue <- do.call(rbind, lapply(p, function(x) x$p.values[pos]))

  idx <- (this.pvalue <= threshold) == T
  #t.matrix=matrix(0.05, nrow=nrow(this.pvalue),ncol=ncol(this.pvalue))
  #t.matrix[,(n.flank+1)] <-threshold
  #idx <- (this.pvalue <= t.matrix) == T
  idx.prod <- rowProds(idx)
  
  # all flanking positions have to be significant as well:
  if (sum(idx.prod) > 0) {

    # this pairwise combination:
    this.comb <- as.matrix(comb[,(idx.prod == 1)])
    
    # check direction and change if neccessary:
    this.sign <- do.call(rbind, lapply(p, function(x) x$sign[pos]))
    idx.invert <- this.sign[(idx.prod == 1),(n.flank+1)] < 0
    if (sum(idx.invert) > 0) {
      this.comb[,idx.invert] <- this.comb[c(2,1),idx.invert]
    }

    # create labels from pairwise comparisons:
    label <- apply(this.comb,2, function(x) paste(x, collapse = ","))
    pattern <- pattern_empty
    pattern[,label] <- 1
    
    if(length(pattern) == length(pattern_empty)){
    	return(data.frame(pos[(n.flank+1)], pattern))
    }
  }
}

##################################################################
# function: search for significant peaks (in pairwise comparisons)
##################################################################

get_peakPattern <- function(p, threshold, IDs, cores, n.flank){
  
  # specifiy all possible pairwise combinations:
  comb <- combn(seq(1:length(IDs)),2)
  comb_full <- cbind(comb, comb[c(2,1),])
  comb_full <- comb_full[,order(comb_full[1,], comb_full[2,])]
  comb_full <- comb_full[,-which(comb_full[1,] == length(IDs))]  

  # emtpy pattern vector:
  pattern_empty <- data.frame(t(rep(0,dim(comb_full)[2])))
  colnames(pattern_empty) <- apply(comb_full,2, function(x) paste(x,collapse = ","))
  
  # define positions that are significant between any pair:
  p.values <- lapply(p, function(x) x$p.values)
  p.values.mat <- do.call(cbind, p.values)
  idx.significant <- which((rowMins(p.values.mat) <= threshold) == T)

  # stop script if there is no significant difference:
  if (length(idx.significant) == 0) {
    cat(paste0(skip(), "no significant positions found for threshold ", threshold, ".\n"))
    q();
  }
  
  # generate list of p value pattern for each bin:
  list <- mclapply(idx.significant,  mc.cores = cores, 
                   function(i) get_positionPattern(	p,
							i,
							threshold,
							pattern_empty,
							comb,
							n.flank)
			)
  
  # remove zeros:
  list <- list[lapply(list,length)>0]

  null.idx <- which(sapply(list, is.null))
  if (length(null.idx) > 0) list <- list[-null.idx]
  
  # create matrix from list:
  mat <- do.call(rbind, list)
  rownames(mat) = mat[,1]
  mat <- mat[, -1]

  return(mat)
}


##################################################################
# function: combine regions
##################################################################

combine_peaks <- function(x, flank){
  
  if (length(x) > 1) {
    
    x.flank = GRanges(	seqnames=seqnames(x),
			ranges=IRanges( start=(start(x) - flank),
					end=(end(x) + flank)
					)
			)
    
    red <- reduce(x.flank,with.revmap = T)
    rm(x.flank)

    best.idx <- unlist(lapply(mcols(red)$revmap, function(y) y[which.min(mcols(x)$score[y])]))
    
    red$score <- x$score[best.idx]
    red$index <- x$index[best.idx]
    red$pattern <- x$pattern[best.idx]
    mcols(red) <- mcols(red)[,-1]
    
    start(red) <- start(red) + flank 
    end(red) <- end(red) - flank
    
    return(red)
  }else{
    return(x)
  }
}

##################################################################
# function: combine peaks by significance pattern
##################################################################

get_combinedDiffPeaks <- function(probs, p, pattern, IDs, len, n.flank){
 
  # reduce GRanges()
  idx <- as.numeric(rownames(pattern))
  gr <- probs[idx,]

  # extend by flanking positions (because flanking positions are also significant):
  w <- width(gr[1])
  start(gr) <- start(gr) - ( n.flank*w )
  end(gr) <- end(gr) + ( n.flank*w )
  
  # add minimim pvalue, bin index and pattern to each bin:
  mcols(gr)$score <- rowMins(abs(do.call(cbind, lapply(p, function(x) x$p.values[idx]))))
  mcols(gr)$index <- idx
  mcols(gr)$pattern <- apply(pattern, 1, function(x) paste0(x, collapse = ""))
  
  # split by p.value pattern
  peaks.split = split(gr, mcols(gr)$pattern)
  
  # remove gr:
  rm(gr)

  peaks <- mclapply(	peaks.split, 
			function(x) combine_peaks(x, flank = len), 
			mc.cores = cores
  			)
  
  rm(peaks.split)

  peaks <- c(do.call("c", unname(peaks)))
  peaks <- combine_peaks(peaks, flank = 0)

  colnames(elementMetadata(peaks))[1:3] <- c(   "best.p.value",
						 "index", 
	 					 "significance.pattern"
						)

  # save enhancer probability value from pos with smalles pvalue
  elementMetadata(peaks)[, unlist(IDs[-length(IDs)])] <-  mcols(probs[mcols(peaks)$index])

  return(peaks)
}


##################################################################
# function: creates a 'super'-pattern, length = length(conditions)
# (1=condition one is significant in any pair)
##################################################################

get_super.pattern <- function(pattern, len.IDs){
	res=c()
	for(i in 1:len.IDs){
		res=c(res, sum(sum(as.numeric(unlist(strsplit(pattern,"")))[((len.IDs-1)*(i-1)+1):((len.IDs-1)*(i-1)+(len.IDs-1))])>0) )
	}
	return(paste0(res, collapse=''))
}

##################################################################
# function: sort differential region by pattern
##################################################################


get_sorted_peaks <- function(peaks, IDs){

	# add new columns to object:
	mcols(peaks)$super.pattern = NA
	mcols(peaks)$cluster = NA

	# get super pattern and transform to matrix:
	pattern <- mcols(peaks)$significance.pattern
	pattern.mat <- do.call(rbind, lapply(pattern, function(x) as.numeric(unlist(strsplit(x, "")))))

	# specifiy all possible pairwise combinations:
	comb <- combn(seq(1:length(IDs)),2)
	comb_full <- cbind(comb, comb[c(2,1),])
	comb_full <- comb_full[,order(comb_full[1,], comb_full[2,])]
	comb_full <- comb_full[,-which(comb_full[1,] == length(IDs))]  

	# get ubiquitous pattern:
	pattern.idx.ubi <- which(comb_full[2,]==length(IDs))
	idx.ubi <- which(rowSums(pattern.mat[, pattern.idx.ubi])==length(pattern.idx.ubi))
	mcols(peaks)$cluster[idx.ubi] = "U"
	mcols(peaks)$super.pattern[idx.ubi] = rep(paste(rep(1,length(pattern.idx.ubi)),collapse=""), length(idx.ubi))

	# define pattern without the ubiquitous regions and define zero lines:
	pattern.mat.r <- pattern.mat[,-pattern.idx.ubi]
	pattern.r <- apply(pattern.mat.r, 1, function(x) paste(x, collapse=""))
	idx.zero <- which(apply(pattern.mat.r, 1, function(x) sum(x) == 0)==T)

	# generate super pattern
	super.pattern.r <- unlist(lapply(pattern.r, function(x) get_super.pattern(x, (length(IDs)-1))))
	super.pattern.r.mat <- do.call(rbind, lapply(super.pattern.r, function(x) as.numeric(unlist(strsplit(x, "")))))

	super.cluster <- c(which(rowSums(super.pattern.r.mat)==1), which(super.pattern.r == paste0(paste0(rep(1, (length(IDs)-1)), collapse=""), 0)))
	super.cluster.table <- sort(table(super.pattern.r[super.cluster]), decreasing=T)

	# assign super pattern:	
	for(this in unique(super.pattern.r[super.cluster])){
		mcols(peaks)$super.pattern[which(super.pattern.r==this)] = this
		mcols(peaks)$cluster[which(super.pattern.r==this)] = which(names(super.cluster.table) == this)
	}

	# assign all the other cluster names (more complex clusters, where more than one condition is significantly different):
	for(this in unique(pattern[-c(idx.zero,idx.ubi,super.cluster)])){
		if(grepl("1",this)){
			mcols(peaks)$cluster[which(pattern==this)] = paste0("r",which(names(sort(table(pattern[-c(idx.zero,idx.ubi,super.cluster)]), decreasing=T)) == this))
			mcols(peaks)$super.pattern[which(pattern==this)] = rep(get_super.pattern(this,  (length(IDs)-1)), length(which(pattern==this)))
		}
	}

	# first: order by super pattern:
	order <- with(data.frame(V0=mcols(peaks)$cluster, V1=mcols(peaks)$super.pattern), order(V1, V0)) 
	peaks=peaks[order]

	# remove NA:
	if(length(which(is.na(mcols(peaks)$cluster))) > 0){
		peaks <- peaks[-which(is.na(mcols(peaks)$cluster))]
	}

	return(peaks)
}

##################################################################
# function: adjust y positioning for heatmap (y axis):
##################################################################

adjust <- function(v){
  
  mid = v[1]
  res = c(mid)
  for (i in 2:length(v)) {
    mid = (mid + v[i - 1]/2 + v[i]/2)
    res = c(res,mid)
  }
  return(res)
}

##################################################################
# function: plot heatmap of all differential enhancer regions
##################################################################

plot_heatmap <- function(mat, IDs, color_low, color_mid, color_high, x_axis_label, y_axis_label, legend_label, out.file){
  
  # define width and height for pdf:
  width <- length((IDs))*2
  height <- 10

  # define colors and labels for cluster:
  cluster <- grep("r",mat$GRID.X)
  cluster.unique <- unique(mat$GRID.X)
  if(length(cluster != 0)){
	cluster.unique <- unique(mat$GRID.X[-grep("r",mat$GRID.X)])
  }

  my_colors <- data.frame( col=sample(rainbow(length(cluster.unique), s=0.3),length(cluster.unique)),
			   cluster=cluster.unique,
			   stringsAsFactors=F)

  colors<- rep("NA", length(which(mat$GRID.Y == "cluster")))
  label <- rep("", length(which(mat$GRID.Y == "cluster")))

  for(this in my_colors$cluster){
	colors[which(mat$GRID.X[which(mat$GRID.Y == "cluster")]==this)] = my_colors$col[which(my_colors$cluster==this)]
	label[median(as.numeric(mat$Y[which(mat$GRID.X[which(mat$GRID.Y == "cluster")]==this)]))] = this
  } 

  plot <- ggplot( mat, aes(x = X,y = Y)) + 
                  geom_tile(aes(fill = FILL, height = 1), size = 0) + 
                  facet_grid( . ~ GRID.Y,
                              scales = "free",
                              space = "free"
                  ) +
                  theme( axis.ticks = element_blank(),
			 axis.title = element_text(size = 15),
			 strip.text = element_text(size = 15),
			 legend.title=element_text(size=15),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         legend.position = "bottom",
			 panel.spacing = unit(0, "lines"),
			 strip.background = element_rect(color="white", fill="white", size=0.5, linetype="solid"
     )
                  ) +
                  scale_fill_gradient2( low = color_low,
                                        mid = color_mid,
                                        high = color_high,
                                        midpoint = 0.5, 
                                        name = legend_label
                  ) +
                  scale_y_continuous( breaks = c(),
                                      expand = c(0, 0),
                                      name = y_axis_label
                  ) +
                  scale_x_discrete(	position = "top",
                                    	name = x_axis_label
                  ) +
		  geom_tile(aes(height = 1, width=0.5), fill=colors, size = 0,subset(mat,GRID.Y %in% c("cluster")),show.legend = F)+
		  geom_text(aes(label=label),subset(mat,GRID.Y %in% c("cluster")))

    ggsave(filename = paste0(out.file,".pdf"), plot = plot, width = width, height = height, device = "pdf", units = 'in')
    ggsave(filename = paste0(out.file,".png"), plot = plot, width = width, height = height, device = "png", units = 'in', dpi=300)
}


##################################################################
# function: get summarized counts in defined bins
##################################################################


get_bamOverlaps <- function(files, IDs, txdb, singleEnd){
  
  # gives information about path of bam and bai file (Object)
  bf = BamFileList(unlist(files))
  
  # gives length of different chromosomes (SeqInfo class)
  si = seqinfo(bf)
  
  # get exons and genes
  expr.gr <- genes(txdb)
  seqlevels(expr.gr) = paste0("chr", gsub("chr","",seqlevels(expr.gr)))
  
  exons <- exonsBy(txdb, by = "gene")
  exons <- exons[which((names(exons) %in% mcols(expr.gr)$gene_id) == T)]
  
  # fix chromosome prefix:
  seqlevels(exons) <- gsub("chr","",seqlevels(exons))
  if (grepl("chr",seqlevels(si)[1])) {
    seqlevels(exons) <- paste0("chr", seqlevels(exons))
  }
  
  se <- summarizeOverlaps(exons,
                          bf,
			  mode="Union",
                          singleEnd = singleEnd,
                          fragments = setdiff(c(FALSE,TRUE), singleEnd))
  
  # get counts and summarize per gene:
  counts.per.exon <- data.frame(assays(se)$counts,
                                gene = rownames(assays(se)$counts))
  counts.split <- split(counts.per.exon,
                        counts.per.exon$gene)
  counts.per.gene <- (do.call(rbind,lapply(counts.split, function(x) colSums(x[,-dim(x)[2]]))))

  # create new symmarized experiment object:
  se0 <- SummarizedExperiment( assays = SimpleList(counts = counts.per.gene),
                               colData = names(bf)
  )
  
  # stabilize the variance across the mean:
  # (The transformed data should be approximated variance stabilized
  # and also includes correction for size factors or normalization factors)

  dds <- suppressMessages(DESeqDataSet(se0, ~ 1))
  vsd <- varianceStabilizingTransformation(	dds,
						blind = TRUE)

  counts_vst <- assay(vsd)
  counts_raw <- assay(dds)
  expr_raw.gr <- expr.gr

  for (i in 1:dim(counts_vst)[2]) {
        mcols(expr.gr)[,unlist(IDs)[i]] <- counts_vst[,i]
        mcols(expr_raw.gr)[,unlist(IDs)[i]] <- counts_raw[,i]
  }
  
  expr.gr <- expr.gr[which(rowVars(counts_vst) > 0)]
  expr_raw.gr <- expr_raw.gr[which(rowVars(counts_vst) > 0)]
  
  return(list(raw = expr_raw.gr, vst = expr.gr))
}


##################################################################
# function: create TAD if not defined
##################################################################

check_TAD <- function(t, TAD, this.region, regions){
  
  if (length(t) == 0) {
    
    precedeT <- precede(TAD, this.region)
    followT  <- follow(TAD, this.region)
    
    if(length(which(!is.na(precedeT))) == 0 ){
	start <- 1
    }else{
    	start <- end(TAD[rev(which(!is.na(precedeT)))[1]])
    	if (length(start) == 0) start <- 1
    }

    if(length(which(!is.na(followT))) == 0 ){
	end <- max(end(regions[which(seqnames(regions) == seqnames(this.region))]))
    }else{
    	end <- start(TAD[(which(!is.na(followT)))[1]])
    	if (length(end) == 0) end <- max(end(regions[which(seqnames(regions) == seqnames(this.region))]))
    }

    t <- GRanges(seqnames(this.region), IRanges(start, width = (end - start + 1)))
  }
  return(t)
}

##################################################################
# function: correlate probabilities and gene exression counts 
##################################################################

get_correlation <- function(i, threshold, regions.gr, expr.gr, TAD.gr, IDs){

  interactions <- data.frame(stringsAsFactors = F)
  
  # get region, associated TAD and genes within TAD
  this.region <- regions.gr[i]
  this.TAD <- subsetByOverlaps(TAD.gr, this.region)
  this.TAD <- check_TAD(this.TAD, TAD.gr, this.region, regions.gr)
  this.genes.idx <- expr.gr %within% this.TAD
  
  if (sum(this.genes.idx) > 0) {
    
    this.genes <- expr.gr[this.genes.idx,]
    cor <- apply(as.matrix(mcols(this.genes)[,IDs]), 1, 
                 function(x) cor(x, as.numeric(unlist(mcols(this.region)[,IDs]))))
    
    for (c in 1:length(cor)) {
      if (!is.na(cor[c]) && cor[c] >= threshold) {

      promoter_start <- start(promoters(this.genes[c]))
      promoter_end <- end(promoters(this.genes[c]))

      if(as.character(strand(this.genes[c]))=="-"){
      	promoter_start <- end(promoters(this.genes[c]))
      	promoter_end <- start(promoters(this.genes[c]))
      }

      interactions <- rbind(	interactions, 
                               data.frame(	data.frame(this.region)[,c(GR_header_short, "cluster", IDs)],
                                           	TAD_COORDINATES = paste0(this.TAD),
                                           	CORRELATED_GENE = paste(mcols(this.genes)[c,"gene_id"]),
						CORRELATED_GENE_CHR = seqnames(this.genes[c]),
						CORRELATED_GENE_PROMOTER_START = promoter_start,
						CORRELATED_GENE_PROMOTER_END = promoter_end,
                                           	CORRELATION = cor[c] ))
      }
    }
    if (length(interactions) > 0) return(makeGRangesFromDataFrame(interactions, keep.extra.columns = T))
  }
}


##################################################################
# function: get nearest gene for a region
##################################################################

get_nearest_gene <- function(i, regions.gr, genes, IDs){

  this.region <- regions.gr[i]

  nearest <- genes[nearest(this.region, genes)]
  distance.to.nearest <- distance(this.region, nearest)

  promoter_start <- start(promoters(nearest))
  promoter_end <- end(promoters(nearest))

  if(as.character(strand(nearest))=="-"){
      promoter_start <- end(promoters(nearest))
      promoter_end <- start(promoters(nearest))
  }

  interactions <-	data.frame(	data.frame(this.region)[,c(GR_header_short, "cluster", IDs)],
                                       	NEAREST_GENE = mcols(nearest)[,"gene_id"],
					NEAREST_GENE_CHR = seqnames(nearest),
					NEAREST_GENE_PROMOTER_START = promoter_start,
					NEAREST_GENE_PROMOTER_END = promoter_end,
                                       	DISTANCE_TO_NEAREST = distance.to.nearest
			)

  return(makeGRangesFromDataFrame(interactions, keep.extra.columns = T))
}


##################################################################
# function: correlate probabilities with gene expression values
#(in same TAD; per cluster)
##################################################################

get_units <- function(regions.gr, expr.gr, TAD.gr, IDs, cores, threshold, txdb){
  
  if(!isEmpty(TAD.gr)){
  	# get correlation for each differential region:
 	 list <- mclapply( seq(length(regions.gr)), 
 	                   function(x) get_correlation(x, threshold, regions.gr, expr.gr, TAD.gr, IDs), 
  	                   mc.cores = cores
 	 )
  }else{

	# get nearest gene:
	genes <- genes(txdb)
	seqlevels(genes) = paste0("chr", gsub("chr|Chr","",seqlevels(genes)))
	list <- mclapply(  seq(length(regions.gr)), 
 	                   function(x) get_nearest_gene(x, regions.gr, genes, IDs), 
  	                   mc.cores = cores
 	 )
  }
  units <-  do.call("c", unname(unlist(list)))  
  return(units)
}


