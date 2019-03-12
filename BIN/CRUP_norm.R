##################################################################
# start run time
##################################################################

start_time <- Sys.time()

##################################################################
# get input paramter
##################################################################

if (!exists('opt')) {
  
  # path to functions
  path <- dirname(getwd())
  source(paste0(path,'/BIN/functions.R'))
  
  # library - getopt
  pkgTest("getopt")
  suppressMessages(library(getopt))

  # get input parameter
  opt <- getopt(spectrum)
}

##################################################################
# check input parameter
##################################################################

parameter <- c("file", "sequencing", "genome")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))
parameter.extra <- c("outdir", "cores", "mapq")

if (!is.null(opt$help)) stop_script(c("norm", parameter, parameter.extra))
if (sum(parameter.bool) == 0) stop_script(c("norm", parameter, parameter.extra))

# check sequencing parameter
if (!(opt$sequencing %in% sequencing_values)) {
  cat("Sequencing parameter is not valid.\n
       Choose one of: ", paste(sequencing_values,","),"\n");
  q();
}

# check summary file
check_file(opt$file)

# check genome
if (!is.null(opt$genome)) { 
  if (!(opt$genome %in% genome_values)) { 
    cat("Genome" , opt$genome, " currently not provided. Choose one of: ",paste(genome_values,","),"\n\n");
    q();}}

# check output directory
opt$outdir <- check_outdir(opt$outdir, opt$file)

# set default values
if (is.null(opt$mapq))   opt$mapq <- 10
if (is.null(opt$cores))  opt$cores <- 1

##################################################################
# define input parameter
##################################################################

startPart("List input parameter")

file        <- normalizePath(opt$file)
genome      <- opt$genome
mapq        <- opt$mapq
sequencing  <- opt$sequencing
cores       <- opt$cores
outdir      <- paste0(normalizePath(opt$outdir),"/")

cat(skip(), "file: ",file, "\n")
cat(skip(), "mapq: ",mapq, "\n")
cat(skip(), "genome: ",genome, "\n")
cat(skip(), "sequencing: ",sequencing, "\n")
cat(skip(), "cores: ",cores, "\n")
cat(skip(), "outdir: ",outdir, "\n")

endPart()

##################################################################
# libraries
##################################################################

startPart("Load packages")

if (genome == "mm10") pkg <- "BSgenome.Mmusculus.UCSC.mm10"
if (genome == "mm9") pkg  <- "BSgenome.Mmusculus.UCSC.mm9"
if (genome == "hg19") pkg <- "BSgenome.Hsapiens.UCSC.hg19"
if (genome == "hg38") pkg <- "BSgenome.Hsapiens.UCSC.hg38"

pkgLoad(pkg)                            # for genome object
assign("txdb", eval(parse(text = pkg)))

pkgLoad("bamsignals")                   # for bamProfile()
pkgLoad("Rsamtools")                    # for scanBamHeader()

endPart()

##################################################################
# get summarized information of ChiP-seq experiments
##################################################################

startPart("Get information of ChIP-seq experiments")

# load info file
info <- read.csv(file	<- file,
                 header	<- TRUE,
                 sep	<- "\t"
)

# check info file:
header.valid <- c("feature", "bam_file", "bam_file_input")
if (!identical(colnames(info), header.valid)) {
  cat(paste0("Header in summary text file is not correct.\n
              Header must be tab separated and contain: ",header.valid,"\n")
      ); 
  q();
} else{
  cat(paste0(skip(),"header is valid.\n"))
}

# get and check feature set
features.valid <-  c("H3K27ac","H3K4me1", "H3K4me3")
if (! all(as.character(info$feature) %in% features.valid)) {
  cat(paste0("Features in the summary text file are not correct.\n
              Possible features: ", features.valid,"\n")
      ); 
  q();
} else{
  cat(paste0(skip(),"feature set is valid.\n"))
  features.valid <- as.character(info$feature)
}

# check if all paths in file exist
bam_files_HM <- as.character(info$bam_file)
bam_files_input <- as.character(info$bam_file_input)

for (this in unique(c(bam_files_HM, bam_files_input))) {
  check_file(this)
}
cat(paste0(skip(),"all files in summary file exist:\n"))

endPart()

##################################################################
# get and adjust binned genome
##################################################################

startPart("Prepare the binned genome")

# bin the genome
cat(paste0(skip(), "get binned genome for ", genome))
binned_genome <- get_binned_genome(txdb,tilewidth <- 100)
done()

# get prefix of chromosome names
cat(paste0(skip(), "get prefix of chromosomes from  file ", bam_files_HM[1]))
bam_header <-  scanBamHeader(bam_files_HM[1])[[1]]$targets
prefix <- gsub("chr", "",binned_genome@seqnames)
done()
 
# check if prefix of chromosome names ("chr1" or "1")
cat(paste0(skip(), "adjust prefix of chromosome names in binned genome"))
if (!("chr1" %in% names(bam_header)))  binned_genome@seqnames = prefix
done()

endPart()

##################################################################
# get counts from ChIP-seq experiments
##################################################################

startPart("Get summarized counts from ChIP-seq experiments")

# prepare feature names:
names = c( features.valid, paste0("Input_",features.valid))

if (length(unique(bam_files_input)) == 1) {
  bam_files_input <- bam_files_input[1]
  names = c( features.valid, "Input_All")
} 

# get counts for ChIP-seq HM experiments
counts <- mclapply( c(bam_files_HM, bam_files_input),
                    get_bamProfile,
                    gr = binned_genome,
                    mapq = mapq,
                    sequencing = sequencing,
                    mc.cores = cores
              )
names(counts) <- names

endPart()

##################################################################
# get counts from ChIP-seq experiments
##################################################################

startPart("Normalize histone modifications by Input")

if ("Input_All" %in% names(counts)) {
  normalized_counts <- lapply( features.valid, function(x) log2((counts[[x]] + 1)/(counts[[paste0("Input_All")]] + 1)))
}else {
  normalized_counts <- lapply( features.valid, function(x) log2((counts[[x]] + 1)/(counts[[paste0("Input_",x)]] + 1)))
}

endPart()

##################################################################
# create data matrix for all normalized ChIP-seq experiments
##################################################################

startPart("Create summarized data matrix")

cat(paste0(skip(), "create matrix"))
mcols(binned_genome) <-  matrix( unlist(normalized_counts), 
                                 ncol = length(info$feature),
                                 byrow = FALSE,
                                 dimnames = list(NULL, features.valid))
seqlevels(binned_genome) = paste0('chr', gsub('chr|Chr','',seqlevels(binned_genome)))
done()

# include log2 H3K4me1/H3K4me3 ratio
cat(paste0(skip(), "include log2 H3K4me1/H3K4me3 ratio"))
nominator   <- mcols(binned_genome)[,"H3K4me1"] + abs(min(mcols(binned_genome)[,"H3K4me1"])) + 1
denominator <- mcols(binned_genome)[,"H3K4me3"] + abs(min(mcols(binned_genome)[,"H3K4me3"])) + 1
mcols(binned_genome)[,"ratio"] <- log2(nominator/denominator)
done()

# save data matrix
out <- paste0(outdir, paste0(gsub(".txt|.info","",basename(file)), ".data_matrix.rds"))
cat(paste0(skip(), "save data matrix as: ", out))
saveRDS(binned_genome, gsub("/{1,}","/",out))
done()

endPart()

##################################################################
# print run time
##################################################################

run_time <- Sys.time() - start_time

startPart("Run time")
cat(paste0(skip(), format(run_time), "\n"))
endPart()
