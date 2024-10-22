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
# define fixed parameters;
##################################################################

ID_prefix         <- "cond"

##################################################################
# check input parameter
##################################################################

parameter <- c("regions", "RNA", "expression", "genome", "sequencing")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))
parameter.extra <- c("outdir", "cores", "TAD", "label", "threshold_c", "nearest")

if (!is.null(opt$help)) stop_script(c("targets", parameter, parameter.extra))
if (is.null(opt$regions)) stop_script(c("targets", parameter, parameter.extra))
if (is.null(opt$RNA) & is.null(opt$expression)) stop_script(c("targets", parameter, parameter.extra))
if (!is.null(opt$RNA) & (is.null(opt$genome) | is.null(opt$sequencing))) stop_script(c("targets", parameter, parameter.extra))
if (!is.null(opt$RNA) & !is.null(opt$expression)) {
  cat(skip(), "normalized expression counts are provided. Ignore RNA-seq files.\n")
  opt$RNA <- NULL
}

# check genome
if (!is.null(opt$genome)) { 
 if (grepl('.fai',opt$genome)) { 
    cat(skip(),"Alternative genomes currently not provided for CRUP-ET. Choose one of: ",paste0(genome_values),"\n\n");
    q();
  } 
 if (!(opt$genome %in% genome_values)) { 
    cat(skip(),"Genome" , opt$genome, " currently not provided. Choose one of: ",paste0(genome_values),"\n\n");
    q();
  }
}

# check input file
check_file(opt$regions)

# check gene expression counts file
if (!is.null(opt$expression)) {
  check_file(opt$expression) 
} else if (!is.null(opt$RNA)) {

  # check RNA-seq bam files:
  rna <- as.list(unlist(strsplit(opt$RNA,',')))
  rna <- lapply(rna, function(x) unlist(strsplit(x,':')))
  for (f in unlist(rna)) check_file(f)
  
  # assign IDs: 
  if (is.null(opt$label)) {
    IDs <- lapply(1:length(rna), function(x) paste0(ID_prefix, x, "_",seq(1:length(rna[[x]]))))
  }else{
    labels <- unlist(strsplit(opt$label,','))
    
    if (length(labels) != length(rna)) {
      cat("\n\tWARNING: Number of alternative condition labels ('l') is not valid.
        Names are set to default labels.\n");
        IDs <- lapply(1:length(rna), function(x) paste0(ID_prefix, x, "_", seq(1:length(rna[[x]]))))
    } else{
      IDs <- rna
      for (i in 1:length(IDs)) {
        for (j in 1:length(IDs[[i]])) {
          IDs[[i]][j] <- paste0(labels[i], "_", j)
        }
      }
    }
  }
}

# check output directory
opt$outdir <- check_outdir(opt$outdir, opt$regions)

# set default values
if (is.null(opt$threshold_c)) {
  opt$threshold_c <- 0.9
} else if (!is.null(opt$threshold_c) & (opt$threshold_c < 0.5 | opt$threshold_c > 1) ) {
  cat(paste0(opt$threshold_c, " is not in range [0.5,1].\n"))
  q();
}
if (is.null(opt$cores)) opt$cores <- 1

if (is.null(opt$nearest)){
	if (is.null(opt$TAD) & is.null(opt$nearest) & opt$genome == "mm10") {
	  opt$TAD   <- normalizePath(paste0(path,"/DATA/mESC_mapq30_KR_all_TADs.bed"))
	}
	if (is.null(opt$TAD) & is.null(opt$nearest))  {
	  cat(paste0("You have to provide your own file with TAD domains (fitting to the genome choice).\n"))
	  q();
	}
}

##################################################################
# define input parameter
##################################################################

startPart("List input parameter")

outdir      <- paste0(normalizePath(opt$outdir),"/")
cores       <- opt$cores
regions     <- paste0(normalizePath(opt$regions))
if (is.null(opt$nearest)) TAD <- opt$TAD
threshold_c <- opt$threshold_c

if (!is.null(opt$expression)) {
  expression  <- normalizePath(opt$expression)
  cat(skip(), "expression: ",expression, "\n")
}else{
  rna     <- lapply(rna, function(x) normalizePath(x))
  
  genome      <- opt$genome
  mapq        <- opt$mapq
  sequencing  <- opt$sequencing
  
  cat(skip(), "RNA: \n", 
      paste0("\t\t",unlist(lapply(1:length(rna), 
                                  function(x) paste(IDs[[x]], rna[[x]], sep = "\t-> "))), "\n"))
  
  cat(skip(), "sequencing: ",sequencing, "\n")
  cat(skip(), "genome: ",genome, "\n")
}

cat(skip(), "regions: ",regions, "\n")
if (is.null(opt$nearest)){
	cat(skip(), "threshold_c: ",threshold_c, "\n")
	cat(skip(), "TAD: ",TAD, "\n")
}
if (!is.null(opt$nearest)){
	genome      <- opt$genome
	cat(paste0("Will choose the nearest gene to a each differential region to build a regulatory unit.\n"))
}

cat(skip(), "cores: ",cores, "\n")
cat(skip(), "outdir: ",outdir, "\n")

endPart()

##################################################################
# libraries
##################################################################

startPart("Load packages")

pkgLoad("ggplot2")        # for ggplot()
pkgLoad("GenomicRanges")  # for GRanges()

txdb <- ""
if (is.null(opt$expression) | !is.null(opt$nearest)) {
  
  if (genome == "mm9") pkg <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
  if (genome == "mm10") pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
  if (genome == "hg19") pkg <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
  if (genome == "hg38") pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"

  pkgLoad(pkg)                  # for genome object
  assign("txdb", eval(parse(text = pkg)))
  
  pkgLoad("GenomicFeatures")    # for exonsBy()

  if (is.null(opt$expression))  {
  	pkgLoad("GenomicAlignments")  # for summarizeOverlaps()
  	pkgLoad("Rsamtools")          # for seqinfo()
  	pkgLoad("DESeq2")             # for varianceStabilizingTransformation()
  }
}

endPart()

##################################################################
# get normalized gene expression values
##################################################################

if (is.null(opt$expression)) {
  startPart("Get normalized gene expression values from RNA-seq experiments")

  # get summarized and normalized counts
  cat(paste0(skip(), "get normalized gene expression counts"))
  list <- get_bamOverlaps( unlist(rna),
                              IDs,
                              txdb,
                              singleEnd = (sequencing == "single"))
  expr.gr <- list[['vst']]
  expr_raw.gr <- list[['raw']]
  done()

  out_raw.rds <- paste0(outdir, paste0("gene_expression_raw.rds"))
  cat(paste0(skip(), "save raw gene expression counts to:  ", out_raw.rds))
  saveRDS(expr_raw.gr, out_raw.rds)
  done()
  
  out_vst.rds <- paste0(outdir, paste0("gene_expression_vst.rds"))
  cat(paste0(skip(), "save normalized gene expression counts to:  ", out_vst.rds))
  saveRDS(expr.gr, out_vst.rds)
  done()
  
  # get IDs:
  IDs <- unlist(IDs)
  
  endPart()
} else{
  startPart("Get normalized gene expression values from saved rds file")
  expr.gr <- readRDS(expression, T)
  IDs <- colnames(mcols(expr.gr)[,-1])
  endPart()
}

##################################################################
# get correlation:
##################################################################

startPart("Build Regulatory Units")

# dynamic enhancer regions:
regions.gr <- makeGRangesFromDataFrame(read.table(regions, header = TRUE), keep.extra.columns = T)
seqlevels(regions.gr) = paste0("chr", gsub("chr|Chr","",seqlevels(regions.gr)))
cluster.U <- which(mcols(regions.gr)$cluster == 'U')
if(length(cluster.U) > 0) {regions.gr <- regions.gr[-cluster.U]}

# TAD/domain regions:
TAD.gr <- GRanges()
if (is.null(opt$nearest)){
	TAD.df <- read.table(TAD, col.names = GR_header_short)
	TAD.df <- TAD.df[which((TAD.df$end-TAD.df$start) > 0),]
	TAD.gr <- makeGRangesFromDataFrame(TAD.df)
	seqlevels(TAD.gr) = paste0("chr", gsub("chr|Chr","",seqlevels(TAD.gr)))
}

units <- get_units(	regions.gr,
			expr.gr,
			TAD.gr,
			IDs,
			cores,
			threshold_c,
			txdb
			)

out.txt <- ""
if (is.null(opt$nearest)){
	out.txt <- paste0(outdir, paste0("RegulatoryUnits_CorrTres_", threshold_c,".txt"))
}else{
	out.txt <- paste0(outdir, paste0("RegulatoryUnitsNearestGene.txt"))
}
cat(paste0(skip(), "save dynamic enhancer gene interactions in txt:  ", out.txt, "\n"))
write.table(units, file = out.txt, quote = F, row.names = F, sep = "\t")
done()

out.interaction <- ""
if (is.null(opt$nearest)){
	out.interaction <- paste0(outdir, paste0("RegulatoryUnits_CorrTres_", threshold_c,".interaction"))
}else{
	out.interaction <- paste0(outdir, paste0("RegulatoryUnitsNearestGene.interaction"))
}
cat(paste0(skip(), "save dynamic enhancer gene interactions in UCSC interaction format:  ", out.interaction, "\n"))

header <- "track type=interact name=\"Dynamic Promoter-Enhancer Pairs\" description=\"Dynamic Promoter-Enhancer Pairs\" interactDirectional=true visibility=full"
enhancer <- as.matrix(data.frame(units)[,c("start","end")])

promoter <-c()
score <- 0
if (is.null(opt$nearest)){
	promoter <- as.matrix(data.frame(mcols(units)$CORRELATED_GENE_PROMOTER_START, mcols(units)$CORRELATED_GENE_PROMOTER_END))
	score <- mcols(units)$CORRELATION
}else{
	promoter <- as.matrix(data.frame(mcols(units)$NEAREST_GENE_PROMOTER_START, mcols(units)$NEAREST_GENE_PROMOTER_END))
	score <- mcols(units)$DISTANCE_TO_NEAREST
}

interaction <- cbind( as.character(seqnames(units)),
                      apply(cbind(enhancer,promoter),1,min),
                      apply(cbind(enhancer,promoter),1,max),
                      rep(".", length(units)),
                      rep(0, length(units)),
                      score,
                      rep(".", length(units)),
                      rep("#7A67EE", length(units)),
                      as.character(seqnames(units)),
                      enhancer,
                      rep(".", length(units)),
                      rep(".", length(units)),
                      as.character(seqnames(units)),
                      promoter,
                      rep(".", length(units)),
                      rep(".", length(units)) )
write.table(header, file = out.interaction, quote = F, col.names = F, row.names = F, sep = "\t")
write.table(interaction, file = out.interaction, quote = F, row.names = F, col.names = F, sep = "\t", append = T)
done()

endPart()

##################################################################
# print run time
##################################################################

run_time <- Sys.time() - start_time

startPart("Run time")
cat(paste0(skip(), format(run_time), "\n"))
endPart()


