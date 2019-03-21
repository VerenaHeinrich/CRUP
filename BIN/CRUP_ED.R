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

# colors for heatmap spectrum:
color_low  <- "white"
color_mid  <- "#84B3D7"
color_high <- "#003964"

# axis and legend labels:
legend_label <- "Probabilities"
y_axis_label <- "Differential Enhancer"
x_axis_label <- ""

# prefixes:
condition_prefix  <- "Condition "
ID_prefix         <- "cond"
sample_prefix     <- "Rep"

##################################################################
# check input parameter
##################################################################

parameter <- c("probabilities")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))
parameter.extra <- c("outdir", "cores", "w_0", "len", "threshold", "names")

if (!is.null(opt$help)) stop_script(c("dynamics", parameter, parameter.extra))
if (sum(parameter.bool) == 0) stop_script(c("dynamics", parameter, parameter.extra))

# check probability files
files <- as.list(unlist(strsplit(opt$probabilities,',')))
files <- lapply(files, function(x) unlist(strsplit(x,':')))
for (f in unlist(files)) check_file(f)

# assign IDs: 
if (is.null(opt$names)) {
  IDs <- lapply(1:length(files), function(x) paste0(ID_prefix, x, "_",seq(1:length(files[[x]]))))
}else{
  labels <- unlist(strsplit(opt$names,','))

  if (length(labels) != length(files)) {
    cat("\n\tWARNING: Number of alternative condition names ('n') is not valid.
        Names are set to default labels.\n");
    IDs <- lapply(1:length(files), function(x) paste0(ID_prefix, x, "_", seq(1:length(files[[x]]))))
  } else{
    IDs <- files
    for (i in 1:length(IDs)) {
      for (j in 1:length(IDs[[i]])) {
        IDs[[i]][j] <- paste0(labels[i], "_", j)
      }
    }
  }
}

# check output directory
opt$outdir <- check_outdir(opt$outdir, unlist(files)[1])

# set default values
if (is.null(opt$cores))       opt$cores <- 1
if (is.null(opt$len))         opt$len <- 1000
if (is.null(opt$threshold)) {
  opt$threshold <- 0.01 
} else if (!is.null(opt$threshold) & (opt$threshold < 0 | opt$threshold > 1) ) {
  cat(paste0(opt$threshold, " is not in range [0,1].\n"))
  q();
}
if (is.null(opt$w_0)) {
  opt$w_0 <- 0.5
} else if (!is.null(opt$w_0) & (opt$w_0 < 0 | opt$w_0 > 1) ) {
  cat(paste0(opt$w_0, " is not in range [0,1].\n"))
  q();
}

##################################################################
# define input parameter
##################################################################

startPart("List input parameter")

files     <- lapply(files, function(x) normalizePath(x))
outdir    <- paste0(normalizePath(opt$outdir),"/")
cores     <- opt$cores
w_0       <- opt$w_0
len       <- opt$len
threshold <- opt$threshold

cat(skip(), "files: \n", 
            paste0("\t\t",unlist(lapply(1:length(files), 
                                        function(x) paste(IDs[[x]], files[[x]], sep = "\t-> "))), "\n"))
cat(skip(), "w_0: ",w_0, "\n")
cat(skip(), "len: ",len, "\n")
cat(skip(), "threshold: ",threshold, "\n")
cat(skip(), "cores: ",cores, "\n")
cat(skip(), "outdir: ",outdir, "\n")

endPart()

##################################################################
# libraries
##################################################################

startPart("Load packages")

pkgLoad("ggplot2")        # for ggplot()
pkgLoad("GenomicRanges")  # for GRanges()
pkgLoad("dplyr")          # for %>%
pkgLoad("matrixStats")    # for rowVars()
pkgLoad("parallel")       # for mclapply()

endPart()

##################################################################
# read and combine probabilities:
##################################################################

startPart("Read enhancer probabilites for all samples and conditions")
probs <- get_probabiltyMatrix(files, IDs)
IDs[[length(IDs)+1]] <- 'background'

endPart()

##################################################################
# get empricial p-values
# (list with empirical p-values per pair comparison and
# direction of group mean differences)
##################################################################

startPart("Calculate (pairwise) empirical p-values")
p <- get_pairwisePvalues(probs, IDs, w_0, cores)
endPart()

##################################################################
# Call dynamically changing enhancer peaks:
##################################################################

startPart("Get condition-specific enhancer peaks")
cat(paste0(skip(), "build significance peak pattern"))
pattern <- get_peakPattern(p, threshold, IDs, cores)
done()

if (is.null(pattern)) {
  cat(skip(), "no significant peaks found between any two conditions.\n")
  q()
}

cat(paste0(skip(), "combine peaks by significance pattern"))
peaks <- get_combinedDiffPeaks(probs, p, pattern, IDs, len)
done()

cat(paste0(skip(), "sort regions by significance peak pattern"))
peaks <- get_sorted_peaks(peaks, IDs)
done()


##################################################################
# OUTPUT:
##################################################################

# output 1: a summarized text file:
out.txt <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".txt"))
cat(paste0(skip(), "Output 1 - save condition-specific enhancer regions to:  ", out.txt))
write.table(data.frame(peaks)[, c(	 GR_header,
					 "best.p.value",
 					 "cluster",
					 unlist(IDs[-length(IDs)]))],
				 file = out.txt, quote = F, row.names = F, sep = "\t")
done()

# output 2: all dynamically changing enhancers in .bed format:
out.bed <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".bed"))
cat(paste0(skip(), "Output 2 - a reduced version is saved in .bed file format:  ", out.bed))
write.table(data.frame(peaks)[, GR_header_short], file = out.bed, quote = F, row.names = F, col.names = F, sep = "\t")
done()

# output 3: main clusters in .bed format:
cat(paste0(skip(), "Output 3 - main clusters (just one condition is active) in separated .bed files."))
for(this in unique(mcols(peaks)$cluster)){
	if(! grepl("r", this)){
		out.bed <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold,".cluster_",this,".bed"))
		write.table(data.frame(peaks)[which(mcols(peaks)$cluster == this), GR_header_short], file = out.bed, quote = F, row.names = F, col.names = F, sep = "\t")
	}
}
done()

# output 4: a visual summary of all dynamically changing enhancers:
IDs=IDs[-length(IDs)]

LABEL_COND <- c(gsub("_.*","", unlist(IDs)), "cluster")
LABEL_REP <- c(gsub(".*_", "", unlist(IDs)), "")
CLUSTER <- mcols(peaks)$cluster
VALUES  <- cbind(elementMetadata(peaks)[, unlist(IDs)],rep(0,length(peaks)))

mat <- data.frame(  X = unlist(lapply(LABEL_REP, function(x) rep(x, length(peaks)))),
                    Y = rep(seq(length(peaks)), (length(unlist(IDs))+1)), 
                    FILL = matrix(unlist(VALUES),  ncol = 1), 
                    GRID.X = factor(rep(CLUSTER, (length(unlist(IDs))+1)), levels = unique(CLUSTER)),
                    GRID.Y = factor(unlist(lapply(LABEL_COND, function(x) rep(x,  length(peaks)))), levels=unique(LABEL_COND))
)

mat<-mat[which(mat$GRID.X !="U"),]

out <- paste0(outdir, paste0("dynamicEnh__w0_",w_0,"__threshold_", threshold))
cat(paste0(skip(), "Output 4 - results are visualized as a heatmap (pdf and png):  ", out,".pdf/png"))

plot_heatmap( mat,
              IDs,
              color_low,
              color_mid,
              color_high,
              x_axis_label,
              y_axis_label,
              legend_label,
              out
)
done()

endPart()

##################################################################
# print run time
##################################################################

run_time <- Sys.time() - start_time

startPart("Run time")
cat(paste0(skip(), format(run_time), "\n"))
endPart()



