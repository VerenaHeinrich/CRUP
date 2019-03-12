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

parameter <- c("matrix")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))
parameter.extra <- c("outdir", "cores", "classifier", "cutoff", "distance")

if (!is.null(opt$help)) stop_script(c("prediction", parameter, parameter.extra))
if (sum(parameter.bool) == 0) stop_script(c("prediction", parameter, parameter.extra))

# check data matrix file
check_file(opt$matrix)

# check classifier directory
if (is.null(opt$classifier)) opt$classifier <- paste0(path, "/DATA/CLASSIFIER/")
if (!dir.exists(opt$classifier)) {
  cat(paste0("Classifier directory '",opt$classifier,"'is not a valid directory. \n
              Classifier directory is set to ",paste0(path, "/DATA/CLASSIFIER/")));
  opt$classifier <- paste0(path, "/DATA/CLASSIFIER/")
}

# check output directory
opt$outdir <- check_outdir(opt$outdir, opt$matrix)

# set default values
if (is.null(opt$cores))     opt$cores <- 1
if (is.null(opt$distance))  opt$distance <- 12500
if (is.null(opt$cutoff)) {
  opt$cutoff <- 0.5
}else if (!is.null(opt$cutoff) & (opt$cutoff < 0 | opt$cutoff > 1) ) {
  cat(paste0(opt$cutoff, " is not in range [0,1].\n"))
  q();
}


##################################################################
# libraries
##################################################################

startPart("Load packages")

pkgLoad("preprocessCore") # for normalize.quantiles.use.target()
pkgLoad("randomForest")   # for predict()
pkgLoad("GenomicRanges")  # for GRanges()
pkgLoad("rtracklayer")    # for export()

endPart()


##################################################################
# define input parameter
##################################################################

startPart("List input parameter")

matrix      <- normalizePath(opt$matrix)
classifier  <- paste0(normalizePath(opt$classifier),"/")
outdir      <- paste0(normalizePath(opt$outdir),"/")
distance    <- opt$distance
cores       <- opt$cores
cutoff      <- opt$cutoff

cat(skip(), "matrix: ",matrix, "\n")
cat(skip(), "classifier: ",classifier, "\n")
cat(skip(), "distance: ",distance, "\n")
cat(skip(), "cutoff: ",cutoff, "\n")
cat(skip(), "cores: ",cores, "\n")
cat(skip(), "outdir: ",outdir, "\n")

endPart()

##################################################################
# Read classifier files and ecdf
##################################################################

startPart("Get classifier and empirical distribution function")

classifier_file1 <- file.path(classifier, "active_vs_inactive.rds")
check_file(classifier_file1)
classifier1 <- readRDS(classifier_file1)
features1 <- unique(gsub('_.*','',names(classifier1$forest$xlevels))) 

classifier_file2 <- file.path(classifier, "enhancer_vs_active_promoter.rds")
check_file(classifier_file2)
classifier2 <- readRDS(classifier_file2)
features2 <- unique(gsub('_.*','',names(classifier2$forest$xlevels))) 

features_all <- unique(c(features1, features2)) 

# ecdf: will be used for quantile normalization
ecdf_file <- paste0(path, "/DATA/ecdf.rds")
check_file(ecdf_file)
ecdf <- readRDS(ecdf_file)

endPart()

##################################################################
# Quantile normalization
##################################################################

startPart("Quantile normalize summarized data matrix")

# load data file
data_matrix <- makeGRangesFromDataFrame( suppressMessages(readRDS(matrix)), keep.extra.columns=T)

# normalize for each feature in data matrix
data_matrix_norm <- data_matrix
for (feature in features_all) {
  cat(paste0(skip(), "quantile normalize counts for feature ", feature))
  
  feature_norm <- get_targetQuantileNorm(ecdf[[feature]])
  mcols(data_matrix_norm)[,feature] <- normalize.quantiles.use.target(matrix(mcols(data_matrix_norm)[,feature]),
                                                                      feature_norm)
  done()
}

endPart()


##################################################################
# Extend data matrix
##################################################################

startPart("Extend data matrix (+/- 5 bins)")
  
# create extended data matrix
cat(paste0(skip(), "create extended data matrix"))

# original matrix
data_matrix_ext      <- extend_dataMatrix(N = 5, df = data.frame(data_matrix), f = features1)
zero.idx <- which(rowSums(data_matrix_ext[,-c(1:3)]) == 0)
rm(data_matrix_ext)

# normalized matrix
data_matrix_norm_ext <- extend_dataMatrix(N = 5, df = data.frame(data_matrix_norm), f = features_all)
data_matrix_norm_ext <- data_matrix_norm_ext[-zero.idx,]
done()

endPart()

##################################################################
# make predictions
##################################################################

startPart("Get enhancer probabilities for each bin")

mid=round(nrow(data_matrix_norm_ext)/2,0)
prediction1 <- predict(classifier1, data_matrix_norm_ext[1:mid,], type = "prob")[,2]
prediction1 <- c(prediction1, predict(classifier1, data_matrix_norm_ext[(mid+1):nrow(data_matrix_norm_ext),], type = "prob")[,2])

prediction2 <- predict(classifier2, data_matrix_norm_ext[1:mid,], type = "prob")[,2]
prediction2 <- c(prediction2, predict(classifier2, data_matrix_norm_ext[(mid+1):nrow(data_matrix_norm_ext),], type = "prob")[,2])

rm(data_matrix_norm_ext)

elementMetadata(data_matrix) <- NULL
mcols(data_matrix)["prob"] <- 0
mcols(data_matrix[-zero.idx])["prob"] <- prediction1 * prediction2

endPart()

##################################################################
# save predicitons
##################################################################

startPart("Save predicted values")

# rds format
out.rds <- paste0(outdir, "prediction.rds")
cat(paste0(skip(), "save predictions in rds format: ", out.rds))
saveRDS(data_matrix, out.rds)
done()

# bw format
out.bw <- paste0(outdir, "prediction.bw")
cat(paste0(skip(), "save predictions in bw format: ", out.bw))
colnames(elementMetadata(data_matrix)) <- "score"
seqlengths(data_matrix) <- end(reduce(data_matrix))
export(data_matrix, out.bw)
done()

##################################################################
# Call enhancer peaks and cluster
##################################################################

startPart("Call enhancer peaks and cluster")

# single enhancer peaks
cat(paste0(skip(), "define single enhancers"))
peaks <- get_enhancerPeaks(data_matrix, cutoff, cores)
cat(paste0(skip(), length(peaks), " single enhancer peak(s) found"))
done()

if (length(peaks) > 0) {
  out.bedgraph <- paste0(outdir, "singleEnh.bedGraph")
  cat(paste0(skip(), "save single enhancers in bedGraph format: ", out.bedgraph))
  export(peaks, out.bedgraph)
  done()
}

# enhancer cluster
cat(paste0(skip(), "define cluster of enhancer peaks"))
cluster <- get_enhancerCluster(peaks, distance, cores)
cat(paste0(skip(), length(cluster), " enhancer cluster found"))
done()

if (length(cluster) > 0) {
  out.bed <- paste0(outdir, "clusterEnh.bed")
  cat(paste0(skip(), "save cluster of enhancer peaks  in bed format: ", out.bed))
  export(cluster, out.bed)
  done()
}

endPart()

##################################################################
# print run time
##################################################################

run_time <- Sys.time() - start_time

startPart("Run time")
cat(paste0(skip(), format(run_time), "\n"))
endPart()



