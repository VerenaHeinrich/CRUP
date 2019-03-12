
##################################################################
# get script path
##################################################################

getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl = TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if (length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if (length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(paste0(normalizePath(script.dir),"/"))
}

path <- getScriptPath()

##################################################################
# path to functions
##################################################################

source(paste0(path,'/BIN/functions.R'))

##################################################################
# libraries
##################################################################

pkgTest("getopt")
suppressMessages(library(getopt))
opt <- getopt(spectrum)

##################################################################
# check input parameter
##################################################################

parameter <- c("norm", "prediction", "dynamics", "targets")
parameter.bool <- sapply(parameter, function(x) !is.null(opt[[x]]))

if ((!is.null(opt$help) & sum(parameter.bool) == 0) | (is.null(opt$help) & sum(parameter.bool) != 1)) {
  cat(headline)
  stop_script(parameter)
}

##################################################################
# create data matrix with normalized read counts
##################################################################

if (!is.null(opt$norm)) source(paste0(path,'/BIN/CRUP_norm.R'))


##################################################################
# run CRUP - (E)nhancer (P)rediction
##################################################################

if (!is.null(opt$prediction)) source(paste0(path,'/BIN/CRUP_EP.R'))


##################################################################
# run CRUP - (E)nhancer (D)ynamics
##################################################################

if (!is.null(opt$dynamics)) source(paste0(path,'/BIN/CRUP_ED.R'))


##################################################################
# run CRUP - (E)nhancer (T)argets
##################################################################

if (!is.null(opt$targets)) source(paste0(path,'/BIN/CRUP_ET.R'))


