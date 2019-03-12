

# **CRUP - Short Description**   

CRUP ([C]ondition-specific [R]egulatory [U]nits [P]rediction) is a workflow consisting of three 
main steps (CRUP - EP, CRUP - ED, CRUP - ED) and an additonal pre-preparing step (CRUP - normalize),
whereas each step build upon one another.
CRUP collapses different layers of epigenetic information into a single list of regulatory units
consisting of dynamically changing enhancers and target genes.

###### General workflow

CRUP-normalize -> CRUP-EP -> CRUP-ED -> CRUP-ET

###### *Contact*

heinrich@molgen.mpg.de

###### *Citation*

I you are using CRUP, please cite the following publication:\
This is still unpublished work.

Preprint:\
https://www.biorxiv.org/content/biorxiv/early/2018/12/19/501601.full.pdf


# **Run CRUP** 

## CRUP - normalize

This function (input) normalizes and summarizes read counts from ChIP-seq experiments in a
(100 bp) binned genome.


## CRUP - EP ([E]nhancer [P]rediction)

The random forest-based enhancer classifier CRUP-EP (Enhancer Prediction) was developed so that 
it can be applied across different cell types and species without the need of being re-trained.
To guarantee a good transferability, binned ChIP-seq counts are quantile normalized to the 
sample that was used to train the classifier.


## CRUP - ED ([E]nhancer [D]ynamics)

CRUP-ED (Enhancer Dynamics) is based on enhancer probabilities and identifies
condition-specific ('dynamic') enhancer regions by applying a permutation test.
Using empricial p-values, pattern of pairwise significance are build to cluster
adjacent regions.


## CRUP - ET ([E]nhancer ([T]argets)

The method CRUP - ET (Enhancer Targets) was developed to correlate condition-specific enhancers
to normalized RNA-seq experiments.


## Get the help message

Run 'Rscript CRUP.R' or 'Rscript CRUP.R -h' to see all possible functions of CRUP.

Usage: CRUP.R [-[-norm|N]] [-[-prediction|P]] [-[-dynamics|D]] [-[-targets|T]] [-[-help|h]]\
    -N|--norm          computes normalized count values for .bam files\
    -P|--prediction    runs CRUP - EP: (E)nhancer (P)rediction from histone modification\
    -D|--dynamics      runs CRUP - ED: assigns (E)nhancer to (D)ynamic conditions\
    -T|--targets       runs CRUP - ET: correlates (E)nhancer to (T)arget genes\
    -h|--help          this help message

    
###### R PACKAGES

All required R packages are installed at the first run.\
Nothing needs to be installed manually.

> Overview of all packages:\
> 
> getopt\
> bamsignals\
> Rsamtools\
> BSgenome.Mmusculus.UCSC.mm9 # when using genome mm9\
> BSgenome.Mmusculus.UCSC.mm10 # when using genome mm10\
> BSgenome.Hsapiens.UCSC.hg19 # when using genome hg19\
> BSgenome.Hsapiens.UCSC.hg38 # when using genome hg38\
> preprocessCore\
> randomForest\
> GenomicRanges\
> rtracklayer\
> ggplot2\
> dplyr\
> matrixStats\
> parallel\
> GenomicFeatures\
> GenomicAlignments\
> DESeq2\
> TxDb.Mmusculus.UCSC.mm10.knownGene


# Preparation:

The only preparation that has to be done is to create a tab delimited info file that lists
the location of all ChIP-seq experiments in bam file format.
All bam files have to be indexed.
The following histone modifications must be present:\
'H3K27ac', 'H3K4me1', 'H3K4me3' (and 'Input')

The required column names are: 'feature', 'bam_file', 'bam_file_input'

feature         -> lists histone modifications\
bam_file        -> location of the alignments in bam format, e.g.: 'TEST/DATA/ChIPseq/condition1.H3K27ac.bam'\
bam_file_input  -> location of the Input experiments that are associated with each bam_file entry. 

Example info files:

TEST/condition1.info.txt\
TEST/condition2.info.txt



## Run CRUP - norm

Run 'Rscript CRUP.R -N ' to see all possible input parameters:

Usage: CRUP.R [-[-norm|N]] [-[-cores|x] <integer>] [-[-file|f] <character>] [-[-genome|g] <character>] [-[-sequencing|s] <character>] [-[-outdir|o] <character>] [-[-mapq|q] <integer>] [-[-help|h]]
    -N|--norm          computes normalized count values for .bam files
    -x|--cores         number of cores to use (DEFAULT:1)
    -f|--file          summary text file for ChIP-seq experiments
    -g|--genome        genome used in the .bam files ('hg19', 'mm10', 'mm9' or 'hg38')
    -s|--sequencing    type of sequencing ('paired' or 'single')
    -o|--outdir        output directory (DEFAULT: same as 'file' directory)
    -q|--mapq          minimum mapping quality (DEFAULT:10)
    -h|--help          this help message


###### Example run:

**Note: the presented example doesn't have any biological meaning.**

This function (input) normalizes and summarizes read counts from ChIP-seq experiments in a (100 bp) binned genome.
The ChIP-seq experiments are listed in the required info file 'condition1.info.txt'.

Run 'Rscript CRUP.R -N -f TEST/condition1.info.txt -g mm10 -s paired -o TEST/RESULTS/0_NORMALIZED_DATA/'

Output:

> Input normalized and summarized counts in rds file format:\
'TEST/RESULTS/0_NORMALIZED_DATA/condition1.data_matrix.rds'


## Run CRUP - EP       

Run 'Rscript CRUP.R -P' to see all possible input parameters:

Usage: CRUP.R [-[-prediction|P]] [-[-cores|x] <integer>] [-[-matrix|m] <character>] [-[-classifier|c] <character>] [-[-cutoff|u] <double>] [-[-distance|d] <integer>] [-[-outdir|o] <character>] [-[-help|h]]
    -P|--prediction    runs CRUP - EP: (E)nhancer (P)rediction from histone modification
    -x|--cores         number of cores to use (DEFAULT:1)
    -m|--matrix        normalized data matrix (rds file format)
    -c|--classifier    directory of enhancer classifier (DEFAULT: DATA/CLASSIFIER/)
    -u|--cutoff        cutoff for probabilities [0,1] (DEFAULT: 0.5)
    -d|--distance      maximum distance (bp) for peak clustering (DEFAULT: 12500)
    -o|--outdir        output directory (DEFAULT: same as 'file' directory)
    -h|--help          this help message
    

###### Example run:

Run 'Rscript CRUP.R -P -m TEST/RESULTS/0_NORMALIZED_DATA/condition1.data_matrix.rds -o TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/'

Output:

> 1. enhancer probabilities for each 100 bp bin in the genome (bigwig (.bw) and .rds fileformat):\
> 'TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/prediction.bw'\
> 'TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/prediction.rds'\
> 2. enhancer peak calls (bedgraph format):\
> 'TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/singleEnh.bedGraph'\
> 3. cluster of peak (bedgraph format):\
> 'TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/clusterEnh.bed'


## Run CRUP - ED 

Run 'Rscript CRUP.R -D' to see all possible input parameters:

Usage: CRUP.R [-[-dynamics|D]] [-[-cores|x] <integer>] [-[-outdir|o] <character>] [-[-probabilities|p] <character>] [-[-names|n] <character>] [-[-w_0|w] <double>] [-[-threshold|t] <double>] [-[-len|l] <integer>] [-[-help|h]]
    -D|--dynamics         runs CRUP - ED: assigns (E)nhancer to (D)ynamic conditions
    -x|--cores            number of cores to use (DEFAULT:1)
    -o|--outdir           output directory (DEFAULT: same as 'file' directory)
    -p|--probabilities    probabilities in rds format. list: delimiter samples: ':', delimiter conditions: ','
    -n|--names            aternative labels for conditions (DEFAULT: cond1,cond2, ..)
    -w|--w_0              minimum difference between group means [0,1]. (DEFAULT: 0.5)
    -t|--threshold        threshold for p-values in [0,1]. (DEFAULT: 0.01)
    -l|--len              length of flanking region for summarizing. (DEFAULT: 1000)
    -h|--help             this help message


###### Example run:

First: \
predict enhancers in another condition (e.g. 'condition2').

Run 'Rscript CRUP.R -N -f TEST/condition2.info.txt -g mm10 -s paired -o TEST/RESULTS/0_NORMALIZED_DATA/'

Run 'Rscript CRUP.R -P -m TEST/RESULTS/0_NORMALIZED_DATA/condition2.data_matrix.rds -o TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_2/'

Then:\
identify condition-specific enhancer regions.

Run 'Rscript CRUP.R -D -p TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_1/prediction.rds,TEST/RESULTS/1_RF_PREDICTIONS/CONDITION_2/prediction.rds -o TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/ -n cond1,cond2'

Output:

> 1. summarized condition-specific enhancer regions, visualized as a heatmap (colored bt lowest p-value in region):\
>'TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/dynamicEnh__w0_0.1__threshold_0.01.pdf'\
>  - main clusters (1,2, ..) are additionally highlighted on the right border. These clusters describe dynamic enhancer regions that are active on just one condition.
>  - Cluster 'U' describes enhancer regions that are active in all conditions.

> 2. summarized condition-specific enhancer regions:\
> 'TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/dynamicEnh__w0_0.1__threshold_0.01.txt'\
> - 'best.p.value'	-> lowest empirical pvalue in condition-specfic enhancer region\
> - 'cluster'	-> cluster obtained from significance pattern\dynamic enhancer cluster that do not belong to the main clusters (1,2, ..) or cluster 'U' start with an 'r' ('remaining')\
> - 'significance.pattern'  -> significance pattern\
>    - e.g. '01' means that condition 1 and condition 2 are significantly different and the enhancer probability  of condition 2 > condition 1\
> - highest enhancer probability values for each region per sample
> 3. all condition-specific enhancer regions in bed file format:\
> 'TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/dynamicEnh__w0_0.1__threshold_0.01.bed'\
> 4. the main clusters (1,2, ..) and cluster 'U' ('ubiquitous') are additionally exported as separated .bed files


## Run CRUP - ET

Run 'Rscript CRUP.R -T' to see all possible input parameters:

CRUP.R [-[-targets|T]] [-[-cores|x] <integer>] [-[-genome|g] <character>] [-[-sequencing|s] <character>] [-[-outdir|o] <character>] [-[-names|n] <character>] [-[-threshold_c|C] <double>] [-[-regions|r] <character>] [-[-RNA|E] <character>] [-[-expression|e] <character>] [-[-TAD|b] <character>] [-[-help|h]]
    -T|--targets        runs CRUP - ET: correlates (E)nhancer to (T)arget genes
    -x|--cores          number of cores to use (DEFAULT:1)
    -g|--genome         genome used in the .bam files ('hg19', 'mm10' or 'mm9')
    -s|--sequencing     type of sequencing ('paired' or 'single')
    -o|--outdir         output directory (DEFAULT: same as 'file' directory)
    -n|--names          aternative labels for conditions (DEFAULT: cond1,cond2, ..)
    -C|--threshold_c    threshold for correlation in [0.5,1]. (DEFAULT: 0.9)
    -r|--regions        text file with condition-specific regions in txt format
    -E|--RNA            RNA-seq experiments in bam format. list: delimiter samples: ':', delimiter conditions: ','
    -e|--expression     gene expression counts for all samples and conditions
    -b|--TAD            .bed file with TADs (DEFAULT: DATA/mESC_mapq30_KR_all_TADs.bed)
    -h|--help           this help message


###### A) Example run with RNA-seq experiments in bam file format:\

Run 'Rscript CRUP.R -T -r TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/dynamicEnh__w0_0.1__threshold_0.01.txt -g mm10 -s paired -E TEST/DATA/RNAseq/Condition1.bam,TEST/DATA/RNAseq/Condition2.bam -o TEST/RESULTS/3_REGULATORY_REGIONS/'

Output:

> 1. normalized gene expression counts:\
> 'gene_expression.rds'\
> 2. dynamic regulatory units in txt format:\
> 'RegulatoryUnits.txt'\
> - 'seqnames'            -> chr of dynamic enhancer region\
> - 'start'               -> start of dynamic enhancer region\
> - 'end'                 -> end of dynamic enhancer region\
> - 'width'               -> width of dynamic enhancer region\
> - 'strand'              -> strand of dynamic enhancer region\
> - 'cluster'             -> associated cluster of dynamic enhancer region\
> - 'TAD_COORDINATES'     -> coordinates of topologically associated domain around dynamic enhancer region\
> - 'CORRELATED_GENE'     -> ID of the gene that is correlated with dynamic enhancer region\
> - 'CORRELATION'         -> correlation value\
> - best probability values for each region per sample\
> 3. dynamic regulatory units in (ucsc) interaction format:\
> 'RegulatoryUnits.interaction'


###### B) Example run with already summarized RNA-seq experiments:

Run 'Rscript CRUP.R -T -r TEST/RESULTS/2_DIFFERENTIAL_ENHANCERS/dynamicEnh__w0_0.1__threshold_0.01.txt -e TEST/RESULTS/3_REGULATORY_REGIONS/gene_expression.rds -o TEST/RESULTS/3_REGULATORY_REGIONS/'

Output:

> Same output as in A), just without 'gene_expression.rds'.


