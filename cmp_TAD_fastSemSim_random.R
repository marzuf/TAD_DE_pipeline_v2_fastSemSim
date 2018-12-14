


startTime <- Sys.time()
cat(paste0("> Rscript cmp_dTAD_fastSemSim_random.R\n"))

on.exit(cat(paste0(startTime, "\n", Sys.time(), "\n")))

# Rscript cmp_TAD_fastSemSim_random.R
# Rscript cmp_TAD_fastSemSim_random.R TCGAcoad_msi_mss chr6_TAD58
# Rscript cmp_TAD_fastSemSim_random.R TCGAcoad_msi_mss chr6_TAD58 SimGIC

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- T

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

vdHeight <- 7
vdWidth <- 7

curr_ds <- "TCGAcoad_msi_mss"
curr_tad <- "chr6_TAD58"
fss_metric <- "SimGIC"
args <- commandArgs(trailingOnly = TRUE)
curr_ds <- args[1]
curr_tad <- args[2]
fss_metric <- args[3]


nRandom <- 1000

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")
source("utils.R")
source("run_fastSemSim.R")


# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

# HARD CODED
gffDT_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(gffDT_file, header=T, stringsAsFactors = F)
gffDT$entrezID <- as.character(gffDT$entrezID)

outFold <- file.path("CMP_TAD_RANDOM", fss_metric, curr_ds, curr_tad)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("cmp_tad_", nRandom, "random_logFile.txt"))
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

# settings for FastSemSim
fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/fastsemsim/mz_data/goa_human_entrezID.gaf")
stopifnot(file.exists(fss_acFile))

IC_table_inFile <- file.path("ALL_ENTREZ_GENES_IC_TABLE/output",
                             paste0("all_entrez_GeneOntology_biological_process_", fss_metric, "_max_IC_table.txt"))
stopifnot(file.exists(IC_table_inFile))

otherFSSparams <- list("--ac_file" = fss_acFile)

txt <- paste0("... fss_acFile:\t", length(fss_acFile), "\n")
printAndLog(txt, logFile)

txt <- paste0("... fss_metric:\t", length(fss_metric), "\n")
printAndLog(txt, logFile)

txt <- paste0("... IC_table_inFile:\t", length(IC_table_inFile), "\n")
printAndLog(txt, logFile)

##########################################################################################
##########################################################################################

# retrieve the genes for curr_tad

### RETRIEVE GENES THAT ARE IN PIPELINE  
cat("... retrieve available genes for this dataset\n")
step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
stopifnot(file.exists(step0_fold))
pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
stopifnot(file.exists(pipelinegeneFile))
pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))

tad_genes <- gene2tad_DT$entrezID[gene2tad_DT$region == curr_tad & 
                                    gene2tad_DT$entrezID %in% pipelineGenes]
stopifnot(length(tad_genes) > 0)

curr_TAD_genes_fssDT <- run_fastSemSim(
  geneList = tad_genes,
  runName = paste0(curr_ds, "_selectTADs_", curr_tad, "_genes"),
  simMetric = fss_metric,
  gafFile = fss_acFile,
  icData =  IC_table_inFile,
  otherFastSemSimParams = otherFSSparams
)
curr_TAD_genes_fssDT$region <- curr_tad

nGenes <- length(tad_genes)

curr_chromo <- unique(gene2tad_DT$chromo[gene2tad_DT$region == curr_tad])
stopifnot(length(curr_chromo) == 1)

av_genes_DT <- gffDT[gffDT$chromo == curr_chromo,]
av_genes_DT <- av_genes_DT[order(av_genes_DT$start, av_genes_DT$end),]

maxPos <- nrow(av_genes_DT) - nGenes + 1
stopifnot(maxPos > 0)

stopifnot(nRandom > 0)

random_fastSemSim_DT <- foreach(i = 1:nRandom, .combine='rbind') %do% {
    # all_ds_DT <- foreach(curr_ds = topDS, .combine='rbind') %dopar% {
    txt <- paste0("*** START random\t", i, "\n")
    printAndLog(txt, logFile)
    
    # select random set of genes
    randomIdx <- sample(x = 1:maxPos, size = 1)
    random_genes <- as.character(av_genes_DT$entrezID[randomIdx:(randomIdx+nGenes-1)])
    stopifnot(length(random_genes) == nGenes)
    stopifnot(!is.na(random_genes))
    
    random_genes_fssDT <- run_fastSemSim(
      geneList = random_genes,
      runName = paste0(curr_ds, "_selectTADs_", curr_tad, "_genes"),
      simMetric = fss_metric,
      gafFile = fss_acFile,
      icData =  IC_table_inFile,
      otherFastSemSimParams = otherFSSparams
    )
    random_genes_fssDT$region <- paste0("random_", i)
    random_genes_fssDT
}

cmp_TAD_random_DT <- rbind(curr_TAD_genes_fssDT, random_fastSemSim_DT)


write.table(cmp_TAD_random_DT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
stopifnot("obj_1" %in% colnames(cmp_TAD_random_DT))
stopifnot("obj_2" %in% colnames(cmp_TAD_random_DT))
cmp_TAD_random_DT$gene1_entrezID <- cmp_TAD_random_DT$obj_1
cmp_TAD_random_DT$gene2_entrezID <- cmp_TAD_random_DT$obj_2
cmp_TAD_random_DT$gene1_symbol <- sapply(cmp_TAD_random_DT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
cmp_TAD_random_DT$gene2_symbol <- sapply(cmp_TAD_random_DT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
cmp_TAD_random_DT$obj_1 <- NULL
cmp_TAD_random_DT$obj_2 <- NULL
newOrder <- c(which(!colnames(cmp_TAD_random_DT) %in% c("ss")), which(colnames(cmp_TAD_random_DT) %in% c("ss")))
cmp_TAD_random_DT <- cmp_TAD_random_DT[, newOrder]
write.table(cmp_TAD_random_DT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
stopifnot(is.numeric(cmp_TAD_random_DT$ss))

cmp_TAD_random_DT$dataset <- curr_ds

cmp_TAD_random_DT <- cmp_TAD_random_DT[,c("dataset", "region", "gene1_entrezID","gene2_entrezID","gene1_symbol","gene2_symbol","ss")]

outFile <-     file.path(outFold, paste0("cmp_TAD_", nRandom, "random_DT.Rdata"))
save(cmp_TAD_random_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# outFile <- "CMP_TAD_RANDOM/SimGIC/TCGAcoad_msi_mss/chr6_TAD58/cmp_TAD_random_DT.Rdata"
# cmp_TAD_random_DT <- eval(parse(text = load(outFile)))
cmp_TAD_random_DT$regionType <- ifelse(cmp_TAD_random_DT$region==curr_tad, curr_tad, "random")

outFile <- file.path(outFold, paste0(curr_tad, "_vs_", nRandom, "random_cmpSS_multiDens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth*1.5))
plot_multiDens(
list(ss_in_TADs = na.omit(cmp_TAD_random_DT$ss[grepl("_TAD", cmp_TAD_random_DT$region)]),
ss_in_random = na.omit(cmp_TAD_random_DT$ss[grepl("random", cmp_TAD_random_DT$region)])),
plotTit = paste0("All SS values distribution"),
my_xlab = paste0("SS (", fss_metric,")")
)
mtext(text = paste0(fss_metric, " - nRandom = ", nRandom), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0(curr_tad, "_vs_", nRandom, "random_cmpSS_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth))
boxplot(ss ~ regionType, data = cmp_TAD_random_DT,
        main = paste0("Pairwise SS"),
        ylab = paste0("SS (", fss_metric,")"))
mtext(text = paste0(fss_metric, " - nRandom = ", nRandom), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


ss_mean_DT <- aggregate(ss ~ region, data = cmp_TAD_random_DT, FUN = mean, na.rm=T)

outFile <- file.path(outFold, paste0(curr_tad, "_vs_", nRandom, "random_cmpMean_dens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width = myWidth*1.5))
plot(density(ss_mean_DT$ss[grepl("^random_", ss_mean_DT$region)]),
     xlab=paste0("mean SS (", fss_metric, ")"),
     main = paste0("Mean SS: distribution random and value observed (", curr_tad,")"))
abline(v = ss_mean_DT$ss[!grepl("^random_", ss_mean_DT$region)], 
       col = "red")
mtext(text = paste0(fss_metric, " - nRandom = ", nRandom), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")



