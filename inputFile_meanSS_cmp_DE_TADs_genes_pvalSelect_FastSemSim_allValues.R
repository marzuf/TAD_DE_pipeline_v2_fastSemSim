startTime <- Sys.time()
cat(paste0("> Rscript inputFile_meanSS_cmp_DE_TADs_genes_pvalSelect_FastSemSim_allValues.R\n"))

# Rscript inputFile_meanSS_cmp_DE_TADs_genes_pvalSelect_FastSemSim_allValues.R

# allValues -> do not aggregate by the mean !

options(scipen=100)
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

buildTable <- FALSE

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

vdHeight <- 7
vdWidth <- 7

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")
source("utils.R")

nTopDS <- 5
rankVar <- "FCC"

args <- commandArgs(trailingOnly = T)
if(length(args) >= 1) {
  if(!is.na(as.numeric(args[1])))
    nTopDS <- as.numeric(args[1])
}
if(length(args) == 2) {
  rankVar <- args[2]
}
stopifnot(rankVar %in% c("FCC", "coexpr", "avg"))

nTopResultsToSave <- 10

#topPercent <- 0.10
# cmp the top topPercent TADs and topPercent genes

# ! IN THIS VERSION, SELECTION OF TAD GENES AND GENES BASED ON PVAL
pvalSelect <- 0.05

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

# HARD CODED
# gffDT_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
# gffDT <- read.delim(gffDT_file, header=T, stringsAsFactors = F)
# gffDT$entrezID <- as.character(gffDT$entrezID)

outFold <- file.path("INPUTFILE_MEANSS_CMP_DE_TADs_GENES_FASTSEMSIM_pvalSelect_ALL_VALUES",nTopDS,rankVar)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS)

tcga_ds <- all_ds[grepl("^TCGA", all_ds)]

#all_ds <- all_ds[all_ds %in% run1_DS]

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

# settings for FastSemSim
fss_metric <- "SimGIC"
fss_combType <- "max"

if(buildTable) {
  
  SS_inFile <- file.path("ALL_GENES_FASTSEMSIM/output", paste0("all_datasets_GeneOntology_biological_process_", fss_metric, "_", fss_combType, "_result_file.txt"))
  stopifnot(file.exists(SS_inFile))
  
  txt <- paste0("... SS_inFile:\t", length(SS_inFile), "\n")
  printAndLog(txt, logFile)
  
  cat("... load all SS input file: ", SS_inFile, "\n")
  ssDT <- fread(SS_inFile) 
  stopifnot(colnames(ssDT) == c("obj_1", "obj_2","ss"))
  
  cat("... done - nrow (# gene pairs) = ", nrow(ssDT), "\n")
  # ssDT <- na.omit(ssDT)
  # cat("... done - nrow (# gene pairs) after na.omit = ", nrow(ssDT), "\n")
  cat("... convert to character\n")
  ssDT$obj_1 <- as.character(ssDT$obj_1)
  ssDT$obj_2 <- as.character(ssDT$obj_2)
  
}

##########################################################################################
##########################################################################################


stopifnot(!is.na(all_ds))

all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  # if(!file.exists(aucFCC_file)) 
  #   cat("aucFCC_file = ", aucFCC_file, "\n")
  stopifnot(file.exists(aucFCC_file))
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}
names(all_ds_aucFCC) <- all_ds

all_ds_aucCoexprDist <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucCoexprDist_file <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", "TopDom"),
                                  "AUC_COEXPRDIST_SORTNODUP", curr_ds, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDist
}
names(all_ds_aucCoexprDist) <- all_ds

maxNds <- min(nTopDS, length(all_ds_aucFCC))

if(rankVar == "FCC"){
  topDS <- names(sort(all_ds_aucFCC, decreasing=TRUE)[1:maxNds])  
} else if(rankVar == "coexpr"){
  topDS <- names(sort(all_ds_aucCoexprDist, decreasing=TRUE)[1:maxNds])
} else if(rankVar == "avg"){
  tmpFCC <- sort(rank(-all_ds_aucFCC, ties = "first" ))
  tmpCoexpr <- sort(rank(-all_ds_aucCoexprDist, ties = "first" ))
  stopifnot(names(tmpFCC) %in% names(tmpCoexpr))
  stopifnot(names(tmpCoexpr) %in% names(tmpFCC))
  tmpFCC <- tmpFCC[names(tmpFCC)]
  tmpCoexpr <- tmpCoexpr[names(tmpFCC)]
  stopifnot(!is.na(tmpFCC))
  stopifnot(!is.na(tmpCoexpr))
  meanScore <- 0.5*(tmpFCC + tmpCoexpr)
  stopifnot(names(meanScore) == names(tmpFCC))
  topDS <- names(sort(meanScore, decreasing=F)[1:maxNds])
} else {
  stop("-- error rankVar -- \n")
}

stopifnot(!is.na(topDS))

##########################################################################################
##########################################################################################

topDS <-  all_ds

if(buildTable){
  all_ds_meanSS_from_file_DT <- foreach(curr_ds = topDS, .combine='rbind') %do% {
    # all_ds_DT <- foreach(curr_ds = topDS, .combine='rbind') %dopar% {
    txt <- paste0("*** START:\t", curr_ds, "\n")
    printAndLog(txt, logFile)
    
    ### RETRIEVE GENES THAT ARE IN PIPELINE  
    cat("... retrieve available genes for this dataset\n")
    step0_fold <- file.path(dsFold, curr_ds, "0_prepGeneData")
    stopifnot(file.exists(step0_fold))
    pipelinegeneFile <- file.path(step0_fold, "pipeline_geneList.Rdata")
    stopifnot(file.exists(pipelinegeneFile))
    pipelineGenes <- eval(parse(text = load(pipelinegeneFile)))
    
    
    ### RETRIEVE TAD GENES AND PVALUES
    cat("... retrieve top TADs genes\n")
    step11_fold <- file.path(dsFold, curr_ds, "11_runEmpPvalCombined")
    stopifnot(file.exists(step11_fold))
    tadpvalFile <- file.path(step11_fold, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tadpvalFile))
    tad_pval <- eval(parse(text = load(tadpvalFile)))
    tad_pval <- p.adjust(tad_pval, method = "BH")
    # -> change here pval selection
    selectTADs <- names(tad_pval[tad_pval <= pvalSelect])
    nSelectTADs <- length(selectTADs)
    
    nTADs <- length(tad_pval)
    
    if(nSelectTADs > 0){
      stopifnot(selectTADs %in% gene2tad_DT$region)
      selectTADs_genes <- gene2tad_DT$entrezID[gene2tad_DT$region %in% selectTADs]  
      selectTADs_genes <- selectTADs_genes[selectTADs_genes %in% pipelineGenes]
      stopifnot(selectTADs_genes %in% gene2tad_DT$entrezID)
      stopifnot(selectTADs_genes %in% pipelineGenes)
      nSelectTADs_genes <- length(selectTADs_genes)
      
    } else {
      nSelectTADs_genes <- 0
    }
    
    ### RETRIEVE GENES PVALUES
    cat("... retrieve signif DE genes\n")
    step1_fold <- file.path(dsFold, curr_ds, "1_runGeneDE")
    stopifnot(file.exists(step1_fold))
    toptableFile <- file.path(step1_fold, "DE_topTable.Rdata")
    stopifnot(file.exists(toptableFile))
    topTable_DT <- eval(parse(text = load(toptableFile)))
    stopifnot(!any(duplicated(topTable_DT$genes)))
    stopifnot(topTable_DT$genes == rownames(topTable_DT))
    gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    stopifnot(names(pipelineGenes) %in% topTable_DT$genes)
    
    topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipelineGenes),]
    genes_entrez <- sapply(topTable_DT$genes, function(x) pipelineGenes[x])
    
    stopifnot(!is.na(genes_entrez))
    stopifnot(length(genes_entrez) == nrow(topTable_DT) )
    stopifnot(names(pipelineGenes) %in% names(gene_pval))
    # gene_pval <- setNames(topTable_DT$adj.P.Val, topTable_DT$genes)
    entrez_pval <- setNames(topTable_DT$adj.P.Val, genes_entrez)
    
    nGenes <- length(entrez_pval)
    
    selectGenes <- names(entrez_pval[entrez_pval <= pvalSelect])
    stopifnot(selectGenes %in% pipelineGenes)
    nSelectGenes <- length(selectGenes)
    
    stopifnot(selectTADs %in% gene2tad_DT$region)
    
    if(nSelectGenes > 0){
      stopifnot(selectGenes %in% gene2tad_DT$entrezID)
      stopifnot(selectGenes %in% pipelineGenes)
    }
    
    if(length(selectGenes) > 1) {
      # selectGenes <- selectGenes[1:10]
      cat("... retrieve SS for selectGenes\n")
      cat("...... check that all selectGenes in ssDT\n")
      stopifnot(selectGenes %in% ssDT$obj_1 | selectGenes %in% ssDT$obj_2)
      cat("...... subset ssDT for selectGenes\n")
      selectGenes_ssDT <- ssDT[ssDT$obj_1 %in% selectGenes & ssDT$obj_2 %in% selectGenes,]
      cat("...... compute meanSS for selectGenes\n")
      selectGenes_meanSS <- mean(selectGenes_ssDT$ss, na.rm=T)
    } else {
      selectGenes_meanSS <- NA
    }
    selectGenes_ssDT <- data.frame(
      dataset = curr_ds,
      geneType = "selectGenes",
      meanSS = selectGenes_meanSS,
      stringsAsFactors = FALSE
    )
    
    if(nSelectTADs_genes > 1) {
      cat("... retrieve SS for selectTADs_Genes\n")
      cat("...... check that all selectTADs_genes in ssDT\n")
      stopifnot(selectTADs_genes %in% ssDT$obj_1 | selectTADs_genes %in% ssDT$obj_2)
      cat("...... subset ssDT for selectGenes\n")
      selectTADs_genes_ssDT <- ssDT[ssDT$obj_1 %in% selectTADs_genes & ssDT$obj_2 %in% selectTADs_genes,]
      cat("...... compute meanSS for selectTADs_genes\n")
      selectTADs_genes_meanSS <- mean(selectTADs_genes_ssDT$ss, na.rm=T)
    }else {
      selectTADs_genes_meanSS <- NA
    }
    selectTADs_genes_ssDT <- data.frame(
      dataset = curr_ds,
      geneType = "selectTADs_genes",
      meanSS = selectTADs_genes_meanSS,
      stringsAsFactors = FALSE
    )    
    rbind(selectGenes_ssDT, selectTADs_genes_ssDT)
  } # end iterating over curr_ds
  outFile <- file.path(outFold, "all_ds_meanSS_from_file_DT.Rdata")
  save(all_ds_meanSS_from_file_DT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {# end if buildTable
  outFile <- file.path(outFold, "all_ds_meanSS_from_file_DT.Rdata")
  all_ds_meanSS_from_file_DT <- eval(parse(text = load(outFile)))
}

head(all_ds_meanSS_from_file_DT)

nSelectGenes <- length(all_ds_meanSS_from_file_DT$dataset[all_ds_meanSS_from_file_DT$geneType=="selectGenes" & !is.na(all_ds_meanSS_from_file_DT$meanSS)])
nSelectTADs_genes <- length(all_ds_meanSS_from_file_DT$dataset[all_ds_meanSS_from_file_DT$geneType=="selectTADs_genes" & !is.na(all_ds_meanSS_from_file_DT$meanSS)])

myylab <- paste0("Avg. SS (", fss_metric, " - ", fss_combType, ")")
myxlab <- paste0("")
mySub <- paste0("# DS: selectGenes=", nSelectGenes, "; selectTADs_genes=", nSelectTADs_genes)
# boxplot(meanSS ~ geneType, data = all_ds_meanSS_from_file_DT)
# without outliers:

outFile <- file.path(outFold, paste0("boxplot_selectGenes_selectTADsgenes_meanSS_no_out.",plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
boxplot(meanSS ~ geneType, data = all_ds_meanSS_from_file_DT, outline=F,
        main=paste0("all DS (n=", length(unique(all_ds_meanSS_from_file_DT$dataset)), ")"),
          xlab = myxlab, ylab = myylab
        )
mtext(side=3, text = mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# plot_multiDens(
#  list(
#    selectGenes = all_ds_meanSS_from_file_DT$meanSS[all_ds_meanSS_from_file_DT$geneType=="selectGenes" & !is.na(all_ds_meanSS_from_file_DT$meanSS)],
#    selectTADs_genes = all_ds_meanSS_from_file_DT$meanSS[all_ds_meanSS_from_file_DT$geneType=="selectTADs_genes" & !is.na(all_ds_meanSS_from_file_DT$meanSS)]
#  ) 
# )

outFile <- file.path(outFold, paste0("multiDens_selectGenes_selectTADsgenes_meanSS_cut.",plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))
plot_multiDens(
  list(
    selectGenes = all_ds_meanSS_from_file_DT$meanSS[all_ds_meanSS_from_file_DT$geneType=="selectGenes" & !is.na(all_ds_meanSS_from_file_DT$meanSS) & all_ds_meanSS_from_file_DT$meanSS < 0.1],
    selectTADs_genes = all_ds_meanSS_from_file_DT$meanSS[all_ds_meanSS_from_file_DT$geneType=="selectTADs_genes" & !is.na(all_ds_meanSS_from_file_DT$meanSS) & all_ds_meanSS_from_file_DT$meanSS < 0.1]
  ) ,
  plotTit = paste0("(cut at 0.1)"),
  my_xlab = "mean SS"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

