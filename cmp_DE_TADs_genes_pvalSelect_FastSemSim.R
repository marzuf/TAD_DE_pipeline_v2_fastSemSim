startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_pvalSelect_FastSemSim.R\n"))

# Rscript cmp_DE_TADs_genes_pvalSelect_FastSemSim.R

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

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

registerDoMC(ifelse(SSHFS, 2, 30))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")
source("utils.R")
source("run_fastSemSim.R")

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
gffDT_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(gffDT_file, header=T, stringsAsFactors = F)
gffDT$entrezID <- as.character(gffDT$entrezID)

outFold <- file.path("CMP_DE_TADs_GENES_FASTSEMSIM_pvalSelect_1412",nTopDS,rankVar)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) >= nTopDS)

#all_ds <- all_ds[all_ds %in% run1_DS]

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

# settings for FastSemSim
fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/fastsemsim/mz_data/goa_human_entrezID.gaf")
stopifnot(file.exists(fss_acFile))
fss_metric <- "SimGIC"

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

# cat("all_ds[60] = ", all_ds[60], "\n")

stopifnot(!is.na(all_ds))


all_ds_aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, curr_ds, "170_score_auc_pval_withShuffle")
  aucFCC_file <- file.path(step17_fold, "allratio_auc_pval.Rdata")
  if(!file.exists(aucFCC_file)) cat("aucFCC_file = ", aucFCC_file, "\n")
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

#all_ds <- all_ds[1:3]
curr_ds="TCGAcoad_msi_mss"

#topDS <- topDS[1:2]

if(buildTable){
  all_ds_semSim_DT <- foreach(curr_ds = topDS, .combine='rbind') %do% {
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
      ### RUN FAST_SEM_SIM FOR THE TADs GENES
      cat("... start run_fastSemSim\n")
      selectGenes_fssDT <- run_fastSemSim(
        geneList = selectGenes[1:10],
        runName = paste0(curr_ds, "_selectGenes"),
        simMetric = fss_metric,
        gafFile = fss_acFile,
        icData =  IC_table_inFile,
        otherFastSemSimParams = otherFSSparams
      )
      write.table(selectGenes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
      # stopifnot("obj_1" %in% colnames(selectGenes_fssDT))
      # stopifnot("obj_2" %in% colnames(selectGenes_fssDT))
      # selectGenes_fssDT$gene1_entrezID <- selectGenes_fssDT$obj_1
      # selectGenes_fssDT$gene2_entrezID <- selectGenes_fssDT$obj_2
      # selectGenes_fssDT$gene1_symbol <- sapply(selectGenes_fssDT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
      # selectGenes_fssDT$gene2_symbol <- sapply(selectGenes_fssDT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
      # newOrder <- c(which(!colnames(selectGenes_fssDT) %in% c("ss")), which(colnames(selectGenes_fssDT) %in% c("ss")))
      # selectGenes_fssDT <- selectGenes_fssDT[, newOrder]
      # check_outFile <- file.path(outFold, paste0(curr_ds, "selectGenes_fssDT.Rdata"))
      # save(selectGenes_fssDT, file = check_outFile)
      # cat("written = ", check_outFile, "\n")
      # # stop("--ok\n")
      # stopifnot(is.numeric(selectGenes_fssDT$ss)) # not true if NA
      selectGenes_nTestedPairs <- nrow(selectTADs_genes_fssDT)
      selectGenes_fssDT <- na.omit(selectGenes_fssDT)
      selectGenes_nSS <- nrow(selectTADs_genes_fssDT)
      if(nrow(selectGenes_fssDT) == 0) {
        selectGenes_meanSS <- NA
      } else{
        selectGenes_meanSS <- mean(selectGenes_fssDT$ss)  
      }
    } else {
      selectGenes_nTestedPairs <- NA
      selectGenes_nSS <- NA
      selectGenes_meanSS <- NA
    }
        
    if(nSelectTADs_genes > 1) {
        
      # selectTADs_genes <- selectTADs_genes[1:10]
      ### RUN FAST_SEM_SIM FOR THE TADs GENES
      cat("... start run_fastSemSim\n")
      selectTADs_genes_fssDT <- run_fastSemSim(
        geneList = selectTADs_genes[1:10],
        runName = paste0(curr_ds, "_selectTADs_genes"),
        simMetric = fss_metric,
        gafFile = fss_acFile,
        icData =  IC_table_inFile,
        otherFastSemSimParams = otherFSSparams
      )
      write.table(selectTADs_genes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
      # stopifnot("obj_1" %in% colnames(selectTADs_genes_fssDT))
      # stopifnot("obj_2" %in% colnames(selectTADs_genes_fssDT))
      # selectTADs_genes_fssDT$gene1_entrezID <- selectTADs_genes_fssDT$obj_1
      # selectTADs_genes_fssDT$gene2_entrezID <- selectTADs_genes_fssDT$obj_2
      # selectTADs_genes_fssDT$gene1_symbol <- sapply(selectTADs_genes_fssDT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
      # selectTADs_genes_fssDT$gene2_symbol <- sapply(selectTADs_genes_fssDT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
      # newOrder <- c(which(!colnames(selectTADs_genes_fssDT) %in% c("ss")), which(colnames(selectTADs_genes_fssDT) %in% c("ss")))
      # selectTADs_genes_fssDT <- selectTADs_genes_fssDT[, newOrder]
      # check_outFile <- file.path(outFold, paste0(curr_ds, "selectTADs_genes_fssDT.Rdata"))
      # save(selectTADs_genes_fssDT, file = check_outFile)
      # cat("written = ", check_outFile, "\n")
      # # stop("--ok\n")
      # stopifnot(is.numeric(selectTADs_genes_fssDT$ss))# not TRUE if NA
      selectTADs_genes_nTestedPairs <- nrow(selectTADs_genes_fssDT)
      selectTADs_genes_fssDT <- na.omit(selectTADs_genes_fssDT)
      selectTADs_genes_nSS <- nrow(selectTADs_genes_fssDT)
      if(nrow(selectTADs_genes_fssDT) == 0) {
        selectTADs_genes_meanSS <- NA
      } else{
        selectTADs_genes_meanSS <- mean(selectTADs_genes_fssDT$ss)  
      }
      
    }else {
      selectTADs_genes_nTestedPairs <- NA
      selectTADs_genes_nSS <- NA
      selectTADs_genes_meanSS <- NA
    }
      
    data.frame(dataset=curr_ds,
      selectTADs_genes_nTestedPairs = selectTADs_genes_nTestedPairs,
      selectTADs_genes_nSS = selectTADs_genes_nSS,
      selectTADs_genes_meanSS = selectTADs_genes_meanSS,
      selectGenes_nTestedPairs = selectGenes_nTestedPairs,
      selectGenes_nSS = selectGenes_nSS,
      selectGenes_meanSS = selectGenes_meanSS,
stringsAsFactors=FALSE
    )    


  } # end iterating over curr_ds

outFile <- file.path(outFold, "all_ds_semSim_DT.Rdata")
save(all_ds_semSim_DT, file=outFile)
cat(paste0("... written: ", outFile, "\n"))
} # end if buildTable
    
    
    
    
    
    


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

