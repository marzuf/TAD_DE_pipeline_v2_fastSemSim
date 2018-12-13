startTime <- Sys.time()
cat(paste0("> Rscript cmp_dataset_DA_TADgenes_pvalSelect_FastSemSim.R\n"))

# Rscript cmp_oneDataset_DA_TADgenes_pvalSelect_FastSemSim.R
# Rscript cmp_oneDataset_DA_TADgenes_pvalSelect_FastSemSim.R TCGAcoad_msi_mss
# Rscript cmp_oneDataset_DA_TADgenes_pvalSelect_FastSemSim.R TCGAcoad_msi_mss SimGIC

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
fss_metric <- "SimGIC"
args <- commandArgs(trailingOnly = TRUE)
curr_ds <- args[1]
fss_metric <- args[2]

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

outFold <- file.path("CMP_ONEDATASET_DA_TADGENES_FASTSEMSIM_pvalSelect_1312", fss_metric, curr_ds)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_tads_genes_go_logFile.txt")
system(paste0("rm -f ", logFile))

dsFold <- "../TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER"
stopifnot(file.exists(dsFold))

all_ds <- list.files(dsFold)
stopifnot(length(all_ds) > 0)

#all_ds <- all_ds[all_ds %in% run1_DS]

txt <- paste0("... found # datasets:\t", length(all_ds), "\n")
printAndLog(txt, logFile)

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


curr_ds="TCGAcoad_msi_mss"

topDS <- "TCGAcoad_msi_mss"

topDS <- all_ds[1:3]
topDS <- all_ds

topDS <- curr_ds

stopifnot(length(topDS) == 1)

if(buildTable){
  all_ds_DA_TADs_fastSemSim_DT <- foreach(curr_ds = topDS, .combine='rbind') %do% {
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
    nTADs <- length(tad_pval)
    # -> change here pval selection
    selectTADs <- names(tad_pval[tad_pval <= pvalSelect])
    nSelectTADs <- length(selectTADs)
    
    
#    selectTADs <- selectTADs[1:2]
    
    if(nSelectTADs > 0){
      
      
      selectTADs_fastSemSim_DT <- foreach(curr_tad = selectTADs, .combine='rbind') %dopar% {
        
        curr_tad_pval <- tad_pval[curr_tad]
        stopifnot(is.numeric(curr_tad_pval))
        
        stopifnot(curr_tad %in% gene2tad_DT$region)
        curr_TAD_genes <- gene2tad_DT$entrezID[gene2tad_DT$region == curr_tad]  
        stopifnot(length(curr_TAD_genes) > 0)
        
        curr_TAD_genes <- curr_TAD_genes[curr_TAD_genes %in% pipelineGenes]
        stopifnot(length(curr_TAD_genes) > 0)
        
        stopifnot(curr_TAD_genes %in% gene2tad_DT$entrezID)
        stopifnot(curr_TAD_genes %in% pipelineGenes)
  
        ### RUN FAST_SEM_SIM FOR THE SIGNIF TADs GENES
        cat(paste0("... ", curr_ds, " - ", curr_tad, " - start run_fastSemSim\n"))
        
        curr_TAD_genes_fssDT <- run_fastSemSim(
          geneList = curr_TAD_genes,
          runName = paste0(curr_ds, "_selectTADs_", curr_tad, "_genes"),
          simMetric = fss_metric,
          gafFile = fss_acFile,
          icData =  IC_table_inFile,
          otherFastSemSimParams = otherFSSparams
        )
        write.table(curr_TAD_genes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
        stopifnot("obj_1" %in% colnames(curr_TAD_genes_fssDT))
        stopifnot("obj_2" %in% colnames(curr_TAD_genes_fssDT))
        curr_TAD_genes_fssDT$gene1_entrezID <- curr_TAD_genes_fssDT$obj_1
        curr_TAD_genes_fssDT$gene2_entrezID <- curr_TAD_genes_fssDT$obj_2
        curr_TAD_genes_fssDT$gene1_symbol <- sapply(curr_TAD_genes_fssDT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
        curr_TAD_genes_fssDT$gene2_symbol <- sapply(curr_TAD_genes_fssDT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
        curr_TAD_genes_fssDT$obj_1 <- NULL
        curr_TAD_genes_fssDT$obj_2 <- NULL
        newOrder <- c(which(!colnames(curr_TAD_genes_fssDT) %in% c("ss")), which(colnames(curr_TAD_genes_fssDT) %in% c("ss")))
        curr_TAD_genes_fssDT <- curr_TAD_genes_fssDT[, newOrder]
        write.table(curr_TAD_genes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
        # stopifnot(is.numeric(curr_TAD_genes_fssDT$ss)) # not always true because of NA
        
        curr_TAD_genes_fssDT$dataset <- curr_ds
        curr_TAD_genes_fssDT$region <- curr_tad
        curr_TAD_genes_fssDT$region_adjPval <- curr_tad_pval
        curr_TAD_genes_fssDT$region_type <- "selectTADs"
        
        
      curr_TAD_genes_fssDT[,c("dataset", "region", "region_adjPval", "region_type", "gene1_entrezID","gene2_entrezID","gene1_symbol","gene2_symbol","ss")]
        
      }
    } else {
      selectTADs_fastSemSim_DT <- data.frame(
        dataset= curr_ds,
        region= NA,
        region_adjPval = NA,
        gene1_entrezID =NA, 
        gene2_entrezID =NA,
        gene1_symbol=NA,
        gene2_symbol=NA,
        ss=NA,
        stringsAsFactors = FALSE
      )
    }
    
    
    
    
    
    
    # the non signif ones
    nsTADs <- names(tad_pval[tad_pval > pvalSelect])
    nNsTADs <- length(nsTADs)
    stopifnot(nNsTADs + nSelectTADs == nTADs)

        
#    nsTADs <- nsTADs[1:2]
    
    if(nNsTADs > 0){
      
      
      nsTADs_fastSemSim_DT <- foreach(curr_tad = nsTADs, .combine='rbind') %dopar% {
        
        stopifnot(curr_tad %in% gene2tad_DT$region)
        curr_TAD_genes <- gene2tad_DT$entrezID[gene2tad_DT$region == curr_tad]  
        stopifnot(length(curr_TAD_genes) > 0)
        
        curr_TAD_genes <- curr_TAD_genes[curr_TAD_genes %in% pipelineGenes]
        stopifnot(length(curr_TAD_genes) > 0)
        
        stopifnot(curr_TAD_genes %in% gene2tad_DT$entrezID)
        stopifnot(curr_TAD_genes %in% pipelineGenes)
        
        ### RUN FAST_SEM_SIM FOR THE SIGNIF TADs GENES
        cat(paste0("... ", curr_ds, " - ", curr_tad, " - start run_fastSemSim\n"))
        
        curr_TAD_genes_fssDT <- run_fastSemSim(
          geneList = curr_TAD_genes,
          runName = paste0(curr_ds, "_nsTADs_", curr_tad, "_genes"),
          simMetric = fss_metric,
          gafFile = fss_acFile,
          icData =  IC_table_inFile,
          otherFastSemSimParams = otherFSSparams
        )
        write.table(curr_TAD_genes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
        stopifnot("obj_1" %in% colnames(curr_TAD_genes_fssDT))
        stopifnot("obj_2" %in% colnames(curr_TAD_genes_fssDT))
        curr_TAD_genes_fssDT$gene1_entrezID <- curr_TAD_genes_fssDT$obj_1
        curr_TAD_genes_fssDT$gene2_entrezID <- curr_TAD_genes_fssDT$obj_2
        curr_TAD_genes_fssDT$gene1_symbol <- sapply(curr_TAD_genes_fssDT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
        curr_TAD_genes_fssDT$gene2_symbol <- sapply(curr_TAD_genes_fssDT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
        curr_TAD_genes_fssDT$obj_1 <- NULL
        curr_TAD_genes_fssDT$obj_2 <- NULL
        newOrder <- c(which(!colnames(curr_TAD_genes_fssDT) %in% c("ss")), which(colnames(curr_TAD_genes_fssDT) %in% c("ss")))
        curr_TAD_genes_fssDT <- curr_TAD_genes_fssDT[, newOrder]
        write.table(curr_TAD_genes_fssDT[1:5,], col.names=T, row.names=F, sep="\t", quote=F, file="")
        # stopifnot(is.numeric(curr_TAD_genes_fssDT$ss)) # not TRUE because of NA
        
        curr_TAD_genes_fssDT$dataset <- curr_ds
        curr_TAD_genes_fssDT$region <- curr_tad
        curr_TAD_genes_fssDT$region_type <- "nsTADs"
        
        
        curr_TAD_genes_fssDT[,c("dataset", "region", "region_type", "gene1_entrezID","gene2_entrezID","gene1_symbol","gene2_symbol","ss")]
        
      }
    } else {
      nsTADs_fastSemSim_DT <- data.frame(
        dataset= curr_ds,
        region= NA,
        region_adjPval = NA,
        gene1_entrezID =NA, 
        gene2_entrezID =NA,
        gene1_symbol=NA,
        gene2_symbol=NA,
        ss=NA,
        stringsAsFactors = FALSE
      )
    }
    
    
    ds_fastSemSim_DT <- rbind(selectTADs_fastSemSim_DT, nsTADs_fastSemSim_DT)
    
    ds_fastSemSim_DT
    
  } # end iterating datasets
    

} # end if buildtable
    
    
    
outFile <-     file.path(outFold, "all_ds_DA_TADs_fastSemSim_DT.Rdata")
save(all_ds_DA_TADs_fastSemSim_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

