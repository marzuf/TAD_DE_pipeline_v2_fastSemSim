# retrieve all the genes from all datasets to create the IC table

library(foreach)
library(doMC)

startTime <- Sys.time()
cat(paste0("> Rscript all_pipeline_genes_prepICtable.R\n"))

# Rscript all_pipeline_genes_prepICtable.R [<metric>]
# Rscript all_pipeline_genes_prepICtable.R Resnik
args <- commandArgs(trailingOnly = TRUE)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

source("utils.R")

pipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER")

outFold <- file.path("ALL_GENES_IC_TABLE")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "allGenesICtable_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))

#************************************************************************************* FAST SEM SIM SETTINGS
fastSemSim_resultFile <- NULL

fss_root <- "biological_process"
fss_species <- "human"
#fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf")
fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/fastsemsim/mz_data/goa_human_entrezID.gaf")
stopifnot(file.exists(fss_acFile))
fss_ontType <- "GeneOntology"
fss_queryType <- "obj"
fss_metric <- ifelse(length(args)==1, args[1], "Resnik")
fss_mixStrategy <- "max"
fss_exec <- "fastsemsim"
fss_queryIn <- "file"
fss_queryMode <- "list"


fastSemSim_queryFile <- file.path(outFold, "input", paste0("all_datasets_",fss_ontType, "_", fss_root, "_", fss_metric, "_all_entrezGenes.txt"))
system(paste0("mkdir -p ", dirname(fastSemSim_queryFile)))


# if fss_ICtable == NULL => run to output the ICtable only
fss_ICtable <- NULL
fss_ICoutFile <- file.path(outFold, "output", paste0("all_datasets_",fss_ontType, "_", fss_root, "_", fss_metric, "_", fss_mixStrategy, "_IC_table.txt"))
system(paste0("mkdir -p ", dirname(fss_ICoutFile)))


txt <- paste0("... ", "fss_root" , "\t=\t", fss_root , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_species" , "\t=\t", fss_species , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_acFile" , "\t=\t", fss_acFile , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_ontType" , "\t=\t", fss_ontType , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_queryType" , "\t=\t", fss_queryType , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_metric" , "\t=\t", fss_metric , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_mixStrategy" , "\t=\t", fss_mixStrategy , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_exec" , "\t=\t", fss_exec , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_queryIn" , "\t=\t", fss_queryIn , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_queryMode" , "\t=\t", fss_queryMode , "\n")
printAndLog(txt, logFile)
#*************************************************************************************

### RETRIEVE ALL THE GENES FROM ALL THE DATASETS

all_gene_files <- system(paste0("ls ", pipOutDir, "/*/0_prepGeneData/pipeline_geneList.Rdata"), intern=TRUE)


all_gene_entrez <- foreach(geneF = all_gene_files, .combine='c') %dopar% {
  as.character(eval(parse(text = load(geneF))))
}
all_gene_entrez <- unique(all_gene_entrez)

cat("# of genes retrieved from all datasets\t=\t", length(all_gene_entrez), "\n")

#****************************************************************************************
# prepare gene list for fast sem sim

fileConn <- file(fastSemSim_queryFile)    
#writeLines(as.character(tad_uniprot), fileConn)  
writeLines(as.character(all_gene_entrez), fileConn)    
close(fileConn)

cat(paste0("... written:\t", fastSemSim_queryFile, "\n"))

#****************************************************************************************
# run fast sem sim to produce the IC table

fastSemSim_cmd <- paste(fss_exec, 
                        "--ontology_type", fss_ontType, 
                        "--ac_species", fss_species,
                        "--ac_file", fss_acFile, 
                        "--query_ss_type", fss_queryType, 
                        "--tss", fss_metric,
                        "--tmix", fss_mixStrategy,
                        "--query_input", fss_queryIn, 
                        "--query_file", fastSemSim_queryFile, 
                        
                        ifelse(is.null(fss_ICtable),
                               # export the IC table:
                               paste("--stats_IC", "--task stats","--output_file",  fss_ICoutFile),
                               # in case to import the IC table:  #   --inject_IC inject_IC, --inject_IC_form_file inject_IC
                               paste("--inject_IC", fss_ICtable, "--task SS", "--output_file",  fastSemSim_resultFile)
                        ),
                        "-vv",  "--query_mode", fss_queryMode, "--root", fss_root)


cat(paste0("... start fastSemSim:\t", Sys.time(), "\n"))
cat(fastSemSim_cmd, "\n")
system(fastSemSim_cmd)
cat(paste0("... end fastSemSim:\t", Sys.time(), "\n"))

if(is.null(fss_ICtable)) {
  stopifnot(file.exists(fss_ICoutFile))
  cat(paste0("... IC table file:\t", fss_ICoutFile, "\n"))
} else {
  stopifnot(file.exists(fastSemSim_resultFile))
  cat(paste0("... fastSemSim result file:\t", fastSemSim_resultFile, "\n"))
}

  
######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
