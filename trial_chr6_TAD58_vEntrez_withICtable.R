startTime <- Sys.time()
cat(paste0("> Rscript trial_chr6_TAD58_vEntrez_withICtable.R\n"))

# Rscript trial_chr6_TAD58_vEntrez_withICtable.R [<metric>]

args <- commandArgs(trailingOnly = TRUE)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

source("utils.R")

pipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER")

curr_dataset <- "TCGAcoad_msi_mss"

curr_tad <- "chr6_TAD58"

outFold <- file.path("FAST_SEM_SIM_vEntrez_withICtable", curr_dataset)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "fastSemSim_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))

#************************************************************************************* FAST SEM SIM SETTINGS
fss_root <- "biological_process"
fss_species <- "human"
#fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf")
fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/fastsemsim/mz_data/goa_human_entrezID.gaf")
stopifnot(file.exists(fss_acFile))
fss_ontType <- "GeneOntology"
fss_queryType <- "obj"
fss_metric <- ifelse(length(args) == 1, args[1], "Resnik")
fss_mixStrategy <- "max"
fss_exec <- "fastsemsim"
fss_queryIn <- "file"
fss_queryMode <- "list"

fastSemSim_resultFile <- file.path(outFold, "output", paste0(curr_tad, "_", fss_ontType, "_", fss_root, "_", fss_metric, "_", fss_mixStrategy, "_result_file.txt"))
system(paste0("mkdir -p ", dirname(fastSemSim_resultFile)))


# if fss_ICtable == NULL => run to output the ICtable only
fss_ICtable <- file.path("ALL_GENES_IC_TABLE/output", paste0("all_datasets",fss_ontType, "_", fss_root, "_", fss_metric, "_", fss_mixStrategy, "_IC_table.txt"))
stopifnot(file.exists(fss_ICtable))
fss_ICoutFile <- NULL

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
txt <- paste0("... ", "fss_exec" , "\t=\t", fss_exec , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_queryIn" , "\t=\t", fss_queryIn , "\n")
printAndLog(txt, logFile)
txt <- paste0("... ", "fss_queryMode" , "\t=\t", fss_queryMode , "\n")
printAndLog(txt, logFile)
#*************************************************************************************


gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

gffDT_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gffDT <- read.delim(gffDT_file, header=T, stringsAsFactors = F)
gffDT$entrezID <- as.character(gffDT$entrezID)

script0_name <- "0_prepGeneData"

geneFile <- file.path(pipOutDir, curr_dataset, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneFile))
geneList <- as.character(eval(parse(text = load(geneFile))))

stopifnot(geneList %in% gene2tad_DT$entrezID)

tad_genes <- gene2tad_DT$entrezID[gene2tad_DT$entrezID %in% geneList & gene2tad_DT$region == curr_tad]
paste0("# genes = ", length(tad_genes), "\n")

stopifnot(length(tad_genes) > 0)


#*****************************************************************************************************
#***************************************************************************************************** PREPARE FAST SEM SIM INPUT FILE
#*****************************************************************************************************

fastSemSim_queryFile <- file.path(outFold, "input", paste0(curr_tad, "_query_file.txt"))
system(paste0("mkdir -p ", dirname(fastSemSim_queryFile)))

# sink(fastSemSim_queryFile)
# paste0(as.character(tad_uniprot), collapse="\n")
# sink()

fileConn <- file(fastSemSim_queryFile)    
#writeLines(as.character(tad_uniprot), fileConn)  
writeLines(as.character(tad_genes), fileConn)    
close(fileConn)

cat(paste0("... written:\t", fastSemSim_queryFile, "\n"))

#*****************************************************************************************************
#***************************************************************************************************** RUN FAST SEM SIM 
#*****************************************************************************************************

# fastsemsim --ontology_type GeneOntology --ac_species human --ac_file /mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf --query_ss_type obj --tss SimGIC --query_input file --query_file query_trial4.txt  -vv --output_file result_trial4.txt --query_mode list --root biological_process
# fastsemsim --ontology_type GeneOntology --ac_species human --ac_file //mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf --query_ss_type obj --tss SimGIC --query_input file --query_file FAST_SEM_SIM/TCGAcrc_msi_mss/input/chr6_TAD58_query_file.txt -vv --output_file FAST_SEM_SIM/TCGAcrc_msi_mss/output/chr6_TAD58_result_file.txt --query_mode list --root biological_process 

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
                        "-vv",  
                        "--query_mode", fss_queryMode,
                        "--root", fss_root)


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
  #*****************************************************************************************************
  #***************************************************************************************************** REFORMAT FAST SEM SIM OUTPUT
  #*****************************************************************************************************
  
  formated_fastSemSim_resultFile <- file.path(outFold, "output", paste0(curr_tad, "_", fss_ontType, "_", fss_root, "_", fss_metric, "_", fss_mixStrategy, "_result_file_formated.txt"))
  
  system(paste0("mkdir -p ", dirname(formated_fastSemSim_resultFile)))
  
  semsim_DT <- read.delim(fastSemSim_resultFile, header=T, stringsAsFactors = FALSE)
  
  semsim_DT$gene1_entrezID <- semsim_DT$obj_1
  semsim_DT$gene2_entrezID <- semsim_DT$obj_2
  
  semsim_DT$gene1_symbol <- sapply(semsim_DT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
  semsim_DT$gene2_symbol <- sapply(semsim_DT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
  
  newOrder <- c(which(!colnames(semsim_DT) %in% c("ss")), which(colnames(semsim_DT) %in% c("ss")))
  semsim_DT <- semsim_DT[, newOrder]
  
  write.table(semsim_DT, col.names=T, row.names=F, sep="\t", quote=F, file = formated_fastSemSim_resultFile)
  
  cat(paste0("... written formated result file:\t", formated_fastSemSim_resultFile, "\n"))
}

###############################################################
cat("***DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


