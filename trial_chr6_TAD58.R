startTime <- Sys.time()
cat(paste0("> Rscript trial_chr6_TAD58.R\n"))

# Rscript trial_chr6_TAD58.R

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

pipOutDir <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom", "OUTPUT_FOLDER")

curr_dataset <- "TCGAcrc_msi_mss"

curr_tad <- "chr6_TAD58"

outFold <- file.path("FAST_SEM_SIM", curr_dataset)
system(paste0("mkdir -p ", outFold))

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
gene2tad_DT$entrezID <- as.character(gene2tad_DT$entrezID)

entrez2uniprot_file <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/fastsemsim/mz_trials/entrezID_to_uniprotID_mapping.txt")
entrez2uniprot_DT <- read.delim(entrez2uniprot_file, stringsAsFactors = F, header=F, col.names=c("entrezID", "uniprotID"))
entrez2uniprot_DT$entrezID <- as.character(entrez2uniprot_DT$entrezID)

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

paste0("available uniprot mapping = ", sum(tad_genes %in% entrez2uniprot_DT$entrezID), "/", length(tad_genes), "\n")

tad_genes <- tad_genes[tad_genes %in% entrez2uniprot_DT$entrezID]

tad_uniprot <- sapply(tad_genes, function(x) entrez2uniprot_DT$uniprotID[entrez2uniprot_DT$entrezID == x][1])

stopifnot(length(tad_genes) == length(tad_uniprot))

#*****************************************************************************************************
#***************************************************************************************************** PREPARE FAST SEM SIM INPUT FILE
#*****************************************************************************************************

fastSemSim_queryFile <- file.path(outFold, "input", paste0(curr_tad, "_query_file.txt"))
system(paste0("mkdir -p ", dirname(fastSemSim_queryFile)))

# sink(fastSemSim_queryFile)
# paste0(as.character(tad_uniprot), collapse="\n")
# sink()

fileConn <- file(fastSemSim_queryFile)    
writeLines(as.character(tad_uniprot), fileConn)    
close(fileConn)

cat(paste0("... written:\t", fastSemSim_queryFile, "\n"))

#*****************************************************************************************************
#***************************************************************************************************** RUN FAST SEM SIM 
#*****************************************************************************************************
fastSemSim_resultFile <- file.path(outFold, "output", paste0(curr_tad, "_result_file.txt"))
system(paste0("mkdir -p ", dirname(fastSemSim_resultFile)))

fss_root <- "biological_process"
fss_species <- "human"
fss_acFile <- file.path(setDir, "/mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf")
stopifnot(file.exists(fss_acFile))
fss_ontType <- "GeneOntology"
fss_queryType <- "obj"
fss_metric <- "SimGIC"
fss_exec <- "fastsemsim"
fss_queryIn <- "file"
fss_queryMode <- "list"

# fastsemsim --ontology_type GeneOntology --ac_species human --ac_file /mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf --query_ss_type obj --tss SimGIC --query_input file --query_file query_trial4.txt  -vv --output_file result_trial4.txt --query_mode list --root biological_process
# fastsemsim --ontology_type GeneOntology --ac_species human --ac_file //mnt/ed4/marie/software/fastsemsim-0.9.4/ebi_files/goa_human.gaf --query_ss_type obj --tss SimGIC --query_input file --query_file FAST_SEM_SIM/TCGAcrc_msi_mss/input/chr6_TAD58_query_file.txt -vv --output_file FAST_SEM_SIM/TCGAcrc_msi_mss/output/chr6_TAD58_result_file.txt --query_mode list --root biological_process 

fastSemSim_cmd <- paste(fss_exec, "--ontology_type", fss_ontType, "--ac_species", fss_species,
                        "--ac_file", fss_acFile, "--query_ss_type", fss_queryType, "--tss", fss_metric,
                        "--query_input", fss_queryIn, "--query_file", fastSemSim_queryFile, 
                        "-vv", "--output_file",  fastSemSim_resultFile, "--query_mode", fss_queryMode, "--root", fss_root)


cat(paste0("... start fastSemSim:\t", Sys.time(), "\n"))
cat(fastSemSim_cmd, "\n")
system(fastSemSim_cmd)
cat(paste0("... end fastSemSim:\t", Sys.time(), "\n"))
cat(paste0("... fastSemSim result file:\t", fastSemSim_resultFile, "\n"))

#*****************************************************************************************************
#***************************************************************************************************** REFORMAT FAST SEM SIM OUTPUT
#*****************************************************************************************************

formated_fastSemSim_resultFile <- file.path(outFold, "output", paste0(curr_tad, "_result_file_formated.txt"))
system(paste0("mkdir -p ", dirname(formated_fastSemSim_resultFile)))

semsim_DT <- read.delim(fastSemSim_resultFile, header=T, stringsAsFactors = FALSE)

semsim_DT$gene1_entrezID <- sapply(semsim_DT$obj_1, function(x) names(tad_uniprot[tad_uniprot == x]))
semsim_DT$gene2_entrezID <- sapply(semsim_DT$obj_2, function(x) names(tad_uniprot[tad_uniprot == x]))

semsim_DT$gene1_symbol <- sapply(semsim_DT$gene1_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])
semsim_DT$gene2_symbol <- sapply(semsim_DT$gene2_entrezID, function(x) gffDT$symbol[gffDT$entrezID == x])

# newOrder <- c(which(!colnames(semsim_DT) %in% c("obj_1", "obj_2")), which(colnames(semsim_DT) %in% c("obj_1", "obj_2")))
newOrder <- c(which(!colnames(semsim_DT) %in% c("ss")), which(colnames(semsim_DT) %in% c("ss")))
semsim_DT <- semsim_DT[, newOrder]

colnames(semsim_DT)[colnames(semsim_DT) == "obj_1"] <- "gene1_uniprotID"
colnames(semsim_DT)[colnames(semsim_DT) == "obj_2"] <- "gene2_uniprotID"

write.table(semsim_DT, col.names=T, row.names=F, sep="\t", quote=F, file = formated_fastSemSim_resultFile)

cat(paste0("... written formated result file:\t", formated_fastSemSim_resultFile, "\n"))

###############################################################
cat("***DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


