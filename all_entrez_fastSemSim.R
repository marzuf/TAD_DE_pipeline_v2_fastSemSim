#  Rscript all_entrez_genes_fastSemSim.R Resnik
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_Resnik_max_result_file.txt

# Rscript all_entrez_genes_fastSemSim.R SimGIC
# 2018-12-10 18:28:02
# 2018-12-29 20:01:06
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt

# Rscript all_entrez_fastSemSim.R


startTime <- Sys.time()

SSHFS <- FALSE

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

setDir <- ifelse(SSHFS, "/media/electron", "")

library(data.table)
library(dplyr)
source("utils.R")

outFold <- "ALL_ENTREZ_FASTSEMSIM"
dir.create(outFold, recursive=T)

simgic_dt <- fread("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt",
                   stringsAsFactors = FALSE)
# Error in fread("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt",  : 
# Opened 13.36GB (14346459291 bytes) file ok but could not memory map it. This is a 64bit process. There is probably not enough contiguous virtual memory available.
               
# naomit_simgic_dt <- eval(parse(text =
#       load("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file_naOmit.Rdata")))

cat("... na.omit the DT\n")
naomit_simgic_dt <- na.omit(simgic_dt)
cat("... as.character obj_1\n")
naomit_simgic_dt$obj_1 <- as.character(naomit_simgic_dt$obj_1)
cat("... as.character obj_2\n")
naomit_simgic_dt$obj_2 <- as.character(naomit_simgic_dt$obj_2)

# TADpos_file <- paste0(setDir, 
#                       "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")
# TADpos_DT <- read.delim(TADpos_file, header=F, stringsAsFactors = FALSE, col.names = c("chromo", "region", "start", "end"))
# TADpos_DT <- TADpos_DT[grep("_TAD", TADpos_DT$region),]

gene2tadDT_file <- paste0(setDir, 
                          "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

g2tDT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2tDT$entrezID <- as.character(g2tDT$entrezID)
g2tDT <- g2tDT[grep("_TAD", g2tDT$region),]

cat("left_join obj_1\n")
merge1 <- left_join(naomit_simgic_dt, g2tDT[,c("entrezID", "region")], by=c("obj_1" = "entrezID"))
colnames(merge1)[colnames(merge1) == "region"] <- "region_obj_1"
cat("left_join obj_2\n")
merge2 <- left_join(merge1, g2tDT[,c("entrezID", "region")], by=c("obj_2" = "entrezID"))
colnames(merge2)[colnames(merge2) == "region"] <- "region_obj_2"

simgic_regions_DT <- merge2
head(simgic_regions_DT)

cat("nrow(simgic_regions_DT) = ", nrow(simgic_regions_DT), "\n")
cat("na.omit(simgic_regions_DT)\n")
simgic_regions_DT <- na.omit(simgic_regions_DT)
cat("nrow(simgic_regions_DT) = ", nrow(simgic_regions_DT), "\n")

simgic_regions_DT$sameTAD <- as.numeric(simgic_regions_DT$region_obj_1 == simgic_regions_DT$region_obj_1 )

aggDT <- aggregate(ss ~ sameTAD, data = simgic_regions_DT, FUN=mean)
cat("aggDT:\n")
write.table(aggDT, file="", sep="\t", quote=F, row.names=F, col.names=T)

outFile <- file.path(outFold, "ss_density_sameTAD_diffTAD.png")
png(outFile, height=300, width=500)
plot_multiDens(list(
  sameTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 1],
  diffTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 0]
))
foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))
    
    
cat(paste0("***DONE\n", startTime,"\n", Sys.time(), "\n"))    










