#  Rscript all_entrez_genes_fastSemSim.R Resnik
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_Resnik_max_result_file.txt

# Rscript all_entrez_genes_fastSemSim.R SimGIC
# 2018-12-10 18:28:02
# 2018-12-29 20:01:06
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt

### USING PIPELINE TISSUE SPECIFIC G2T FILE

# Rscript all_entrez_fastSemSim_otherTADfile.R GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401

cat("> START all_entrez_fastSemSim_otherTADfile.R\n")

startTime <- Sys.time()

SSHFS <- FALSE

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

setDir <- ifelse(SSHFS, "/media/electron", "")

library(data.table)
library(dplyr)
source("utils.R")

dataset="GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
dataset <- args[1]

outFold <- file.path("ALL_ENTREZ_FASTSEMSIM", dataset)
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

# gene2tadDT_file <- paste0(setDir, 
#                           "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

# /mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401/genes2tad/all_genes_positions.txt
gene2tadDT_file <- paste0(setDir, 
                          "/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/", dataset, "/genes2tad/all_genes_positions.txt") 
stopifnot(file.exists(gene2tadDT_file))
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

cat("create sameTAD column\n")
simgic_regions_DT$sameTAD <- as.numeric(simgic_regions_DT$region_obj_1 == simgic_regions_DT$region_obj_2 )

cat("simgic_regions_DT:\n")
write.table(simgic_regions_DT[1:10,], file="", sep="\t", quote=F, row.names=F, col.names=T)

cat("draw multidens\n")
outFile <- file.path(outFold, "ss_density_sameTAD_diffTAD.png")
png(outFile, height=300, width=500)
plot_multiDens(list(
  sameTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 1],
  diffTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 0]
))
foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))

cat("aggregate mean SS by sameTAD column \n")
mean_aggDT <- aggregate(ss ~ sameTAD, data = simgic_regions_DT, FUN=mean, na.rm=TRUE)
colnames(mean_aggDT)[colnames(mean_aggDT) == "ss"] <- "meanSS"
cat("mean_aggDT:\n")
write.table(mean_aggDT, file="", sep="\t", quote=F, row.names=F, col.names=T)

cat("aggregate # SS by sameTAD column \n")
nbr_aggDT <- aggregate(ss ~ sameTAD, data = simgic_regions_DT, FUN=function(x) length(na.omit(x)))
colnames(nbr_aggDT)[colnames(nbr_aggDT) == "ss"] <- "nbrSS"
cat("nbr_aggDT:\n")
write.table(nbr_aggDT, file="", sep="\t", quote=F, row.names=F, col.names=T)

aggDT <- merge(nbr_aggDT, mean_aggDT, by="sameTAD")

cat("aggDT:\n")
write.table(aggDT, file="", sep="\t", quote=F, row.names=F, col.names=T)

outFile <- file.path(outFold, "aggDT.txt")
write.table(aggDT, file=outFile, sep="\t", quote=F, row.names=F, col.names=T)
cat(paste0("written: ", outFile, "\n"))

# outFile <- file.path(outFold, "meanSS_density_sameTAD_diffTAD.png")
# png(outFile, height=300, width=500)
# plot_multiDens(list(
#   sameTAD = aggDT$ss[aggDT$sameTAD == 1],
#   diffTAD = aggDT$ss[aggDT$sameTAD == 0]
# ))
# foo <- dev.off()
# cat(paste0("written: ", outFile, "\n"))
    
cat(paste0("***DONE\n", startTime,"\n", Sys.time(), "\n"))    


dataDT <-	read.table(textConnection("
dataset sameTAD nbrPairs meanSS FC
Pipeline_consensus	0	104770672	0.04941	7.03
Pipeline_consensus	1	71768	0.34720	7.03
Astro	0	62095858	0.04920	5.25
Astro	1	70967	0.25854	5.25
A549_NCIH460	0	44997448	0.04997	5.05
A549_NCIH460	1	65823	0.25251	5.05
Caki2_G401	0	43766581	0.05088	4.42
Caki2_G401	1	70985	0.22508	4.42
DLD1	0	90003529	0.04975	5.49
DLD1	1	78224	0.27315	5.49
MCF7	0	56598882	0.04945	6.84
MCF7	1	53808	0.33839	6.84
Panc1	0	89772452	0.04981	4.8
Panc1	1	94669	0.23894	4.8
RPMI7951_SKMEL5	0	44893549	0.05035	4.55
RPMI7951_SKMEL5	1	65354	0.22925	4.55
                                    "),
                     header=TRUE)

dataDT$nbrPairs_log10 <- log10(dataDT$nbrPairs)

plot(meanSS ~ nbrPairs, 
     data = dataDT[dataDT$sameTAD == 0,],
     pch=16,
     cex=1.2,
     main="mean SS vs. # gene pairs (diffTAD)",
     xlab="# pairs", ylab="meanSS",
     cex.axis=1.2,cex.lab=1.2)

plot(meanSS ~ nbrPairs_log10, 
     data = dataDT[dataDT$sameTAD == 0,],
     pch=16,
     cex=2,
     main="mean SS vs. # gene pairs (diffTAD)",
     xlab="# pairs [log10]", ylab="meanSS",
     cex.axis=2,cex.lab=2,cex.main=2)



plot(meanSS ~ nbrPairs, 
     data = dataDT[dataDT$sameTAD == 1,],
     pch=16,
     cex=1.2,
     main="mean SS vs. # gene pairs (sameTAD)",
     xlab="# pairs", ylab="meanSS",
     cex.axis=1.2,cex.lab=1.2)


plot(meanSS ~ nbrPairs_log10, 
     data = dataDT[dataDT$sameTAD == 1,],
     pch=16,
     cex=2,
     main="mean SS vs. # gene pairs (sameTAD)",
     xlab="# pairs [log10]", ylab="meanSS",
     cex.axis=2,cex.lab=2,cex.main=2)








