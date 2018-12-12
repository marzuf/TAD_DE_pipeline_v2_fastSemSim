startTime <- Sys.time()
cat(paste0("> Rscript cmp_dataset_DA_FastSemSim_results.R\n"))

# Rscript cmp_dataset_DA_FastSemSim_results.R
# Rscript cmp_dataset_DA_FastSemSim_results.R TCGAcoad_msi_mss
# Rscript cmp_dataset_DA_FastSemSim_results.R TCGAcoad_msi_mss SimGIC

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


outFold <- file.path("CMP_DATASET_FASTSEMSIM", fss_metric, curr_ds)
system(paste0("mkdir -p ", outFold))

resultFile <- file.path("CMP_DATASET_DA_TADGENES_FASTSEMSIM_pvalSelect",
                        curr_ds, fss_metric, "all_ds_DA_TADs_fastSemSim_DT.Rdata")
all_ds_DA_TADs_fastSemSim_DT <- eval(parse(text = load(resultFile)))

head(all_ds_DA_TADs_fastSemSim_DT)


cat(paste0("... discard NA:\t", sum(is.na(all_ds_DA_TADs_fastSemSim_DT$ss)), "\n"))
all_ds_DA_TADs_fastSemSim_DT <- all_ds_DA_TADs_fastSemSim_DT[!is.na(all_ds_DA_TADs_fastSemSim_DT$ss),]

all_ds_DA_TADs_fastSemSim_DT$region_typeLabel <- ifelse(all_ds_DA_TADs_fastSemSim_DT$region_type == "selectTADs", "signif. TADs",
                                                        ifelse(all_ds_DA_TADs_fastSemSim_DT$region_type == "nsTADs", "non signif. TADs",
                                                               NA))
stopifnot(!is.na(all_ds_DA_TADs_fastSemSim_DT$region_typeLabel))

boxplot(ss ~ region_typeLabel, data= all_ds_DA_TADs_fastSemSim_DT,
        main = paste0(curr_ds, " - ", fss_metric),
        ylab = paste0("sem. sim. (", fss_metric, ")")
)

regionTypeDT <- all_ds_DA_TADs_fastSemSim_DT[,c("region", "region_type")]
regionTypeDT <- unique(regionTypeDT)
legend("topleft", bty="n",
       legend = c(paste0("# signif. TADs = ", sum(regionTypeDT$region_type=="selectTADs")),
                  paste0("# gene pairs from signif. TADs = ", sum(all_ds_DA_TADs_fastSemSim_DT$region_type=="selectTADs")),
                  paste0("# ns. TADs = ", sum(regionTypeDT$region_type=="nsTADs")),
                  paste0("# gene pairs from ns. TADs = ", sum(all_ds_DA_TADs_fastSemSim_DT$region_type=="nsTADs")))
         )