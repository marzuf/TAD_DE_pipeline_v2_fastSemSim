#  Rscript all_entrez_genes_fastSemSim.R Resnik
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_Resnik_max_result_file.txt

# Rscript all_entrez_genes_fastSemSim.R SimGIC
# 2018-12-10 18:28:02
# 2018-12-29 20:01:06
# => ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt

### USING PIPELINE TISSUE SPECIFIC G2T FILE

# Rscript all_entrez_fastSemSim_otherTADfile_correctDist.R GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401

cat("> START all_entrez_fastSemSim_otherTADfile_correctDist.R\n")

startTime <- Sys.time()

buildDT <- FALSE
buildModel <- FALSE

cat(paste0(" !!! ... buildDT = ", as.character(buildDT), "\n"))
cat(paste0(" !!! ... buildModel = ", as.character(buildModel), "\n"))

SSHFS <- FALSE

if(SSHFS) setwd("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_fastSemSim")

setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

library(data.table)
library(dplyr)
library(flux)
source("utils.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.2

distLimit <- 500*10^3
nbrLoessPoints <- 1000
distVect <- seq(from=0, to = distLimit, length.out = nbrLoessPoints)
nSamp_shortModel <- 20000 # FOR PLOTTING CI

dataset="GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
dataset <- args[1]

outFold <- file.path("ALL_ENTREZ_FASTSEMSIM", dataset)
dir.create(outFold, recursive=T)


sameTADcol <- "darkorange1"
diffTADcol <- "darkslateblue"

if(buildDT) {

  distFile <- file.path(setDir, "mnt/etemp/marie/Dixon2018_integrative_data/CREATE_DIST_SORTNODUP", dataset, "all_dist_pairs.Rdata")
  stopifnot(file.exists(distFile))
  cat("... load inter-gene distance data\n")
  distDT <- eval(parse(text=load(distFile)))
  head(distDT)
  #       gene1     gene2 chromo     dist
  # 1 100009667 100038246  chr10 33099818
  # 2 100009667     10006  chr10 42677402
  # 3 100009667 100118954  chr10 63950905
  distDT$gene1 <- as.character(distDT$gene1)
  distDT$gene2 <- as.character(distDT$gene2)

                              #
  cat("... load SS data\n")
  simgic_dt <- fread("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt",
                     stringsAsFactors = FALSE)
  # Error in fread("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file.txt",  :
  # Opened 13.36GB (14346459291 bytes) file ok but could not memory map it. This is a 64bit process. There is probably not enough contiguous virtual memory available.
  
  # naomit_simgic_dt <- eval(parse(text =
  #       load("ALL_ENTREZ_GENES_FASTSEMSIM/output/all_entrez_GeneOntology_biological_process_SimGIC_max_result_file_naOmit.Rdata")))
  cat("... na.omit the DT\n")
  naomit_simgic_dt <- na.omit(simgic_dt)
  cat("... remove gene1==gene2\n")
  naomit_simgic_dt <- naomit_simgic_dt[naomit_simgic_dt$obj_1 != naomit_simgic_dt$obj_2,]
  cat("... as.character obj_1\n")
  naomit_simgic_dt$obj_1 <- as.character(naomit_simgic_dt$obj_1)
  cat("... as.character obj_2\n")
  naomit_simgic_dt$obj_2 <- as.character(naomit_simgic_dt$obj_2)
  head(naomit_simgic_dt)
  # obj_1  obj_2         ss
  # 1: 79501  79501 1.00000000
  # 2: 79501  26155 0.01938086
  
  naomit_simgic_dt$gene1 <- pmin(naomit_simgic_dt$obj_1, naomit_simgic_dt$obj_2)
  naomit_simgic_dt$gene2 <- pmax(naomit_simgic_dt$obj_1, naomit_simgic_dt$obj_2)
  # obj_1  obj_2         ss gene1  gene2
  # # 1: 79501  79501 1.00000000 79501  79501
  # # 2: 79501  26155 0.01938086 26155  79501
  
  cat("... merge simgic_dt with distDT\n")
  naomit_simgic_dt_withDist <- left_join(naomit_simgic_dt, distDT, by = c("gene1", "gene2"))

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
  merge1 <- left_join(naomit_simgic_dt_withDist, g2tDT[,c("entrezID", "region")], by=c("obj_1" = "entrezID"))
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

  ###################################################################
  ###################################################################   SS ~ distance (-> to correct for the # of genes by distance)
  ###################################################################

  simgic_regions_DT$distKb <- simgic_regions_DT$dist/1000

  sameTAD_DT <- simgic_regions_DT[simgic_regions_DT$sameTAD == 1,]
  diffTAD_DT <- simgic_regions_DT[simgic_regions_DT$sameTAD == 0,]

  outFile = file.path(outFold, "sameTAD_DT.Rdata")
  save(sameTAD_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))

  outFile = file.path(outFold, "diffTAD_DT.Rdata")
  save(diffTAD_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  cat(paste0("... load sameTAD_DT\t", Sys.time(), "\t"))
  sameTAD_DT <- eval(parse(text = load( file.path(outFold, "sameTAD_DT.Rdata"))))
  cat(paste0(Sys.time(), "\n"))
      
  cat(paste0("... load diffTAD_DT\t", Sys.time(), "\t"))
  diffTAD_DT <- eval(parse(text = load(file.path(outFold, "diffTAD_DT.Rdata"))))
  cat(paste0(Sys.time(), "\n"))
  
  simgic_regions_DT <- rbind(sameTAD_DT, diffTAD_DT)
  
}




# rather subset random data instead
# binSizeBp <- 100
# binSameTAD_DT <- sameTAD_DT
# binSameTAD_DT$distBin <- binSameTAD_DT$dist %/% binSizeBp
# length(unique(binSameTAD_DT$distBin))
# meanBinSS_sameTAD_DT <- aggregate(ss ~ distBin, data=binSameTAD_DT, FUN=mean, na.rm=T)
# colnames(meanBinSS_sameTAD_DT)[colnames(meanBinSS_sameTAD_DT) == "ss"] <- "meanSS"
# head(meanBinSS_sameTAD_DT)

my_ylab <- paste0("Semantic similarity")
my_xlab <- paste0("Distance between the 2 genes (kb)")
my_sub <- paste0(dataset)

# PREDICT WITH ORIGINAL DISTANCE VALUES
my_xlab <- paste0("Distance between the 2 genes (bp)")

                ### HERE I AM JUST INTERESTED WITH THE DIST VECT
                # # smooth_vals_sameTAD <- predict(loess(ss ~ dist, data = sameTAD_DT), sort(sameTAD_DT$dist))
                # # smooth_vals_diffTAD <- predict(loess(ss ~ dist, data = diffTAD_DT), sort(diffTAD_DT$dist))
                # smooth_vals_sameTAD <- predict(loess(ss ~ dist, data = sameTAD_DT), se =TRUE)
                # smooth_vals_diffTAD <- predict(loess(ss ~ dist, data = diffTAD_DT), se =TRUE)
                # 
                # 
                # auc_diffTAD_obsDist <- auc(x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD)
                # auc_sameTAD_obsDist <- auc(x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD)
                # 
                # 
                # outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_originalDist", ".", plotType))
                # do.call(plotType, list(outFile, height = myHeight, width = myWidth))
                # plot(NULL,
                #      xlim = range(simgic_regions_DT$dist), 
                #      ylim = range(c(smooth_vals_sameTAD, smooth_vals_diffTAD)),
                #      # xlab="", 
                #      # ylab="",
                #      xlab=my_xlab, 
                #      ylab=my_ylab,
                #      main=paste0(dataset, ": ss ~ dist loess fit"))
                # mtext(text = "observed distance values", side = 3)
                # # lines( x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD, col = sameTADcol)
                # # lines( x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD, col = diffTADcol)
                # 
                # lines(sort(sameTAD_DT$dist),smooth_vals_sameTAD$fit, col = sameTADcol)
                # lines(sameTAD_DT$dist,smooth_vals_sameTAD$fit+2*smooth_vals_sameTAD$s, lty=2, col = sameTADcol) #rough & ready CI
                # lines(sameTAD_DT$dist,smooth_vals_sameTAD$fit-2*smooth_vals_sameTAD$s, lty=2, col = sameTADcol)
                # 
                # 
                # lines(sort(diffTAD_DT$dist),smooth_vals_diffTAD$fit)
                # lines(diffTAD_DT$dist,smooth_vals_diffTAD$fit+2*smooth_vals_diffTAD$s, lty=2, col = diffTADcol) #rough & ready CI
                # lines(diffTAD_DT$dist,smooth_vals_diffTAD$fit-2*smooth_vals_diffTAD$s, lty=2, col = diffTADcol)
                # 
                # 
                # legend("topright", 
                #        legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_obsDist, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_obsDist, 2))),  
                #        col = c(sameTADcol, diffTADcol),
                #        lty=1,
                #        bty = "n")
                # 
                # foo <- dev.off()
                # cat(paste0("... written: ", outFile, "\n"))
if(buildModel){
sameTAD_loess <- loess(ss ~ dist, data = sameTAD_DT)
outFile <- file.path(outFold, "sameTAD_loess.Rdata")
save(sameTAD_loess, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

diffTAD_loess <- loess(ss ~ dist, data = diffTAD_DT)
outFile <- file.path(outFold, "diffTAD_loess.Rdata")
save(diffTAD_loess, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


} else {
outFile <- file.path(outFold, "sameTAD_loess.Rdata")
cat("...load sameTAD_loess\n")
sameTAD_loess = eval(parse(text = load(outFile)))

outFile <- file.path(outFold, "diffTAD_loess.Rdata")
cat("...load diffTAD_loess\n")
diffTAD_loess = eval(parse(text = load(outFile)))
}

# REDO THE MODELS WITH A SUBSET OF DATA
#stopifnot( nrow(sameTAD_DT) == nrow(diffTAD_DT) )
sampIdxSame <- sample(x = 1:nrow(sameTAD_DT), size=nSamp_shortModel, replace = FALSE)
sampIdxDiff <- sample(x = 1:nrow(diffTAD_DT), size=nSamp_shortModel, replace = FALSE)


sameTAD_DT_short <- sameTAD_DT[sampIdxSame,]
sameTAD_loess_short <- loess(ss ~ dist, data = sameTAD_DT_short)

diffTAD_DT_short <- diffTAD_DT[sampIdxDiff,]
diffTAD_loess_short <- loess(ss ~ dist, data = diffTAD_DT_short)


### WITH THE FULL DATA
# smooth_vals_sameTAD_distVect <- predict(object=sameTAD_loess, 
#                                         newdata=data.frame(dist=distVect), 
#                                         se =TRUE)
# smooth_vals_diffTAD_distVect <- predict(object=diffTAD_loess, 
#                                         newdata=data.frame(dist=distVect), 
#                                         se =TRUE)
# => this does not work because memory issue ! (even if distVect is a small vector)
# see ?predict.loess => use predict on a subset of data!
smooth_vals_sameTAD_distVect <- predict(object=sameTAD_loess,
                                        newdata=data.frame(dist=distVect),
                                        se =FALSE)
smooth_vals_diffTAD_distVect <- predict(object=diffTAD_loess,
                                        newdata=data.frame(dist=distVect),
                                        se =FALSE)

### WITH THE SUBSETTED DATA
smooth_vals_sameTAD_short_distVect <- predict(object=sameTAD_loess_short,
                                                newdata=data.frame(dist=distVect),
                                                se =TRUE)
sameTAD_short_fit <- smooth_vals_sameTAD_short_distVect$fit
sameTAD_short_fit_CIlow <- sameTAD_short_fit - smooth_vals_sameTAD_short_distVect$se.fit*qnorm(1-.05/2)
sameTAD_short_fit_CIhigh <- sameTAD_short_fit + smooth_vals_sameTAD_short_distVect$se.fit*qnorm(1-.05/2)

smooth_vals_diffTAD_short_distVect <- predict(object=diffTAD_loess_short,
                                              newdata=data.frame(dist=distVect),
                                              se =TRUE)
diffTAD_short_fit <- smooth_vals_diffTAD_short_distVect$fit
diffTAD_short_fit_CIlow <- diffTAD_short_fit - smooth_vals_diffTAD_short_distVect$se.fit*qnorm(1-.05/2)
diffTAD_short_fit_CIhigh <- diffTAD_short_fit + smooth_vals_diffTAD_short_distVect$se.fit*qnorm(1-.05/2)

# plot  1) the fit, without CI with all data points
#       2) the fit with CI with the subset of data points

outFile <- file.path(outFold, "smooth_vals_sameTAD_distVect.Rdata")
save(smooth_vals_sameTAD_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "smooth_vals_diffTAD_distVect.Rdata")
save(smooth_vals_diffTAD_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# smooth_vals_sameTAD_short_distVect <- smooth_vals_sameTAD_short_distVect$fit  # because was computed with se=TRUE
outFile <- file.path(outFold, "smooth_vals_sameTAD_short_distVect.Rdata")
save(smooth_vals_sameTAD_short_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# smooth_vals_diffTAD_short_distVect <- smooth_vals_diffTAD_short_distVect$fit
outFile <- file.path(outFold, "smooth_vals_diffTAD_short_distVect.Rdata")
save(smooth_vals_diffTAD_short_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "distVect.Rdata")
save(distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

### compute and save AUC to compare the loess with data subset
cat("auc1\n")
auc_diffTAD_distVect <- auc(x = distVect, y = smooth_vals_diffTAD_distVect)
cat("auc2\n")
auc_sameTAD_distVect <- auc(x = distVect, y = smooth_vals_sameTAD_distVect)
cat("auc3\n")
auc_diffTAD_short_distVect <- auc(x = distVect, y = smooth_vals_diffTAD_short_distVect$fit)
cat("auc4\n")
auc_sameTAD_short_distVect <- auc(x = distVect, y = smooth_vals_sameTAD_short_distVect$fit)

outFile <- file.path(outFold, "auc_diffTAD_distVect.Rdata")
save(auc_diffTAD_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "auc_sameTAD_distVect.Rdata")
save(auc_sameTAD_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "auc_diffTAD_short_distVect.Rdata")
save(auc_diffTAD_short_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, "auc_sameTAD_short_distVect.Rdata")
save(auc_sameTAD_short_distVect, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


###################################################################
### PLOT ONLY THE FITTED DATA WITH THE FULL DATA                        - LOESS - FULL DATA -NO CI
###################################################################

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_vectDist.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(NULL,
     xlim = range(distVect), 
     ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), na.omit(smooth_vals_diffTAD_distVect))),
     xlab=my_xlab,
     ylab=my_ylab,
     main=paste0(dataset, ": SS ~ dist loess fit"))
mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)

lines(distVect,smooth_vals_sameTAD_distVect, col = sameTADcol)

lines(distVect,smooth_vals_diffTAD_distVect, col = diffTADcol)
# lines(diffTAD_DT$dist,smooth_vals_diffTAD_distVect$fit+2*smooth_vals_diffTAD_distVect$s, lty=2, col = diffTADcol) #rough & ready CI
# lines(diffTAD_DT$dist,smooth_vals_diffTAD_distVect$fit-2*smooth_vals_diffTAD_distVect$s, lty=2, col = diffTADcol)

legend("topright", 
       legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2))), 
       col = c(sameTADcol, diffTADcol),
       lty=1,
       bty = "n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

###################################################################
### PLOT ONLY THE FITTED DATA WITH THE FULL DATA                        - LOESS - DATA SUBSET -WITH CI
###################################################################

mytit <- paste0(dataset, ": SS ~ dist loess fit (subset: ", nSamp_shortModel, "/", nrow(sameTAD_DT), ")")

outFile <- file.path(outFold, paste0("sameTAD_diffTAD_subset_loessFit_withCI_vectDist.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(NULL,
     xlim = range(distVect), 
     ylim = range(c(diffTAD_short_fit_CIlow,diffTAD_short_fit_CIhigh,
                    sameTAD_short_fit_CIlow,sameTAD_short_fit_CIhigh ), na.rm=TRUE),
     xlab=my_xlab,
     ylab=my_ylab,
     main=mytit)
mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)

lines(distVect, sameTAD_short_fit, col = sameTADcol)
lines(distVect, sameTAD_short_fit_CIlow, col = sameTADcol, lty=2)
lines(distVect, sameTAD_short_fit_CIhigh, col = sameTADcol,lty=2)

lines(distVect, diffTAD_short_fit, col = sameTADcol)
lines(distVect, diffTAD_short_fit_CIlow, col = sameTADcol,lty=2)
lines(distVect, diffTAD_short_fit_CIhigh, col = sameTADcol,lty=2)

legend("topright", 
       legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_short_distVect, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_short_distVect, 2))), 
       col = c(sameTADcol, diffTADcol),
       lty=1,
       bty = "n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("*** Done\n", startTime, "\n", Sys.time(), "\n"))

stop("-- ok \n")

###################################################################
###################################################################   density
###################################################################

cat("draw multidens\n")
outFile <- file.path(outFold, paste0("ss_density_sameTAD_diffTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  sameTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 1],
  diffTAD = simgic_regions_DT$ss[simgic_regions_DT$sameTAD == 0]
))
foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))

###################################################################
###################################################################   aggregate meanSS
###################################################################

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


###################################################################
###################################################################   meanSS ~ total # of gene pairs
###################################################################


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

# plot(meanSS ~ nbrPairs, 
#      data = dataDT[dataDT$sameTAD == 0,],
#      pch=16,
#      cex=1.2,
#      main="mean SS vs. # gene pairs (diffTAD)",
#      xlab="# pairs", ylab="meanSS",
#      cex.axis=1.2,cex.lab=1.2)
outFile <- file.path(outFold, paste0("meanSS_nbrGenePairs_log10_diffTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(meanSS ~ nbrPairs_log10, 
     data = dataDT[dataDT$sameTAD == 0,],
     pch=16,
     cex=2,
     main="mean SS vs. # gene pairs (diffTAD)",
     xlab="# pairs [log10]", ylab="meanSS",
     cex.axis=2,cex.lab=2,cex.main=2)
foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))

# plot(meanSS ~ nbrPairs, 
#      data = dataDT[dataDT$sameTAD == 1,],
#      pch=16,
#      cex=1.2,
#      main="mean SS vs. # gene pairs (sameTAD)",
#      xlab="# pairs", ylab="meanSS",
#      cex.axis=1.2,cex.lab=1.2)

outFile <- file.path(outFold, paste0("meanSS_nbrGenePairs_log10_sameTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(meanSS ~ nbrPairs_log10, 
     data = dataDT[dataDT$sameTAD == 1,],
     pch=16,
     cex=2,
     main="mean SS vs. # gene pairs (sameTAD)",
     xlab="# pairs [log10]", ylab="meanSS",
     cex.axis=2,cex.lab=2,cex.main=2)
foo <- dev.off()
cat(paste0("written: ", outFile, "\n"))


# very different number of genes !!! 
# need to correct for distance





