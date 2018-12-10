# gene_list => vector of genes query for the SemSim
# runName => name used for naming file
# icData => possibly a path to existing table
# fastSemSimParams => list of parameters for  FastSemSim


fss_exec <- "fastsemsim"

run_fastSemSim <- function(
geneList,
runName,
simMetric,
gafFile,
otherFastSemSimParams=NULL,
icData=NULL, 
runFolder = file.path(".", tempdir()),
rmTmp = FALSE
) {

  ################################ some default params - overwritten if provided in otherFastSemSimParams
default_root <- "biological_process"
default_species <- "human"
default_ontType <- "GeneOntology"
default_queryType <- "obj"
default_mixStrategy <- "max"
default_queryIn <- "file"
default_queryMode <- "list"
  ################################

  cat("... start run_fastSemSim\n")

  fastSemSim_queryFile <- file.path(runFolder, "input", paste0(runName, "_query_file.txt") )

  fastSemSim_resultFile <- file.path(runFolder, "output", paste0(runName, "_result_file.txt") )

  fastSemSimParams_default <- list(
                          "--ontology_type" =  default_ontType, 
                          "--ac_species" = default_species,
                          "--query_ss_type" = default_queryType, 
                          "--tss" = simMetric,
                          "--tmix" = default_mixStrategy,
                          "--query_input" = default_queryIn, 
                          "--query_file"= fastSemSim_queryFile, 
                          "--output_file" =  fastSemSim_resultFile,
                          "--task" = "SS", 
                          "-vv",  
                          "--query_mode"= default_queryMode,
                          "--root" = default_root
  )
  
  if(!is.null(icData)) fastSemSimParams_default[["--inject_IC"]] <- icData
  
  # update if provide other settings:
  if(!is.null(otherFastSemSimParams)) {
    keep_default <- fastSemSimParams_default[!names(fastSemSimParams_default) %in% names(otherFastSemSimParams)]
    fastSemSimParams <- c(keep_default, otherFastSemSimParams)
  } else {
    fastSemSimParams <- fastSemSimParams_default 
  }
  rm(fastSemSimParams_default)
  argsTxt <- paste0(paste(names(fastSemSimParams), "\t=\t", unlist(fastSemSimParams)), collapse="\n")
  cat("... found following parameters for FastSemSim:\n")
  cat(argsTxt, "\n")
  
  argsCmd <- paste0(paste(names(fastSemSimParams), unlist(fastSemSimParams)), collapse=" ")
  
  # in case of reset by function call
  fastSemSim_queryFile <- fastSemSimParams[["--query_file"]]
  fastSemSim_resultFile <- fastSemSimParams[["--output_file"]]
  
  system(paste0("mkdir -p ", dirname(fastSemSim_queryFile)))
  system(paste0("mkdir -p ", dirname(fastSemSim_resultFile)))
  
  
  fileConn <- file(fastSemSim_queryFile)    
  writeLines(as.character(geneList), fileConn)    
  close(fileConn)
  cat(paste0("... written:\t", fastSemSim_queryFile, "\n"))
  
  #************* RUN FAST SEM SIM 
  cat(paste0("... start fastSemSim:\t", Sys.time(), "\n"))
  fastSemSim_cmd <- paste(fss_exec, argsCmd)
  cat(fastSemSim_cmd, "\n")
  system(fastSemSim_cmd)
  cat(paste0("... end fastSemSim:\t", Sys.time(), "\n"))

  stopifnot(file.exists(fastSemSim_resultFile))
  cat(paste0("... fastSemSim result file:\t", fastSemSim_resultFile, "\n"))

  
  #*************REFORMAT FAST SEM SIM OUTPUT
  semsim_DT <- read.delim(fastSemSim_resultFile, header=T, stringsAsFactors = FALSE)
  
  if(rmTmp) {
    cat("... make clean: ", paste0("rm -rf  ", runFolder), "\n")
    system(paste0("rm -rf  ", runFolder))
  } else {
    cat("... temp files stored in: ", paste0(runFolder), "\n")
  }
  
  return(semsim_DT)
  
}






