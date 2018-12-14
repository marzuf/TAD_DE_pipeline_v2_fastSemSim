specific_IC_table_file <- file.path("DATASETS_TADs_IC_TABLES",
"TCGAcoad_msi_mss", "chr6_TAD58_GeneOntology_biological_process_SimGIC_max_IC_table.txt")
stopifnot(file.exists(specific_IC_table_file))

specific_results_file <- file.path("cmp_IC_table",
                              "TCGAcoad_msi_mss/output",
                              "chr6_TAD58_GeneOntology_biological_process_SimGIC_max_result_file_dataset_tad_specific_IC_table.txt")
stopifnot(file.exists(specific_results_file))

all_IC_table_file <- file.path("ALL_GENES_IC_TABLE",
                          "output",
                          "all_datasets_GeneOntology_biological_process_SimGIC_max_IC_table.txt")
stopifnot(file.exists(all_IC_table_file))

all_results_file <- file.path("cmp_IC_table",
                         "TCGAcoad_msi_mss",
                         "output",
                         "chr6_TAD58_GeneOntology_biological_process_SimGIC_max_result_file_all_datasets_IC_table.txt")
stopifnot(file.exists(all_results_file))


specific_IC_table_DT <- read.delim(specific_IC_table_file, header=F, stringsAsFactors = FALSE, col.names=c("GO", "SS"))
head(specific_IC_table_DT)
specific_IC_table_DT$SS <- ifelse(as.character(specific_IC_table_DT$SS) == "None", NA, as.numeric(as.character(specific_IC_table_DT$SS)))

specific_results_DT <- read.delim(specific_results_file, header=T, stringsAsFactors = FALSE)
head(specific_results_DT)

stopifnot(is.numeric(na.omit(specific_results_DT$ss)))

all_IC_table_DT <- read.delim(all_IC_table_file, header=F, stringsAsFactors = FALSE, col.names=c("GO", "SS"))
head(all_IC_table_DT)
all_IC_table_DT$SS <- ifelse(as.character(all_IC_table_DT$SS) == "None", NA, as.numeric(as.character(all_IC_table_DT$SS)))

all_results_DT <- read.delim(all_results_file, header=T, stringsAsFactors = FALSE)
head(all_results_DT)

stopifnot(is.numeric(na.omit(all_results_DT$ss)))

stopifnot(all_results_DT$obj_1 == specific_results_DT$obj_1)
stopifnot(all_results_DT$obj_2 == specific_results_DT$obj_2)

plot(x=all_results_DT$ss,
     y= specific_results_DT$ss)
cor.test(x=all_results_DT$ss,
         y= specific_results_DT$ss)


stopifnot(specific_IC_table_DT$GO == all_IC_table_DT$GO)
stopifnot(na.omit(specific_IC_table_DT$SS) == na.omit(all_IC_table_DT$SS))

plot(x=specific_IC_table_DT$SS,
     y= all_IC_table_DT$SS)
cor.test(x=specific_IC_table_DT$SS,
         y= all_IC_table_DT$SS)