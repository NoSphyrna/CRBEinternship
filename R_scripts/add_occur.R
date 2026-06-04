# ============================================================================ #
# ===================== Script to add occurences by country ================== #
# ============================================================================ #

# ============== Set Workspace accordingly ======== #
# Get input directory and output directory :

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage : Rscript add_occur.R INPUT_PQ OUTPUT_PQ"
  )
}
INPUT_PQ <- args[1]
OUTPUT_PQ <- args[2]

if (!dir.exists(INPUT_PQ)) {
  stop("[R] : Input OTU folder doesn't exist")
}

if (!dir.exists(OUTPUT_PQ)) {
  dir.create(OUTPUT_PQ, recursive = TRUE)
}
# ============== Libraries ============ #
library(taxinfo)
library(MiscMetabar)
library(stringr)


# ============= Utility Functoins ==================#

# Get a phyloseq object from registered files with write_pq
read_pq_corr <- function(path) {
  OTU_table <- read.csv(
    paste0(path, "/otu_table.csv"),
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )
  TAX_table <- read.csv(
    paste0(
      path,
      "/tax_table.csv"
    ),
    header = TRUE,
    row.names = 1,
    sep = "\t"
  ) |>
    mutate_all(as.character) # taxa needs to be char

  ## Get the OTU table (as a matrix for the phyloseq object)
  OTU <- otu_table(as.matrix(OTU_table), taxa_are_rows = TRUE)

  ## Get the TAX table (as a matrix for the phyloseq object)
  TAX <- tax_table(
    TAX_table |>
      mutate_all(as.character) |>
      as.matrix()
  )

  PHYLOSEQ <- phyloseq(OTU, TAX)

  return(PHYLOSEQ)
}
# =============== Workflow for verifying names on OTUs from samples ================ #

# INPUT_PQ <- "Database/cleaned_OTU_tables"
# INPUT_TRAITS <- "Database/traitsTable"
# OUTPUT_PLOT <- "Database/plot_taxInfo_traits"
# OUTPUT_PQ <- "Database/phyloseq_traits"
methods <- c("Illumina", "Tedersoo", "Getplage")
sections_Getplage <- c("ITS1", "ITS_full", "ITS_none")
sections_Tedersoo <- c("ITS1", "ITS_full")
sections_Illumina <- c("ITS1")
clusters <- c("sintax", "vsearch")
# method <- "Tedersoo"
# section <- "ITS1"
# cluster <- "vsearch"
#

for (method in methods) {
  switch(
    method,
    "Getplage" = sections <- sections_Getplage,
    "Tedersoo" = sections <- sections_Tedersoo,
    "Illumina" = sections <- sections_Illumina
  )
  for (section in sections) {
    for (cluster in clusters) {
      cat("\n# ========================== #\n")
      cat("# Method | Section | Cluster #\n")
      cat("#", method, "|", section, "|", cluster, "#\n")
      cat("# ========================== #\n")

      # Save the phyloseq object
      data_clean <- read_pq_corr(
        paste0(
          INPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_verified"
        )
      )

      data_clean_occ <- tax_gbif_occur_pq(data_clean, by_country = TRUE)

      write_pq(
        data_clean_occ,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_gbif_occ_by_country"
        )
      )
    }
  }
}
