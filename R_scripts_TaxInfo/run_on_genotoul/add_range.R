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
# Expand taxa columns
# Corrections :
# -> k or d for the kingdom
# -> exclude (Fungi) at the end of Genus name
expand_taxnames <- function(OTU_table, taxonomy) {
  OTU_table <- OTU_table |>
    mutate(
      Kingdom = str_match(taxonomy, "[kd]:([^, \n]+)")[, 2],
      Phylum = str_match(taxonomy, "p:([^, \n]+)")[, 2],
      Class = str_match(taxonomy, "c:([^, \n]+)")[, 2],
      Order = str_match(taxonomy, "o:([^, \n]+)")[, 2],
      Family = str_match(taxonomy, "f:([^, \n]+)")[, 2],
      Genus = str_match(taxonomy, "g:([^, \n(]+)")[, 2],
      Species = str_match(taxonomy, "s:([^, \n]+)")[, 2]
    )
  return(OTU_table)
}

# Get a phyloseq object from the cean tables with the taxa
phyloseq_from_table <- function(OTU_table) {
  ## Get the OTU table (as a matrix for the phyloseq object)
  otu_matrix <- OTU_table |>
    select(OTU, where(is.numeric)) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix()

  ## here OTU names are the row so we have to set taxa_are_rows to TRUE
  OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)

  ## Get the Taxa matrix
  tax_matrix <- OTU_table |>
    select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix()

  TAX <- tax_table(tax_matrix)

  PHYLOSEQ <- phyloseq(OTU, TAX)
  return(PHYLOSEQ)
}
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

# getwd()
# setwd("Database")
# INPUT_PQ <- "Database/cleaned_OTU_tables"
# INPUT_TRAITS <- "Database/traitsTable"
# OUTPUT_PLOT <- "Database/plot_taxInfo_traits"
# OUTPUT_PQ <- "Database/phyloseq_traits"
methods <- c("Tedersoo", "Illumina", "Getplage")
sections_Getplage <- c("ITS1", "ITS_full", "ITS_none")
sections_Tedersoo <- c("ITS_full", "ITS1")
sections_Illumina <- c("ITS1")
clusters <- c("vsearch", "sintax")
# method <- "Tedersoo"
# section <- "ITS1"
# cluster <- "vsearch"
#
# OTU_table <- read.csv(
#   "cleaned_OTU_tables/OTU_table_Tedersoo_ITS_full_vsearch_fungi_clean.csv",
#   header = TRUE,
#   sep = ";"
# )
#
# OTU_table <- expand_taxnames(OTU_table, taxonomy)
#
# # Then get a phyloseq object to use the taxInfo functions
# data <- phyloseq_from_table(OTU_table)
#
# # Check names according to Taxref (210) conventions
# data_clean <- gna_verifier_pq(data, data_sources = 210)
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

      sample_df <- data.frame(
        row.names = sample_names(data_clean),
        site = substr(sample_names(data_clean), 1, 1)
      )
      #   latitude = case_when(
      #     substr(sample_names(data_clean), 1, 1) == "C" ~ 44.63095,
      #     substr(sample_names(data_clean), 1, 1) == "M" ~ 42.48102,
      #     substr(sample_names(data_clean), 1, 1) == "V" ~ 47.95948
      #   ),
      #   longitude = case_when(
      #     substr(sample_names(data_clean), 1, 1) == "C" ~ 6.00120,
      #     substr(sample_names(data_clean), 1, 1) == "M" ~ 3.02871,
      #     substr(sample_names(data_clean), 1, 1) == "V" ~ 6.92549
      #   )
      # )
      sample_data(data_clean) <- sample_df
      ps_C <- subset_samples(data_clean, site == "C")
      ps_M <- subset_samples(data_clean, site == "M")
      ps_V <- subset_samples(data_clean, site == "V")

      # test_merge <- merge_phyloseq(ps_C, ps_M, ps_V)

      # sample_data(data_clean)
      ps_C_check <- tax_occur_check_pq(
        ps_C,
        longitude = 6.00120,
        latitude = 44.63095,
        radius_km = 100,
        n_occur = 1000,
      )
      ps_M_check <- tax_occur_check_pq(
        ps_M,
        longitude = 3.02871,
        latitude = 42.48102,
        radius_km = 100,
        n_occur = 1000,
      )
      ps_V_check <- tax_occur_check_pq(
        ps_V,
        longitude = 6.92549,
        latitude = 47.95948,
        radius_km = 100,
        n_occur = 1000,
      )
      write_pq(
        ps_C_check,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_gbif_check/C"
        )
      )
      write_pq(
        ps_M_check,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_gbif_check/M"
        )
      )
      write_pq(
        ps_V_check,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_gbif_check/V"
        )
      )
    }
  }
}
