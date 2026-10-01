# ============================================================================ #
# ===================== Script to verify names with taxref =================== #
# ============================================================================ #

# ============== Installation of packages =================== #

# devtools::install_github is deprecated use pak::pak instead

# TaxInfo
# pak::pak("adrientaudiere/MiscMetabar")
# pak::pak("adrientaudiere/taxinfo")
#
# # FugalTraits
# pak::pak("ropenscilabs/datastorr")
# pak::pak("traitecoevo/fungaltraits")

# ============== Set Workspace accordingly ======== #
# Get input directory and output directory :

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage : Rscript verify_gna_from_clean_data.R INPUT_OTU_DIRECTORY OUTPUT_PHYLOSEQ_DIRECTORY"
  )
}
INPUT_OTU <- args[1]
OUTPUT_PQ <- args[2]

if (!dir.exists(INPUT_OTU)) {
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

# To treat Illumina tables as the pacbio tables
standardise_illumina <- function(OTU_table) {
  OTU_table <- OTU_table |>
    select(
      !c(
        abundance,
        length,
        chimera,
        spread,
        identity,
        V5,
        Genus_species,
        proba_last
      )
    ) |>
    rename(OTU = amplicon)
  return(OTU_table)
}

# Expand taxa columns
# Corrections :
# -> k or d for the kingdom
# -> exclude (Fungi) at the end of Genus name
expand_taxnames <- function(OTU_table) {
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
# =============== Workflow for verifying names on OTUs from samples ================ #

# INPUT_OTU <- "Database/cleaned_OTU_tables"
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
      # Get the OTU_table
      OTU_table <- read.csv(
        paste0(
          INPUT_OTU,
          "/OTU_table_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_fungi_clean.csv"
        ),
        header = TRUE,
        sep = ";"
      )

      if (method == "Illumina") {
        OTU_table <- standardise_illumina(OTU_table)
      }

      # First expand the taxnames to match phyloseq conventions
      OTU_table <- expand_taxnames(OTU_table)

      # Then get a phyloseq object to use the taxInfo functions
      data <- phyloseq_from_table(OTU_table)

      # Check names according to Taxref (210) conventions
      data_clean <- gna_verifier_pq(data, data_sources = 210)

      # Save the phyloseq object
      write_pq(
        data_clean,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_verified"
        )
      )
    }
  }
}
