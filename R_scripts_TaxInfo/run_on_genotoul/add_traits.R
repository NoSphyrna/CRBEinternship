# ============================================================================ #
# ======================= Script to plot trophic abundance =================== #
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

if (length(args) != 3) {
  stop(
    "Usage : Rscript add_traits.R INPUT_PQ OUTPUT_PQ TRAITS_TABLE"
  )
}
# INPUT_DIR <- "~/Dext2/BioInfoMatser/HAU802_BILL/ONT_sequencing/Analyses/SNP/Isec/"
INPUT_PQ <- args[1]
OUTPUT_PQ <- args[2]
TRAITS_TABLE <- args[3]

if (!dir.exists(INPUT_PQ)) {
  stop("[R] : Input PQ folder doesn't exist")
}

if (!file.exists(TRAITS_TABLE)) {
  stop("[R] : Traits tabel doesn't exist")
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
# Normalise FungalTRaits trophicMode (adapted from TaxInfo)
# This reduces the traits to three categories to match the trophicMode from FunGuild

ft_to_trophic_mode <- function(x) {
  dplyr::case_when(
    x %in%
      c(
        "algal_decomposer", # new
        "animal_decomposer", # new
        "dung_saprotroph",
        "fungal_decomposer", # new
        "litter_saprotroph",
        "myxomycete_decomposer", # new
        "nectar/tap_saprotroph", # new
        "pollen_saprotroph",
        "resin_saprotroph", # new
        "rock-inhabiting", # new
        "soil_saprotroph",
        "unsepcified_saprotroph", # new
        "unspecified_saprotroph",
        "wood_saprotroph"
      ) ~ "Saprotroph",
    x %in%
      c(
        "algal_parasite",
        "algivorous/protistivorous", # new
        "animal_parasite",
        "arthropod_parasite", # new
        "bacterivorous", # new
        "bryophilous", # new
        "fish_parasite", # new
        "invertebrate_parasite", # new
        "lichen_parasite",
        "moss_parasite", # new
        "mycoparasite",
        "nematophagous", # new
        "plant_pathogen",
        "protistan_parasite",
        "sooty_mold",
        "unspecified_pathotroph"
      ) ~ "Pathotroph",
    x %in%
      c(
        "algal_ectosymbiont", # new
        "algal_symbiont", # new
        "animal-associated",
        "animal_endosymbiont",
        "arbuscular_mycorrhizal",
        "arthropod-associated",
        "coral-associated", # new
        "ectomycorrhizal",
        "epiphyte",
        "ericoid_mycorrhizal", # new
        "foliar_endophyte",
        "insect-associated", # new
        "invertebrate-associated", # new
        "lichenized",
        "liverwort-associated", # new
        "moss_symbiont",
        "root-associated", # new
        "root_endophyte",
        "root_endophyte_dark_septate", # new
        "termite_symbiont", # new
        "unspecified_symbiotroph",
        "vertebrate-associated" # new
      ) ~ "Symbiotroph",
    is.na(x) |
      x %in% c("unspecified", "", "0", "fatty_acid_producer") ~ NA_character_, # new
    .default = "Other"
  )
}

# Add column trophicMode to a table with traits assigned by fungaltraits enhanced table
add_trophicMode_ft <- function(table_traits) {
  table_traits <- table_traits |>
    mutate(
      tmp_primary_troph = ft_to_trophic_mode(ft_primary_lifestyle),
      tmp_secondary_troph = ft_to_trophic_mode(ft_Secondary_lifestyle),
      ft_trophicMode = ifelse(
        is.na(tmp_primary_troph) |
          is.na(tmp_secondary_troph) |
          tmp_primary_troph == tmp_secondary_troph,
        coalesce(tmp_primary_troph, tmp_secondary_troph),
        paste(tmp_primary_troph, tmp_secondary_troph, sep = "-")
      )
    ) |>
    select(-c(tmp_primary_troph, tmp_secondary_troph))
  return(table_traits)
}


# =============== Workflow for assigning traits on OTUs from samples ================ #

# fungal_traits_file <- read.table(
#   paste0(INPUT_TRAITS, "/FUNGALT_DB_MROY041125.csv"),
#   header = TRUE,
#   sep = ";"
# )
#
# fungal_traits_file <- fungal_traits_file |>
#   mutate(
#     tmp_primary_troph = ft_to_trophic_mode(primary_lifestyle),
#     tmp_secondary_troph = ft_to_trophic_mode(Secondary_lifestyle),
#     trophicMode = ifelse(
#       is.na(tmp_primary_troph) |
#         is.na(tmp_secondary_troph) |
#         tmp_primary_troph == tmp_secondary_troph,
#       coalesce(tmp_primary_troph, tmp_secondary_troph),
#       paste(tmp_primary_troph, tmp_secondary_troph, sep = "-")
#     )
#   ) |>
#   select(-c(tmp_primary_troph, tmp_secondary_troph))
#
# unique(fungal_traits_file$trophicMode)
#
# other_values <- unique(fungal_traits_file$primary_lifestyle[
#   ft_to_trophic_mode(fungal_traits_file$primary_lifestyle) == "Other"
# ])
# print(sort(other_values))

#INTPUT_PQ <- "Database/phyloseq_traits"
# INPUT_TRAITS <- "Database/traitsTable"
# OUTPUT_PLOT <- "Database/plot_taxInfo_traits"
#OUTPUT_PQ <- "Database/phyloseq_traits"
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

      data_traits <- fungal_traits_guilds(
        data_clean,
        fungal_traits_file = paste0(TRAITS_TABLE),
        ft_taxonomic_rank = "genusEpithet",
        ft_csv_rank = "GENUS",
        ft_sep = ";",
        ft_col_prefix = "ft_",
        fg_tax_levels = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        ),
        fg_col_prefix = "fg_",
        db_url = "http://www.stbates.org/funguild_db_2.php",
        add_consensus = TRUE,
        consensus_col_prefix = "cons_",
        add_to_phyloseq = TRUE,
        verbose = TRUE
      )

      # Save the phyloseq object
      write_pq(
        data_traits,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_FungalTraits_enhanced_FunGuild"
        )
      )
    }
  }
}
