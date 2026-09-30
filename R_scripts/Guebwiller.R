# ================================================================================== #
# ====== Script to add traits to table OTUs with FunGuild and Fungal Traits ======== #
# ================================================================================== #

# ============================ Installation of packages ============================ #

# <NOTE> devtools::install_github is deprecated use pak::pak instead

## To install the package taxinfo uncomment and execute the following:
# pak::pak("adrientaudiere/MiscMetabar")
pak::pak("adrientaudiere/taxinfo")

# <NOTE> You might need to install pak before with the command:
# install.packages("pak")

## Installation of stringr (to manipulate strings)
# install.packages("stringr")

# ==================================== Settings ==================================== #
# Here are the default settings with your files and directory locations

# Uncomment and adapt the following to make sure you are in the right directory
# if you run this cript not from command lines
getwd()
setwd("../../Database/Guebwiller/")

INPUT_TABLE <- "taxonomy_MITI2_OTU97_ITS_full_vsearch - Guebwiller_species.csv"
FUNGAL_TRAITS_TABLE <- "../traitsTable/FUNGALT_DB_MROY041125.csv"

# =================================== Libraries ==================================== #
library(taxinfo)
library(stringr)

# Check loaded libraries:
# (.packages())
# =============================== Utility Functions ================================ #

# Normalise FungalTraits trophicMode (adapted from TaxInfo)
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
add_trophicMode_ft <- function(data_traits) {
  table_traits <- as.data.frame(tax_table(data_traits))
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
  tax_table(data_traits) <- tax_table(as.matrix(table_traits))
  return(data_traits)
}


# Allows to clean the column guild from FunGuild traits
# directly from a phyloseq object
# e.g. "Endophyte-|Plant Pathogen|" -> "endophyte-plant_pathogen"
clean_guild <- function(data_traits) {
  table_traits <- as.data.frame(tax_table(data_traits))
  new_table <- table_traits |>
    mutate(
      fg_guild = fg_guild |>
        str_squish() |>
        str_to_lower() |>
        str_replace_all(" ", "_") |>
        str_replace_all("\\|", "")
    )
  tax_table(data_traits) <- tax_table(as.matrix(new_table))
  return(data_traits)
}

# Here with this function we remove unecessary spaces in troph column
clean_troph <- function(data_traits) {
  table_traits <- as.data.frame(tax_table(data_traits))
  new_table <- table_traits |>
    mutate(
      fg_trophicMode = fg_trophicMode |>
        str_squish()
    )
  tax_table(data_traits) <- tax_table(as.matrix(new_table))
  return(data_traits)
}

# =============== Workflow for assigning traits on OTUs from samples ================ #

cat("Getting table\n")
# Get the OTU_table
TAX_table <- read.csv(
  INPUT_TABLE,
  header = TRUE,
  sep = ","
)

# Clean the taxo columns
# We also remove the (Fungi) from Lactarius(Fungi)
TAX_table <- TAX_table |>
  mutate(
    Kingdom = str_match(Kingdom, "[kd]__([^; \n]+)")[, 2],
    Phylum = str_match(Phylum, "p__([^; \n]+)")[, 2],
    Class = str_match(Class, "c__([^; \n]+)")[, 2],
    Order = str_match(Order, "o__([^; \n]+)")[, 2],
    Family = str_match(Family, "f__([^; \n]+)")[, 2],
    Genus = str_match(Genus, "g__([^; \n(]+)")[, 2],
    Species = str_match(Species, "s__([^; \n]+)")[, 2]
  )


# Create the phyloseq object
TAX_pq <- TAX_table |>
  tibble::column_to_rownames("otu97") |>
  as.matrix() |>
  tax_table()

OTU_dummy <- TAX_table |>
  select(otu97) |>
  mutate(
    dummy_sample = as.integer(!is.na(otu97))
  ) |>
  tibble::column_to_rownames("otu97") |>
  as.matrix() |>
  otu_table(taxa_are_rows = TRUE)

data_pq <- phyloseq(OTU_dummy, TAX_pq)


cat("Check names\n")
# Check names according to Taxref (210) conventions
data_clean <- gna_verifier_pq(data_pq, data_sources = 210)
# Warning: program compiled against libxml 210 using older 209
# ℹ Some GNA `data_sources` are older than 365 days; name resolution may miss
#   recent taxa:
#   TAXREF (id 210): last updated 2025-04-02
# ℹ Compare update dates at <https://verifier.globalnames.org/data_sources>.
# ✔ GNA verification summary:
# • Total taxa in phyloseq: 4151
# • Taxa submitted for verification: 2665
# • Genus-level only taxa: 0
# • Total matches found: 1156
# • Synonyms: 105 (including 5 uninomial)
# • Accepted names: 1051 (including 346 uninomial)
# ℹ 346 uninomial accepted name(s) have `currentCanonicalSimple` set to "NA"
#   (`species_only` = TRUE)
cat("Get traits from FunGuild and FungalTraits\n")

# Assign traits from FunGuild and FungalTraits
data_traits <- fungal_traits_guilds(
  data_clean,
  fungal_traits_file = FUNGAL_TRAITS_TABLE,
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
    "genusEpithet",
    "specificEpithet"
  ),
  fg_col_prefix = "fg_",
  db_url = "http://www.stbates.org/funguild_db_2.php",
  add_consensus = TRUE,
  consensus_col_prefix = "cons_",
  add_to_phyloseq = TRUE,
  verbose = TRUE
)
# ✔ Added 17 columns from ../traitsTable/FUNGALT_DB_MROY041125.csv with information for 21
# 79 taxa in the tax_table slot of the phyloseq object
# ✔ Added 11 FUNGuild columns to tax_table.
# ✔ Added 2 consensus columns:
# "cons_trophicMode" and "cons_trophicMode_agreement".
# Messages d'avis :
# 1: Dans grepl(ft_norm, fg_norm, fixed = TRUE) :
#   l'argument pattern a une longueur > 1 et seul le premier élément est utilisé
# 2: Dans grepl(ft_norm2, fg_norm2, fixed = TRUE) :
#   l'argument pattern a une longueur > 1 et seul le premier élément est utilisé
#
# == Clean and add personnalized columns for further analysis == #

data_traits_no_check <- fungal_traits_guilds(
  data_clean,
  fungal_traits_file = FUNGAL_TRAITS_TABLE,
  ft_taxonomic_rank = "Genus",
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
cat("Add personnalised columns and cleaning of fg_guild\n")

# Clean and standardize the "fg_guild" column
# e.g. "Endophyte-|Plant Pathogen|" -> "endophyte-plant_pathogen"
data_traits <- clean_guild(data_traits)
data_traits <- clean_troph(data_traits)

# Add a ft_trophicMode column to compare FunGuild and FungalTraits
data_traits <- add_trophicMode_ft(data_traits)
data_traits_no_check <- clean_guild(data_traits_no_check)
data_traits_no_check <- clean_troph(data_traits_no_check)

# Add a ft_trophicMode column to compare FunGuild and FungalTraits
data_traits_no_check <- add_trophicMode_ft(data_traits_no_check)


tax_table_traits <- as.data.frame(tax_table(data_traits)) |>
  tibble::rownames_to_column("otu97")
tax_table_traits_no_check <- as.data.frame(tax_table(data_traits_no_check)) |>
  tibble::rownames_to_column("otu97")

sum(!is.na(tax_table_traits$fg_trophicMode))
sum(!is.na(tax_table_traits_no_check$fg_trophicMode))

sum(!is.na(tax_table_traits$ft_primary_lifestyle))
sum(!is.na(tax_table_traits_no_check$ft_primary_lifestyle))


write.csv(
  tax_table_traits,
  paste0(
    "pq_",
    stringr::str_replace_all(
      tools::file_path_sans_ext(basename(INPUT_TABLE)),
      " ",
      ""
    ),
    "_FungalTraits_FunGuild.csv"
  ),
  row.names = FALSE
)
