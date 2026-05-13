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

# ============== Libraries ============ #
library(taxinfo)
library(MiscMetabar)
library(stringr)

# ============== Set Workspace accordingly ======== #
# Get input directory and output directory :

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage : Rscript plot_trophic_abundance.R INPUT_OTU_DIRECTORY INPUT_FUNTRAITS_DIRECTORY OUTPUT_DIRECTORY"
  )
}
# INPUT_DIR <- "~/Dext2/BioInfoMatser/HAU802_BILL/ONT_sequencing/Analyses/SNP/Isec/"
INPUT_OTU <- args[1]
INPUT_TRAITS <- args[2]
OUTPUT_DIR <- args[3]

if (!dir.exists(INPUT_OTU)) {
  stop("[R] : Input OTU folder doesn't exist")
}

if (!dir.exists(INPUT_TRAITS)) {
  stop("[R] : Input traits folder doesn't exist")
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============= Utility Functoins ==================#

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
        "nectar/tap_saprotroph", # new
        "pollen_saprotroph",
        "rock-inhabiting", # new
        "soil_saprotroph",
        "unspecified_saprotroph",
        "wood_saprotroph"
      ) ~ "Saprotroph",
    x %in%
      c(
        "algal_parasite",
        "animal_parasite",
        "arthropod_parasite", # new
        "bryophilous", # new
        "lichen_parasite",
        "mycoparasite",
        "nematophagous", # new
        "plant_pathogen",
        "protistan_parasite",
        "sooty_mold",
        "unspecified_pathotroph"
      ) ~ "Pathotroph",
    x %in%
      c(
        "animal-associated",
        "animal_endosymbiont",
        "arbuscular_mycorrhizal",
        "arthropod-associated",
        "ectomycorrhizal",
        "epiphyte",
        "ericoid_mycorrhizal", # new
        "foliar_endophyte",
        "lichenized",
        "moss_symbiont",
        "root-associated", # new
        "root_endophyte",
        "root_endophyte_dark_septate", # new
        "unspecified_symbiotroph",
        "vertebrate-associated" # new
      ) ~ "Symbiotroph",
    is.na(x) | x %in% c("unspecified", "") ~ NA_character_,
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

# This function plots the abundance of each trophicMode in each samples of a melted phyloseq object
plot_trophic_abundance_per_sample <- function(
  table_traits,
  troph_col,
  title_end
) {
  cat("Replacing NA by \"Unknown\"\n")
  table_traits[[troph_col]][is.na(
    table_traits[[troph_col]]
  )] <- "Unknown"

  cat("Trimming spaces\n")
  # Trim spaces around and inside the column
  table_traits[[troph_col]] <- str_squish(table_traits[[troph_col]])

  # Sum the relative abundance by trophic mode and sample
  cat("Creating temporary grouped table\n")
  trophic_long <- table_traits |>
    group_by(Sample, .data[[troph_col]]) |>
    summarise(abundance = sum(Abundance), .groups = "drop")

  cat("Set the color of unknown to grey\n")
  trophs <- unique(trophic_long[[troph_col]])
  others <- trophs[trophs != "Unknown"]

  cols <- setNames(
    scales::hue_pal()(length(others)),
    others
  )

  # Change color of Unknown to grey
  cols <- c(cols, "Unknown" = "grey70")
  cat("Plotting graph\n")
  p <- ggplot(
    trophic_long,
    aes(x = Sample, y = abundance, fill = .data[[troph_col]])
  ) +
    # here it's important to put identity so the bars correspond to the relative abundance
    # in each samples and not the number of different abundance found in each sample
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = cols
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = paste("Trophic diversity per sample : ", title_end),
      x = "Sample",
      y = "Relative abundance",
      fill = ""
    ) +
    theme_idest() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )
  return(p)
}

# =============== Workflow for assigning traits on OTUs from samples ================ #

INPUT_OTU <- "cleaned_OTU_tables"
INPUT_TRAITS <- "traitsTable"
OUTPUT_DIR <- "plot_taxInfo_traits"
method <- "Tedersoo"
section <- "ITS_full"
cluster <- "vsearch"


# Get the OTU_table
OTU_table <- read.csv(
  paste0(
    INPUT_DIR,
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

# First expand the taxnames to match phyloseq conventions
OTU_table <- expand_taxnames(OTU_table)

# Then get a phyloseq object to use the taxInfo functions
data <- phyloseq_from_table(OTU_table)

# Check names according to Taxref (210) conventions
data_clean <- gna_verifier_pq(data, data_sources = 210)

# Get traits from the enhanced fungalTraits and funGuild
data_traits <- fungal_traits_guilds(
  data_clean,
  fungal_traits_file = paste0(INPUT_TRAITS, "/FUNGALT_DB_MROY041125.csv"),
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

# Get a melted dataframe from the phyloseq object obtain (melted datafram correpsond to reapeated lines for each sample and each OTU)
# Example :

# df
#    OTU Sample1 Sample2 Kingdom             Phylum
# 1 OTU1     0.1     0.0   Fungi         Ascomycota
# 2 OTU2     0.0     0.6   Fungi      Basidiomycota
# 3 OTU3     0.3     0.0   Fungi      Glomeromycota
# 4 OTU4     0.6     0.4   Fungi Blastocladiomycota

# ps : phyloseq object from df
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4 taxa and 2 samples ]
# tax_table()   Taxonomy Table:    [ 4 taxa by 2 taxonomic ranks ]

# df_melt <- psmelt(ps)
# df_melt
#    OTU  Sample Abundance Kingdom             Phylum
# 4 OTU2 Sample2       0.6   Fungi      Basidiomycota
# 7 OTU4 Sample1       0.6   Fungi Blastocladiomycota
# 8 OTU4 Sample2       0.4   Fungi Blastocladiomycota
# 5 OTU3 Sample1       0.3   Fungi      Glomeromycota
# 1 OTU1 Sample1       0.1   Fungi         Ascomycota
# 2 OTU1 Sample2       0.0   Fungi         Ascomycota
# 3 OTU2 Sample1       0.0   Fungi      Basidiomycota
# 6 OTU3 Sample2       0.0   Fungi      Glomeromycota

table_traits <- psmelt(data_traits)

# Reduce lifestyle to the 3 categories of trophic mode (Saprotroph, Pathotroph, Symbiotroph)
# To match the match FungalTraits with FunGuild

table_traits <- add_trophicMode_ft(table_traits)

# If there are "Other" in categories, you can check which lifestyles didn't match with the following :

# other_values <- unique(table_traits$ft_Secondary_lifestyle[
#   ft_to_trophic_mode(table_traits$ft_Secondary_lifestyle) == "Other"
# ])
# print(sort(other_values))
#
# Plot trophic abundance graph with fungatraits assignments
plot_trophic_abundance_per_sample(
  table_traits,
  "ft_trophicMode",
  "FungalTraits enhanced"
)
ggsave(
  paste0(
    OUTPUT_DIR,
    "/trophic_abundance_",
    method,
    "_",
    section,
    "_",
    cluster,
    "_FungalTraits_enhanced.png"
  ),
  width = 20,
  height = 10
)


# Plot trophic abundance for funGuild assignments
plot_trophic_abundance_per_sample(
  table_traits,
  "fg_trophicMode",
  "FunGuild"
)
ggsave(
  paste0(
    OUTPUT_DIR,
    "/trophic_abundance_",
    method,
    "_",
    section,
    "_",
    cluster,
    "_FunGuild.png"
  ),
  width = 20,
  height = 10
)
