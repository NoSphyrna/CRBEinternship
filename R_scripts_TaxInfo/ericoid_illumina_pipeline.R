# ================================================================================== #
# ==== Script to add traits to Illumina OTUs with FunGuild and ericoid database ==== #
# ================================================================================== #

# ============================ Installation of packages ============================ #

# <NOTE> devtools::install_github is deprecated use pak::pak instead

## To install the package taxinfo uncomment and execute the following:
# pak::pak("adrientaudiere/MiscMetabar")
# pak::pak("adrientaudiere/taxinfo")

# <NOTE> You might need to install pak before with the command:
# install.packages("pak")

## Installation of stringr (to manipulate strings)
# install.packages("stringr")

## Installation of tidyr (to manipulate tables in the plot functions)
# install.packages("tidyr")

# ==================================== Settings ==================================== #
# Here are the default settings with your files and directory locations

# Uncomment and adapt the following to make sure you are in the right directory
# if you run this cript not from command lines
# getwd()
# setwd("../../Database/Puymorens/")

INPUT_TABLE <- "OTU_table_PUY_ITS1_OTU97_VSEARCH.csv"
OUTPUT_DIR <- "."
ERICOID_TABLE <- "ericoid_table.csv"
FUNGAL_TRAITS_TABLE <- "../traitsTable/FUNGALT_DB_MROY041125.csv"

# ======================== Run the script from command line ======================== #

# This part allows this script to be run with a command line (the command is given by
# the "Usage :")
# If you want to execute it directly from an R editor as Rstudio skip this part and
# enter your file locations in the "settings" part above

# Get input directory and output directory :

args <- commandArgs(trailingOnly = TRUE)

# Usage thrown there is not the correct number of arguments
if (length(args) != 4) {
  stop(
    "Usage : Rscript ericoid_illumina_pipeline.R INPUT_TABLE OUTPUT_DIR ERICOID_TABLE FUNGAL_TRAITS_TABLE"
  )
}
INPUT_TABLE <- args[1]
OUTPUT_DIR <- args[2]
ERICOID_TABLE <- args[3]
FUNGAL_TRAITS_TABLE <- args[4]

# Checks if the gievn files exists
if (!file.exists(INPUT_TABLE)) {
  stop("[R] : Input table doesn't exist")
}

if (!file.exists(ERICOID_TABLE)) {
  stop("[R] : Ericoid table doesn't exist")
}

if (!file.exists(FUNGAL_TRAITS_TABLE)) {
  stop("[R] : FungalTraits table doesn't exist")
}


# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}


# =================================== Libraries ==================================== #
library(taxinfo)
library(stringr)
library(tidyr) # this one is used only in the plot functions

# Check loaded libraries:
# (.packages())
# =============================== Utility Functions ================================ #

# Get a phyloseq object from registered files with the "write_pq" function
# Note this a corrected version of the already existing read_pq function since it
# didn't work properly
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


# Discard unused columns and keep only the OTUs, abundances and taxonomy
standardise_illumina <- function(OTU_table) {
  OTU_table <- OTU_table |>
    select(
      !c(
        abundance,
        length,
        chimera,
        spread,
        identity
      )
    ) |>
    rename(OTU = amplicon)
  return(OTU_table)
}

# Expand taxa columns
# Corrections :
# -> k or d for the kingdom
# -> exclude (Fungi) at the end of Genus name e.g. Lactarius(Fungi) -> Lactarius
expand_taxnames <- function(OTU_table) {
  OTU_table <- OTU_table |>
    mutate(
      Kingdom = str_match(taxonomy, "[kd]:([^, \n]+)")[, 2],
      Phylum = str_match(taxonomy, "p:([^, \n]+)")[, 2],
      Class = str_match(taxonomy, "c:([^, \n]+)")[, 2],
      Order = str_match(taxonomy, "o:([^, \n]+)")[, 2],
      Family = str_match(taxonomy, "f:([^, \n]+)")[, 2],
      Genus = str_match(taxonomy, "g:([^, \n(]+)")[, 2], #here for the exclusion of (Fungi)
      Species = str_match(taxonomy, "s:([^, \n]+)")[, 2]
    )
  return(OTU_table)
}

# Get a phyloseq object from the table with the taxa
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

# This function adds a trophic column based on the following priority:
# Mycorrhizal > Pathogen > Endophyte > Saprotroph > Other
# NA is captured as "Unknown"
# For the myco_type column, the ericoid table is prioritised else it's all
# the mycotypes given in fg_guild
fg_to_trophic_plant <- function(data_traits) {
  table_traits <- as.data.frame(tax_table(data_traits))
  table_traits <- table_traits |>
    mutate(
      plant_trophic = case_when(
        !is.na(er_guild) | str_detect(fg_guild, "mycorrhizal") ~ "Mycorrhizal",
        is.na(fg_guild) ~ "Unknown", #placed her to avoid unnecessary comparisons
        str_detect(fg_guild, "plant_pathogen|plant_parasite") ~ "Pathogen",
        str_detect(fg_guild, "endophyte") ~ "Endophyte",
        str_detect(fg_guild, "saprotroph") ~ "Saprotroph",
        TRUE ~ "Other"
      ),
      myco_type = case_when(
        !is.na(er_guild) &
          str_detect(er_guild, "endophyte") ~ "ericoid_mycorrhizal",
        !is.na(er_guild) ~ er_guild,
        TRUE ~ fg_guild |>
          str_split("-") |>
          map_chr(
            ~ {
              myco <- str_subset(.x, "mycorrhizal")
              if (length(myco) == 0) {
                NA_character_
              } else {
                paste(myco, collapse = "-")
              }
            }
          )
      )
    )
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

# =================================== plot function ================================= #

plot_trophic_abundance <- function(data_rel, troph_col, title) {
  otu <- as.data.frame(otu_table(data_rel))
  tax <- as.data.frame(tax_table(data_rel))

  troph <- select(tax, all_of(troph_col))
  # all_of() is used because troph_col is a variable

  otu <- merge(otu, troph, by = 0) |> # merge is joining otu table and troph col
    # based on rownames here rownames are OTUs (by=0)
    group_by(across(all_of(troph_col))) |> # the across allows troph_col to be a variable
    summarise(
      across(where(is.numeric), sum) # sum the abundances by group
    )
  trophs <- unique(otu[[troph_col]])
  others <- trophs[trophs != "Unknown"]

  cols <- setNames(
    scales::hue_pal()(length(others)),
    others
  )

  # Change color of Unknown to grey
  cols <- c(cols, "Unknown" = "grey70")

  #Pivot the table to long fromat
  otu <- otu |>
    pivot_longer(
      cols = -all_of(troph_col),
      names_to = "Sample",
      values_to = "abundance"
    )

  cat("Plotting graph\n")
  p <- ggplot(
    otu,
    # .data[[troph_col]] because troph col is a variable
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
      title = title,
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


plot_pie_for_a_sample <- function(data_rel, troph_col, sample, title) {
  otu <- as.data.frame(otu_table(data_rel))
  tax <- as.data.frame(tax_table(data_rel))

  troph <- select(tax, all_of(troph_col)) # select only the column of interest
  # all_of() is used because troph_col is a variable

  otu <- merge(otu, troph, by = 0) |> # merge is joining otu table and troph col
    # based on rownames here rownames are OTUs (by=0)
    filter(!is.na(.data[[troph_col]])) |> # for filter we have to use
    # .data[[troph_col]] because troph col is a variable
    group_by(across(all_of(troph_col))) |> # the across(all_of()) allows troph_col to be a variable
    summarise(
      across(where(is.numeric), sum) # sum the abundances by group
    )

  trophs <- unique(otu[[troph_col]])

  cols <- setNames(
    scales::hue_pal()(length(trophs)),
    trophs
  ) # get the colour palette

  cat("Plotting graph\n")
  p <- ggplot(
    otu,
    aes(x = "", y = .data[[sample]], fill = .data[[troph_col]])
  ) +
    # here it's important to put identity so the bars correspond to the relative abundance
    # in each samples and not the number of different abundance found in each sample
    geom_bar(stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(
      values = cols
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = title,
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

cat("Getting table\n")
# Get the OTU_table
OTU_table <- read.csv(
  INPUT_TABLE,
  header = TRUE,
  sep = ";"
)

# Remove unused columns
OTU_table <- standardise_illumina(OTU_table)

# First expand the taxnames to match phyloseq conventions
OTU_table <- expand_taxnames(OTU_table)

# Then get a phyloseq object to use the taxInfo functions
data <- phyloseq_from_table(OTU_table)
cat("Check names\n")
# Check names according to Taxref (210) conventions
data_clean <- gna_verifier_pq(data, data_sources = 210)

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


cat("Get traits from ericoid table\n")

# Add the erricoid traits from the ericoid table
data_traits <- tax_info_pq(
  data_traits,
  taxonomic_rank = "genusEpithet",
  file_name = ERICOID_TABLE,
  csv_taxonomic_rank = "Genus",
  col_prefix = "er_",
  sep = ";",
  verbose = TRUE,
)

# == Clean and add personnalized columns for further analysis == #

cat("Add personnalised columns and cleaning of fg_guild\n")

# Clean and standardize the "fg_guild" column
# e.g. "Endophyte-|Plant Pathogen|" -> "endophyte-plant_pathogen"
data_traits <- clean_guild(data_traits)

# Add a ft_trophicMode column to compare FunGuild and FungalTraits
data_traits <- add_trophicMode_ft(data_traits)

# Add a trophic mode on a plant point of view and mycorrhizal type column
data_traits <- fg_to_trophic_plant(data_traits)


cat("Saving phyloseq object\n")
# Save the phyloseq object in the OUTPUT_PQ directory
write_pq(
  data_traits,
  path = file.path(
    OUTPUT_DIR,
    paste0(
      "/pq_",
      tools::file_path_sans_ext(basename(INPUT_TABLE)),
      "_FungalTraits_FunGuild_er_table"
    )
  )
)

# To get the phyloseq object from the save uncomment the following:
data_traits <- read_pq_corr(
  path = file.path(
    OUTPUT_DIR,
    paste0(
      "/pq_",
      tools::file_path_sans_ext(basename(INPUT_TABLE)),
      "_FungalTraits_FunGuild_er_table"
    )
  )
)

# ============================ Graphs exemples ==================================== #

# First lets use relative abundance instead of read number
data_rel <- transform_sample_counts(data_traits, function(x) x / sum(x))

# Plot the trophic abundance of all samples
p <- plot_trophic_abundance(
  data_rel,
  "plant_trophic",
  "Trophic diversity per sample : FunGuild and ericoid table"
)

# Print p if you want to see the plot in Rstudio: (uncomment the following)
# p

cat("Saving plot\n")
# And we save it in the plot folder
ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    paste0(
      "trophic_abundance",
      tools::file_path_sans_ext(basename(INPUT_TABLE)),
      ".png"
    )
  ),
  plot = p,
  width = 20,
  height = 10,
  create.dir = TRUE
)

# Plot a pie chart of the abundance of mycorrhizal types in sample P001
p <- plot_pie_for_a_sample(
  data_rel,
  "myco_type",
  "P001",
  "Relative abundance of mycrrhizal types in sample P001"
)
# Print p if you want to see the plot in Rstudio: (uncomment the following)
# p

cat("Saving plot\n")
# And we save it in the plot folder
ggsave(
  file.path(
    OUTPUT_DIR,
    "plots",
    paste0(
      "myco_type_pie_P001",
      tools::file_path_sans_ext(basename(INPUT_TABLE)),
      ".png"
    )
  ),
  plot = p,
  width = 15,
  height = 10,
  create.dir = TRUE
)
