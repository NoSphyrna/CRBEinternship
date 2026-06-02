# ============================================================================ #
# =========================== Test file for taxinfo ========================== #
# ============================================================================ #

# ============== Installation of packages =================== #

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.22")
# devtools::install_github("adrientaudiere/MiscMetabar")
# devtools::install_github("adrientaudiere/taxinfo")
# pak::pak("adrientaudiere/taxinfo")

# devtools::install_github is deprecated use pak::pak instead

# TaxInfo
# pak::pak("adrientaudiere/MiscMetabar")
# pak::pak("adrientaudiere/taxinfo")
#
#
# # FugalTraits
# pak::pak("ropenscilabs/datastorr")
# pak::pak("traitecoevo/fungaltraits")

# ============== Libraries ============ #
library(taxinfo)
library(MiscMetabar)
library(fungaltraits)
library(stringr)
library(tibble)

# ============== Set Workspace accordingly ======== #
# Modify here where tou want to work
getwd()
setwd("Database/")

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

clean_guild <- function(table_traits) {
  table_traits |>
    mutate(
      fg_guild = fg_guild |>
        str_squish() |>
        str_to_lower() |>
        str_replace_all(" ", "_") |>
        str_replace_all("\\|", "")
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

estimate_assignation_per_sample <- function(
  table_traits,
  troph_col,
  title_end
) {
  cat("Replacing NA by \"unknown\"\n")
  table_traits[[troph_col]][is.na(
    table_traits[[troph_col]]
  )] <- "unknown"

  cat("Trimming spaces\n")
  # Trim spaces around and inside the column
  table_traits[[troph_col]] <- str_squish(table_traits[[troph_col]])

  # Sum the relative abundance by trophic mode and sample
  cat("Creating temporary grouped table\n")
  trophic_long <- table_traits |>
    group_by(Sample, .data[[troph_col]]) |>
    summarise(abundance = sum(Abundance), .groups = "drop") |>
    group_by(.data[[troph_col]]) |>
    summarise(mean_abundance = mean(abundance))

  return(trophic_long)
}

# This function plots the abundance of each trophicMode in each samples of a melted phyloseq object
plot_trophic_abundance_per_sample <- function(
  table_traits,
  troph_col,
  title_end
) {
  cat("Replacing NA by \"unknown\"\n")
  table_traits[[troph_col]][is.na(
    table_traits[[troph_col]]
  )] <- "unknown"

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
  others <- trophs[trophs != "unknown"]

  cols <- setNames(
    scales::hue_pal()(length(others)),
    others
  )

  # Change color of Unknown to grey
  cols <- c(cols, "unknown" = "grey70")
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

# ================ Workflow for assigning traits on OTUs from samples ================= #

# A test on Tedersoo full sintax
OTU_table <- read.csv(
  "cleaned_OTU_tables/OTU_table_Tedersoo_ITS_full_vsearch_fungi_clean.csv",
  header = TRUE,
  sep = ";"
)
# Test for illumina
OTU_table_il <- read.csv(
  "cleaned_OTU_tables/OTU_table_Illumina_ITS1_vsearch_fungi_clean.csv",
  header = TRUE,
  sep = ";"
)

rm(OTU_table_il)
rm(data)
rm(data_clean)

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

OTU_table_il <- standardise_illumina(OTU_table_il)
colnames(OTU_table)[!colnames(OTU_table) %in% colnames(OTU_table_il)]

head(OTU_table_il$amplicon)

head(OTU_table_il$proba_last)

OTU_table <- expand_taxnames(OTU_table_il, taxonomy)

# First expand the taxnames to match phyloseq conventions
OTU_table <- expand_taxnames(OTU_table, taxonomy)

# Then get a phyloseq object to use the taxInfo functions
data <- phyloseq_from_table(OTU_table)

# Check names according to Taxref (210) conventions
data_clean <- gna_verifier_pq(data, data_sources = 210)

# Tests for the modified taxa by the verifier
taxa_ps <- as.data.frame(tax_table(data_clean))
mismatches <- taxa_ps |>
  filter(!is.na(Genus) & !is.na(genusEpithet)) |>
  filter(Genus != genusEpithet) |>
  group_by(Genus, genusEpithet) |> #, Species, currentCanonicalSimple) |>
  summarise(number_mismatch = n(), .groups = "drop") |>
  arrange(desc(number_mismatch))

print.data.frame(mismatches)
sum(mismatches$number_mismatch)

# Get traits from the enhanced fungalTraits and funGuild
data_traits <- fungal_traits_guilds(
  data_clean,
  fungal_traits_file = "traitsTable/FUNGALT_DB_MROY041125.csv",
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

colnamed(as.data.frame(tax_table(data_clean)))

# Occurence from GBIF
data_clean_occ <- tax_gbif_occur_pq(data_clean, by_country = TRUE)

taxa <- read.table(
  file = "tax_table_gbif_tedersoo_ITS_full_vsearch_clean_names.csv",
  header = TRUE,
  sep = ";",
  row.names = "OTU"
)
taxa <- as.data.frame(tax_table(data_clean_occ))
otu_mat <- as.data.frame(otu_table(data_clean_occ))

sample_mat <- t(otu_mat)
str(sample_mat)

# Select only country datas and get the relative occurance found
occ_mat <- taxa |>
  select(
    -c(
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      taxa_name,
      currentName,
      currentCanonicalSimple,
      genusEpithet,
      specificEpithet,
      namePublishedInYear,
      authorship,
      bracketauthorship,
      scientificNameAuthorship
    )
  ) |>
  mutate_all(~ replace(., is.na(.), 0)) |>
  mutate(
    row_sum = rowSums(across(where(is.numeric)))
  ) |>
  mutate(row_sum = if_else(row_sum == 0, 1, row_sum)) |>
  mutate_at(vars(-("row_sum")), ~ . / row_sum) |>
  select(-row_sum)

head(occ_mat)
str(occ_mat)

# Get row names
taxa <- taxa |>
  rownames_to_column(var = "OTU")

write.table(
  taxa,
  file = "tax_table_gbif_tedersoo_ITS_full_vsearch_clean_names.csv",
  sep = ";",
  row.names = FALSE
)

data_clean_not_verified <- data_clean


tax_gbif_occur_pq(
  data_clean,
  taxonomic_rank = c("Genus", "Species"),
  by_country = TRUE
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

# ============= If previous phyloseq regstered ============= #

table_traits <- read_pq(
  path = "phyloseq_traits/pq_Tedersoo_ITS_full_vsearch_FungalTraits_enhanced_FunGuild",
  taxa_are_rows = TRUE
)

cols_to_delete <- c(
  "scientificNameAuthorship",
  "bracketauthorship",
  "authorship",
  "namePublishedInYear",
  "fg_citationSource",
  "fg_notes",
  "fg_growthForm"
)
for (col in cols_to_delete) {
  table_traits[[col]] <- NULL
  gc()
}

# Plot trophic abundance graph with fungatraits assignments
plot_trophic_abundance_per_sample(
  table_traits,
  "ft_trophicMode",
  "FungalTraits enhanced"
)
ggsave(
  "plot_tests/trophic_abundance_Illumina_ITS1_vsearch_FungalTraits_enhanced.png",
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
  "plot_tests/trophic_abundance_Illumina_ITS1_vsearch_FungalTraits_enhanced.png",
  width = 20,
  height = 10
)

# Get mean assignation
estimate_assignation_per_sample(
  table_traits,
  "fg_trophicMode",
  "FunGuild"
)
# ==================== Test of unasigned taxa ============================== #

df_traits <- as.data.frame(tax_table(data_traits))
df_traits <- add_trophicMode_ft(df_traits)
no_trophic <- df_traits |>
  filter(is.na(ft_trophicMode) & !is.na(fg_trophicMode)) |>
  group_by(Phylum, Family, Genus) |> #, Species, currentCanonicalSimple) |>
  summarise(number_mismatch = n(), .groups = "drop") |>
  arrange(desc(number_mismatch)) |>
  filter(number_mismatch > 15)

print.data.frame(no_trophic)
sum(no_trophic$number_mismatch)

no_trophic <- df_traits |>
  filter(is.na(ft_trophicMode)) |>
  group_by(Phylum) |> #, Species, currentCanonicalSimple) |>
  summarise(number_mismatch = n(), .groups = "drop") |>
  arrange(desc(number_mismatch)) |>
  filter(number_mismatch > 15)

print.data.frame(no_trophic)
sum(no_trophic$number_mismatch)


# ======== Venn Diagrams ============= #
head(df_traits$cons_trophicMode_agreement)

library(VennDiagram)

X <- list(
  FunGuild = rownames(df_traits |> filter(!is.na(fg_trophicMode))),
  FungalTraits = rownames(df_traits |> filter(!is.na(ft_trophicMode)))
)

getwd()

grid.newpage()

venn_plot <- venn.diagram(
  x = X,
  category.names = c("FunGuild", "Fungal Traits"),
  filename = NULL,
  output = TRUE,
  imagetype = "png",

  lwd = 1,
  col = c("orange", 'orchid'),
  fill = c(alpha("orange", 0.5), alpha('orchid', 0.5)),

  cex = 0.5,
  fontfamily = "sans",
  ext.text = FALSE,

  cat.cex = 0.3,
  cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("orange", 'orchid')
)


grid.draw(venn_plot)

png(
  'Venn_diagrams/trophic_funguild_vs_fungaltraits_Tedersoo_ITS_full_vsearch.png',
  height = 480,
  width = 480,
  res = 300
)

grid.draw(venn_plot)
dev.off()

X <- list(
  FunGuild = rownames(df_traits |> filter(!is.na(fg_trophicMode))),
  FungalTraits = rownames(df_traits |> filter(!is.na(ft_trophicMode))),
  consensus = rownames(df_traits |> filter(!is.na(cons_trophicMode)))
)


venn.diagram(
  x = X,
  category.names = c("FunGuild", "Fungal Traits", "consensus"),
  filename = 'Venn_diagrams/trophic_funguild_vs_fungaltraits_Tedersoo_ITS_full_vsearch.png',
  output = TRUE,
  imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = c("orange", 'orchid', 'slateblue'),
  fill = c(alpha("orange", 0.5), alpha('orchid', 0.5), alpha('slateblue', 0.5)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("orange", 'orchid', 'slateblue')
)


print("yeay")
# ==================== Tests and unused functions ========================== #

# == Test for the table obtain from the fungaltrait r package == #
# Get traits from the package fungalTraits and funGuild
# Get the Fungal traits csv file

# fun_traits <- fungal_traits()
# ft_col_names <- colnames(fun_traits)
# write.csv(fun_traits, file = "traitsTable/fungalTraits.csv")

# data_traits <- fungal_traits_guilds(
#   data_clean,
#   fungal_traits_file = "traitsTable/fungalTraits.csv",
#   ft_taxonomic_rank = "genusEpithet",
#   ft_csv_rank = "Genus",
#   ft_sep = ",",
#   ft_col_prefix = "ft_",
#   ft_csv_cols_select = ft_col_names,
#   fg_tax_levels = c(
#     "Kingdom",
#     "Phylum",
#     "Class",
#     "Order",
#     "Family",
#     "Genus",
#     "Species"
#   ),
#   fg_col_prefix = "fg_",
#   db_url = "http://www.stbates.org/funguild_db_2.php",
#   add_consensus = TRUE,
#   consensus_col_prefix = "cons_",
#   add_to_phyloseq = TRUE,
#   verbose = TRUE
# )

# table_traits <- psmelt(data_traits)

# # Plot trophic abundance graph with fungatraits assignments
# plot_trophic_abundance_per_sample(
#   table_traits,
#   "ft_trophic_mode_fg",
#   "FungalTraits R Package"
# )
# ggsave(
#   "plot_tests/trophic_abundance_GetPlage_ITS1_sintax_FungalTraits_package.png",
#   width = 20,
#   height = 10
# )
# test_tax <- tax_table(data_traits)
# tax_table_traits <- as.data.frame(tax_table(data_traits))
# unique(tax_table_traits$ft_trait_fg)
# unique(tax_table_traits$ft_guild_fg)
# unique(tax_table_traits$ft_trophic_mode_fg)

unique(table_traits$ft_primary_lifestyle)

unique(table_traits$fg_trait)
unique(table_traits$fg_guild)
table_traits <- clean_guild(table_traits)

plot_trophic_abundance_per_sample(
  table_traits,
  "fg_guild",
  "FunGuild guilds"
)

unique(table_traits$fg_guild)
unique(table_traits$fg_trophicMode)

#
# Test the tax bar function (not useful for trophic per sample)
# tax_bar_pq(data_traits, taxa = "ft_trophic_mode_fg")

# ======================== Examples =========================== #
# df <- data.frame(
#   OTU = c("OTU1", "OTU2", "OTU3", "OTU4"),
#   Sample1 = c(0.1, 0, 0.3, 0.6),
#   Sample2 = c(0, 0.6, 0, 0.4),
#   Kingdom = c("Fungi", "Fungi", "Fungi", "Fungi"),
#   Phylum = c(
#     "Ascomycota",
#     "Basidiomycota",
#     "Glomeromycota",
#     "Blastocladiomycota"
#   )
# )
# df
# #    OTU Sample1 Sample2 Kingdom             Phylum
# # 1 OTU1     0.1     0.0   Fungi         Ascomycota
# # 2 OTU2     0.0     0.6   Fungi      Basidiomycota
# # 3 OTU3     0.3     0.0   Fungi      Glomeromycota
# # 4 OTU4     0.6     0.4   Fungi Blastocladiomycota
# otu_matrix <- df |>
#   select(OTU, where(is.numeric)) |>
#   tibble::column_to_rownames("OTU") |>
#   as.matrix()
#
# ## here OTU names are the row so we have to set taxa_are_rows to TRUE
# OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
#
# ## Get the Taxa matrix
# tax_matrix <- df |>
#   select(OTU, Kingdom, Phylum) |>
#   tibble::column_to_rownames("OTU") |>
#   as.matrix()
#
# TAX <- tax_table(tax_matrix)
#
# ps <- phyloseq(OTU, TAX)
# ps
# # phyloseq-class experiment-level object
# # otu_table()   OTU Table:         [ 4 taxa and 2 samples ]
# # tax_table()   Taxonomy Table:    [ 4 taxa by 2 taxonomic ranks ]
# df_melt <- psmelt(ps)
# df_melt
# #    OTU  Sample Abundance Kingdom             Phylum
# # 4 OTU2 Sample2       0.6   Fungi      Basidiomycota
# # 7 OTU4 Sample1       0.6   Fungi Blastocladiomycota
# # 8 OTU4 Sample2       0.4   Fungi Blastocladiomycota
# # 5 OTU3 Sample1       0.3   Fungi      Glomeromycota
# # 1 OTU1 Sample1       0.1   Fungi         Ascomycota
# # 2 OTU1 Sample2       0.0   Fungi         Ascomycota
# # 3 OTU2 Sample1       0.0   Fungi      Basidiomycota
# # 6 OTU3 Sample2       0.0   Fungi      Glomeromycota
#
