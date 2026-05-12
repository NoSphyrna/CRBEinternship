# ============================================================================ #
# =========================== Test file for taxinfo ========================== #
# ============================================================================ #

# ============== Installation of packages =================== #

# devtools::install_github is deprecated use pak::pak instead

# TaxInfo
pak::pak("adrientaudiere/MiscMetabar")
pak::pak("adrientaudiere/taxinfo")

# FugalTraits
pak::pak("ropenscilabs/datastorr")
pak::pak("traitecoevo/fungaltraits")

# ============== Libraries ============ #
library(taxinfo)
library(MiscMetabar)
library(fungaltraits)
library(stringr)
# ============== Set Workspace accordingly ======== #

getwd()

setwd("Database/")

# Get the Fungal traits csv file

fun_traits <- fungal_traits()
ft_col_names <- colnames(fun_traits)
summarise(fun_traits)
str(fun_traits)
ft_col_names
write.csv(fun_traits, file = "traitsTable/fungalTraits.csv")

unique(fun_traits$guild_fg)
unique(fun_traits$trait_fg)
unique(fun_traits$trophic_mode_fg)

# test sintax
OTU_table <- read.csv(
  "OTU_table_Getplage_ITS1_sintax_fungi_clean.csv",
  header = TRUE,
  sep = ";"
)


head(OTU_table)

# test vsearch
OTU_table_vs <- read.csv(
  "OTU_table_Getplage_ITS1_vsearch_fungi_clean.csv",
  header = TRUE,
  sep = ";"
)

head(OTU_table_vs$taxonomy)

# Expand taxa columns
expand_taxnames <- function(OTU_table) {
  OTU_table <- OTU_table |>
    mutate(
      Kingdom = str_match(taxonomy, "[kd]:([^, \n]+)")[, 2],
      Phylum = str_match(taxonomy, "p:([^, \n]+)")[, 2],
      Class = str_match(taxonomy, "c:([^, \n]+)")[, 2],
      Order = str_match(taxonomy, "o:([^, \n]+)")[, 2],
      Family = str_match(taxonomy, "f:([^, \n]+)")[, 2],
      Genus = str_match(taxonomy, "g:([^, \n]+)")[, 2],
      Species = str_match(taxonomy, "s:([^, \n]+)")[, 2],
    )
  OTU_table
}
OTU_table <- expand_taxnames(OTU_table)

head(OTU_table$Kingdom)


# Get a phyloseq object from the table
phyloseq_from_table <- function(OTU_table) {
  ## Get the OTU table (as a matrix for the phyloseq object)

  otu_matrix <- OTU_table |>
    select(OTU, where(is.numeric)) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix()

  # here OTU names are the row so we have to set taxa_are_rows to TRUE
  OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)

  ## Get the Taxo matrix

  tax_matrix <- OTU_table |>
    select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix()

  TAX <- tax_table(tax_matrix)

  PHYLOSEQ <- phyloseq(OTU, TAX)
  PHYLOSEQ
}

data_sintax <- phyloseq_from_table(OTU_table)

# Check names

data_sintax_clean <- gna_verifier_pq(data_sintax, data_sources = 210)

# Get traits from fungalTraits and funGuild

data_sintax_traits <- fungal_traits_guilds(
  data_sintax_clean,
  fungal_traits_file = "traitsTable/fungalTraits.csv",
  ft_taxonomic_rank = "genusEpithet",
  ft_csv_rank = "Genus",
  ft_sep = ",",
  ft_col_prefix = "ft_",
  ft_csv_cols_select = ft_col_names,
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

test_tax <- tax_table(data_sintax_traits)
tax_table_sintax_traits <- as.data.frame(tax_table(data_sintax_traits))
unique(tax_table_sintax_traits$ft_trait_fg)
unique(tax_table_sintax_traits$ft_guild_fg)
unique(tax_table_sintax_traits$ft_trophic_mode_fg)

unique(tax_table_sintax_traits$fg_trait)
unique(tax_table_sintax_traits$fg_guild)
unique(tax_table_sintax_traits$fg_trophicMode)

tax_bar_pq(data_sintax_traits, taxa = "ft_trophic_mode_fg")

str(data_sintax_traits)

table_sintax_traits <- psmelt(data_sintax_traits)

table_sintax_traits$ft_trophic_mode_fg[is.na(
  table_sintax_traits$ft_trophic_mode_fg
)] <- "Unknown"

trophic_long <- table_sintax_traits %>%
  group_by(Sample, ft_trophic_mode_fg) %>%
  summarise(abundance = sum(Abundance))

# Set the color of unknown to grey
trophs <- unique(trophic_long$ft_trophic_mode_fg)
others <- trophs[trophs != "Unknown"]

cols <- setNames(
  scales::hue_pal()(length(others)),
  others
)

# Ajouter Unknown en gris
cols <- c(cols, "Unknown" = "grey70")
ggplot(
  trophic_long,
  aes(x = Sample, y = abundance, fill = ft_trophic_mode_fg)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cols
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Trophic diversity per sample",
    x = "Sample",
    y = "Relative abundance",
    fill = ""
  ) +
  theme_idest() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

ggsave("plot_tests/trophic_abundance_ITS1_sintax.png", width = 20, height = 10)
