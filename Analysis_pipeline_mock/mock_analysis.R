# ============== Libraries ============ #
library(stringi)
library(vegan)
library(compositions)
library(zCompositions)
library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::select)
library(fs)
library(gtools)
library(taxinfo)
library(iNEXT)

# ====================== Utility Functions =================== #

# Expand taxa columns
# Corrections :
# -> k or d for the kingdom
# -> exclude (Fungi) at the end of Genus name
expand_taxnames <- function(OTU_table, taxonomy, bdd) {
  if (bdd == "euk") {
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
  } else {
    OTU_table <- OTU_table |>
      mutate(
        Kingdom = str_match(taxonomy, "[kd]:([^, \n]+)")[, 2],
        Phylum = str_match(taxonomy, "p:([^, \n]+)")[, 2],
        Class = str_match(taxonomy, "c:([^, \n]+)")[, 2],
        Order = str_match(taxonomy, "o:([^, \n]+)")[, 2],
        Family = str_match(taxonomy, "f:([^, \n]+)")[, 2],
        Genus = str_match(taxonomy, "g:([^, \n(]+)")[, 2],
        Species = str_match(taxonomy, "s:[^_]+_([^, \n]+)")[, 2]
      )
  }
  return(OTU_table)
}


### get_otu
# Description:
# This function allows to get the otu_table with the total abundances column
# It also returns the vector of the sample_columns and the info columns
# Paramters:
# - file_otu_table: The path to the otu table as tsv or csv
# - techno: The sequencing techno thus the pipeline used either "illumina", "pacbio" or "nanopore"
# - format: The format of the file given in input
# - onfile: True if the taxonomy has already been added False otherwise (default: True)
# Returns:
# - res: list of [otu_table, sample_cols, info_cols]
# -- otu_table: the corresponding otu table with sample in mixeddsort order and the total_abundances column
# -- sample_cols: the vector of the sample columns names
# -- info_cols: the vector with the info columns previously attached to sample cols
get_otu <- function(
  file_otu_table,
  techno,
  sep = ";",
  sample_cols = NULL,
  onefile = TRUE,
  keep_info = TRUE
) {
  # We read the otu files
  otu = read.csv(
    file_otu_table,
    header = TRUE,
    sep = sep
  )

  # for pacbio
  if (techno == "pacbio" | techno == "nanopore") {
    info_cols <- c()
  }
  if (techno == "illumina") {
    info_cols <- c(
      "abundance",
      "length",
      "chimera",
      "spread",
      "identity"
    )
    colnames(otu) <- str_replace(colnames(otu), "amplicon", "OTU")
  }
  if (onefile) {
    info_cols <- append(info_cols, "taxonomy")
  }
  # Get the sample columns sorted
  if (is.null(sample_cols)) {
    sample_cols <- mixedsort(setdiff(colnames(otu), append(info_cols, "OTU")))
  }

  infos <- otu[, c("OTU", info_cols)]
  # Sort and get only sample and OTU cols
  otu <- otu[, c("OTU", sample_cols)]

  # sum of abundances
  otu <- otu |>
    mutate(total_abundances = rowSums(across(all_of(sample_cols))))

  # remove otus with no abundances in all samples
  otu <- otu[
    otu$total_abundances > 0,
  ]

  cat("test")
  if (keep_info) {
    otu <- merge(otu, infos, by = "OTU")
  }

  res <- list(otu, sample_cols, info_cols)
  return(res)
}


# VSEARCH (with best match)
# file_otu_table <- "Nanopore/OTU/OTU_table_mumu.tsv"
# file_taxo_table <- "Nanopore/taxo/taxonomy_OTU_vsearch_unite.tsv"
# techno <- "nanopore"
add_taxo_vsearch <- function(otu_table, file_taxo_table, techno, bdd) {
  # read taxo
  taxo = read.csv(
    file_taxo_table,
    header = FALSE,
    sep = "\t"
  ) |>
    rename("OTU" = "V1", "ID_PERCT" = "V2", "TAXO" = "V3") |>
    group_by(OTU) |>
    slice_max(order_by = ID_PERCT, n = 1, with_ties = FALSE) |>
    ungroup()

  if (techno == "illumina") {
    taxo <- taxo |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  otu_taxo <- merge(otu_table, taxo, by = "OTU")

  otu_taxo <- expand_taxnames(
    otu_taxo,
    taxonomy = otu_taxo$TAXO,
    bdd
  )

  # # we can test if there are no OTUs with no sample attached :
  # unique(nano_otu_sintax_unite$total_abundances == 0)
  #
  # sort(unique(nano_otu_sintax_unite$TAXO[
  #   nano_otu_sintax_unite$Kingdom == "unidentified"
  # ]))

  tax_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )

  # replacing "NA by unidentified"

  otu_taxo <- otu_taxo |>
    mutate(across(
      all_of(tax_cols),
      ~ if_else(is.na(.x) | (.x == "unidentified"), "unclassified", .x)
    ))

  return(otu_taxo)
}


# ============================= Workflow ============================== #

getwd()
setwd("../../Database/Mock/")


# Get otu tables

# Illumina

res <- get_otu(
  "Illumina/OTU_table_ITS1_Mock_OTU97_VSEARCH_clean.csv",
  "illumina",
  sep = ";"
)

otu_table_ill <- res[[1]]
sample_cols <- res[[2]]
info_cols <- res[[3]]

otu_table_ill <- expand_taxnames(otu_table_ill, "taxonomy", "euk")

# PacBio
res <- get_otu(
  "PacBio/OTU_table_LULU_clean0006.csv",
  "pacbio",
  sep = ",",
  sample_cols = sample_cols[sample_cols != "Tpos3_L214"],
  onefile = FALSE,
  keep_info = FALSE
)

otu_table_pacbio <- res[[1]]
sample_cols <- res[[2]]
info_cols_pacbio <- res[[3]]

otu_table_pacbio_vs <- add_taxo_vsearch(
  otu_table_pacbio,
  "PacBio/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "pacbio",
  "euk"
)

unique(otu_table_pacbio_vs$Genus)
unique(otu_table_pacbio_vs$Species)

# Reading of the original composition
compo <- read.csv(
  "Samples_plan_plaque_MITI_run2_2026 - SELECTION MOCK.tsv",
  sep = "\t"
)

unique(compo[1:46, ]$Genre.espèce)

colnames(compo) <- str_replace(colnames(compo), "Mock.([0-9]+)", "Mock\\1")

compo <- compo[1:46, ] |>
  mutate(
    Mock1 = as.numeric(Mock1),
    Genus = str_match(Genre.espèce, "([^ \n]+) ?([^\n]+)?")[, 2],
    Species = str_match(Genre.espèce, "([^ \n]+) ?([^\n]+)?")[, 3]
  ) |>
  mutate(
    Species = if_else(is.na(Species), "unclassified", Species)
  )

unique(compo$Genus)
unique(compo$Species)

# Get genus and species tables

# Illumina
genus_species_ill <- otu_table_ill |>
  select(c("Genus", "Species", sample_cols)) |>
  group_by(Genus, Species) |>
  summarise(across(all_of(sample_cols), sum), n_otu = n(), .groups = 'drop') |>
  mutate(across(all_of(sample_cols), \(x) x / sum(x)))

write.csv(genus_species_ill, "Genus_secies_ill.csv")

# PacBio
genus_species_pacbio <- otu_table_pacbio_vs |>
  select(c("Genus", "Species", sample_cols)) |>
  group_by(Genus, Species) |>
  summarise(across(all_of(sample_cols), sum), n_otu = n(), .groups = 'drop') |>
  mutate(across(all_of(sample_cols), \(x) x / sum(x)))

write.csv(genus_species_pacbio, "Genus_secies_pacbio.csv")


# Original Composition
genus_species_compo <- compo |>
  select(c("Genus", "Species", sample_cols[sample_cols != "Tpos3_L214"])) |>
  mutate(across(all_of(sample_cols[sample_cols != "Tpos3_L214"]), function(x) {
    if_else(is.na(x), 0, x)
  })) |>
  group_by(Genus, Species) |>
  summarise(
    across(all_of(sample_cols[sample_cols != "Tpos3_L214"]), sum),
    .groups = 'drop'
  ) |>
  mutate(across(all_of(sample_cols[sample_cols != "Tpos3_L214"]), \(x) {
    x / sum(x)
  }))


colnames(genus_species_compo) <- str_replace(
  colnames(genus_species_compo),
  "Mock([0-9]+)",
  "Mock\\1_gt"
)
write.csv(genus_species_compo, "Genus_secies_compo.csv")

# merge it

# Illumina
genus_species_merged <- full_join(
  genus_species_ill,
  genus_species_compo,
  by = c("Genus", "Species")
)

genus_species_merged_filled <- genus_species_merged |>
  mutate(across(where(is.numeric), .fns = function(x) if_else(is.na(x), 0, x)))

cols <- setdiff(colnames(genus_species_merged_filled), c("Genus", "Species"))

genus_species_merged_filled <- genus_species_merged_filled |>
  select(Genus, Species, all_of(mixedsort(cols))) |>
  arrange(Genus, Species)


write.csv(genus_species_merged_filled, "Genus_species_merged_filled.csv")

# PacBio

genus_species_merged_pacbio <- full_join(
  genus_species_pacbio,
  genus_species_compo,
  by = c("Genus", "Species")
)

genus_species_merged_filled_pacbio <- genus_species_merged_pacbio |>
  mutate(across(where(is.numeric), .fns = function(x) if_else(is.na(x), 0, x)))

cols <- setdiff(
  colnames(genus_species_merged_filled_pacbio),
  c("Genus", "Species")
)

genus_species_merged_filled_pacbio <- genus_species_merged_filled_pacbio |>
  select(Genus, Species, all_of(mixedsort(cols))) |>
  arrange(Genus, Species)


write.csv(
  genus_species_merged_filled_pacbio,
  "Genus_species_merged_filled_pacbio.csv"
)

# Full join 3 tables
colnames(genus_species_pacbio) <- str_replace(
  colnames(genus_species_pacbio),
  "Mock([0-9]+)",
  "Mock\\1_pac"
)
colnames(genus_species_ill) <- str_replace(
  colnames(genus_species_ill),
  "Mock([0-9]+)",
  "Mock\\1_ill"
)
colnames(genus_species_pacbio) <- str_replace(
  colnames(genus_species_pacbio),
  "n_otu",
  "n_otu_pac"
)
colnames(genus_species_ill) <- str_replace(
  colnames(genus_species_ill),
  "n_otu",
  "n_otu_ill"
)
genus_species_merged_pacbio <- full_join(
  genus_species_pacbio,
  genus_species_compo,
  by = c("Genus", "Species")
)
genus_species_merged_all <- full_join(
  genus_species_merged_pacbio,
  genus_species_ill,
  by = c("Genus", "Species")
)

genus_species_merged_filled_all <- genus_species_merged_all |>
  mutate(across(where(is.numeric), .fns = function(x) if_else(is.na(x), 0, x)))

cols <- setdiff(
  colnames(genus_species_merged_filled_all),
  c("Genus", "Species")
)

genus_species_merged_filled_all <- genus_species_merged_filled_all |>
  select(Genus, Species, all_of(mixedsort(cols))) |>
  arrange(Genus, Species)


write.csv(
  genus_species_merged_filled_all,
  "Genus_species_merged_filled_all.csv"
)


# Pie charts

# Create folder :
if (!dir.exists("pie_charts")) {
  dir.create("pie_charts")
}

table <- genus_species_merged_filled_all

# Merge Taxonomy
table$taxonomy <- paste(table$Genus, table$Species)

# Get the list of samples
list_samples <- grep("_(gt|ill|pac)$", colnames(table), value = TRUE)

# Get the list of unique samples
unique_samples <- unique(sub("_(gt|ill|pac)$", "", list_samples))

# Do the pie charts for each mock
for (sample in unique_samples) {
  merged_plot <- c()
  # Here a little adaptation just so Mock1_gt get separated from Mock10_gt for example
  reps <- list_samples[startsWith(list_samples, paste0(sample, "_"))]

  for (rep in reps) {
    extract <- table[, c("taxonomy", rep)]
    extract$method <- sub(paste0("^", sample, "_"), "", rep) # gt / ill / pac
    colnames(extract) <- c("taxonomy", "plot", "method")
    extract <- extract[extract$plot > 0, ]
    merged_plot <- rbind(merged_plot, extract)
  }

  merged_plot$method <- factor(
    merged_plot$method,
    levels = c("gt", "ill", "pac"),
    labels = c("Mock", "Illumina", "PacBio")
  )

  pdf(
    paste0("pie_charts/fungal_composition_", sample, ".pdf"),
    width = 12,
    height = 5
  )
  print(
    ggplot(merged_plot, aes(x = factor(1), y = plot, fill = factor(taxonomy))) +
      geom_bar(stat = "identity", width = 1) +
      facet_wrap(. ~ method) +
      coord_polar(theta = "y") +
      theme_classic() +
      ylab("") +
      xlab("") +
      labs(fill = "") +
      theme(
        strip.text = element_text(face = "bold"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
  )
  dev.off()
}
