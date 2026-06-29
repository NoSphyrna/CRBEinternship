# ============== Libraries ============ #
library(taxinfo)
library(MiscMetabar)
library(stringr)
library(tibble)
library(stringi)
library(ggsci)
library(tidyr)
library(sf)
library(rnaturalearth)
library(CoordinateCleaner)
library(vegan)
library(compositions)
library(tidyverse)
# ============== Set Workspace accordingly ======== #
# Modify here where tou want to work
getwd()
setwd("Database/data_taxinfo/")

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
plot_mean_trophic_abundance <- function(
  phyloseq_list,
  troph_col = "fg_trophicMode",
  title = "Mean trophic abundance accross samples by method"
) {
  # here impa_dfr allows to get the indexes of the phyloseq objects in the list, then map applies the function to all elements and _dfr allows to use rbind at the end to combine all dataframes at the end
  combined <- purrr::imap_dfr(phyloseq_list, function(pq, group_pq) {
    otu <- as.matrix(otu_table(pq))
    tax <- as.data.frame(tax_table(pq))

    troph <- tax[[troph_col]]
    rm(tax)
    gc()

    cat("Replacing NA by \"unknown\"\n")
    troph[is.na(
      troph
    )] <- "Unknown"

    cat("Trimming spaces\n")
    # Trim spaces around and inside the column
    troph <- str_squish(troph)

    #
    troph <- troph |>
      stri_split_fixed("-") |> # liste de vecteurs
      purrr::map_chr(~ paste(sort(unique(stri_trim_both(.x))), collapse = "-"))

    troph_mat <- rowsum(otu, group = troph)

    rm(otu)
    gc()
    rm(troph)
    gc()

    trophic_mean <- data.frame(
      fg_trophicMode = rownames(troph_mat),
      MeanRelAbundance = rowMeans(troph_mat, na.rm = TRUE),
      Group = group_pq
    )

    return(trophic_mean)
  })

  cat("Set the color of unknown to grey\n")
  trophs <- unique(combined[[troph_col]])
  others <- trophs[trophs != "Unknown"]

  cols <- setNames(
    scales::pal_viridis(option = "turbo")(length(others)),
    others
  )

  # Change color of Unknown to grey
  cols <- c(cols, "Unknown" = "grey70")

  cat("Plotting graph\n")
  p <- ggplot(
    combined,
    aes(x = Group, y = MeanRelAbundance, fill = .data[[troph_col]])
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
      x = NULL,
      y = "Abondance relative",
      fill = ""
    ) +
    theme_idest() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 13), # Taille du texte
      legend.title = element_text(size = 14, face = "bold"), # Titre légende
      legend.key.size = unit(1.2, "cm") # Taille des carrés de couleur
    )

  return(p)
}

getwd()
pq_Tedersoo <- read_pq_corr(
  "traits_pq/pq_Tedersoo_ITS_full_vsearch_FungalTraits_enhanced_FunGuild"
)
pq_GetPlage <- read_pq_corr(
  "traits_pq/pq_Getplage_ITS_none_vsearch_FungalTraits_enhanced_FunGuild"
)
pq_Illumina <- read_pq_corr(
  "traits_pq/pq_Illumina_ITS1_vsearch_FungalTraits_enhanced_FunGuild"
)

pq_list <- list(
  # "GetPlage" = pq_GetPlage,
  "Tedersoo" = pq_Tedersoo,
  "Illumina" = pq_Illumina
)

plot_mean_trophic_abundance(pq_list)

# ramp <- scales::colour_ramp(c("darkred","darkgreen","darkblue"))
# scales::show_col(c(ramp(seq(0,1, length = 7)),"grey70"))
ggsave(
  "trophic_abundance_Illumina_tedersoo.png",
  width = 10,
  height = 8
)


tax <- as.data.frame(tax_table(pq_Tedersoo))

tax <- as.data.frame(tax_table(pq_Illumina))

head(tax)
unique(tax$Genus)

tax <- add_trophicMode_ft(tax)

summary_table <- tax %>%
  summarise(
    # Nombre total d'OTU
    nb_OTU = n(),

    # % OTU assignés au niveau du Genre
    pct_Genus = round(
      sum(Genus != "unclassified", na.rm = TRUE) / n() * 100,
      2
    ),

    # % OTU assignés au niveau de l'Espèce
    pct_Species = round(
      sum(Species != "unclassified", na.rm = TRUE) / n() * 100,
      2
    ),

    # % OTU avec niche trophique FunGuild
    pct_fg_trophicMode = round(sum(!is.na(fg_trophicMode)) / n() * 100, 2),

    # % OTU avec trait fongique
    pct_ft_trophicMode = round(sum(!is.na(ft_trophicMode)) / n() * 100, 2)
  )

print(summary_table)


write.csv(tax, "taxa.csv")
# test
ted_col <- pq_Tedersoo |>
  otu_table() |>
  as.data.frame() |>
  colnames() |>
  unique()

ill_col <- pq_Illumina |>
  otu_table() |>
  as.data.frame() |>
  colnames() |>
  unique() |>
  str_replace("_", ".")

shared_col <- intersect(ted_col, ill_col)
shared_col


get_mat <- function(pq, troph_col) {
  otu <- as.matrix(otu_table(pq))
  tax <- as.data.frame(tax_table(pq))

  troph <- tax[[troph_col]]
  rm(tax)
  gc()

  cat("Replacing NA by \"unknown\"\n")
  troph[is.na(
    troph
  )] <- "Unknown"

  cat("Trimming spaces\n")
  # Trim spaces around and inside the column
  troph <- str_squish(troph)

  #
  troph <- troph |>
    stri_split_fixed("-") |> # liste de vecteurs
    purrr::map_chr(~ paste(sort(unique(stri_trim_both(.x))), collapse = "-"))

  troph_mat <- rowsum(otu, group = troph)
  return(troph_mat)
}

troph_tedersoo <- get_mat(pq_Tedersoo, "fg_trophicMode")
troph_illumina <- get_mat(pq_Illumina, "fg_trophicMode")

colnames(troph_illumina) <- str_replace(colnames(troph_illumina), "_", ".")

mat_ted <- troph_tedersoo[, shared_col]
mat_ill <- troph_illumina[, shared_col]

troph_tedersoo
troph_illumina

clr_ill <- as.matrix(clr(mat_ill))
clr_ted <- as.matrix(clr(mat_ted))

mat_combined <- rbind(clr_ill, clr_ted)

shared_samples <- intersect(rownames(mat_ill), rownames(mat_ted))
meta <- data.frame(
  methode = factor(c(
    rep("Illumina", length(shared_samples)),
    rep("Tedersoo", length(shared_samples))
  )),
  echantillon = factor(rep(shared_samples, 2)), # <- paired values
  row.names = rownames(mat_combined)
)


## range

pq_C <- read_pq_corr(
  "range_pq/pq_Tedersoo_ITS_full_vsearch_gbif_check/C"
)

tax <- as.data.frame(tax_table(pq_C))

tax |>
  filter(!is.na(taxa_name)) |>
  mutate(
    count_in_radius = replace_na(as.numeric(count_in_radius), 0),
    total_count_in_world = replace_na(as.numeric(total_count_in_world), 0)
  ) |>
  summarise(
    count_in_radius = max(count_in_radius),
    total_count_in_world = max(total_count_in_world),
    .by = taxa_name
  ) |>
  mutate(
    taxa_name = forcats::fct_reorder(taxa_name, count_in_radius, .fun = max)
  ) |>
  ggplot(aes(x = count_in_radius, y = taxa_name, fill = total_count_in_world)) +
  geom_col()

tax |>
  filter(!is.na(taxa_name)) |>
  mutate(
    count_in_radius = replace_na(as.numeric(count_in_radius), 0),
    total_count_in_world = replace_na(as.numeric(total_count_in_world), 0)
  ) |>
  summarise(
    count_in_radius = max(count_in_radius),
    total_count_in_world = max(total_count_in_world),
    .by = taxa_name
  ) |>
  mutate(
    occurence_pattern = case_when(
      count_in_radius == 0 & total_count_in_world == 0 ~ "Non Observée",
      count_in_radius > 10 &
        total_count_in_world > 10000 ~ "Cosmopolite étendue",
      count_in_radius > 10 &
        total_count_in_world <= 10000 ~ "Endémique régionale",
      count_in_radius <= 10 &
        total_count_in_world > 10000 ~ "Cosmopolite, rare localement",
      count_in_radius <= 10 &
        total_count_in_world <= 10000 ~ "Intéressante (rare) ou erreur",
    )
  ) |>
  summarise(
    n = n(),
    .by = occurence_pattern
  ) |>
  mutate(prop_pattern = n / sum(n)) |>
  ggplot(aes(
    x = "",
    fill = occurence_pattern,
    y = prop_pattern
  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d(option = "turbo") + # Palette turbo (discrète)
  scale_y_continuous(labels = scales::percent) + # Axe Y en %
  labs(
    x = NULL, # Supprime label axe X
    y = NULL,
    fill = "Patrons d'occurrence"
  ) +
  theme_idest() +
  theme(
    legend.text = element_text(size = 13), # Taille du texte
    legend.title = element_text(size = 14, face = "bold"), # Titre légende
    legend.key.size = unit(1.2, "cm") # Taille des carrés de couleur
  )

ggsave(
  "occurrence_pattern_tedersoo.png",
  width = 8,
  height = 8
)

## map

pq_counrty <- read_pq_corr(
  "occur_pq/pq_Tedersoo_ITS_full_vsearch_gbif_occ_by_country"
)
taxa <- as.data.frame(tax_table(pq_counrty))
otu_mat <- as.matrix(otu_table(pq_counrty))

otu_mat <- otu_mat[, "C1"]
head(otu_mat)

sample_mat <- t(otu_mat)
head(sample_mat)
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
  mutate_all(as.numeric) |>
  mutate(
    rowsum = rowSums(across(where(is.numeric)))
  ) |>
  mutate(
    rowsum = ifelse(rowsum == 0, 1, rowsum)
  ) |>
  as.matrix()

rel_occ_mat <- sweep(occ_mat, 1, occ_mat[, "rowsum"], "/")
rel_occ_mat <- rel_occ_mat[, colnames(rel_occ_mat) != "rowsum"]

head(rel_occ_mat)

rel_occ_sample <- sample_mat %*% rel_occ_mat

rel_occ_sample <- as.data.frame(t(rel_occ_sample))
head(rel_occ_sample)

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

data(countryref)
head(countryref)
centroides <- countryref |>
  select(iso2, lon = centroid.lon, lat = centroid.lat) |>
  filter(!is.na(iso2), !is.na(lon), !is.na(lat)) |>
  distinct(iso2, .keep_all = TRUE)

centroides

df <- data.frame(iso2 = row.names(rel_occ_sample), valeur = rel_occ_sample) |>
  left_join(centroides, by = "iso2")

head(df)
ggplot() +
  geom_sf(
    data = world_sf,
    fill = "#ECEFF4",
    colour = "#B0BAC6",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(x = lon, y = lat, size = C1, colour = C1),
    alpha = 0.75
  ) +
  scale_size_continuous(
    name = "",
    range = c(4, 20)
  ) +
  scale_colour_viridis_c(name = "", option = "turbo", begin = 0, end = 0.5) +
  coord_sf(xlim = c(-160, 160), ylim = c(-55, 80)) +
  theme_idest() +
  labs(
    title = "Occurrence relative pondérée par la présence dans l'échantillon"
  ) +
  theme(
    legend.text = element_text(size = 13), # Taille labels légende
    legend.title = element_text(size = 14), # Taille titre légende
    legend.key.size = unit(1.5, "cm"), # Taille des clés
    legend.spacing.y = unit(0.4, "cm") # Espacement entre entrées
  ) +
  guides(
    size = guide_legend(reverse = TRUE)
  ) # Valeur haute en haut

ggsave(
  "occurrence_map_C1_tedersoo.png",
  width = 20,
  height = 10
)

## Assignations taxonomiques

pq_Tedersoo <- read_pq_corr(
  "traits_pq/pq_Tedersoo_ITS_full_vsearch_FungalTraits_enhanced_FunGuild"
)
pq_GetPlage <- read_pq_corr(
  "traits_pq/pq_Getplage_ITS_none_vsearch_FungalTraits_enhanced_FunGuild"
)
pq_Illumina <- read_pq_corr(
  "traits_pq/pq_Illumina_ITS1_vsearch_FungalTraits_enhanced_FunGuild"
)

pq_list <- list(
  "GetPlage" = pq_GetPlage,
  "Tedersoo" = pq_Tedersoo,
  "Illumina" = pq_Illumina
)


taxa_pq <- as.data.frame(tax_table(pq_Illumina))
head(taxa)
mismatches <- taxa_pq |>
  filter(!is.na(Genus) & !is.na(genusEpithet)) |>
  mutate(original_name = paste(Genus, Species)) |>
  filter(original_name != currentCanonicalSimple) |>
  group_by(Family, original_name, currentCanonicalSimple) |> #, Species, currentCanonicalSimple) |>
  summarise(number_mismatch = n(), .groups = "drop") |>
  arrange(desc(number_mismatch))

print.data.frame(mismatches)
sum(mismatches$number_mismatch)
