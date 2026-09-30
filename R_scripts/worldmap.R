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
setwd("../../Database/Taaf/")

## map

# Read the file
data_taaf <- read.csv(
  "ASV_Crozet_Peat_GNA_GBIF - ASV_Crozet_Peat_GNA_GBIF.tsv",
  header = TRUE,
  sep = "\t"
)


# Get the phyloseq object

# First get the otu_matrix
otu_matrix <- data_taaf |>
  select(c(
    ASV_ID,
    X16,
    X18,
    X19,
    X6,
    somme
  )) |>
  tibble::column_to_rownames("ASV_ID") |>
  as.matrix()

## here OTU names are the row so we have to set taxa_are_rows to TRUE
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)

## Get the Taxa matrix
tax_matrix <- data_taaf |>
  select(
    c(
      ASV_ID,
      weight,
      cluster,
      score,
      center,
      bestID,
      ID,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      GenusFG,
      Trophic,
      Guild,
      Sequence,
      taxa_name
    )
  ) |>
  tibble::column_to_rownames("ASV_ID") |>
  as.matrix()

TAX <- tax_table(tax_matrix)

data_taaf_pq <- phyloseq(OTU, TAX)

# Verify names with GBIF Backbone (11):
data_taaf_verified_pq <- gna_verifier_pq(
  data_taaf_pq,
  data_sources = 11
)

data_taaf_occ_pq <- tax_gbif_occur_pq(data_taaf_verified_pq, by_country = TRUE)

#Select only ASV where currentCanonicalSimple is not NA
not_na <- taxa_names(data_taaf_occ_pq)[
  !is.na(data_taaf_occ_pq@tax_table[, "currentCanonicalSimple"])
]

data_taaf_occ_pq_trim <- prune_taxa(taxa = not_na, x = data_taaf_occ_pq)


# First get the otu_table and tax table from the phyloseq object
taxa <- as.data.frame(tax_table(data_taaf_occ_pq_trim))

# Here we add a relative sum column based on the column "somme"
otu_mat <- otu_table(data_taaf_occ_pq_trim) |>
  as.data.frame() |>
  mutate(
    rel_sum = somme / sum(somme)
  ) |>
  as.matrix()

# Here we select the column of the sample we want to apply the function
# WARING: here the column selected is supposed to be relative occurence
otu_mat <- otu_mat[, "rel_sum", drop = FALSE]

head(otu_mat)

# We transpose the otu_mat
sample_mat <- t(otu_mat)
head(sample_mat)
str(sample_mat)

# Select only country datas and get the relative occurance found
occ_mat <- taxa |>
  select(
    -c(
      weight,
      cluster,
      score,
      center,
      bestID,
      ID,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      GenusFG,
      Trophic,
      Guild,
      Sequence,
      taxa_name,
      currentName,
      currentCanonicalSimple,
      genusEpithet,
      genusSpeciesEpithet,
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

# We calculate relative occurences by country for each taxa
rel_occ_mat <- sweep(occ_mat, 1, occ_mat[, "rowsum"], "/")
rel_occ_mat <- rel_occ_mat[, colnames(rel_occ_mat) != "rowsum"]

head(rel_occ_mat)
head(sample_mat)

# we do a mtrice product to get the sum relative abundance of each taxa mutiplied by relative occurences of each taxa by country for the chosen sample (this is why we need to transpose before)
# WARNING, here sample_mat is supposed to be relative abundances
rel_occ_sample <- sample_mat %*% rel_occ_mat
head(rel_occ_sample)

# We transpose it to have countries as rows and the column being the name of the sample
rel_occ_sample <- as.data.frame(t(rel_occ_sample))
head(rel_occ_sample)

# Now we use sf package to get the coordinate of each countries
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

data(countryref)
head(countryref)
centroides <- countryref |>
  select(iso2, lon = centroid.lon, lat = centroid.lat) |>
  filter(!is.na(iso2), !is.na(lon), !is.na(lat)) |>
  distinct(iso2, .keep_all = TRUE)

centroides

# We then join the centroides to each countries (rows)
df <- data.frame(iso2 = row.names(rel_occ_sample), valeur = rel_occ_sample) |>
  left_join(centroides, by = "iso2")

head(df)
# Now we need
ggplot() +
  geom_sf(
    data = world_sf,
    fill = "#ECEFF4",
    colour = "#B0BAC6",
    linewidth = 0.2
  ) +
  geom_point(
    data = df,
    aes(x = lon, y = lat, size = rel_sum, colour = rel_sum),
    alpha = 0.75
  ) + # here we set size and coulour of points with the indice with the name beaing the sample
  scale_size_continuous(
    name = "",
    range = c(4, 20)
  ) +
  scale_colour_viridis_c(name = "", option = "turbo", begin = 0, end = 0.5) +
  coord_sf(xlim = c(-160, 160), ylim = c(-55, 80)) +
  theme_idest() +
  labs(
    title = "Occurrence observée relative des taxons dans gbif pondérée par la présence dans les échantillons"
  ) +
  theme(
    legend.text = element_text(size = 13), # Size of legen labels
    legend.title = element_text(size = 14), # Size of title legend
    legend.key.size = unit(1.5, "cm"), # Size of keys
    legend.spacing.y = unit(0.4, "cm") # Spacing between entries
  ) +
  guides(
    size = guide_legend(reverse = TRUE)
  ) # We set the higher value on top

# we then save the fig
ggsave(
  "occurrence_map_taaf_gbif.png",
  width = 20,
  height = 10
)
