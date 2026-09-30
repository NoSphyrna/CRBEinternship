# ============================================================================ #
# ============ Script to add altitudes of taxa gbif observations ============= #
# ============================================================================ #

# ============== Set Workspace accordingly ======== #
# Get input directory and output directory :

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage : Rscript add_occur.R INPUT_PQ OUTPUT_PQ"
  )
}
INPUT_PQ <- args[1]
OUTPUT_PQ <- args[2]

if (!dir.exists(INPUT_PQ)) {
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
# =============== Workflow for verifying names on OTUs from samples ================ #

getwd()
setwd("../../")
INPUT_PQ <- "Database/data_taxinfo/verified_pq"
# INPUT_TRAITS <- "Database/traitsTable"
OUTPUT_PLOT <- "Database/plot_taxInfo_traits"
OUTPUT_PQ <- "Database/data_taxinfo/alt_pq"
methods <- c("Illumina", "Tedersoo", "Getplage")
sections_Getplage <- c("ITS1", "ITS_full", "ITS_none")
sections_Tedersoo <- c("ITS1", "ITS_full")
sections_Illumina <- c("ITS1")
clusters <- c("sintax", "vsearch")
method <- "Tedersoo"
section <- "ITS_full"
cluster <- "vsearch"
#
clean_guild <- function(data_traits) {
  table_traits <- as.data.frame(tax_table(data_traits))
  new_table <- table_traits |>
    mutate(
      fg_guild = fg_guild |>
        str_squish() |>
        str_to_lower() |>
        str_replace_all(" ", "_") |>
        str_replace_all("\\|", "") |>
        str_replace_all("--", "-")
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

      readRenviron("/home/nosphyrna/.Renviron")
      data_clean_alt <- tax_gbif_alt(data_clean)

      write_pq(
        data_clean_alt,
        path = paste0(
          OUTPUT_PQ,
          "/pq_",
          method,
          "_",
          section,
          "_",
          cluster,
          "_gbif_alt"
        )
      )

      FUNGAL_TRAITS_TABLE <- "Database/traitsTable/FUNGALT_DB_MROY041125.csv"
      data_traits_alt <- fungal_traits_guilds(
        data_clean_alt,
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

      data_traits_alt <- clean_guild(data_traits_alt)
      data_traits_alt <- clean_troph(data_traits_alt)

      tax_table_alt <- as.data.frame(tax_table(data_traits_alt))
      data_traits_alt@tax_table |>
        as.data.frame() |>
        tibble() |>
        filter(as.numeric(altitude_n_records) > 100) |>
        distinct(currentCanonicalSimple, .keep_all = TRUE) |>
        ggplot(aes(
          y = as.numeric(altitude_mean),
          x = currentCanonicalSimple,
          fill = ft_primary_lifestyle
        )) +
        geom_col() +
        coord_flip() +
        geom_errorbar(
          aes(ymin = as.numeric(altitude_q05), ymax = as.numeric(altitude_q95)),
          width = 0.2
        ) +
        geom_label(aes(label = paste0("n=", altitude_n_records)), size = 2) +
        labs(
          title = "Mean altitude with 5%-95% quantiles (only taxa with >100 records)",
          subtitle = "Labels depict the number of gbif records with altitude data, \n
    color depict ecological Guild",
          x = "Taxa names",
          y = "Mean altitude",
          fill = "Fungal traits primary lifestyle"
        ) +
        theme(legend.position = "bottom")

      ggsave(
        paste0(OUTPUT_PLOT, "/test_alt.png"),
        width = 100,
        height = 120,
        units = "cm"
      )
    }
  }
}
