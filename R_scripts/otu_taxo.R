# ============== Libraries ============ #
library(stringr)
library(tibble)
library(stringi)
library(vegan)
library(compositions)
library(zCompositions)
library(tidyverse)
library(phyloseq)
library(conflicted)
conflicts_prefer(dplyr::select)
library(fs)
library(gtools)


getwd()
setwd("../../results/")

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

#exapnd with __ and ;
expand_taxnames_gen <- function(OTU_table, taxonomy, bdd) {
  if (bdd == "euk") {
    OTU_table <- OTU_table |>
      mutate(
        Kingdom = str_match(taxonomy, "[kd]__([^; \n]+)")[, 2],
        Phylum = str_match(taxonomy, "p__([^; \n]+)")[, 2],
        Class = str_match(taxonomy, "c__([^; \n]+)")[, 2],
        Order = str_match(taxonomy, "o__([^; \n]+)")[, 2],
        Family = str_match(taxonomy, "f__([^; \n]+)")[, 2],
        Genus = str_match(taxonomy, "g__([^; \n(]+)")[, 2],
        Species = str_match(taxonomy, "s__([^; \n]+)")[, 2]
      )
  } else {
    OTU_table <- OTU_table |>
      mutate(
        Kingdom = str_match(taxonomy, "[kd]__([^; \n]+)")[, 2],
        Phylum = str_match(taxonomy, "p__([^; \n]+)")[, 2],
        Class = str_match(taxonomy, "c__([^; \n]+)")[, 2],
        Order = str_match(taxonomy, "o__([^; \n]+)")[, 2],
        Family = str_match(taxonomy, "f__([^; \n]+)")[, 2],
        Genus = str_match(taxonomy, "g__([^; \n(]+)")[, 2],
        Species = str_match(taxonomy, "s__[^_]+_([^; \n]+)")[, 2]
      )
  }
  return(OTU_table)
}

expand_sintax_bs_proba <- function(OTU_table, bs_taxo) {
  OTU_table <- OTU_table |>
    mutate(
      BS_Kingdom = as.numeric(str_match(
        bs_taxo,
        "[kd]:[^, \n]+\\(([0-9.]+)\\)"
      )[, 2]),
      BS_Phylum = as.numeric(str_match(bs_taxo, "p:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ]),
      BS_Class = as.numeric(str_match(bs_taxo, "c:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ]),
      BS_Order = as.numeric(str_match(bs_taxo, "o:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ]),
      BS_Family = as.numeric(str_match(bs_taxo, "f:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ]),
      BS_Genus = as.numeric(str_match(bs_taxo, "g:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ]),
      BS_Species = as.numeric(str_match(bs_taxo, "s:[^, \n]+\\(([0-9.]+)\\)")[,
        2
      ])
    )
  return(OTU_table)
}


# organise sample names
replace_name <- function(df) {
  corresp = read.csv(
    "~/Dext2/StageCRBE/Database/corres_pool_single.csv",
    header = TRUE,
    sep = ";"
  )
  lookup_vec <- setNames(corresp$Code, corresp$Single)
  new_col <- str_match(colnames(df), "([A-Za-z0-9]+)(.[0-9]+)?")

  new_col <- paste0(
    if_else(
      new_col[, 2] %in% names(lookup_vec),
      lookup_vec[new_col[, 2]],
      new_col[, 2]
    ),
    replace_na(new_col[, 3], "")
  )
  colnames(df) <- new_col
  sample_cols <- setdiff(new_col, "OTU")
  new_order <- c("OTU", mixedsort(sample_cols))
  return(df[, new_order])
}

# Standardise illumina columns
group_illumina <- function(ill_otu) {
  ill_control <- ill_otu |>
    select(matches("^[PTNB][A-Za-z0-9_]+"))
  ill_otu <- ill_otu |>
    select(-matches("^[PTNB][A-Za-z0-9_]+"))
  colnames(ill_otu) <- str_replace_all(colnames(ill_otu), "_", ".")
  colnames(ill_otu) <- str_replace(colnames(ill_otu), "amplicon", "OTU")

  #sum a and b

  cols_num <- setdiff(colnames(ill_otu), "OTU")
  base_names <- sub("\\.[ab]$", "", cols_num)

  res <- ill_otu |>
    select(all_of(cols_num)) |>
    as.matrix() |>
    t() |>
    rowsum(base_names) |>
    t() |>
    as.data.frame()

  ill_otu <- bind_cols(OTU = ill_otu$OTU, res[, unique(base_names)])

  ill_otu <- replace_name(ill_otu)

  cols_num <- setdiff(colnames(ill_otu), "OTU")
  base_names <- sub("\\..*$", "", cols_num)

  res <- ill_otu |>
    select(all_of(cols_num)) |>
    as.matrix() |>
    t() |>
    rowsum(base_names) |>
    t() |>
    as.data.frame()

  ill_otu <- bind_cols(ill_otu, res[, unique(base_names)])

  keep_cols <- c(
    "OTU",
    "C1",
    "C1.6",
    "C2",
    "C2.12",
    "C3",
    "C3.7",
    "C4",
    "C4.15",
    "C5",
    "C5.5",
    "C6",
    "C6.1",
    "C7",
    "C7.3",
    "C8",
    "C8.6",
    "C9",
    "C9.1",
    "C10",
    "C10.14",
    "C11",
    "C11.16",
    "C12",
    "C12.16",
    "M1",
    "M1.15",
    "M2",
    "M2.10",
    "M3",
    "M3.3",
    "M4",
    "M4.14",
    "M5",
    "M5.11",
    "M6",
    "M6.1",
    "M7",
    "M7.5",
    "M8",
    "M8.5",
    "M9",
    "M9.6",
    "M10",
    "M10.12",
    "M11",
    "M11.16",
    "M12",
    "M12.3",
    "V1",
    "V1.3",
    "V2",
    "V2.5",
    "V3",
    "V3.16",
    "V4",
    "V5",
    "V5.12",
    "V6",
    "V6.5",
    "V7",
    "V7.15",
    "V8",
    "V8.15",
    "V9",
    "V9.4",
    "V10",
    "V10.11",
    "V11",
    "V11.4",
    "V12",
    "V12.2"
  )

  ill_otu <- ill_otu[, keep_cols]
  return(list(ill_otu, ill_control))
}

# file_otu_table = "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt"
# file_taxo_table = "Illumina/taxo/taxonomy_OTU_sintax_unite.tsv"
# techno = "illumina"

get_otu <- function(file_otu_table, techno) {
  # We read the otu files
  nano_otu = read.csv(
    file_otu_table,
    header = TRUE,
    sep = "\t"
  )

  # for pacbio
  if (techno == "pacbio") {
    colnames(nano_otu) <- sub(".ITS9.ITS4", "", colnames(nano_otu))
    ## => little issue, the V4.6 is missing in nanopore (V1440-6 not present in linked_adapters)
    nano_otu <- nano_otu |>
      select(-c("V1440.6"))
  }
  if (techno == "illumina") {
    res <- group_illumina(nano_otu)
    nano_otu <- res[[1]]
    nano_control <- res[[2]]
  } else {
    # Sort and adjust names
    nano_otu <- replace_name(nano_otu)

    # Removal of control samples
    nano_control = nano_otu |>
      select(c(
        "OTU",
        "B1",
        "B2",
        "B3",
        "B4",
        "B5",
        "P1",
        "P2",
        "P3",
        "P4",
        "P5",
        "CLT251",
        "T152",
        "T153",
        "T154",
        "T155",
        "T251",
        "T252",
        "T253",
        "T254",
        "T255",
        "Tneg",
        "Tnegextract.1",
        "Tnegextract.2",
        "Tnegextract.3"
      ))

    nano_otu = nano_otu |>
      select(
        -c(
          "B1",
          "B2",
          "B3",
          "B4",
          "B5",
          "P1",
          "P2",
          "P3",
          "P4",
          "P5",
          "CLT251",
          "T152",
          "T153",
          "T154",
          "T155",
          "T251",
          "T252",
          "T253",
          "T254",
          "T255",
          "Tneg",
          "Tnegextract.1",
          "Tnegextract.2",
          "Tnegextract.3"
        )
      )
  }
  # sum of abundances

  sample_cols <- setdiff(colnames(nano_otu), c("OTU"))
  nano_otu <- nano_otu |>
    mutate(total_abundances = rowSums(across(all_of(sample_cols))))

  sample_cols <- setdiff(colnames(nano_control), c("OTU"))
  nano_control <- nano_control |>
    mutate(total_abundances = rowSums(across(all_of(sample_cols))))

  # remove otus with no abundaces in all samples
  nano_otu <- nano_otu[
    nano_otu$total_abundances > 0,
  ]
  nano_control <- nano_control[
    nano_control$total_abundances > 0,
  ]

  res <- list(nano_otu, nano_control)
  return(res)
}

# VSEARCH (with best match)
# file_otu_table <- "Nanopore/OTU/OTU_table_mumu.tsv"
# file_taxo_table <- "Nanopore/taxo/taxonomy_OTU_vsearch_unite.tsv"
# techno <- "nanopore"
get_otu_taxo_vsearch <- function(file_otu_table, file_taxo_table, techno, bdd) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = FALSE,
    sep = "\t"
  ) |>
    rename("OTU" = "V1", "ID_PERCT" = "V2", "TAXO" = "V3") |>
    group_by(OTU) |>
    slice_max(order_by = ID_PERCT, n = 1, with_ties = FALSE) |>
    ungroup()

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_otu, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$TAXO,
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}

get_control_taxo_vsearch <- function(
  file_otu_table,
  file_taxo_table,
  techno,
  bdd
) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = FALSE,
    sep = "\t"
  ) |>
    rename("OTU" = "V1", "ID_PERCT" = "V2", "TAXO" = "V3") |>
    group_by(OTU) |>
    slice_max(order_by = ID_PERCT, n = 1, with_ties = FALSE) |>
    ungroup()

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_control, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$TAXO,
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}
# SINTAX
get_otu_taxo_sintax <- function(file_otu_table, file_taxo_table, techno, bdd) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = FALSE,
    sep = "\t"
  ) |>
    select(-"V3") |>
    rename("OTU" = "V1", "BS_TAXO" = "V2", "TAXO" = "V4")

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_otu, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$TAXO,
    bdd
  )
  nano_otu_sintax_unite <- expand_sintax_bs_proba(
    nano_otu_sintax_unite,
    nano_otu_sintax_unite$BS_TAXO
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}

# control samples
get_control_taxo_sintax <- function(
  file_otu_table,
  file_taxo_table,
  techno,
  bdd
) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = FALSE,
    sep = "\t"
  ) |>
    select(-"V3") |>
    rename("OTU" = "V1", "BS_TAXO" = "V2", "TAXO" = "V4")

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_control, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$TAXO,
    bdd
  )
  nano_otu_sintax_unite <- expand_sintax_bs_proba(
    nano_otu_sintax_unite,
    nano_otu_sintax_unite$BS_TAXO
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}
# DNABARCODER
# file_otu_table <- "Nanopore/OTU/OTU_table_mumu.tsv"
# file_taxo_table <- "Nanopore/taxo/rep_seqs.unite2025ITS_BLAST.classified"
# techno <- "nanopore"
get_otu_taxo_dnabarcoder <- function(
  file_otu_table,
  file_taxo_table,
  techno,
  bdd
) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = TRUE,
    sep = "\t"
  ) |>
    rename("OTU" = "ID")

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_otu, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames_gen(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$Full.classification,
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}

get_control_taxo_dnabarcoder <- function(
  file_otu_table,
  file_taxo_table,
  techno,
  bdd
) {
  res <- get_otu(file_otu_table, techno)
  nano_otu <- res[[1]]
  nano_control <- res[[2]]

  # read taxo
  nano_sintax_unite = read.csv(
    file_taxo_table,
    header = TRUE,
    sep = "\t"
  ) |>
    rename("OTU" = "ID")

  if (techno == "illumina") {
    nano_sintax_unite <- nano_sintax_unite |>
      mutate(OTU = sub(";.*$", "", OTU))
  }
  # merging and expanding of taxa
  nano_otu_sintax_unite <- merge(nano_control, nano_sintax_unite, by = "OTU")

  nano_otu_sintax_unite <- expand_taxnames_gen(
    nano_otu_sintax_unite,
    taxonomy = nano_otu_sintax_unite$Full.classification,
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

  # replacing "unidentified by NA"

  nano_otu_sintax_unite <- nano_otu_sintax_unite |>
    mutate(across(all_of(tax_cols), ~ na_if(.x, "unidentified")))

  return(nano_otu_sintax_unite)
}

get_vsearch_summary <- function(
  nano_otu_sintax_unite,
  summary_path,
  summary_name
) {
  tax_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  meta_cols <- c(
    "OTU",
    "ID_PERCT",
    "TAXO",
    "total_abundances",
    tax_cols
  )

  sample_cols <- setdiff(colnames(nano_otu_sintax_unite), meta_cols)

  # str(nano_otu_sintax_unite)

  # percentage of OTU assigned at each ranks
  OTU_assigned <- sapply(nano_otu_sintax_unite[tax_cols], function(x) {
    mean(!is.na(x)) * 100
  })
  # OTU_assigned

  #percentage of OTU assigned ponderatede by abundances
  total_reads <- sum(nano_otu_sintax_unite$total_abundances)

  OTU_abundances_assigned <- sapply(
    tax_cols,
    function(col) {
      assigned <- !is.na(nano_otu_sintax_unite[[col]])
      sum(nano_otu_sintax_unite$total_abundances[assigned]) / total_reads * 100
    }
  )

  # OTU_abundances_assigned

  # mean bs rate of all (include below 0.51)
  mean_id_percent <- mean(nano_otu_sintax_unite$ID_PERCT, na.rm = TRUE)
  # mean_bs_OTU

  # mean bs rate of all ponderd by abundances

  mean_id_perct_abundances <- weighted.mean(
    nano_otu_sintax_unite$ID_PERCT,
    nano_otu_sintax_unite$total_abundances,
    na.rm = TRUE
  )

  # mean_bs_OTU_abundances

  # summary table

  summary_table <- data.frame(
    percent_OTU_assigned = OTU_assigned,
    precent_read_assigned = OTU_abundances_assigned,
    id_mean_OTU = mean_id_percent,
    id_mean_reads = mean_id_perct_abundances
  )
  # summary_table

  #Save summary table
  dir_create(summary_path)
  write.csv(
    summary_table,
    paste0(summary_path, summary_name),
    row.names = TRUE
  )
}

get_sintax_summary <- function(
  nano_otu_sintax_unite,
  summary_path,
  summary_name
) {
  tax_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  bs_cols <- c(
    "BS_Kingdom",
    "BS_Phylum",
    "BS_Class",
    "BS_Order",
    "BS_Family",
    "BS_Genus",
    "BS_Species"
  )
  meta_cols <- c(
    "OTU",
    "BS_TAXO",
    "TAXO",
    "total_abundances",
    tax_cols,
    bs_cols
  )

  sample_cols <- setdiff(colnames(nano_otu_sintax_unite), meta_cols)

  # str(nano_otu_sintax_unite)

  # percentage of OTU assigned at each ranks
  OTU_assigned <- sapply(nano_otu_sintax_unite[tax_cols], function(x) {
    mean(!is.na(x)) * 100
  })
  # OTU_assigned

  #percentage of OTU assigned ponderatede by abundances
  total_reads <- sum(nano_otu_sintax_unite$total_abundances)

  OTU_abundances_assigned <- sapply(
    tax_cols,
    function(col) {
      assigned <- !is.na(nano_otu_sintax_unite[[col]])
      sum(nano_otu_sintax_unite$total_abundances[assigned]) / total_reads * 100
    }
  )

  # OTU_abundances_assigned

  # mean bs rate of all (include below 0.51)
  mean_bs_OTU <- sapply(nano_otu_sintax_unite[bs_cols], mean, na.rm = TRUE)
  # mean_bs_OTU

  # mean bs rate of all ponderd by abundances

  mean_bs_OTU_abundances <- sapply(bs_cols, function(col) {
    weighted.mean(
      nano_otu_sintax_unite[[col]],
      nano_otu_sintax_unite$total_abundances,
      na.rm = TRUE
    )
  })
  # mean_bs_OTU_abundances

  # mean bs rate but only for assigned
  mean_bs_OTU_assigned <- sapply(seq_along(bs_cols), function(i) {
    mean(
      nano_otu_sintax_unite[[bs_cols[i]]][
        !is.na(nano_otu_sintax_unite[[tax_cols[i]]])
      ],
      na.rm = TRUE
    )
  })

  names(mean_bs_OTU_assigned) <- tax_cols
  # mean_bs_OTU_assigned

  # mean bs rate but only for assigned
  mean_bs_OTU_assigned_abundances <- sapply(seq_along(bs_cols), function(i) {
    weighted.mean(
      nano_otu_sintax_unite[[bs_cols[i]]][
        !is.na(nano_otu_sintax_unite[[tax_cols[i]]])
      ],
      nano_otu_sintax_unite$total_abundances[
        !is.na(nano_otu_sintax_unite[[tax_cols[i]]])
      ],
      na.rm = TRUE
    )
  })

  names(mean_bs_OTU_assigned_abundances) <- tax_cols
  # mean_bs_OTU_assigned_abundances

  # summary table

  summary_table <- data.frame(
    percent_OTU_assigned = OTU_assigned,
    precent_read_assigned = OTU_abundances_assigned,
    bs_mean_OTU = mean_bs_OTU,
    bs_mean_reads = mean_bs_OTU_abundances,
    bs_mean_OTU_above_threshold = mean_bs_OTU_assigned,
    bs_mean_reads_above_threshold = mean_bs_OTU_assigned_abundances
  )
  # summary_table

  #Save summary table
  dir_create(summary_path)
  write.csv(
    summary_table,
    paste0(summary_path, summary_name),
    row.names = TRUE
  )
}

get_dnabarcoder_summary <- function(
  nano_otu_sintax_unite,
  summary_path,
  summary_name
) {
  tax_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )
  meta_cols <- c(
    "OTU",
    "total_abundances",
    "Given.label",
    "Prediction",
    "Full.classification",
    "Rank",
    "Cut.off",
    "Confidence",
    "ReferenceID",
    "BLAST.score",
    "BLAST.sim",
    "BLAST.coverage",
    tax_cols
  )

  sample_cols <- setdiff(colnames(nano_otu_sintax_unite), meta_cols)

  # str(nano_otu_sintax_unite)

  # percentage of OTU assigned at each ranks
  OTU_assigned <- sapply(nano_otu_sintax_unite[tax_cols], function(x) {
    mean(!is.na(x)) * 100
  })
  # OTU_assigned

  #percentage of OTU assigned ponderatede by abundances
  total_reads <- sum(nano_otu_sintax_unite$total_abundances)

  OTU_abundances_assigned <- sapply(
    tax_cols,
    function(col) {
      assigned <- !is.na(nano_otu_sintax_unite[[col]])
      sum(nano_otu_sintax_unite$total_abundances[assigned]) / total_reads * 100
    }
  )

  # OTU_abundances_assigned

  # mean blast score
  mean_blastscore <- mean(nano_otu_sintax_unite$BLAST.score, na.rm = TRUE)
  # mean_bs_OTU

  # mean blast score ponderd by abundances
  mean_blastscore_abundances <- weighted.mean(
    nano_otu_sintax_unite$BLAST.score,
    nano_otu_sintax_unite$total_abundances,
    na.rm = TRUE
  )

  # mean blast score
  mean_confidence <- mean(
    as.numeric(nano_otu_sintax_unite$Confidence),
    na.rm = TRUE
  )
  # mean_bs_OTU

  # mean blast score ponderd by abundances
  mean_confidence_abundances <- weighted.mean(
    as.numeric(nano_otu_sintax_unite$Confidence),
    nano_otu_sintax_unite$total_abundances,
    na.rm = TRUE
  )

  # mean blast score
  mean_blastcoverage <- mean(nano_otu_sintax_unite$BLAST.coverage, na.rm = TRUE)
  # mean_bs_OTU

  # mean blast score ponderd by abundances
  mean_blastcoverage_abundances <- weighted.mean(
    nano_otu_sintax_unite$BLAST.coverage,
    nano_otu_sintax_unite$total_abundances,
    na.rm = TRUE
  )

  # summary table
  summary_table <- data.frame(
    percent_OTU_assigned = OTU_assigned,
    precent_read_assigned = OTU_abundances_assigned,
    mean_blastscore = mean_blastscore,
    weighted_mean_blastscore = mean_blastscore_abundances,
    mean_confidence = mean_confidence,
    weighted_mean_confidence = mean_confidence_abundances,
    mean_blastcoverage = mean_blastcoverage,
    weighted_mean_blastcoverage = mean_blastcoverage_abundances
  )
  # summary_table

  #Save summary table
  dir_create(summary_path)
  write.csv(
    summary_table,
    paste0(summary_path, summary_name),
    row.names = TRUE
  )
}
## ============== workflow ===================== ##

## ========== Nanopore Sintax ================ ##
nano_otu_sintax_unite <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_unite.tsv",
  "nanopore",
  "unite"
)

nano_otu_sintax_euk <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_euk.tsv",
  "nanopore",
  "euk"
)

nb_OTUs <- nrow(nano_otu_sintax_euk)
nb_OTUs_nanopore
total_reads <- sum(nano_otu_sintax_unite$total_abundances)
total_reads_nanopore

get_sintax_summary(
  nano_otu_sintax_unite,
  "Nanopore/stats/",
  "summary_sintax_unite.csv"
)
get_sintax_summary(
  nano_otu_sintax_euk,
  "Nanopore/stats/",
  "summary_sintax_euk.csv"
)

## =========== Nanopore vsearch ============ ##
nano_otu_vsearch_unite <- get_otu_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "nanopore",
  "unite"
)
nano_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "nanopore",
  "euk"
)
get_vsearch_summary(
  nano_otu_vsearch_unite,
  "Nanopore/stats/",
  "summary_vsearch_unite.csv"
)
get_vsearch_summary(
  nano_otu_vsearch_euk,
  "Nanopore/stats/",
  "summary_vsearch_euk.csv"
)

## =========== Nanopore dnabarcoder ============ ##
nano_otu_dnabarcoder <- get_otu_taxo_dnabarcoder(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/rep_seqs.unite2025ITS_BLAST.classified",
  "nanopore",
  "unite"
)
get_dnabarcoder_summary(
  nano_otu_dnabarcoder,
  "Nanopore/stats/",
  "summary_dnabarcoder.csv"
)

## ========== PacBio Sintax ================ ##
pacbio_otu_sintax_unite <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_otu_sintax_euk <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_euk.tsv",
  "pacbio",
  "euk"
)

pacbio_nb_OTUs <- nrow(pacbio_otu_sintax_unite)
pacbio_nb_OTUs
pacbio_total_reads <- sum(pacbio_otu_sintax_unite$total_abundances)
pacbio_total_reads

get_sintax_summary(
  pacbio_otu_sintax_unite,
  "PacBio/stats/",
  "summary_sintax_unite.csv"
)
get_sintax_summary(
  pacbio_otu_sintax_euk,
  "PacBio/stats/",
  "summary_sintax_euk.csv"
)

## =========== PacBio vsearch ============ ##
pacbio_otu_vsearch_unite <- get_otu_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "pacbio",
  "euk"
)

get_vsearch_summary(
  pacbio_otu_vsearch_unite,
  "PacBio/stats/",
  "summary_vsearch_unite.csv"
)
get_vsearch_summary(
  pacbio_otu_vsearch_euk,
  "PacBio/stats/",
  "summary_vsearch_euk.csv"
)
## =========== PacBio dnabarcoder ============ ##
pacbio_otu_dnabarcoder <- get_otu_taxo_dnabarcoder(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/OTUs_LULU.unite2025ITS_BLAST.classified",
  "pacbio",
  "unite"
)
get_dnabarcoder_summary(
  pacbio_otu_dnabarcoder,
  "PacBio/stats/",
  "summary_dnabarcoder.csv"
)
## ========= Illumina sintax =============== ##
ill_otu_sintax_unite <- get_otu_taxo_sintax(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/taxonomy_OTU_sintax_unite.tsv",
  "illumina",
  "unite"
)
ill_otu_sintax_euk <- get_otu_taxo_sintax(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/taxonomy_OTU_sintax_euk.tsv",
  "illumina",
  "euk"
)

ill_nb_OTUs <- nrow(ill_otu_sintax_unite)
ill_nb_OTUs
ill_total_reads <- sum(ill_otu_sintax_unite$total_abundances)
ill_total_reads

get_sintax_summary(
  ill_otu_sintax_unite,
  "Illumina/stats/",
  "summary_sintax_unite.csv"
)
get_sintax_summary(
  ill_otu_sintax_euk,
  "Illumina/stats/",
  "summary_sintax_euk.csv"
)
## =========== Illumina dnabarcoder ============ ##
ill_otu_dnabarcoder <- get_otu_taxo_dnabarcoder(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/reads_OTU_nonchimeras.unite2025ITS1_BLAST.classified",
  "illumina",
  "unite"
)
get_dnabarcoder_summary(
  ill_otu_dnabarcoder,
  "Illumina/stats/",
  "summary_dnabarcoder.csv"
)

## ================ rarefaction ================ ##

pacbio_otu <- get_otu(
  "PacBio/OTU/OTU_table_LULU.txt",
  "pacbio"
)[[1]]
nano_otu <- get_otu(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "nanopore"
)[[1]]
ill_otu <- get_otu(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "illumina"
)[[1]]

pacbio_otu_matrix <- pacbio_otu |>
  remove_rownames() |>
  column_to_rownames("OTU") |>
  select(-c("total_abundances")) |>
  as.matrix()

get_rarecurve <- function(otu_table, techname, step = 200) {
  otu_matrix <- otu_table |>
    remove_rownames() |>
    column_to_rownames("OTU") |>
    select(-c("total_abundances")) |>
    as.matrix() |>
    t()
  raw <- rarecurve(otu_matrix, step = step, tidy = TRUE)
  raw$Techno <- techname
  return(raw)
}

ill_curve <- get_rarecurve(ill_otu, "Illumina")
pacbio_curve <- get_rarecurve(pacbio_otu, "PacBio")
nano_curve <- get_rarecurve(nano_otu, "Nanopore")

all_curve <- bind_rows(ill_curve, pacbio_curve, nano_curve)
all_curve <- all_curve |>
  mutate(
    SampleType = if_else(str_detect(Site, "\\."), "Single", "Pool"),
    Techno_Type = paste(Techno, SampleType, sep = " - ")
  )
ggplot(
  all_curve,
  aes(
    x = Sample,
    y = Species,
    group = interaction(Site, Techno),
    color = Techno_Type
  )
) +
  geom_line(linewidth = 1, alpha = 0.8) +
  facet_wrap(~Techno, nrow = 1) +
  scale_color_manual(
    values = c(
      "Illumina - Pool" = "darkgoldenrod",
      "Illumina - Single" = "darkgoldenrod1",
      "PacBio - Pool" = "orchid4",
      "PacBio - Single" = "orchid1",
      "Nanopore - Pool" = "steelblue4",
      "Nanopore - Single" = "steelblue1"
    )
  ) +
  labs(
    x = "Number of reads",
    y = "Number of observed OTU",
    color = "Techno / Type"
  ) +
  xlim(0, 1e5) +
  ylim(0, 5000) +
  theme_minimal()


ggsave(
  "~/Dext2/StageCRBE/results/fig/rarefaction_lim.png",
  width = 20,
  height = 10
)

## ====== Venn diagrams ======= #

library(VennDiagram)

nano_otu_sintax_unite <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_unite.tsv",
  "nanopore",
  "unite"
)

nano_otu_sintax_euk <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_euk.tsv",
  "nanopore",
  "euk"
)

nano_otu_vsearch_unite <- get_otu_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "nanopore",
  "unite"
)
nano_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "nanopore",
  "euk"
)

nano_otu_dnabarcoder <- get_otu_taxo_dnabarcoder(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/rep_seqs.unite2025ITS_BLAST.classified",
  "nanopore",
  "unite"
)
X <- list(
  SintaxEuk = nano_otu_sintax_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  SintaxUnite = nano_otu_sintax_unite |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  VsearchEuk = nano_otu_vsearch_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  VsearchUnite = nano_otu_vsearch_unite |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  DnabarcoderUnite = nano_otu_dnabarcoder |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique()
)

pacbio_otu_sintax_unite <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_otu_sintax_euk <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_euk.tsv",
  "pacbio",
  "euk"
)
pacbio_otu_vsearch_unite <- get_otu_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "pacbio",
  "euk"
)
pacbio_otu_dnabarcoder <- get_otu_taxo_dnabarcoder(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/OTUs_LULU.unite2025ITS_BLAST.classified",
  "pacbio",
  "unite"
)
X <- list(
  SintaxEuk = pacbio_otu_sintax_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  SintaxUnite = pacbio_otu_sintax_unite |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  VsearchEuk = pacbio_otu_vsearch_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  VsearchUnite = pacbio_otu_vsearch_unite |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  DnabarcoderUnite = pacbio_otu_dnabarcoder |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique()
)
grid.newpage()

venn_plot <- venn.diagram(
  x = X,
  category.names = c(
    "Sintax Euk",
    "Sintax Unite",
    "Vsearch Euk",
    "Vsearch Unite",
    "Dnabarcoder Unite"
  ),
  filename = NULL,
  output = TRUE,
  imagetype = "png",

  lwd = 1,
  col = c("orange", 'orangered', 'orchid', 'mediumpurple', 'palegreen'),
  fill = c(
    alpha("orange", 0.5),
    alpha("orangered", 0.5),
    alpha('orchid', 0.5),
    alpha("mediumpurple", 0.5),
    alpha("palegreen", 0.5)
  ),

  cex = 0.5,
  fontfamily = "sans",
  ext.text = FALSE,

  cat.cex = 0.3,
  # cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("orange", 'orangered', 'orchid', 'mediumpurple', 'palegreen')
)


grid.draw(venn_plot)

png(
  'fig/venn_pacbio_genus.png',
  height = 1080,
  width = 1920,
  res = 300
)

grid.draw(venn_plot)
dev.off()

## Venn diagrams techno

nano_otu_sintax_euk <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_euk.tsv",
  "nanopore",
  "euk"
)
pacbio_otu_sintax_euk <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_euk.tsv",
  "pacbio",
  "euk"
)
ill_otu_sintax_euk <- get_otu_taxo_sintax(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/taxonomy_OTU_sintax_euk.tsv",
  "illumina",
  "euk"
)

X <- list(
  Illumina = ill_otu_sintax_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  Nanopore = nano_otu_sintax_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  PacBio = pacbio_otu_sintax_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique()
)
grid.newpage()

venn_plot <- venn.diagram(
  x = X,
  category.names = c(
    "Illumina",
    "Nanopore",
    "PacBio"
  ),
  filename = NULL,
  output = TRUE,
  imagetype = "png",

  lwd = 1,
  col = c('darkgoldenrod1', 'steelblue', 'orchid'),
  fill = c(
    alpha("darkgoldenrod", 0.5),
    alpha("steelblue", 0.5),
    alpha("orchid", 0.5)
  ),

  cex = 0.5,
  fontfamily = "sans",
  ext.text = FALSE,

  cat.cex = 0.3,
  # cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c('darkgoldenrod', 'steelblue', 'orchid')
)


grid.draw(venn_plot)

png(
  'fig/venn_techno_genus.png',
  height = 1080,
  width = 1920,
  res = 300
)

grid.draw(venn_plot)
dev.off()
# vsearch
nano_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "nanopore",
  "euk"
)
pacbio_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "pacbio",
  "euk"
)
ill_otu_vsearch_euk <- get_otu_taxo_vsearch(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "illumina",
  "euk"
)
X <- list(
  Illumina = ill_otu_vsearch_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  Nanopore = nano_otu_vsearch_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique(),
  PacBio = pacbio_otu_vsearch_euk |>
    dplyr::filter(!is.na(Genus)) |>
    pull(Genus) |>
    unique()
)
grid.newpage()

venn_plot <- venn.diagram(
  x = X,
  category.names = c(
    "Illumina",
    "Nanopore",
    "PacBio"
  ),
  filename = NULL,
  output = TRUE,
  imagetype = "png",

  lwd = 1,
  col = c('darkgoldenrod1', 'steelblue', 'orchid'),
  fill = c(
    alpha("darkgoldenrod", 0.5),
    alpha("steelblue", 0.5),
    alpha("orchid", 0.5)
  ),

  cex = 0.5,
  fontfamily = "sans",
  ext.text = FALSE,

  cat.cex = 0.3,
  # cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c('darkgoldenrod', 'steelblue', 'orchid')
)


grid.draw(venn_plot)

png(
  'fig/venn_techno_genus.png',
  height = 1080,
  width = 1920,
  res = 300
)

grid.draw(venn_plot)
dev.off()


## ============ control samples ============== #
pacbio_control_sintax_unite <- get_control_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_control_sintax_euk <- get_control_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_euk.tsv",
  "pacbio",
  "euk"
)
pacbio_control_vsearch_unite <- get_control_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "pacbio",
  "unite"
)

pacbio_control_vsearch_euk <- get_control_taxo_vsearch(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "pacbio",
  "euk"
)
pacbio_control_dnabarcoder <- get_control_taxo_dnabarcoder(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/OTUs_LULU.unite2025ITS_BLAST.classified",
  "pacbio",
  "unite"
)

nano_control_sintax_unite <- get_control_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_unite.tsv",
  "nanopore",
  "unite"
)

nano_control_sintax_euk <- get_control_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_euk.tsv",
  "nanopore",
  "euk"
)

nano_control_vsearch_unite <- get_control_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_unite.tsv",
  "nanopore",
  "unite"
)
nano_control_vsearch_euk <- get_control_taxo_vsearch(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_vsearch_euk.tsv",
  "nanopore",
  "euk"
)

nano_control_dnabarcoder <- get_control_taxo_dnabarcoder(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/rep_seqs.unite2025ITS_BLAST.classified",
  "nanopore",
  "unite"
)

plot_control_abundance_per_sample <- function(
  otu_table,
  title_end
) {
  otu_table <- otu_table |>
    dplyr::mutate(
      control = dplyr::case_when(
        Genus == "Pleurotus" & Species == "eryngii" ~ "Pleurotus eryngii",
        Genus == "Pleurotus" ~ "Pleurotus",
        is.na(Genus) ~ "Unidentified",
        TRUE ~ "Assigned"
      )
    ) |>
    select(c(
      "control",
      "B1",
      "B2",
      "B3",
      "B4",
      "B5",
      "P1",
      "P2",
      "P3",
      "P4",
      "P5",
      "CLT251",
      "T152",
      "T153",
      "T154",
      "T155",
      "T251",
      "T252",
      "T253",
      "T254",
      "T255",
      "Tneg",
      "Tnegextract.1",
      "Tnegextract.2",
      "Tnegextract.3"
    )) |>
    group_by(control) |>
    summarise(across(everything(), sum), .groups = 'drop') |>
    pivot_longer(!control, names_to = "Sample", values_to = "abundance")

  # Change color of Unknown to grey
  cols <- c(
    "Assigned" = "royalblue",
    "Pleurotus" = "chartreuse",
    "Pleurotus eryngii" = "darkgoldenrod1",
    "Unidentified" = "grey70"
  )
  cat("Plotting graph\n")
  p <- ggplot(
    otu_table,
    aes(x = Sample, y = abundance, fill = control)
  ) +
    # here it's important to put identity so the bars correspond to the relative abundance
    # in each samples and not the number of different abundance found in each sample
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = cols
    ) +
    # scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = paste("Genus and Pleurotus assignation", title_end),
      x = "Sample",
      y = "Abundance (number of reads)",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )
  return(p)
}

plot_control_abundance_per_sample(
  pacbio_control_sintax_euk,
  "Pacbio Sintax Eukaryome"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_pacbio_sintax_euk.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  pacbio_control_sintax_unite,
  "Pacbio Sintax Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_pacbio_sintax_unite.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  pacbio_control_vsearch_euk,
  "Pacbio Vsearch Eukaryome"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_pacbio_vsearch_euk.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  pacbio_control_vsearch_unite,
  "Pacbio Vsearch Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_pacbio_vsearch_unite.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  pacbio_control_dnabarcoder,
  "Pacbio Dnabarcoder Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_pacbio_dnabarcoder_unite.png",
  width = 20,
  height = 10
)
# Nanopore
plot_control_abundance_per_sample(
  nano_control_sintax_euk,
  "Nanopore Sintax Eukaryome"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_nano_sintax_euk.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  nano_control_sintax_unite,
  "Nanopore Sintax Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_nano_sintax_unite.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  nano_control_vsearch_euk,
  "Nanopore Vsearch Eukaryome"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_nano_vsearch_euk.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  nano_control_vsearch_unite,
  "Nanopore Vsearch Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_nano_vsearch_unite.png",
  width = 20,
  height = 10
)
plot_control_abundance_per_sample(
  nano_control_dnabarcoder,
  "Nanopore Dnabarcoder Unite"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/control_abundance_nano_dnabarcoder_unite.png",
  width = 20,
  height = 10
)

## ====== PCoA Aitchison, Bray-Curtis, permANOVA ==== ##

## On prend sintax euk = assignation la plus prudente
## On pourrait aussi prendre vsearch euk (plus consensuelle)
nano_otu_sintax_euk <- get_otu_taxo_sintax(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "Nanopore/taxo/taxonomy_OTU_sintax_euk.tsv",
  "nanopore",
  "euk"
)
pacbio_otu_sintax_euk <- get_otu_taxo_sintax(
  "PacBio/OTU/OTU_table_LULU.txt",
  "PacBio/taxo/taxonomy_OTU_sintax_euk.tsv",
  "pacbio",
  "euk"
)
ill_otu_sintax_euk <- get_otu_taxo_sintax(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "Illumina/taxo/taxonomy_OTU_sintax_euk.tsv",
  "illumina",
  "euk"
)

# test_nano <- nano_otu_sintax_euk |>
#   dplyr::filter(Kingdom == "Fungi")
# rm(test_nano)
# unique(test_nano$Kingdom)
# unique(nano_otu_sintax_euk$Kingdom)
# Here we prepare the tables with only Genus and abundances
# We also filter to get only fungi
prep_table <- function(otu_table, techno) {
  otu_table <- otu_table |>
    dplyr::filter(Kingdom == "Fungi") |>
    dplyr::select(c(
      "Genus",
      "C1",
      "C1.6",
      "C2",
      "C2.12",
      "C3",
      "C3.7",
      "C4",
      "C4.15",
      "C5",
      "C5.5",
      "C6",
      "C6.1",
      "C7",
      "C7.3",
      "C8",
      "C8.6",
      "C9",
      "C9.1",
      "C10",
      "C10.14",
      "C11",
      "C11.16",
      "C12",
      "C12.16",
      "M1",
      "M1.15",
      "M2",
      "M2.10",
      "M3",
      "M3.3",
      "M4",
      "M4.14",
      "M5",
      "M5.11",
      "M6",
      "M6.1",
      "M7",
      "M7.5",
      "M8",
      "M8.5",
      "M9",
      "M9.6",
      "M10",
      "M10.12",
      "M11",
      "M11.16",
      "M12",
      "M12.3",
      "V1",
      "V1.3",
      "V2",
      "V2.5",
      "V3",
      "V3.16",
      "V4",
      "V5",
      "V5.12",
      "V6",
      "V6.5",
      "V7",
      "V7.15",
      "V8",
      "V8.15",
      "V9",
      "V9.4",
      "V10",
      "V10.11",
      "V11",
      "V11.4",
      "V12",
      "V12.2"
    )) |>
    group_by(Genus) |>
    summarise(across(everything(), sum), .groups = 'drop') |>
    dplyr::filter(!is.na(Genus)) |>
    pivot_longer(!Genus, names_to = "Sample", values_to = "abundance") |>
    mutate(Technology = techno)
  return(otu_table)
}

# We then binds the data in a single dataframe (long)
techno_long <- bind_rows(
  prep_table(ill_otu_sintax_euk, "Illumina"),
  prep_table(nano_otu_sintax_euk, "Nanopore"),
  prep_table(pacbio_otu_sintax_euk, "PacBio")
)

# we move Genus to the columns and make a row id as Techno_Sample
techno_wide <- techno_long |>
  dplyr::mutate(RowID = paste(Technology, Sample, sep = "_")) |>
  dplyr::select(RowID, Genus, abundance) |>
  pivot_wider(names_from = Genus, values_from = abundance, values_fill = 0) |>
  column_to_rownames("RowID")

# We keep a metadata tha binds Technology, Sample and RowID for paired permANOVA

metadata <- techno_long |>
  dplyr::distinct(Technology, Sample) |>
  dplyr::mutate(RowID = paste(Technology, Sample, sep = "_")) |>
  column_to_rownames("RowID")

# order to be sure it's the same
metadata <- metadata[rownames(techno_wide), ]

#clr transform + zero handleing with geometric bayesian multiplicative (GBM)
# We use the the pseudo count output to keep abundances and not relative abundances so the clr is not modified after
# We set z.delete to false to keep all data even if it results in a biased imputation

# techno_wide_filt <- techno_wide[, colSums(techno_wide) > 0]
# rm(techno_wide_filt)

# we filter to get at least to lines that are not 0 for each genera kept so that the bayesian inference can be at least done
min_prevalence <- 2
techno_wide_filt <- techno_wide[, colSums(techno_wide > 0) >= min_prevalence]

techno_wide_nozero <- cmultRepl(
  techno_wide_filt,
  method = "GBM",
  output = "p-counts",
  z.warning = 0.95,
  z.delete = FALSE
)

clr_techno <- clr(techno_wide_nozero)
clr_mat <- as.matrix(clr_techno)

# Ensuite on réalise la méthode d'Aitchison qui est une distance euclidienne sur la matrice clr :
aitch_dist <- vegdist(clr_mat, method = "euclidian")

# We want to test also the Bray-Curtis meric but we then have to normalise our matrix to have relative abundance using vegan decostand function
techno_wide_rel <- decostand(techno_wide_filt, method = "total")

#Then the Bray Curtis dist
bray_dist <- vegdist(techno_wide_rel, method = "bray")

# ========== PCoA Aitchison  ======== #
## We prepare the PCoA (k = 2 (dimensions), eig= TRUE(valeurs propres retournées)
pcoa_aitch <- cmdscale(aitch_dist, k = 2, eig = TRUE)
# We calculate the percentage of explained variance with the eigen values of both axes of the PCoA
var_expl_aitch <- round(
  100 * pcoa_aitch$eig[1:2] / sum(pcoa_aitch$eig[pcoa_aitch$eig > 0]),
  digits = 1
)

# we prepare the graph and add the sample and Techno data
graph_aitch <- as.data.frame(pcoa_aitch$points) |>
  setNames(c("Axis1", "Axis2")) |>
  rownames_to_column("RowID") |>
  left_join(metadata |> rownames_to_column("RowID"), by = "RowID")

ggplot(graph_aitch, aes(Axis1, Axis2, color = Technology)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_manual(
    values = c(
      "Illumina" = "darkgoldenrod1",
      "PacBio" = "orchid1",
      "Nanopore" = "steelblue1"
    )
  ) +
  theme_minimal() +
  labs(
    title = "PCoA Aitchison",
    x = "Axis1",
    y = "Axis2",
    subtitle = paste0(
      "Axis1=",
      var_expl_aitch[1],
      "%, Axis2=",
      var_expl_aitch[2],
      "%"
    )
  )

ggsave(
  "~/Dext2/StageCRBE/results/fig/PCoA_aitchison_techno_sintax_euk.png",
  width = 20,
  height = 10
)

# ========= PCoA Bray-Curtis ============= #

## We prepare the PCoA (k = 2 (dimensions), eig= TRUE(valeurs propres retournées)
pcoa_bray <- cmdscale(bray_dist, k = 2, eig = TRUE)
# We calculate the percentage of explained variance with the eigen values of both axes of the PCoA
var_expl_bray <- round(
  100 * pcoa_bray$eig[1:2] / sum(pcoa_bray$eig[pcoa_bray$eig > 0]),
  digits = 1
)

# we prepare the graph and add the sample and Techno data
graph_bray <- as.data.frame(pcoa_bray$points) |>
  setNames(c("Axis1", "Axis2")) |>
  rownames_to_column("RowID") |>
  left_join(metadata |> rownames_to_column("RowID"), by = "RowID")

ggplot(graph_bray, aes(Axis1, Axis2, color = Technology)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_manual(
    values = c(
      "Illumina" = "darkgoldenrod1",
      "PacBio" = "orchid1",
      "Nanopore" = "steelblue1"
    )
  ) +
  theme_minimal() +
  labs(
    title = "PCoA Bray-Curtis",
    x = "Axis1",
    y = "Axis2",
    subtitle = paste0(
      "Axis1=",
      var_expl_bray[1],
      "%, Axis2=",
      var_expl_bray[2],
      "%"
    )
  )

ggsave(
  "~/Dext2/StageCRBE/results/fig/PCoA_bray_techno_sintax_euk.png",
  width = 20,
  height = 10
)

# ============= PermANOVA =============== #

# Aitchison

# Here we will be cautiaus that data is paired so we need to restrict the permutations to only the Technologies (not permute the samples)
param_perm <- how(nperm = 999, blocks = metadata$Sample)

permanova_aitch <- adonis2(
  aitch_dist ~ Technology,
  data = metadata,
  permutations = param_perm
)

permanova_aitch
# Permutation test for adonis under reduced model
# Blocks:  metadata$Sample
# Permutation: free
# Number of permutations: 999
#
# adonis2(formula = aitch_dist ~ Technology, data = metadata, permutations = param_perm)
#           Df SumOfSqs      R2      F Pr(>F)
# Model      2    72688 0.11908 14.193  0.001 ***
# Residual 210   537729 0.88092
# Total    212   610417 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# To interpret Correctly, we need to do a homogeneity dispersion test

disp_test_aitch <- betadisper(aitch_dist, metadata$Technology)
permutest(disp_test_aitch, permutations = param_perm)
#
# Permutation test for homogeneity of multivariate dispersions
# Blocks:  metadata$Sample
# Permutation: free
# Number of permutations: 999
#
# Response: Distances
#            Df Sum Sq Mean Sq     F N.Perm Pr(>F)
# Groups      2  123.1  61.562 2.032    999  0.001 ***
# Residuals 210 6362.2  30.296
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# So we can say that permANOVA Technology explain 12% of variance in fungal community composition (R2= 0.119, P=0.001) but this result as to be taken cautiously because the result of the betadisper test can be interpreted as dispersion is significantly heterogenous (P=0.001)  So permANOVA migth just detect heterogeneity.

# Bray-Curtis

# Here we will be cautiaus that data is paired so we need to restrict the permutations to only the Technologies (not permute the samples)
param_perm <- how(nperm = 999, blocks = metadata$Sample)

permanova_bray <- adonis2(
  bray_dist ~ Technology,
  data = metadata,
  permutations = param_perm
)

permanova_bray
# Permutation test for adonis under reduced model
# Blocks:  metadata$Sample
# Permutation: free
# Number of permutations: 999
#
# adonis2(formula = bray_dist ~ Technology, data = metadata, permutations = param_perm)
#           Df SumOfSqs      R2      F Pr(>F)
# Model      2    3.363 0.06617 7.4397  0.001 ***
# Residual 210   47.467 0.93383
# Total    212   50.831 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# To interpret Correctly, we need to do a homogeneity dispersion test

disp_test_bray <- betadisper(bray_dist, metadata$Technology)
permutest(disp_test_bray, permutations = param_perm)
#
# Permutation test for homogeneity of multivariate dispersions
# Blocks:  metadata$Sample
# Permutation: free
# Number of permutations: 999
#
# Response: Distances
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.59608 0.298041 57.614    999  0.001 ***
# Residuals 210 1.08634 0.005173
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# So we can say that permANOVA Technology explain 6.6% of variance in fungal community composition (R2= 0.066, P=0.001) but this result as to be taken cautiously because the result of the betadisper test can be interpreted as dispersion is significantly heterogenous (P=0.001)  So permANOVA migth just detect heterogeneity.

# ============ Box plot beta disper =========== #

png(filename = "~/Dext2/StageCRBE/results/fig/boxplot_disp_aitch_techno.png")
svg(filename = "~/Dext2/StageCRBE/results/fig/boxplot_disp_aitch_techno.svg")
boxplot(
  disp_test_aitch,
  main = "Dispersion par technologie - Aitchison",
  ylab = "Distance d'Aitchison au centroïde",
  xlab = "",
  col = c("darkgoldenrod1", "steelblue1", "orchid1"),
  notch = FALSE
)
dev.off()

png(filename = "~/Dext2/StageCRBE/results/fig/boxplot_disp_bray_techno.png")
svg(filename = "~/Dext2/StageCRBE/results/fig/boxplot_disp_bray_techno.svg")
boxplot(
  disp_test_bray,
  main = "Dispersion par technologie - Bray-Curtis",
  ylab = "Distance de Bray-Curtis au centroïde",
  xlab = "",
  col = c("darkgoldenrod1", "steelblue1", "orchid1"),
  notch = FALSE
)
dev.off()

disp_graph_aitch <- data.frame(
  Distance = disp_test_aitch$distances,
  Technology = metadata$Technology
)

## ====== indice de Shanon =========== #

pacbio_otu <- get_otu(
  "PacBio/OTU/OTU_table_LULU.txt",
  "pacbio"
)[[1]]
nano_otu <- get_otu(
  "Nanopore/OTU/OTU_table_mumu.tsv",
  "nanopore"
)[[1]]
ill_otu <- get_otu(
  "Illumina/OTU/OTU_table_ITS1_OTU97_VSEARCH_LULU.txt",
  "illumina"
)[[1]]

# For Vegan, we need the samples to be the row and the OTUs to be columns

pacbio_otu_matrix <- pacbio_otu |>
  remove_rownames() |>
  column_to_rownames("OTU") |>
  dplyr::select(-c("total_abundances")) |>
  as.matrix() |>
  t()

metada_otu_pacbio <- data.frame(
  Sample = colnames(
    pacbio_otu |>
      dplyr::select(-c("OTU", "total_abundances"))
  )
) |>
  mutate(
    Site = case_when(
      str_detect(Sample, "^C") ~ "Bois du Chapitre",
      str_detect(Sample, "^M") ~ "Forêt de la Massane",
      str_detect(Sample, "^V") ~ "Grand Ventron"
    ),
    SampleType = if_else(str_detect(Sample, "\\."), "Single", "Pool"),
    Site_Type = paste(Site, SampleType, sep = " - ")
  )

shannon <- diversity(pacbio_otu_matrix, index = "shannon")

df_shannon <- data.frame(
  Sample = names(shannon),
  Shannon = shannon
)

df_shannon <- merge(df_shannon, metada_otu_pacbio, by = "Sample")

ggplot(df_shannon, aes(x = Site, y = Shannon, fill = Site)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  labs(
    title = "Indice de Shannon par site d'échantillonnage - PacBio",
    x = "",
    y = "Indice de Shannon"
  )
ggsave(
  "~/Dext2/StageCRBE/results/fig/boxplot_shannon_pacbio_sites.png",
  width = 20,
  height = 10
)
ggplot(df_shannon, aes(x = Site_Type, y = Shannon, fill = Site_Type)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  labs(
    title = "Indice de Shannon par site et type d'échantillonnage - PacBio",
    x = "",
    y = "Indice de Shannon"
  )
ggsave(
  "~/Dext2/StageCRBE/results/fig/boxplot_shannon_pacbio_site_type.png",
  width = 20,
  height = 10
)
