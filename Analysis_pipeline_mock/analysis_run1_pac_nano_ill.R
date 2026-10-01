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
library(scales)
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


get_shanono_rarecurve <- function(otu_table) {
  otu_matrix <- otu_table |>
    remove_rownames() |>
    column_to_rownames("OTU") |>
    select(-c("total_abundances"))
  otu_rare <- iNEXT(
    otu_matrix,
    q = 1,
    datatype = "abundance",
    se = FALSE,
    knots = 15
  )
  return(otu_rare)
}


pacbio_shan_rare <- get_shanono_rarecurve(pacbio_otu)
nano_shan_rare <- get_shanono_rarecurve(nano_otu)
ill_shan_rare <- get_shanono_rarecurve(ill_otu)

ggiNEXT(
  pacbio_shan_rare,
  type = 1,
  se = FALSE,
  facet.var = "None",
  color.var = "Assemblage",
  grey = FALSE
)
png(
  'figs/pacbio_shannon_rare.png',
  height = 1080,
  width = 1920,
  res = 300
)

plot(pacbio_shan_rare)
dev.off()
png(
  'figs/nano_shannon_rare.png',
  height = 1080,
  width = 1920,
  res = 300
)

plot(nano_shan_rare)
dev.off()


png(
  'figs/ill_shannon_rare.png',
  height = 1080,
  width = 1920,
  res = 300
)

plot(ill_shan_rare)
dev.off()


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
    x = "Nombre de reads",
    y = "Nombre d'OTU observés",
    color = "Techno - type"
  ) +
  scale_x_continuous(
    limits = c(0, 1e5),
    labels = label_number(scale = 1e-3, suffix = "k") # 20 000 -> 20k
  ) +
  scale_y_continuous(limits = c(0, 5000)) +
  theme_minimal(base_size = 16) + # augmente la taille de base de tout le texte
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    strip.text = element_text(size = 16, face = "bold"), # titres des facettes
  )


ggsave(
  "~/Dext2/StageCRBE/results/fig/rarefaction_lim_2.png",
  width = 24,
  height = 12,
  units = "cm",
  dpi = 600
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

# Test of matching taxonomies

matching_unite <- merge(
  pacbio_otu_sintax_unite |>
    select(c("OTU", "Genus")) |>
    mutate(Genus_sintax = Genus) |>
    select(-c("Genus")),
  pacbio_otu_vsearch_unite |>
    select(c("OTU", "Genus")) |>
    mutate(Genus_vsearch = Genus) |>
    select(-c("Genus")),
  by = "OTU"
) |>
  mutate(
    match = (ifelse(
      !is.na(Genus_sintax) &
        !is.na(Genus_vsearch) &
        Genus_sintax == Genus_vsearch,
      1,
      0
    )),
    not_match = (ifelse(
      !is.na(Genus_sintax) &
        !is.na(Genus_vsearch) &
        Genus_sintax != Genus_vsearch,
      1,
      0
    )),
    na = (ifelse(
      is.na(Genus_sintax) |
        is.na(Genus_vsearch),
      1,
      0
    )),
    xor_na = (ifelse(
      (is.na(Genus_sintax) |
        is.na(Genus_vsearch)) &
        !(is.na(Genus_sintax) &
          is.na(Genus_vsearch)),
      1,
      0
    ))
  )

percent_matching <- sum(matching_unite$match) /
  nrow(matching_unite) *
  100

percent_matching
percent_not_matching <- sum(matching_unite$not_match) /
  nrow(matching_unite) *
  100
percent_not_matching
percent_na <- sum(matching_unite$na) /
  nrow(matching_unite) *
  100
percent_na
percent_xor_na <- sum(matching_unite$xor_na) /
  nrow(matching_unite) *
  100
percent_xor_na

matching_euk <- merge(
  pacbio_otu_sintax_euk |>
    select(c("OTU", "Genus")) |>
    mutate(Genus_sintax = Genus) |>
    select(-c("Genus")),
  pacbio_otu_vsearch_euk |>
    select(c("OTU", "Genus")) |>
    mutate(Genus_vsearch = Genus) |>
    select(-c("Genus")),
  by = "OTU"
) |>
  mutate(
    match = (ifelse(
      !is.na(Genus_sintax) &
        !is.na(Genus_vsearch) &
        Genus_sintax == Genus_vsearch,
      1,
      0
    )),
    not_match = (ifelse(
      !is.na(Genus_sintax) &
        !is.na(Genus_vsearch) &
        Genus_sintax != Genus_vsearch,
      1,
      0
    )),
    na = (ifelse(
      is.na(Genus_sintax) |
        is.na(Genus_vsearch),
      1,
      0
    )),
    xor_na = (ifelse(
      (is.na(Genus_sintax) |
        is.na(Genus_vsearch)) &
        !(is.na(Genus_sintax) &
          is.na(Genus_vsearch)),
      1,
      0
    ))
  )

percent_matching <- sum(matching_unite$match) /
  nrow(matching_unite) *
  100

percent_matching
percent_not_matching <- sum(matching_unite$not_match) /
  nrow(matching_unite) *
  100
percent_not_matching
percent_na <- sum(matching_unite$na) /
  nrow(matching_unite) *
  100
percent_na
percent_xor_na <- sum(matching_unite$xor_na) /
  nrow(matching_unite) *
  100
percent_xor_na


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
    unique()
  # DnabarcoderUnite = pacbio_otu_dnabarcoder |>
  #   dplyr::filter(!is.na(Genus)) |>
  #   pull(Genus) |>
  #   unique()
)
grid.newpage()

venn_plot <- venn.diagram(
  x = X,
  category.names = c(
    "SINTAX - EUK",
    "SINTAX - Unite",
    "Vsearch - EUK",
    "Vsearch - Unite" #,
    # "Dnabarcoder Unite"
  ),
  filename = NULL,
  output = TRUE,
  imagetype = "png",

  lwd = 1,
  col = c("orange", 'orangered', 'orchid', 'mediumpurple'), #, 'palegreen'),
  fill = c(
    alpha("orange", 0.5),
    alpha("orangered", 0.5),
    alpha('orchid', 0.5),
    alpha("mediumpurple", 0.5) #,
    # alpha("palegreen", 0.5)
  ),

  cex = 0.5,
  fontfamily = "sans",
  ext.text = FALSE,

  cat.cex = 1,
  cat.dist = 0.22,
  # cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("orange", 'orangered', 'orchid', 'mediumpurple'), #, 'palegreen'),
  margin = 0.1
)


grid.draw(venn_plot)

png(
  'figs/venn_pacbio_nodna_genus.png',
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

  cat.cex = 1,
  cat.dist = 0.06,
  # cat.pos = c(-27, 27),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c('darkgoldenrod', 'steelblue', 'orchid'),
  margin = 0.05
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
  # stat_ellipse() +
  scale_color_manual(
    values = c(
      "Illumina" = "darkgoldenrod1",
      "PacBio" = "orchid1",
      "Nanopore" = "steelblue1"
    )
  ) +
  labs(
    title = "PCoA Aitchison",
    x = "Axe 1",
    y = "Axe 2",
    color = "Technologie",
    subtitle = paste0(
      "Axe 1 = ",
      var_expl_aitch[1],
      "%, Axe 2 = ",
      var_expl_aitch[2],
      "%"
    )
  ) +

  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm")
  )

ggsave(
  "~/Dext2/StageCRBE/results/fig/PCoA_aitchison_techno_sintax_euk_2.pdf",
  width = 30,
  height = 16,
  units = "cm",
  dpi = 600
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
  # stat_ellipse() +
  scale_color_manual(
    values = c(
      "Illumina" = "darkgoldenrod1",
      "PacBio" = "orchid1",
      "Nanopore" = "steelblue1"
    )
  ) +
  labs(
    title = "PCoA Bray-Curtis",
    x = "Axe 1",
    y = "Axe 2",
    color = "Technologie",
    subtitle = paste0(
      "Axe 1 = ",
      var_expl_bray[1],
      "%, Axe 2 = ",
      var_expl_bray[2],
      "%"
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm")
  )


ggsave(
  "~/Dext2/StageCRBE/results/fig/PCoA_bray_techno_sintax_euk.pdf",
  width = 30,
  height = 16,
  units = "cm",
  dpi = 600
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

## get otu and tax for phyloseq

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

# =============== functions for taxinfo ================= #

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


# test_nano <- nano_otu_sintax_euk |>
#   dplyr::filter(Kingdom == "Fungi")
# rm(test_nano)
# unique(test_nano$Kingdom)
# unique(nano_otu_sintax_euk$Kingdom)
# Here we prepare the tables with only Genus and abundances
# We also filter to get only fungi
get_pq <- function(otu_table) {
  otu_mat <- otu_table |>
    dplyr::filter(Kingdom == "Fungi") |>
    dplyr::select(c(
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
      "V12.2",
    )) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix() |>
    t() |>
    decostand(method = "total") |>
    t()
  ## here OTU names are the row so we have to set taxa_are_rows to TRUE
  OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

  ## Get the Taxa matrix
  tax_matrix <- otu_table |>
    select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
    tibble::column_to_rownames("OTU") |>
    as.matrix()

  TAX <- tax_table(tax_matrix)

  PHYLOSEQ <- phyloseq(OTU, TAX)
  return(PHYLOSEQ)
}

nano_pq_sintax_euk <- get_pq(nano_otu_sintax_euk)
pacbio_pq_sintax_euk <- get_pq(pacbio_otu_sintax_euk)
ill_pq_sintax_euk <- get_pq(ill_otu_sintax_euk)


FUNGAL_TRAITS_TABLE <- "~/Dext2/StageCRBE/Database/traitsTable/FUNGALT_DB_MROY041125.csv"

nano_pq_sintax_euk_clean <- gna_verifier_pq(
  nano_pq_sintax_euk,
  data_sources = 210
)
# ✔ GNA verification summary:
# • Total taxa in phyloseq: 18876
# • Taxa submitted for verification: 13661
# • Genus-level only taxa: 11267
# • Total matches found: 1225
# • Synonyms: 93 (including 6 uninomial)
# • Accepted names: 1132 (including 671 uninomial)
# ℹ 671 uninomial accepted name(s) have `currentCanonicalSimple` set to
#   "NA" (`species_only` = TRUE)
pacbio_pq_sintax_euk_clean <- gna_verifier_pq(
  pacbio_pq_sintax_euk,
  data_sources = 210
)
# ✔ GNA verification summary:
# • Total taxa in phyloseq: 12604
# • Taxa submitted for verification: 8059
# • Genus-level only taxa: 5539
# • Total matches found: 1705
# • Synonyms: 121 (including 8 uninomial)
# • Accepted names: 1584 (including 856 uninomial)
# ℹ 856 uninomial accepted name(s) have `currentCanonicalSimple` set to "NA"
#   (`species_only` = TRUE)
ill_pq_sintax_euk_clean <- gna_verifier_pq(
  ill_pq_sintax_euk,
  data_sources = 210
)
# ✔ GNA verification summary:
# • Total taxa in phyloseq: 39698
# • Taxa submitted for verification: 26923
# • Genus-level only taxa: 20981
# • Total matches found: 2729
# • Synonyms: 222 (including 14 uninomial)
# • Accepted names: 2507 (including 1309 uninomial)
# ℹ 1309 uninomial accepted name(s) have `currentCanonicalSimple` set to "NA"
#   (`species_only` = TRUE)

nano_traits_sintax_euk <- fungal_traits_guilds(
  nano_pq_sintax_euk_clean,
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
pacbio_traits_sintax_euk <- fungal_traits_guilds(
  pacbio_pq_sintax_euk_clean,
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
ill_traits_sintax_euk <- fungal_traits_guilds(
  ill_pq_sintax_euk_clean,
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

nano_traits_sintax_euk <- add_trophicMode_ft(nano_traits_sintax_euk)

pacbio_traits_sintax_euk <- add_trophicMode_ft(pacbio_traits_sintax_euk)
ill_traits_sintax_euk <- add_trophicMode_ft(ill_traits_sintax_euk)
# Save the phyloseq object
write_pq(
  nano_traits_sintax_euk,
  path = "~/Dext2/StageCRBE/results/Nanopore/PQ/nano_sintax_euk_FungalTraits_enhanced_FunGuild"
)

write_pq(
  pacbio_traits_sintax_euk,
  path = "~/Dext2/StageCRBE/results/PacBio/PQ/pacbio_sintax_euk_FungalTraits_enhanced_FunGuild"
)

write_pq(
  ill_traits_sintax_euk,
  path = "~/Dext2/StageCRBE/results/Illumina/PQ/ill_sintax_euk_FungalTraits_enhanced_FunGuild"
)

write_pq(
  nano_traits_sintax_euk,
  path = "~/Dext2/StageCRBE/results/Nanopore/PQ/nano_sintax_euk_FungalTraits_enhanced_FunGuild"
)

pacbio_traits_sintax_euk <- read_pq_corr(
  path = "~/Dext2/StageCRBE/results/PacBio/PQ/pacbio_sintax_euk_FungalTraits_enhanced_FunGuild"
)

write_pq(
  ill_traits_sintax_euk,
  path = "~/Dext2/StageCRBE/results/Illumina/PQ/ill_sintax_euk_FungalTraits_enhanced_FunGuild"
)


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
      MeanRelAbundance = rowMeans(troph_mat, na.rm = TRUE),
      Group = group_pq
    )
    trophic_mean[[troph_col]] <- rownames(troph_mat)

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

plot_mean_all_trophic_abundance <- function(
  col_list,
  pq,
  title = "Mean trophic abundance accross samples by database"
) {
  # here impa_dfr allows to get the indexes of the phyloseq objects in the list, then map applies the function to all elements and _dfr allows to use rbind at the end to combine all dataframes at the end
  combined <- purrr::imap_dfr(col_list, function(troph_col, name_col) {
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
      MeanRelAbundance = rowMeans(troph_mat, na.rm = TRUE),
      TrophicMode = rownames(troph_mat),
      Group = name_col
    )

    return(trophic_mean)
  })

  cat("Set the color of unknown to grey\n")
  trophs <- unique(combined$TrophicMode)
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
    aes(x = Group, y = MeanRelAbundance, fill = TrophicMode)
  ) +
    # here it's important to put identity so the bars correspond to the relative abundance
    # in each samples and not the number of different abundance found in each sample
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = cols
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = str_wrap(title, width = 50),
      x = NULL,
      y = "Abondance relative",
      fill = ""
    ) +
    theme_idest(base_size = 16) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 13), # Taille du texte
      legend.title = element_text(size = 14, face = "bold"), # Titre légende
      legend.key.size = unit(1.2, "cm") # Taille des carrés de couleur
    )

  return(p)
}
col_list <- list(
  "FunGuild" = "fg_trophicMode",
  "FungalTraits" = "ft_trophicMode",
  "Consensus" = "cons_trophicMode"
)
plot_mean_all_trophic_abundance(
  col_list,
  pacbio_traits_sintax_euk,
  title = "Abondance relative moyenne des modes trophiques par base de données - PacBio - SINTAX - EUKARYOME"
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/troph_rel_col_pacbio_sintax_euk.png",
  width = 10,
  height = 8
)
ggsave(
  "~/Dext2/StageCRBE/results/fig/troph_rel_col_pacbio_sintax_euk.pdf",
  width = 10,
  height = 8
)
test <- as.data.frame(tax_table(nano_traits_sintax_euk))

plot_mean_trophic_abundance(
  pq_list,
  title = "Abondance relative moyenne des modes trophiques par technologies - FunGuild"
)

# ramp <- scales::colour_ramp(c("darkred","darkgreen","darkblue"))
# scales::show_col(c(ramp(seq(0,1, length = 7)),"grey70"))
ggsave(
  "~/Dext2/StageCRBE/results/fig/troph_rel_techno_fg_sintax_euk.png",
  width = 10,
  height = 8
)


plot_mean_trophic_abundance(
  pq_list,
  "ft_trophicMode",
  title = "Abondance relative moyenne des modes trophiques par technologies - FungalTrait"
)

# ramp <- scales::colour_ramp(c("darkred","darkgreen","darkblue"))
# scales::show_col(c(ramp(seq(0,1, length = 7)),"grey70"))
ggsave(
  "~/Dext2/StageCRBE/results/fig/troph_rel_techno_ft_sintax_euk.png",
  width = 10,
  height = 8
)

plot_mean_trophic_abundance(
  pq_list,
  "cons_trophicMode",
  title = "Abondance relative moyenne des modes trophiques par technologies - Consensus"
)

# ramp <- scales::colour_ramp(c("darkred","darkgreen","darkblue"))
# scales::show_col(c(ramp(seq(0,1, length = 7)),"grey70"))
ggsave(
  "~/Dext2/StageCRBE/results/fig/troph_rel_techno_cons_sintax_euk.png",
  width = 10,
  height = 8
)
