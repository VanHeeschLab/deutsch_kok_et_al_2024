# Author: Leron Kok
# This scripts loads the supplementary tables needed for the analyses

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)
library(readxl)
library(here)

#-------------------------------------------------------------------------------
# Define general functions
make_biotype_levels <- function(df) {
  # Transforms the orf_biotype column in df to a factor
  #
  # Args:
  #   df: dataframe which contains at least the column orf_biotype
  #
  # Returns:
  #   The input dataframe with orf_biotype as a factor
  df %>%
    mutate(
      orf_biotype = str_replace(
        orf_biotype, "Processed transcript ORFs",
        "Processed\ntranscript ORFs"
      ),
      orf_biotype =
        factor(orf_biotype, levels = c(
          "uORFs", "uoORFs", "intORFs", "doORFs",
          "dORFs", "lncRNA ORFs",
          "Processed\ntranscript ORFs"
        )),
    )
}

#-------------------------------------------------------------------------------
# Define required file locations
raw_dir <- here("raw")

# List of 7,264 ncORFs from doi: 10.1038/s41587-022-01369-0
orf_file <- file.path(raw_dir, "ncorf_list.xlsx")

# Supplementary tables about detected peptides and ncORFs
table_s2_file <- file.path(raw_dir, "20250801_VanHeesch_EDtable2.xlsx")
table_s3_file <- file.path(raw_dir, "20250801_VanHeesch_EDtable3.xlsx")
table_s6_file <- file.path(raw_dir, "20250801_VanHeesch_EDtable6.xlsx")
table_s7_file <- file.path(raw_dir, "20250801_VanHeesch_EDtable7.xlsx")
table_s8_file <- file.path(raw_dir, "20250801_VanHeesch_EDtable8.xlsx")

#-------------------------------------------------------------------------------
# Load tables

# Data of 7,624 ncORFs
orf_sequences <- read_excel(orf_file, sheet = 3)[
  , c(
    "orf_name", "orf_biotype", "orf_sequence", "gene_name", "gene_id",
    "transcript", "strand", "chrm", "starts", "ends", "PhyloCSF (120 mammals)"
  )
] %>%
  bind_rows(read_excel(orf_file, sheet = 4)[
    , c(
      "orf_name", "orf_biotype", "orf_sequence", "gene_name", "gene_id",
      "transcript", "strand", "chrm", "starts", "ends", "PhyloCSF (120 mammals)"
    )
  ]) %>%
  mutate(
    orf_sequence = str_replace(orf_sequence, "[*]", ""),
    orf_biotype = ifelse(orf_biotype == "processed_transcript",
      "Processed transcript ORF",
      orf_biotype
    ),
    orf_biotype = ifelse(orf_biotype == "lncRNA", "lncRNA ORF", orf_biotype),
    orf_biotype = str_replace(orf_biotype, "ORF", "ORFs"),
    protein_accession = paste0("CONTRIB_GENCODE_", orf_name),
  ) %>%
  dplyr::rename(protein = orf_sequence, phylocsf = `PhyloCSF (120 mammals)`) %>%
  mutate(phylocsf = as.numeric(phylocsf))

# Annotation data of the MS-runs
annotations <- read_excel(table_s8_file) %>%
  rename(dataset_id = dataset, ms_name = ms_run_name, hla_class = HLA_class) %>%
  mutate(
    is_cancer = case_when(
      grepl("TUMORAL", disease_state) ~ "Cancer",
      TRUE ~ "Non-cancer"
    ),
    is_cellline = ifelse(str_detect(
      sample_category,
      regex("Cell Line", ignore_case = TRUE)
    ),
    "cell line", "non-cell line"
    ),
    sample_type = paste0(is_cancer, "_", is_cellline)
  )

# Supplementary tables S2-5
nonhla_peptide_table <- read_excel(table_s2_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas identifier",
    orf_name = "Ribo-Seq_ORF"
  )

nonhla_orf_table <- read_excel(table_s3_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF",
    initial_tier =
      "Initial.tier.for.level.of.evidence.prior.to.manual.inspection",
    final_tier = "Final.tier.for.level.of.evidence"
  )

hla_peptide_table <- read_excel(table_s6_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF"
  )

hla_orf_table <- read_excel(table_s7_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF",
    initial_tier =
      "Initial.tier.for.level.of.evidence.prior.to.manual.inspection",
    final_tier = "Final.tier.for.level.of.evidence"
  ) %>%
  rename_all(~ make.names(.))
