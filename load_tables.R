#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)
library(readxl)
library(here)

#-------------------------------------------------------------------------------
# Define general function
make_biotype_levels <- function(df) {
  # Transforms the orf_biotype column in df to a factor
  df %>%
    mutate(
      orf_biotype = str_replace(
        orf_biotype, "Processed transcript ORFs",
        "Processed\ntranscript ORFs"
      ),
      orf_biotype = factor(orf_biotype,
                           levels = c(
                             "uORFs", "uoORFs", "intORFs", "doORFs",
                             "dORFs", "lncRNA ORFs",
                             "Processed\ntranscript ORFs"
                           )
      ),
    )
}

#-------------------------------------------------------------------------------
# Define required file locations
here::i_am("load_peptideatlas_data.R")

raw_dir <- here("raw")

# List of 7,264 ncORFs from doi: 10.1038/s41587-022-01369-0
orf_file <- file.path(raw_dir, "ncorf_list.xlsx")

# Supplementary tables about detected peptides and ncORFs
table_s2_file <- file.path(raw_dir, "060924_Supp_Table_S2.xlsx")
table_s3_file <- file.path(raw_dir, "060924_Supp_Table_S3.xlsx")
table_s4_file <- file.path(raw_dir, "060924_Supp_Table_S4.xlsx")
table_s5_file <- file.path(raw_dir, "060924_Supp_Table_S5.xlsx")

# Suplementary table S6 containing MS-runs and their annotation
annotation_file <- file.path(raw_dir, "060924_Supp_Table_S6.xlsx")

#-------------------------------------------------------------------------------
#Load tables
orf_sequences <- read_excel(orf_file, sheet = 3)[
  , c("orf_name", "orf_biotype", "orf_sequence", "gene_name", "gene_id",
      "transcript", "strand", "chrm", "starts", "ends")
] %>%
  bind_rows(read_excel(orf_file, sheet = 4)[
    , c("orf_name", "orf_biotype", "orf_sequence", "gene_name", "gene_id",
        "transcript", "strand", "chrm", "starts", "ends")
  ]) %>%
  mutate(
    orf_sequence = str_replace(orf_sequence, "[*]", ""),
    orf_biotype = ifelse(orf_biotype == "processed_transcript",
                         "Processed transcript ORF",
                         orf_biotype
    ),
    orf_biotype = ifelse(orf_biotype == "lncRNA", "lncRNA ORF", orf_biotype),
    orf_biotype = str_replace(orf_biotype, "ORF", "ORFs"),
    protein_accession = paste0("CONTRIB_GENCODE_", orf_name)
  ) %>%
  dplyr::rename(protein = orf_sequence)

annotations <- read_excel(annotation_file) %>%
  rename(dataset_id = dataset, ms_name = ms_run_name, hla_class = HLA_class) %>%
  mutate(
    is_cancer = case_when(
      grepl("TUMORAL", disease_state) ~ "Cancer",
      TRUE ~ "Non-cancer"
    ),
    is_cellline = ifelse(str_detect(
      sample_category,
      regex("Cell Line", ignore_case = T)
    ),
    "cell line", "non-cell line"
    ),
    sample_type = paste0(is_cancer, "_", is_cellline)
  )

table_s2 <- read_excel(table_s2_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas identifier",
    orf_name = "Ribo-Seq_ORF"
  )

table_s3 <- read_excel(table_s3_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF",
    initial_tier =
      "Initial.tier.for.level.of.evidence.prior.to.manual.inspection",
    final_tier = "Final.tier.for.level.of.evidence"
  )

table_s4 <- read_excel(table_s4_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF"
  )

table_s5 <- read_excel(table_s5_file) %>%
  dplyr::rename(
    identifier = "PeptideAtlas.identifier",
    orf_name = "Ribo-Seq_ORF",
    initial_tier =
      "Initial.tier.for.level.of.evidence.prior.to.manual.inspection",
    final_tier = "Final.tier.for.level.of.evidence"
  ) %>%
  rename_all(~ make.names(.))
