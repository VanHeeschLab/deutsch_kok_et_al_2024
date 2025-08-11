# Author: Leron Kok
# This scripts loads the PeptideAtlas data needed for the analyses

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)
library(readxl)
library(data.table)
library(seqinr)
library(here)

#-------------------------------------------------------------------------------
# Define required file locations
here::i_am("load_peptideatlas_data.R")

raw_dir <- here("raw")  # Directory where raw data is stored

# File paths for peptide and mapping data
peptide_file <- file.path(raw_dir, "peptide_sample_msrun_counts.tsv")
## Either needs mapping file
mapping_file <- file.path(raw_dir, "peptide_mapping.tsv")
## Or filtered peptides file
filtered_peptides_file <- file.path(raw_dir, "filtered_peptides.Rdata")

# Protein sequence data
## Either needs file containing the sequences of proteins/ncORFs from
##  PeptideAtlas
seq_file <- file.path(raw_dir, "Homo_sapiens.fasta")
## Or can_nonc_seq file
can_nonc_seq_file <- file.path(raw_dir, "can_nonc_seq.RDS")

# Canonical protein list from PeptideAtlas
canonical_file <- file.path(raw_dir, "Core20k.txt")

#-------------------------------------------------------------------------------
# Load initial files

# Load canonical protein IDs
can_id <- read_tsv(canonical_file, col_names = "protein_accession")

# Load peptide data and rename columns for clarity
peptide <- fread(peptide_file)
colnames(peptide) <- c(
  "peptide", "sam_id", "count", "dataset_id",
  "sample_tag", "ms_name"
)

# Load protein sequences from RDS if available, otherwise parse FASTA
if (file.exists(can_nonc_seq_file)) {
  # Load pre-processed sequence data
  can_nonc_seq <- readRDS(can_nonc_seq_file)
} else {
  # Read FASTA file and process into a tidy data frame
  prot_seq_df <- read.fasta(seq_file, seqtype = "AA", as.string = TRUE) %>%
    enframe() %>%
    dplyr::rename(protein_accession = name, protein = value) %>%
    mutate(protein = as.character(protein))

  # Keep only canonical and specific noncanonical entries
  can_nonc_seq <- prot_seq_df %>%
    semi_join(can_id, "protein_accession") %>%
    bind_rows(prot_seq_df %>%
                filter(str_starts(protein_accession, "CONTRIB_GENCODE_c")))
}

# Load or generate filtered peptide sets (canonical and noncanonical)
if (file.exists(filtered_peptides_file)) {
  # Load pre-filtered peptide sets (likely from GitHub or prior run)
  load(filtered_peptides_file)
} else {
  # Load raw mapping file
  mapping <- fread(mapping_file)
  colnames(mapping) <- c("peptide_accession", "peptide", "protein_accession",
                         "start_loc", "end_loc", "prev_aa", "next_aa")

  # Extract peptides that map to canonical proteins
  can_peptides <- mapping %>%
    inner_join(can_id, "protein_accession") %>%
    distinct(peptide)

  # Keep only peptides that exist in the 'peptide' file
  mapping <- mapping %>% semi_join(peptide, "peptide")
  setDT(mapping)  # Convert to data.table for efficiency

  # Define unwanted protein accession prefixes
  prefix_filter <- paste0("DECOY|CONTRIB_smORFs_Cui|CONTRIB_sORFs|",
                          "CONTRIB_Fedor|CONTRIB_Bazz|CONTRIB_HLA")
  mapping <- mapping[!str_starts(protein_accession, prefix_filter)]

  # Process noncanonical peptides
  pep_nonc_long <- mapping[str_starts(protein_accession, "CONTRIB_GENCODE")] %>%
    distinct(peptide) %>%
    inner_join(mapping, "peptide") %>%
    filter(!str_starts(protein_accession,
                       paste0(prefix_filter,
                              "|CONTRIB_GENCODE_nearcognate"))) %>%
    anti_join(can_peptides, "peptide") %>%  # Exclude canonical matches
    group_by(peptide) %>%
    mutate(n_pep = n()) %>%  # Count how many mappings each peptide has
    filter(n_pep <= 10) %>%  # Exclude promiscuous peptides
    inner_join(orf_sequences, "protein_accession") %>%  # Match to ORF sequences
    mutate(pep_len = nchar(peptide)) %>%
    filter(pep_len >= 8) %>%
    ungroup()

  # Process canonical peptides
  pep_can <- mapping %>%
    semi_join(can_peptides, "peptide") %>%
    group_by(peptide) %>%
    mutate(n_pep = n(), pep_len = nchar(peptide)) %>%
    filter(n_pep <= 30, pep_len >= 8)
}

# Prepare condensed noncanonical peptide annotation (one ncORF pep peptide)
pep_nonc <- pep_nonc_long %>%
  inner_join(
    hla_peptide_table %>%
      select(-c(orf_name, transcript, starts, ends, orf_biotype, gene_name,
                gene_id, strand, chrm)) %>%
      dplyr::rename(other_mappings = "Other.mappings"),
    c("peptide" = "sequence", "protein_accession" = "identifier")
  ) %>%
  ungroup()

#-------------------------------------------------------------------------------
# Annotate peptides with HLA dataset info

# Match all peptides with the MS run annotation
peptide_hla <- peptide %>%
  inner_join(annotations, c("ms_name", "dataset_id"))

# Combine canonical and noncanonical peptides with types
pep_can_nonc <- peptide %>%
  distinct(peptide) %>%
  left_join(
    bind_rows(
      pep_can %>%
        distinct(peptide) %>%
        mutate(pep_type = "canonical"),
      pep_nonc %>%
        distinct(peptide) %>%
        mutate(pep_type = "noncanonical")
    ),
    "peptide"
  ) %>%
  mutate(pep_type = ifelse(is.na(pep_type), "other", pep_type))

# Convert to data.table for fast merging
setDT(pep_can_nonc)
setDT(peptide_hla)
# Filter peptide_hla to only include filtered peptides
peptide_hla <- peptide_hla[pep_can_nonc, on = "peptide"]

# Filter annotations for HLA class I peptides and clean up types
hla_df <- annotations %>%
  select(
    dataset_id, ms_name, sample_type, is_cancer, is_cellline, hla_class,
    hla_typing
  ) %>%
  filter(hla_class == "I", !is.na(hla_typing)) %>%
  separate_rows(hla_typing, sep = ",") %>%
  dplyr::rename(hla_type = hla_typing) %>%
  distinct() %>%
  group_by(dataset_id, ms_name) %>%
  mutate(
    n_hla = n(), hla_len = nchar(hla_type), min_len = min(hla_len),
    # Exclude HLA outliers not recognized ny netMHCpan
    hla_outlier = ifelse(str_detect(hla_type, "A2401|B4301|C1201"),
                         TRUE, FALSE),
    has_hla_outlier = max(hla_outlier)
  ) %>%
  ungroup() %>%
  # Remove datasets with outlier HLA types or incomplete (no 4 digit) HLA types
  filter(min_len >= 5, has_hla_outlier == 0) %>%
  select(-min_len, -hla_len, has_hla_outlier, hla_outlier) %>%
  group_by(dataset_id, ms_name) %>%
  mutate(hla_typing = paste0(hla_type, collapse = ",")) %>%
  inner_join(
    peptide %>%
      mutate(pep_len = nchar(peptide)) %>%
      filter(pep_len >= 8, pep_len <= 12) %>%
      select(peptide, dataset_id, ms_name, pep_len),
    c("ms_name", "dataset_id"),
    relationship = "many-to-many"
  )

# Filter pep_can_nonc for those peptides with at least one known HLA typing
pep_can_nonc_hladf <- pep_can_nonc %>% semi_join(
  hla_df %>% distinct(peptide),
  "peptide"
)
setDT(hla_df)
# Filter hla_df for the filtered peptides
hla_df <- hla_df[pep_can_nonc_hladf, on = "peptide"]
