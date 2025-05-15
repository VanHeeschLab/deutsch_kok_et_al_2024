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

raw_dir <- here("raw")
processed_dir <- here("processed")
plot_dir <- here("plots")

# List of 7,264 ncORFs from doi: 10.1038/s41587-022-01369-0
orf_file <- file.path(raw_dir, "ncorf_list.xlsx")

# Peptide files and file mapping peptides to proteins/ncORFs from PeptideAtlas
peptide_file <- file.path(raw_dir, "peptide_sample_msrun_counts.tsv")
## Either needs mapping file
mapping_file <- file.path(raw_dir, "peptide_mapping.tsv")
## Or filtered peptides file
filtered_peptides_file <- file.path(raw_dir, "filtered_peptides.Rdata")

# Protein sequences
## Either needs file containing the sequences of proteins/ncORFs from
##  PeptideAtlas
seq_file <- file.path(raw_dir, "Homo_sapiens.fasta")
## Or prot_seq_df file
can_nonc_seq_file <- file.path(raw_dir, "can_nonc_seq.RDS")

# List of canonical CDSs from PeptideAtlas
canonical_file <- file.path(raw_dir, "Core20k.txt")

# Supplementary tables about detected peptides and ncORFs
table_s2_file <- file.path(raw_dir, "060924_Supp_Table_S2.xlsx")
table_s3_file <- file.path(raw_dir, "060924_Supp_Table_S3.xlsx")
table_s4_file <- file.path(raw_dir, "060924_Supp_Table_S4.xlsx")
table_s5_file <- file.path(raw_dir, "060924_Supp_Table_S5.xlsx")

# Suplementary table S6 containing MS-runs and their annotation
annotation_file <- file.path(raw_dir, "060924_Supp_Table_S6.xlsx")

# GTF file
gtf_file <- file.path("Homo_sapiens.GRCh38.102.gtf")

# GTEX FPKM mean expression data (excluding testis)
gtex_fpkm_file <- file.path(raw_dir, "GTEX_FPKMmean_expression.txt")

#-------------------------------------------------------------------------------
#Load initial files
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

gtex_data <- read.csv(gtex_fpkm_file)

can_id <- read_tsv(canonical_file, col_names = "protein_accession")

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

peptide <- fread(peptide_file)
colnames(peptide) <- c(
  "peptide", "sam_id", "count", "dataset_id",
  "sample_tag", "ms_name"
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

#-------------------------------------------------------------------------------
# Load protein sequences either from pre-generated file (from github) if
# available, or from PeptideAtlas fasta file

if (file.exists(can_nonc_seq_file)) {
  can_nonc_seq <- readRDS(can_nonc_seq_file)
} else {
  prot_seq_df <- read.fasta(seq_file, seqtype = "AA", as.string = T) %>%
    enframe() %>%
    dplyr::rename(protein_accession = name, protein = value) %>%
    mutate(protein = as.character(protein))
  
  can_nonc_seq <- prot_seq_df %>%
    semi_join(can_id, "protein_accession") %>%
    bind_rows(prot_seq_df %>%
                filter(str_starts(protein_accession, "CONTRIB_GENCODE_c")))
}

#-------------------------------------------------------------------------------
### Filter canonical and non-canonical peptides either using earlier generated 
# data (from Github) or using the mapping file from PeptideAtlas

if (file.exists(filtered_peptides_file)) {
  load(filtered_peptides_file)
} else {
  mapping <- fread(mapping_file)
  colnames(mapping) <- c("peptide_accession", "peptide", "protein_accession",
                         "start_loc", "end_loc", "prev_aa", "next_aa")
  
  can_peptides <- mapping %>%
    inner_join(can_id, "protein_accession") %>%
    distinct(peptide)
  
  mapping <- mapping %>% semi_join(peptide, "peptide")
  setDT(mapping)
  prefix_filter <- paste0("DECOY|CONTRIB_smORFs_Cui|CONTRIB_sORFs|",
                          "CONTRIB_Fedor|CONTRIB_Bazz|CONTRIB_HLA")
  mapping <- mapping[!str_starts(protein_accession, prefix_filter)]
  
  # Filter noncanonical peptides
  pep_nonc_long <- mapping[str_starts(protein_accession, "CONTRIB_GENCODE")] %>%
    distinct(peptide) %>%
    inner_join(mapping, "peptide") %>%
    filter(!str_starts(protein_accession, 
                       paste0(prefix_filter, 
                              "|CONTRIB_GENCODE_nearcognate"))) %>%
    anti_join(can_peptides, "peptide") %>%
    group_by(peptide) %>%
    mutate(n_pep = n()) %>%
    filter(n_pep <= 10) %>%
    inner_join(orf_sequences, "protein_accession") %>%
    mutate(pep_len = nchar(peptide)) %>%
    filter(pep_len >= 8) %>%
    ungroup()
  
  # Filter canonical peptides
  pep_can <- mapping %>%
    semi_join(can_peptides, "peptide") %>%
    group_by(peptide) %>%
    mutate(n_pep = n(), pep_len = nchar(peptide)) %>%
    filter(n_pep <= 30, pep_len >= 8)
}

#-------------------------------------------------------------------------------
# Prepare data frames for analysis
pep_nonc <- pep_nonc_long %>%
  inner_join(
    table_s4 %>% 
      select(-c(orf_name, transcript, starts, ends, orf_biotype, gene_name,
                gene_id, strand, chrm)) %>% 
      dplyr::rename(other_mappings = "Other.mappings"),
    c("peptide" = "sequence", "protein_accession" = "identifier")
  ) %>%
  ungroup()

seq_stats <- can_nonc_seq %>%
  filter(!str_detect(protein, "U")) %>%
  mutate(
    pro_cat = ifelse(str_starts(protein_accession, "CONTRIB_GENCODE_"),
                     "Non-canonical", "Canonical"
    ),
    len = nchar(protein)
  ) %>%
  filter(len >= 16) %>%
  left_join(
    bind_rows(
      pep_can %>%
        ungroup() %>%
        distinct(protein_accession) %>%
        semi_join(can_id, "protein_accession"),
      pep_nonc_long %>% distinct(protein_accession)
    ) %>%
      mutate(detected = "Detected"), "protein_accession"
  ) %>%
  left_join(
    orf_sequences %>% distinct(protein_accession, orf_biotype),
    "protein_accession"
  ) %>%
  replace_na(list(detected = "Undetected", orf_biotype = "CDS"))


# Combine peptides and annotation, filter for HLA-I
peptide_hla <- peptide %>%
  inner_join(annotations, c("ms_name", "dataset_id"))

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

setDT(pep_can_nonc)
setDT(peptide_hla)
peptide_hla <- peptide_hla[pep_can_nonc, on = "peptide"]

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
    hla_outlier = ifelse(str_detect(hla_type, "A2401|B4301|C1201"), T, F),
    has_hla_outlier = max(hla_outlier)
  ) %>%
  ungroup() %>%
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

pep_can_nonc_hladf <- pep_can_nonc %>% semi_join(
  hla_df %>% distinct(peptide),
  "peptide"
)
setDT(hla_df)
hla_df <- hla_df[pep_can_nonc_hladf, on = "peptide"]
