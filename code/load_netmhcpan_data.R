# Author: Leron Kok
# This scripts loads the NetMCHpan data (or prepares the input data)

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)
library(data.table)
library(here)

#-------------------------------------------------------------------------------
# Define functions
write_net <- function(sample_num, net_data) {
  # Write peptides to net_input folder for netMHCpan predictions
  sample_name <- paste0("/netmhcpan_input", sample_num, ".txt")
  net_data <- net_data %>%
    dplyr::select(peptide) %>%
    distinct()
  write.table(net_data, paste0(net_input, sample_name),
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

read_netmhcpan <- function(net_path) {
  # Read netMHCpan predictions to dataframe
  net_data <- read.delim(net_path) %>% distinct()
  rank_data <- net_data[
    2:nrow(net_data),
    c(2, which(grepl("EL_Rank", net_data[1, ])))
  ]
  colnames(rank_data) <- c("peptide", "rank")
  rank_data$hla_type <- colnames(net_data)[grepl("HLA",
                                                 colnames(net_data))][1] %>%
    str_replace_all("HLA|[.]", "")
  rank_data <- rank_data %>% mutate(rank = as.numeric(rank))
  return(rank_data)
}

#-------------------------------------------------------------------------------
# Define required file locations
here::i_am("load_netmhcpan_data.R")
processed_dir <- here("processed")

# For running predictions yourself specify
## Directory used for NetMHCpan binding predictions input
net_input <- file.path(processed_dir, "netmhcpan_input")
## Directory containing binding predictions for detected peptides
net_output <- file.path(processed_dir, "netmhcpan_output_annotatedalleles")
# Or use data previously generated (from github)
net_output_rds <- file.path(raw_dir, "netmhcpan.RDS")

#-------------------------------------------------------------------------------
# Load initial files
## Check if pre-generated NetMCHpan data exists
if (file.exists(net_output_rds)) {
  netmhcpan <- readRDS(net_output_rds)
  ## Else try to locate NetMCHpan prediction output
} else {
  net_files <- list.files(
    path = net_output,
    pattern = paste0("net_.*.xls"), full.names = TRUE
  )
  if (length(net_files) > 0) {
    netmhcpan <- lapply(net_files, read_netmhcpan) %>% bind_rows()
  }
}

## If netMHCpan output data does not exist, create input for NetMHCpan
##  predictions and create warning that NetMHCpan should run first
if (!exists("netmhcpan")) {
  netmhcpan_input <- pep_net_orf %>%
    distinct(peptide, hla_type) %>%
    rename("allele" = 2) %>%
    mutate(
      allele = paste0(
        "HLA-", str_extract(allele, "^[A-Z]\\d\\d"), ":",
        str_extract(allele, "(?<=^[A-Z]\\d\\d)\\d+")
      )
    ) %>%
    select(peptide, allele)

  netmhcpan_input_s <- netmhcpan_input %>% split(netmhcpan_input$allele)
  write_net_res <- sapply(
    1:length(netmhcpan_input_s),
    function(x) write_net(x, netmhcpan_input_s[[x]])
  )
  write.table(names(netmhcpan_input_s), paste0(net_input, "/net_alleles.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  # Raise warning
  stop(paste0("No NetMHCpan predictions present.\n",
              "Run NetMHCpan first on generated input data"))
}

# Create dataframe with minimum rank found per sample-peptide combination
setDT(hla_df)
setDT(netmhcpan)

pep_net_orf <- netmhcpan[hla_df, on = c("peptide", "hla_type"),
                         nomatch = NULL] %>%
  group_by(
    peptide, dataset_id, ms_name, pep_len, sample_type, is_cancer,
    is_cellline, pep_type
  ) %>%
  mutate(rank = as.numeric(rank)) %>%
  summarise(
    min_rank = min(rank), hla_type = hla_type[which.min(rank)],
    hla_typing = first(hla_typing),
    is_binding = ifelse(min_rank <= 2, TRUE, FALSE),
    .groups = "drop"
  ) %>%
  left_join(pep_nonc %>% distinct(peptide, orf_biotype), "peptide",
            relationship = "many-to-many"
  ) %>%
  mutate(is_orf = ifelse(is.na(orf_biotype), FALSE, TRUE))
