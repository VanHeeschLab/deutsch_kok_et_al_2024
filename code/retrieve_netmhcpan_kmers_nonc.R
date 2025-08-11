# Script to process binding predictions on ncORF 9mers
library(tidyverse)
library(here)

processed_dir <- here("processed")
# Load Rds file created by 'Figures_BindingPredictions.Rmd'
load(file.path(processed_dir, "retrieve_net_kmers_data.Rds"))

read_netmhcpan_nonc <- function(net_path) {
  net_data <- read.delim(net_path) %>% distinct()
  rank_data <- net_data[2:nrow(net_data),
                        c(2, which(grepl("EL_Rank", net_data[1,])))]
  colnames(rank_data) <- c("peptide", "rank")
  rank_data <- rank_data %>% mutate(rank = as.numeric(rank)) %>% 
    semi_join(nonc_kmers, c("peptide" = "kmer"))
  rank_data$allele <- colnames(net_data)[grepl("HLA", colnames(net_data))][1]
  return(rank_data)
}

# Read netmhcpan output files and combine results
net_kmers_nonc <- lapply(net_files, read_netmhcpan_nonc) %>% bind_rows() %>% 
  mutate(hla_type = str_replace(allele, "HLA", ""),
         hla_type = str_replace_all(hla_type, "[.]", "")) %>% 
  select(-allele)

# Save results as Rds file for usage by 'Figures_BindingPredictions.Rmd'
saveRDS(net_kmers_nonc, file = file.path(processed_dir, "net_kmers_nonc.Rds"))