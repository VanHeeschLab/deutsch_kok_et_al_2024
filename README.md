# Deutsch, Kok, Mudge et al. 2024

This repository contains the code for the analyses and figures of the manuscript: High-quality peptide evidence for annotating non-canonical open reading frames as human proteins, Deutsch, Kok, Mudge et al., 2024 ([link](https://www.biorxiv.org/content/10.1101/2024.09.09.612016v1))

The main analysis script is annotate_analyze_samples.Rmd. The necessary input files are described in this script. The script can be run from start to end (after specifying work_dir under the header "Define directories, files and colors"), and is divided in code blocks based on the figure(panel) that is being created. Some blocks only have to be run once, and for some blocks only one of the two has to be run. This is specified in the script

Nearly all necessary input files for the script are located in the "raw" directory in this repository. A couple of RDS files are used to serve as alternative (smaller) input files compared to the files that were originally used. All files are described below:
Supplementary tables from the manuscript:
- 060924_Supp_Table_S2.xlsx
- 060924_Supp_Table_S3.xlsx
- 060924_Supp_Table_S4.xlsx
- 060924_Supp_Table_S5.xlsx
- 060924_Supp_Table_S6.xlsx

List of cancer genes according to the Cancer Gene Census
- Census_allThu_Jan_4_14_08_58_2024.csv

List of canonical protein IDs from PeptideAtlas
- Core20k.txt

Mean FPKM expression of genes in GTEX (excluding testis)
- GTEX_FPKMmean_expression.txt

Statistics from PeptideAtlas on the human HLA (2023-11) and non-HLA (2023-06) builds
- HLA2023-09_experiment_summary.xlsx
- Non-HLA2023-09_experiment_summary.tsv

R data object containing the canonical protein and ncORF sequences (alternative to Homo_sapiens.fasta from PeptideAtlas)
- can_nonc_seq.RDS

Filtered canonical and non-canonical peptides (alternative to peptide_mapping.tsv from PeptideAtlas)
- filtered_peptides.tar.bz2

List of the 7,264 ncORFs (from doi: 10.1038/s41587-022-01369-0)
- ncorf_list.xlsx

netMHCpan predictions for c17norep146, c5norep142 and all detected peptides
- net_c17_146.RDS
- net_c5_142.RDS
- netmhcpan.RDS

File linking each detected peptide to the MS-run
- peptide_sample_msrun_counts.tar.bz

The following files are **not** provided but are publicly available:
- Homo_sapiens.GRCh38.102.gtf
- GWIPS-viz data: All human files from ‘Initiating Ribosomes (P-site)’ with track ‘Global Aggregate’, and all human files from ‘Elongating Ribosomes (A-site)’ with track ‘Global Aggregate’. ([link](https://gwips.ucc.ie/downloads/index.html))

At several instances, netMHCpan predictions were performed, for which the script run_netmchpan.sh is used. This script uses a netMHCpan docker container. Prediction results are also provided as RDS files in the raw directory. Only for the predictions for figures S4A & B the predictions could not be uploaded due to size limits. 

At one part in the code, the script retrieve_netmhcpan_kmers.R has to be run to process the data. This was done because this step of the data processing takes a while.
