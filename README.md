# Deutsch, Kok, Mudge et al. 2024

This repository contains the code for most of the analyses and figures of the manuscript: High-quality peptide evidence for annotating non-canonical open reading frames as human proteins, Deutsch, Kok, Mudge et al., 2024 ([link](https://www.biorxiv.org/content/10.1101/2024.09.09.612016v1))

## Code
The code directory contains all the necessary scripts. Below, the purpose of each script is briefly described; more detailed information is given in the scripts themselves.

**R files with 'load_' prefix** - These scripts are always called from within an R Markdown file
+ load_tables.R - Loads specific extended data tables
+ load_peptideatlas_data.R - Loads PeptideAtlas data related to the ncORF and canonical peptides detected in immunopeptidomics data
+ load_netmhcpan_data.R - Loads (or prepares) the NetMHCpan binding predictions for the immunopeptidomics data
+ load_plot_aesthetics.R - Defines default theme settings and colors

**R Markdown files** - These scripts are used to analyze the data and generate the plots. The scripts can individually be run from start to end, and are divided in code blocks based on the figure(panel) that is being created.
+ Figures_BindingPredictions.Rmd - For figures that use NetMHCpan binding predictions
+ Figures_nonHLAData.Rmd - For figures based on the non-HLA PeptideAtlas build
+ Figures_ORBL.Rmd - For figures based on the ORBL tool
+ Figures_PeptideAtlas_Immunopeptidomics.Rmd - For figures based on the PeptideAtlas immunopeptidomics data (that do not use NetMHCpan binding predictions)
+ Figures_StructurePredictions.Rmd - For figures related to the structure predictions for ncORFs (e.g. with AlphaFold3)
+ Tier_Peptidein_Counts - For determining number of ncORFs per Tier, and the number of ncORFs with Protein or Peptidein status

**Other scripts** 
+ retrieve_netmhcpan_kmers_nonc.R - Script to process binding predictions for all possible ncORF peptides of length 9. Usage of this script is indicated in the script `Figures_BindingPredictions.Rmd`
+ run_netmhcpan.sh - Script to run netMHCpan using a NetMHCpan docker container. 

## Input files
Nearly all necessary input files for the script are located in the "raw" directory in this repository. A couple of RDS files are used to serve as alternative (smaller) input files compared to the files that were originally used. All files are described below:
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

netMHCpan predictions for c17norep146 and all detected ncORF peptides
- net_c17_146.RDS
- netmhcpan.RDS

File linking each detected peptide to the MS-run
- peptide_sample_msrun_counts.tar.bz

Files containing analysis results about the usage of different proteases, search engines and FDR cutoffs for ncORF detection
+ 20250318_PXD010154_results.xlsx
+ 20250318_R2R_build_different_FDRs.xlsx

ORBL results (will be uploaded upon publication):
+ MSA_Species.xlsx
+ Ribo-Seq_ORF_ORBL.2025-06-03.txt
+ RiboSeqORFs.ORBLq.placental.txt
+ RiboSeqORFs.ORBLq.primate.tsv
+ RiboVsManeOrblv.placental.tsv
+ RiboVsManeOrblv.primate.tsv

The following files are **not** provided but are publicly available:
- Homo_sapiens.GRCh38.102.gtf
- Extended data tables 2, 3, 5, 6, 7, 8, 12, 14 (available upon publication)

At several instances, netMHCpan predictions were performed, for which the script run_netmchpan.sh is used. This script uses a netMHCpan docker container. Prediction results are also provided as RDS files in the raw directory. Only the predictions for all possible length 9 ncORF peptides could not be uploaded due to the size limit. However, these predictions can be performed with the provided scripts.
