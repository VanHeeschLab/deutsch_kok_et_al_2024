# Deutsch, Kok, Mudge et al. 2024

This repository contains the code for the analyses and figures of the manuscript: High-quality peptide evidence for annotating non-canonical open reading frames as human proteins, Deutsch, Kok, Mudge et al., 2024

The main analysis script is annotate_analyze_samples.Rmd. The necessary input files are described in this script. The script can be run from start to end, and is divided in code blocks based on the figure(panel) that is being created.

At several instances, netMHCpan predictions were performed, for which the scripts run_netmchpan.sh, netmhcpan.sh and config.config are used. These files together allow netMHCpan to run predictions in a parallel slurm job. 

At one part in the code, the script retrieve_netmhcpan_kmers.R has to be run to process the data. This was done because this step of the data processing takes a while.
