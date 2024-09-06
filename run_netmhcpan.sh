#!/bin/bash

######################################################################
#
# Authors: l.w.kok-15@prinsesmaximacentrum.nl
# Date: 22-11-2022
#
######################################################################

function usage() {
    cat <<EOF
SYNOPSIS
  app_mainscript.sh ./path/to/config.config [run_identifier] - run NetMHCpan antigen prediction
  app_mainscript.sh help - display this help message
DESCRIPTION
  Run NetMHCpan to predict antigens

AUTHOR
  Leron Kok, MSc
EOF
}

# Source all variables from the config file
CONFIG=$1

# Read or create an id to identify the run and the generated output
runid=${2:-$(date +"%y%m%d_%H%M%S")}

# Show help message if there was no config file location given on the commandline
if [[ -z $1 ]]; then usage; exit; fi

source ${CONFIG}

# Create log directory
mkdir -p ${logdir}

################################################################################
#
# Run the pipeline
#
################################################################################

echo -e "\n`date` Running NetMHCpan ..."
echo -e "====================================================================================== \n"

#MHCFlurry. Predict how well each peptide binds to mHC-I

netmhcpan_jobid=()

netmhcpan_jobid+=($(sbatch --parsable \
    --time=48:00:00 \
    -c 1 \
    --mem 32G \
    --gres=tmpspace:20G \
    --array 1-1%${simul_array_runs} \
    --output ${logdir}/${runid}_slurm-%A_%a.out \
    ${scriptdir}/netmhcpan.sh \
    ${CONFIG} 
))

echo -e "\n ====== `date` Started all jobs! ====== \n"