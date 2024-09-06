#!/bin/bash

######################################################################
#
# Authors: l.w.kok-15@prinsesmaximacentrum.nl
# Date: 22-10-2021
#
######################################################################

# Show start time, program and input file
start=$(date +"%d-%m-%Y %T")
echo "Start : $start"
echo "Program : NetMHCpan"

# Load parameters from main script
source $1
mapfile -t allele_list < ${net_alleles}
net_allele=${allele_list[$((SLURM_ARRAY_TASK_ID-1))]}

netmhcpan_out="${netmhcpan_outdir}/net_annotatedalleles5_${SLURM_ARRAY_TASK_ID}.txt"
netmhcpan_outxls="${netmhcpan_outdir}/net_annotatedalleles5_${SLURM_ARRAY_TASK_ID}.xls"
netmhcpan_input="${netmhcpan_indir}/netmhcpan_input${SLURM_ARRAY_TASK_ID}.txt"

echo "Input read from ${netmhcpan_input}"
# Create output directory
mkdir -p ${netmhcpan_outdir}

# Run NetMHCpan
singularity exec --bind ${wd} --bind ${TMPDIR} --writable-tmpfs --containall --env TMPDIR=${TMPDIR} /hpc/local/Rocky8/pmc_vanheesch/singularity_images/netmhcpan-4.1b.sif /app/package/netMHCpan-4.1/netMHCpan ${netmhcpan_input} -BA -a ${net_allele} -p -xls -xlsfile ${netmhcpan_outxls} > ${netmhcpan_out}

# Show end time
echo "Output written to ${netmhcpan_out}"
end=$(date +"%d-%m-%Y %T")
echo "End : $end"