#!/bin/bash

# Function to print usage
print_usage() {
    echo "Usage: $0 -p <peptide_file(s)> -a <allele_file> [-o <output_dir>]"
    echo "  -p: Path to peptide file(s). Use space-separated list for multiple files."
    echo "      Wildcard characters are supported, e.g., '*.txt' or 'peptide_*.fasta'"
    echo "  -a: Path to allele file."
    echo "  -o: Output directory (optional, default: current directory)"
    echo "  -l: Maximum number of concurrent jobs (optional, default: 20)"
    exit 1
}

# Parse command line arguments
PEPTIDE_FILES=()
JOB_LIMIT=20
while getopts "p:a:o:l:" opt; do
    case $opt in
        p) 
            # Handle wildcard expansion
            for file in $OPTARG; do
                PEPTIDE_FILES+=("$(readlink -f "$file")")
            done
            ;;
        a) ALLELE_FILE=$(readlink -f "$OPTARG") ;;
        o) OUTPUT_DIR=$(readlink -f "$OPTARG") ;;
        l) JOB_LIMIT=$OPTARG ;;
        *) print_usage ;;
    esac
done

# Check if required arguments are provided
if [ ${#PEPTIDE_FILES[@]} -eq 0 ] || [ -z "$ALLELE_FILE" ]; then
    print_usage
fi

# Set default output directory if not provided
OUTPUT_DIR=${OUTPUT_DIR:-.}

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Array to store job IDs
declare -a JOB_IDS

# Function to submit a single netMHCpan job
submit_job() {
    local peptide_file=$1
    local allele=$2
    local job_name="netMHCpan_${peptide_file##*/}_${allele}"
    local output_file="${OUTPUT_DIR}/${job_name}.out"

    # Check if we need to add a dependency
    local dependency=""
    if [ ${#JOB_IDS[@]} -ge $JOB_LIMIT ]; then
        local oldest_job=${JOB_IDS[0]}
        dependency="#SBATCH --dependency=afterany:$oldest_job"
        JOB_IDS=("${JOB_IDS[@]:1}") # Remove the oldest job from the array
    fi

    
printf """#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$output_file
#SBATCH --ntasks=1
#SBATCH --gres=tmpspace:20G
#SBATCH --mem=32G
#SBATCH --time=48:00:00
$dependency


# Ensure TMPDIR is set and exists
TMPDIR=\${TMPDIR:-/tmp}
mkdir -p \$TMPDIR

# Run Singularity container
singularity exec \\
    --bind /hpc \\
    --bind \$TMPDIR \\
    --writable-tmpfs \\
    --containall \\
    --env TMPDIR=\$TMPDIR \\
    /hpc/local/Rocky8/pmc_vanheesch/singularity_images/netmhcpan-4.1b.sif \\
    /app/package/netMHCpan-4.1/netMHCpan \\
    -p "$peptide_file" \\
    -BA \\
    -a $allele \\
    -xls \\
    -xlsfile "${output_file}.xls" \\
    > "${output_file}.txt"
""" > "${output_file}_script.txt"

# Submit the job and capture the job ID
    local job_id=$(sbatch  --parsable ${output_file}_script.txt)

    # Add the new job ID to the array
    JOB_IDS+=($job_id)

    echo "Submitted job $job_id for $peptide_file with allele $allele"
}

# Read alleles from file
mapfile -t ALLELES < "$ALLELE_FILE"

# Case 1 & 2: Single peptide file
if [ ${#PEPTIDE_FILES[@]} -eq 1 ]; then
    if [ ${#ALLELES[@]} -eq 1 ]; then
        # Case 1: Single peptide file, single allele
        submit_job "${PEPTIDE_FILES[0]}" "${ALLELES[0]}"
    else
        # Case 2: Single peptide file, multiple alleles
        for allele in "${ALLELES[@]}"; do
            submit_job "${PEPTIDE_FILES[0]}" "$allele"
        done
    fi
else
    # Case 3 & 4: Multiple peptide files
    if [ ${#ALLELES[@]} -eq 1 ]; then
        # Case 3: Multiple peptide files, single allele
        for peptide_file in "${PEPTIDE_FILES[@]}"; do
            submit_job "$peptide_file" "${ALLELES[0]}"
        done
    else
        # Case 4: Multiple peptide files, multiple alleles
        num_jobs=$(( ${#PEPTIDE_FILES[@]} < ${#ALLELES[@]} ? ${#PEPTIDE_FILES[@]} : ${#ALLELES[@]} ))
        for ((i=0; i<num_jobs; i++)); do
            submit_job "${PEPTIDE_FILES[i]}" "${ALLELES[i]}"
        done
    fi
fi

echo "All jobs submitted successfully."
echo "Number of peptide files processed: ${#PEPTIDE_FILES[@]}"
echo "Number of alleles processed: ${#ALLELES[@]}"
echo "Total number of jobs submitted: $((${#PEPTIDE_FILES[@]} * ${#ALLELES[@]}))"