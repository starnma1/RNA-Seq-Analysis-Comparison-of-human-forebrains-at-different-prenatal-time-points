mkdir -p logs        # Slurm won't create the log folder
N=$(ls ../raw/fastq/*.fastq.gz | wc -l)
sbatch --array=0-$((N-1))%5 ./scripts/slurm-align-fastq.sh