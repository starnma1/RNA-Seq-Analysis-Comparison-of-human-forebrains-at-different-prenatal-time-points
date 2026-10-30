mkdir -p logs        # Slurm won't create the log folder
N=$(ls ../fastq_prenatal_forebrain/fastq/*.fastq.gz | wc -l)
sbatch --array=0-$((N-1))%5 ./scripts/slurm-align-fastq.sh