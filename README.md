# System Genomics

## Set up environment

### Install `micromamba` on Euler, our environment manager:
1. First run:
    ```
    mkdir -p ~/.local/bin
    cd ~/.local/bin
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba --strip-components=1
    ```
1. Set up your shell to configure micromamba when you log in, add this to your `~/.bashrc`:
    ```
    export PATH=$HOME/.local/bin:$PATH
    eval "$(micromamba shell hook -s bash)"
    ```
1. Then run `bash ~/.bashrc`

### Environment setup
#### Normal build from lock
1. To load the environment with all of the packages we need, from the root of the repo, run `bash ./scripts/sync-align-environment.sh`
1. Then copy and paste the contents of `./scripts/activate-align-environment.sh` into your terminal.
#### From scratch (de novo)
1. If you want to build the environment *from scratch*, from the root of the repo, run `bash ./scripts/build-align-environment.sh`

## Gene counting on Euler

The counting job uses R 4.4.1 from the `2025-06` stack and requires the
Bioconductor packages `Rsubread` and `rtracklayer` in an R library visible to
the job. For example, after configuring a writable personal R library, install
missing packages in R with:

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("Rsubread", "rtracklayer"), ask = FALSE, update = FALSE)
```

Start R with `module load stack/2025-06` and `module load r/4.4.1` first,
and use the same library for the job. From the repository root, submit with:

```sh
bash alignment/scripts/submit-feature-counting.sh
```

The job reads BAM files under `alignment/mapping` and the GENCODE annotation
at `index/reference_genome/gencode.v50.basic.annotation.gtf`. It writes
unique and multimapping count outputs to `counts/`, with logs in
`alignment/logs/`. The `64G`, 8-CPU, 24-hour resource request is an initial
estimate and can be adjusted in the Slurm script.