# System Genomics

## Set up environment

### Install `micromamba` on Euler, our environment manage:
1. First run:
    ```
    mkdir -p ~/.local/bin
    cd ~/.local/bin
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba --strip-components=1
    ```
1. Set up your shell to configure micromamba when you log in, add this to your ~/.bashrc:
    ```
    export MAMBA_ROOT_PREFIX=$HOME/micromamba
    export PATH=$HOME/.local/bin:$PATH
    eval "$(micromamba shell hook -s bash)"
    ```

### Environment setup
#### Normal build from lock
1. To load the environment with all of the packages we need, from the root of the repo, run `source ./scripts/build-environment.sh`
1. Then run `micromamba activate sys-genomics-align`
#### From scratch (de novo)
1. If you want to build the environment *from scratch*, from the root of the repo, run `source ./scripts/build-environment.sh`