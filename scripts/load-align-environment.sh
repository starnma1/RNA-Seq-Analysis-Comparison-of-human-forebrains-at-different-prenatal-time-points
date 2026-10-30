# Run from the root of the repo, as `source ./scripts/load-align-environment.sh`

export MAMBA_ROOT_PREFIX=$PWD/micromamba
export MAMBA_NO_RC=true
export MAMBA_CHANNEL_PRIORITY=strict
micromamba clean --all --yes
micromamba create -n sys-genomics-align -f ./environment.align.lock.txt