# Run from the root of the repo, as `source ./scripts/build-align-environment.sh`

export MAMBA_ROOT_PREFIX=$PWD/micromamba
export MAMBA_NO_RC=true
export MAMBA_CHANNEL_PRIORITY=strict
micromamba clean --all --yes
micromamba create -f ./environment.align.yml
micromamba env export -n sys-genomics-align --explicit --md5 > environment.align.lock.txt