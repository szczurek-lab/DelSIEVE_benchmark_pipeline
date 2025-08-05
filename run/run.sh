#!/bin/bash

# USAGE: ./run/run.sh path/to/snakefile


##################################
#           USER SETUP           #
##################################

# conda environment name
ENV_NAME="snake"

# whether to use taskset (1) or not (0)
USETASKSET=1

# if use taskset, set starting and ending cpu indices, both included
STARTCPU=8
ENDCPU=95

# if not use taskset, set the number of cores you want to use
CORENUM=1

# total memory limitation in MB
MEMORY=100000

# other flags of snakemake to attach
ATTACHMENTS="-kp"


###################################
#           DO NOT EDIT           #
###################################

if [[ -z $1 ]]; then
  echo "> Unset snakemake file; use \"Snakefile\" by default."
  SNAKEFILE=Snakefile
else
  SNAKEFILE=$1
fi

# the start of analysis
SECONDS=0

if [[ ${USETASKSET} -eq 1 ]]; then
  echo "> Using ${STARTCPU} to ${ENDCPU} cpus..."
  ((CORENUM = ENDCPU - STARTCPU + 1))

  taskset -c ${STARTCPU}-${ENDCPU} mamba run -n "${ENV_NAME}" snakemake -s "${SNAKEFILE}" --use-conda --resources mem_mb="${MEMORY}" --rerun-triggers mtime --rerun-incomplete --cores "${CORENUM}" ${ATTACHMENTS}
elif [[ ${USETASKSET} -eq 0 ]]; then
  mamba run -n "${ENV_NAME}" snakemake -s "${SNAKEFILE}" --use-conda --resources mem_mb="${MEMORY}" --rerun-triggers mtime --rerun-incomplete --cores "${CORENUM}" ${ATTACHMENTS}
else
  echo "> Set USETASKSET to 1 or 0 in order to use taskset or not!"
fi

# the end of analysis
DURATION=$SECONDS

echo "> Time elapsed: $((DURATION / 3600)) hours, $(((DURATION / 60) % 60)) minutes and $((DURATION % 60)) seconds"
