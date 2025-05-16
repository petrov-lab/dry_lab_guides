#!/usr/bin/env bash

# Modified from https://github.com/jdblischak/smk-simple-slurm/blob/main/extras/status-sacct.sh
# deal with sacct calls constantly failing on Sherlock
# Check status of Slurm job

jobid="$1"

if [[ "$jobid" == Submitted ]]
then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk '{print $1}' || echo "RUNNING"`

# slurm failed job codes from https://slurm.schedmd.com/squeue.html#lbAG

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(BOOT_FAIL|CANCELLED|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|REVOKED|SPECIAL_EXIT|STOPPED|SUSPENDED|TIMEOUT).* ]]
then
  echo failed
else
  echo running
fi
