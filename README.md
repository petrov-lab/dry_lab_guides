### Introduction
https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

Snakemake is a text-based workflow system.

Workflows are defined with rules, which determine how to create output files from input files.

Dependencies between rules are determined, creating a directed-acyclic graph (DAG) of jobs that can be parallelized.

Snakemake allows each rules to be constrained custom resources, and it provides support for distributed computing (ie clusters).

Snakemake integrates with Conda and Singularity so that the required software also becomes part of the workflow itself.

Snakemake extends Python, so a working knowledge of Python syntax is recommended for learning Snakemake.

### When to use Snakemake

Snakemake is really useful for computational pipelines that require multiple steps, especially if those steps require different resources and/or can be completed in parallel.

Snakemake is most useful for pipelines that will be run often for multiple samples, but it can also be useful for general reproducibility purposes (i.e. for publications)

Snakemake is probably overkill for one-step commands or commands you will only ever run once.

### Conceptualizing Snakemake workflows

Snakemake determines which rules to run by traversing the DAG. It starts with nodes that do not have any outgoing edges and works backwards, collecting dependencies.

It can be useful to write a Snakemake pipeline in a similar way, starting with the final output file(s), and working backwards, one step/rule at a time.

In our case, the nodes of the graph will be rules, and the edges connect rules based on a sharing of input/output files.
### Rules

Snakemake rules should be a single, discrete step (at your discretion) in the workflow. A rule will take in some number of `input` files and spit out some number of `output` files. The rule will transform the `input` files to the `output` files using either `shell` or `script` commands. As so much bioinformatics happens with bash, I've only ever used `shell`. Additional `params` necessary for the execution of the rule can also be passed. All rules are defined in a `Snakefile`.

```python
rule align_reads:
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/dmel1.fastq.gz"
    output:
        bamfile="/scratch/users/jahemker/bams/dmel1.reads.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/dmel1.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```

Importantly, the `shell` section does not follow bash syntax for variables, but instead uses Snakemake's variables syntax.

Ultimately our goal is to be submitting these jobs to Sherlock, and for that we will need to define the resources each rule requires. We do this in the `resources` field of the rule. Later, we will see how to set-up default values for the different resources so that we do not have to always list everything out. (For example, the partition is not explicitly stated here, because we are ok with using with the default partition value is)
```python
rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/dmel1.fastq.gz"
    output:
        bamfile="/scratch/users/jahemker/bams/dmel1.reads.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/dmel1.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```

---

If we are starting at the very end of our workflow, how do we write the rule for the final files we want? We can use a `rule all` to anchor our workflow. This rule only has input files, and naturally will come at the end of the DAG.

```python
rule all:
	input:
		"/scratch/users/jahemker/bams/dmel1.reads.sorted.bam"
```

Snakemake does not blindly run every rule written in the workflow. If a rule is not a dependency for something else, it will not be run. It's important to make sure that the `rule all` accurately lists all necessary outputs.

### Scaling rules to multiple samples
The power of Snakemake is it's ability to easily scale a workflow to many samples. We do this by generalizing our rules using wildcards. For every value of a wildcard, snakemake will run that rule, populating all instances of the wildcard with the same value.
```python
rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/{sample}.fastq.gz"
    output:
        bamfile="/scratch/users/jahemker/bams/{sample}.reads.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/{sample}.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```
Note: if you need to use wildcards within the `shell` section of the rule, you have to use `{wildcards.sample}` instead of just `{sample}` like in the `input`/`output`/`params` sections.

How does Snakemake know what values to fill the wildcards with? It again goes back to defining the necessary output files. If we are using wildcards, then it's likely we will have multiple, almost-identical files. Snakemake has a special function `expand` that will take a pattern and automatically populate it with values from a given list. We can use this in our `rule all`.
```python
rule all:
	input:
		expand("/scratch/users/jahemker/bams/{sample}.reads.sorted.bam", sample=SAMPLES)
```
We just need to make sure that the list of `SAMPLES` is defined somewhere, and we can just use basic python for that. In full, our `Snakefile` looks like
```python
SAMPLES = ["dmel1","dmel2","dmel3"]

rule all:
	input:
		expand("/scratch/users/jahemker/bams/{sample}.sorted.bam", sample=SAMPLES)

rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/{sample}.fastq.gz"
    output:
        bamfile="/scratch/users/jahemker/bams/{sample}.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/{sample}.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```
This will give us three output `bam` files. 
```
/scratch/users/jahemker/bams/dmel1.sorted.bam
/scratch/users/jahemker/bams/dmel2.sorted.bam
/scratch/users/jahemker/bams/dmel3.sorted.bam
```
### Running Snakemake on Sherlock

Snakemake is (relatively) easy to set up and run on Sherlock. This lets us take advantage of the SLURM scheduling system and automatically submit jobs as dependencies are fulfilled.

An over-arching Snakemake job can be submitted to run for a long time, and that job will submit each rule as needed. An example `sbatch` submission script:
```bash
#! /bin/bash
#
#SBATCH --job-name=snakemake
#SBATCH --time=7-00:00
#SBATCH --partition=dpetrov,hns
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

. /home/users/jahemker/.bashrc
conda activate snakemake

snakemake --profile profile/
```
Before submitting this job to Sherlock, you should do a dry-run of your Snakemake workflow. This will allow you to make sure that there aren't any syntax errors with the `Snakefile` and other config files. It also allows you to see what rules will be run.
```bash
snakemake --profile profile/ -n
```
You can either install Snakemake manually or you can install it into a conda environment. The following command creates a new environment called snakemake, which contains Snakemake and the extra plugin we need to be able to use Snakemake with Sherlock.
```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake snakemake-executor-plugin-cluster-generic
```
Next, we can look at the actual Snakemake command that we will be running. It is surprisingly simple.
```bash
snakemake --profile profile/
```
It's so simple because of how we are going to set up our Snakemake directory . All of the needed information will be kept in config files instead of in a command with 10+ arguments and flags.
### Snakemake organization

```
workflow_1/
├── logs/
├── profile/
│ ├── config.yaml
│ └── status_sacct.sh
├── Snakefile
├── parameters_config.yaml
├── submit_snakemake.sh
└── README.md
```
For a given workflow, we need a directory that has a `Snakefile`, a `parameters_config.yaml`, and a subdirectory called `profile` that holds a `config.yaml` and a Sherlock-specific script called `status_sacct.sh`.

`Snakefile` holds the rule definitions for the workflow. 

`parameters_config.yaml` holds all of the user-changeable parameters (file locations, mutable command values) that are found in `Snakefile`.

`submit_snakemake.sh` is the sbatch script we will use to send the Snakemake pipeline to Sherlock.

`README.md` is optional, but always good practice. Really useful when there are a lot of user-defined parameters.

Within the `profile/` subdirectory, `config.yaml` holds all of the arguments for the `snakemake` command we will be submitting to Sherlock. `status_sacct.sh` is a command that allows Snakemake to monitor the status of jobs on Sherlock. This lets it know if a job fails or times out, or successfully completes. 

We can see that this subdirectory `profile/` is the only parameter we are passing to the `snakemake` command. It will know to parse the directory for a `config.yaml` and find the relevant settings.

The subdirectory `logs` will be automatically created and will hold the SLURM output files.

Note: `parameters_config.yaml`, `submit_snakemake.sh`, and `profile/` can be renamed to other things if you wish.

---

Let's briefly look at `profile/config.yaml` to see what arguments we will be running the `snakemake` command with.
```bash
executor: cluster-generic
configfile: parameters_config.yaml
cluster-generic-submit-cmd:
  mkdir -p logs/ &&
  sbatch
    --partition={resources.partition}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name={rule}
    --output=logs/{rule}-{wildcards}-%j.out
    --time={resources.time}
    --parsable
cluster-generic-cancel-cmd:
  scancel
cluster-generic-status-cmd:
  profile/status_sacct.sh
default-resources:
  - partition=hns,normal,dpetrov
  - mem_mb=4000
  - time="1:00:00"
  - cpus_per_task=1
rerun-triggers: mtime
restart-times: 0
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 100
keep-going: True
rerun-incomplete: True
printshellcmds: True
```
And for reference, this is what `profile/status_sacct.sh` looks like. It shouldn't ever need to be modified, but you do need to give it executable permissions.
```bash
chmod +x profile/status_sacct.sh
```

```bash
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
```
### Conda and Singularity with Snakemake

One of the most powerful aspects of Snakemake is the ability to have containers and virtual environments built into the pipeline. This is one of the key reasons why Snakemake is so reproducible. For each rule, you can provide a replicable environment or container that will be able to run everything needed in that rule.

When using conda with Snakemake, you need to provide a path to a `yaml` file that describes the environment. To generate a `yaml` file of a given environment use the following command.
```bash
# with the conda env active
conda env export --no-builds > environment.yaml
```
Then we, 1) pass this file location to the rule that needs the conda env.
```python
rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    conda:
		"/scratch/users/jahemker/environment.yaml"
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/{sample}.fastq.gz"
    output:
		bamfile="/scratch/users/jahemker/bams/{sample}.reads.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/{sample}.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```
Next, we need to add some new arguments to our `profile/config.yaml` to tell Snakemake to use conda.
```bash
# Add these to use conda
use-conda: True
conda-prefix: /home/groups/dpetrov/jahemker/miniconda3/envs/
```
`use-conda` simply says to use conda
`conda-prefix` is optional, but it tells Snakemake where to download the conda env to.

TBD actual conda downloading behaviors. I think it follows `conda-prefix` and if there is no `conda-prefix` set then it downloads into the typical conda directories.

---
Singularity is generally easier to use, but requires a lot more effort to make a container if you don't already have it. Here, we need to put the location of the singularity image. Docker works similarly, but Sherlock doesn't allow Docker.
```python
rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    singularity:
	    "/scratch/users/jahemker/aligner.sif"
    input:
        ref="/scratch/users/jahemker/D.melanogaster.fa",
        ref_idx="/scratch/users/jahemker/D.melanogaster.fa.fai",
        reads="/scratch/users/jahemker/reads/{sample}.fastq.gz"
    output:
        bamfile="/scratch/users/jahemker/bams/{sample}.reads.sorted.bam"
    params:
        samfile="/scratch/users/jahemker/bams/{sample}.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```
Again we need to add certain arguments to `profile/config.yaml`.
```bash
# Add these to use singularity
sdm: apptainer
# If using GPU nodes with singularity
singularity-args: "--nv"
```

### Generalizing workflows with user-defined parameters
Currently, our rules have hard-coded paths and values that make it hard for other users to easily use the pipeline. We can remedy this by having all of these values in an easy-to-edit config file (`parameters_config.yaml`) that then populates the main `Snakefile`.
```bash
# contents of parameters_config.yaml
reference:
  location: "/home/groups/dpetrov/jahemker/reference/v6.63/D.melanogaster.fa"
project:
  dir: "/scratch/users/jahemker"
samples:
  "/scratch/users/jahemker/sv_popgen/samples.txt"
envs:
  genometools: "/home/users/jahemker/drosophila_sv_popgen/envs/genometools.yaml"
```
We can then add values from this file into the `Snakefile` rules.
```python
SAMPLES = open(config["samples"], 'r').readlines()
SAMPlES = [x.rstrip() for x in SAMPLES]

rule all:
	input:
		expand(f"{config['project']['dir']}/bams/{{sample}}.sorted.bam", sample=SAMPLES)

rule align_reads:
	resources:
		time = "4:00:00",
        cpus = 16,
        mem_mb = 32000
    conda:
	    config['envs']['genometools']
    input:
        ref=config['reference']['location'],
        ref_idx=f"{config['project']['dir']}/D.melanogaster.fa.fai",
        reads=f"{config['project']['dir']}/reads/{{sample}}.fastq.gz"
    output:
        bamfile=f"{config['project']['dir']}/bams/{sample}.reads.sorted.bam"
    params:
        samfile=f"{config['project']['dir']}/bams/{{sample}}.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```

### Generalizing resource usage for variable-resource rules
If you have rules that can greatly change the amount of resources needed depending on the sample, it can be useful to have auto-changing resource values. Using this setup, if a job fails (with the implication being it didn't have enough resources), then it will resubmit the job with incrementally higher resources. Changes need to be made in `profile/config.yaml` to enable job resubmission as well as the rules themselves.
```bash
# Change the following argument in profile/config.yaml
restart-times: 3
```
New `resources` section:
```python
rule align_reads:
	resources:
		time = lambda wildcards, attempt: 240 + 240 * (attempt-1),
        cpus = lambda wildcards, attempt: 16 + 16 * (attempt-1),
        mem_mb = lambda wildcards, attempt: 16000 + 16000 * (attempt-1)
    conda:
	    config['envs']['genometools']
    input:
        ref=config['reference']['location'],
        ref_idx=f"{config['project']['dir']}/D.melanogaster.fa.fai",
        reads=f"{config['project']['dir']}/reads/{{sample}}.fastq.gz"
    output:
        bamfile=f"{config['project']['dir']}/bams/{sample}.reads.sorted.bam"
    params:
        samfile=f"{config['project']['dir']}/bams/{{sample}}.reads.sam"
    shell:
        '''
        minimap2 -a -x map-ont {input.ref} {input.reads} > {params.samfile}
        samtools sort -o {output.bamfile} {params.samfile}
        samtools index -b {output.bamfile}
        rm {params.samfile}
        '''
```

### Debugging Snakemake on Sherlock
Maybe the worst part of using Snakemake is debugging before you have a working pipeline. To make sure the `Snakefile` and all requisite config files are working correctly, you can do a dry-run of the workflow with the `-n` flag.
```bash
snakemake --profile profile/ -n
```
This will make sure that the syntax and logic of all of the Snakemake rules and files adds up. What it won't do is make sure that the commands that are run within the `shell` portion of the rules work.

In order to debug specific rules themselves, you'll just have to look at the rule-specific slurm logs that are generated and dumped into the `log/` directory.

In isolated cases, jobs will fail in a way that does not notify the overarching Snakemake job. This means that the overarching job continues running even though it won't be submitting anything else. In this case, the overarching job needs to be manually canceled. When a Snakemake workflow is running, it 'locks' the directory to prevent other instances of the same workflow from running and messing up files. Normally, it unlocks the directory when it's finished, but manually cancelling the overarching job keeps the directory locked, so we have to manually unlock it.
```bash
#manually cancel a job by using its job id.
#job id can be found by running 'squeue -u $USER'
scancel [job_id]
snakemake --unlock -s Snakefile --profile profile/
```

### Warnings when using Snakemake
The only big warning is to make sure that your job counts aren't going to be insane when running a workflow. Sherlock limits job submissions to something like 1000 per hour (and labs to like 3000 per hour). Depending on how many samples you have and how many potentially parallel rules can be run this is a feasible limit to hit. You also want to make sure that you aren't asking for copious resources that you aren't using, as low usage of resources can get your account flagged. 

Make sure to condense rules sufficiently. If a command takes just a few seconds to run, its best to lump multiple together, as the workflow bottleneck will come at the time it takes to actually submit the job.

### More info
There are a lot more Snakemake tricks and functions that can be used to help make workflows more efficient. I'd definitely recommend looking through the documentation or googling to see what other people have done. https://snakemake.readthedocs.io/en/stable/index.html

For working Snakemake workflows that use different kinds of scripts, singularity containers, GPU nodes, and other stuff, see https://github.com/jahemker/drosophila_ultralong_sv_calling/workflows/
