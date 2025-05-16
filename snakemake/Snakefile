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
