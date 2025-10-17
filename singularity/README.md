# Using Singularity on Sherlock

## Why use Singularity?

Singularity is a container system similar to Docker. Sherlock uses Singularity instead of Docker. It is analogous to `conda` except it represents an entire operating system. This means that when you access a Singularity container, you are basically accessing a virtual computer completely different from the one you are currently using (as far as software goes. Hardware stays the same). 

Why is this useful for us? If you have ever had to try and install a specific program (ie `bcftools` or `bedtools`), you may have seen that the programs has dependencies (other installed programs required to make it run). Those dependencies can have dependencies of their own, and you can get a whole chain of installations required for the one program you actually care about.

Alternatively, a singularity or docker container can have that specific program and all of its dependencies already installed. You can then access this container (which is just a single file on your computer) specifically when you need to run the program. You no longer have to worry about dependencies or versions of programs.

Using containers also leads to much better reproducibility, as once a container is made, it is fixed and cannot be changed. If you give that same container to someone else, they will be running the exact same version of the program that you did. Of course, programs get updated and you might want to use the newer version. In that case you would create or download a new version of the container.

A container does not have to have just a single program on it (remember it is a virtual operating system, think of how many programs are installed on your computer... you could convert your whole computer to a container, essentially freezing it in time with all of its programs, but that's not a good idea in reality). Similar programs can be grouped together to provide ease of use. For example `samtools` and `bcftools` are closely related and used programs, and they normally come bundled together.

What does this mean in practice?

Let's say I need to call SNPs from some reads. I need an aligner program to align my reads to a reference (ie. bwa-mem). I need a variant caller (ie. GATK) to call SNPs. I want a program to easily allow me to manipulate the resulting VCF file (ie bcftools). There are a number of ways I could install these programs onto sherlock so that I can use them. Sherlock might already have them installed and I can activate them with `module load`. I could install them into a conda environment with `conda install` (assuming they are on conda). I could manually install them from their github repos or websites.

I could also grab the (already existing) singularity containers for these programs. This will leave me with three container files on my computer to deal with. I can run all of my code with it, and then when I'm writing the paper I can say I used bwa-mem (https://hub.docker.com/r/intelliseqngs/bwa-mem/tags tag:3.2.0), gatk (https://hub.docker.com/r/broadinstitute/gatk/tags tag:4.6.2.0), and bcftools (https://hub.docker.com/r/staphb/bcftools/tags tag:1.22). Anyone reading that can instantly access the same version of tools you used.

One of the key words in the last paragraph was **already existing**. Almost every bioinformatics tool I have searched for already has a container made for it. This is why it's so easy for us to use them. Making containers is definitely harder and more time-consuming but you basically will never have to do that.

The place where all these containers are stored is called [dockerhub](https://hub.docker.com/). Technically, all these containers are docker containers, but singularity can read them and build singularity image files out of them.

## Find and install a Singularity image

How do you actually get and use a container? First look for the program you want on [dockerhub](https://hub.docker.com/). For this example, we will look for a container of `bcftools`.

Searching for bcftools on [dockerhub](https://hub.docker.com/) gives you these results: https://hub.docker.com/search?q=bcftools. As you can see there are a lot of options, as anyone can upload a docker container (similar to a github repo). Generally look for well updated containers that have been downloaded a lot. 

We will use the first option that comes up: https://hub.docker.com/r/staphb/bcftools. Clicking on it will give you an overview of the container (if they wrote a readme). There is also a tab called "Tags". Tags are just what dockerhub calls the version number.

If we look at the Tags tab and sort by newest, we see that there is a tag called "latest" and a tag call "1.22" which both were pushed 4 months ago. Pretty much every container has a tag called "latest". This allows you to just download the most recent version without having to actually know the version number. You can then check and see what the most recent non-latest tag is. That is the actual version number that you will want to record for your records.

Even though it doesn't really make a difference to us, let's say we will download the container with the tag "1.22" on to Sherlock. Note how big the container is. For this container its only ~50Mb, but some containers can easily be 5-10Gb.

On Sherlock, find a place to put the containers. I have a directory where I install them all on `$SCRATCH`, as they can take up space and if they get deleted they are super fast to reinstall. When we install these containers on to Sherlock, we are actually converting the docker container to a singularity image.

When downloading from dockerhub, the command is of the generic form:
```
singularity build [your_name_for_the_image].sif docker://[username][/container_name]:[tag]
```
For this bcftools container it is specifically:
```
singularity build bcftools.sif docker://staphb/bcftools:1.22
```

Note: with larger containers, this can be a bit computationally intensive, and if you do it on a login node, sherlock will sometimes kill the command. I have added a simple SLURM script that will build the image.

## Use Singularity images to run programs

We now have a singularity image with our program of choice installed on Sherlock.

We can access the operating system of the container using:
```
singularity shell bcftools.sif
```
This puts us in a new shell that's inside the container. If we run `bcftools`, we will now get the standard help/error message showing that it has been properly installed.

If we have a specific `bcftools` command we want to run, we can also just run it without having to open the shell of the container:
```
singularity exec bcftools.sif bcftools view my_file.vcf
```

Finally, another great reason to use singularity containers is that they have great functionality with Snakemake. You can designate specific containers to Snakemake rules, and those rules will be automatically run inside those containers.
