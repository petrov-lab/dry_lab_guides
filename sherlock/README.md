# Sherlock

The best way to learn about Sherlock is to read the documentation: https://www.sherlock.stanford.edu/docs/

Other good resources are the Slack channels `#sherlock-users` and `#sherlock-announce`

---

### Filesystems

Sherlock has multiple filesystem that serve different purposes. In the Petrov lab, there are 4 (maybe 5) filesystems to be aware of.

`$HOME` aka `/home/users/[stanford_id]/` -- This is your own personal directory that only you can access. Everything stored here is permenant as long as you have an account on Sherlock. There are 15GB of storage, which is not a lot of space. This means you're forced to keep only important, personal files here. I personally keep all of my scripts and github repos here.

`$HOME_GROUP` aka `/home/groups/dpetrov/` -- This is the lab's directory that only lab members can access. Everything stored here is permenant as long as the lab exists. There is 1TB of storage here, which can fill up fast if everyone dumps data into it. Generally it's recommended to have a personal directory within `$GROUP_HOME` and you can install any programs or tools that you may need. I personally keep download tools here as well as certain small files I consistently use such as reference fastas.

`$SCRATCH` aka `/scratch/users/[stanford_id]/` -- This is your own personal directory that only you can access. Everything stored here will be deleted after 3 months. The timer resets when a file is modified. Don't store stuff in scratch long term as it will be lost. There is a 100TB of storeage here, which is a ton. Use scratch for active analyses, for storing data (that you have backed-up elsewhere), and for intermediate files.

`$GROUP_SCRATCH` aka `/scratch/groups/dpetrov/' -- This is the lab's scratch directory. Everything stored here will be deleted after 3 months. The timer resets when a file is modified. There is a 100TB of storage here also. Group scratch can be useful for working on data with other members in lab.

`$OAK` aka `/oak/stanford/groups/dpetrov/' -- This is permenant storage for the lab. You need to be granted access to $OAK. There are 250TB of storage, which the lab pays for. $OAK is really only used by PCG.

Sherlock summarizes all the filesystem space info everytime you login:
```
+---------------------------------------------------------------------------+
| Disk usage for jahemker (group: dpetrov)                                  |
+---------------------------------------------------------------------------+
|    Directory |  volume /   limit [          use%] | inodes /  limit (use%)|
+---------------------------------------------------------------------------+
          HOME |   4.5GB /  15.0GB [|||        30%] |      - /      - (  -%) 
    GROUP_HOME | 781.9GB /   1.0TB [|||||||    76%] |      - /      - (  -%) 
       SCRATCH |  10.9TB / 100.0TB [|          10%] | 252.0K /  20.0M (  1%) 
 GROUP_SCRATCH |   7.2TB / 100.0TB [            7%] | 101.5K /  20.0M (  0%) 
           OAK | 225.4TB / 250.0TB [||||||||   90%] |   6.8M /  37.5M ( 18%) 
+---------------------------------------------------------------------------+
 WARNING: files on SCRATCH and GROUP_SCRATCH are automatically purged
 90 days after their last content modification.
+---------------------------------------------------------------------------+
```


### Nodes

Orthogonal to the Sherlock filesystems, there are nodes. Nodes (as defined in Sherlock) are individual computers that have some number of CPUs, GPUs, RAM, etc. 

There are two basic types of nodes, login nodes and compute nodes. When you first connect to Sherlock via `ssh` you will be on a login node (as denoted by the red `login` next to your username on the command line). Login nodes are not meant for running computationally intensive tasks. If you run a program that uses up too many resources, it will be killed. On a login node you can traverse the filesystems, submit jobs to the cluster, and perform other basic computational tasks. I do a lot of bash file manipulation on login nodes.

Compute nodes are designed for computationally heavy workloads. You have to request a compute node with specific resource counts for a certain amount of time. When you submit and run jobs, these are done on compute nodes. You can also request access to a compute node from the command line, and this would allow you to do more computationally intensive tasks by hand instead of submitting jobs.

Nodes and filesystems are the two halves to the whole Sherlock cluster. You can access any filesystem from any node, but the specific combination of nodes and filesystem that you will want will depend on what you are trying to do.


### Slurm Scheduling System

SLURM documentation: https://slurm.schedmd.com/documentation.html

The real power of using Sherlock is the access to the huge amounts of compute resources. You generally access these resources by submitting a job script explaining what you want to do to the SLURM scheduler, which is in charge of fairly allocating resources to all requesting users.

Job Submission on Sherlock: https://www.sherlock.stanford.edu/docs/getting-started/submitting/?h=job+su#batch-scripts

#### Partitions

All of the nodes in the cluster are divided up into specific partitions. There are a certain number of nodes in each partition, and different partitions have different types of nodes. Each partition has its own specific resource limits with respect to time, CPUs, RAM, GPUs, etc. When you submit a job you want to be smart in choosing the right partitions to minimize your wait times. For example there is a specific partition `gpu` that has GPU resources. There is another partition `bigmem` that has increased amounts of available RAM.

You can see all of the partitions you have access to with the command `sh_part`.
My output is:
```
partition           || nodes         | CPU cores             | GPUs                 || job runtime     | mem/core        | per-node
 name         public ||   idle  total |   idle  total  queued |   idle  total queued || default maximum | default maximum |    cores   mem(GB)  gpus
-----------------------------------------------------------------------------------------------------------------------------------------------------
 normal*      yes    ||      0    218 |    616   5844   12941 |      0      0      0 ||      2h      7d |     6GB     8GB |    20-64   128-384     0
 bigmem       yes    ||      0     11 |    413    824      69 |      0      0      0 ||      2h      1d |     6GB    64GB |   24-256  384-4096     0
 gpu          yes    ||      0     33 |    403   1068     739 |     32    136    123 ||      1h      2d |     8GB    32GB |    20-64  191-2048   3-8
 dev          yes    ||      0      4 |     67    104       0 |     59     64      0 ||      1h      2h |     6GB     8GB |    20-32   128-256  0-32
 service      yes    ||      6      6 |    132    132       0 |      0      0      0 ||      1h      2h |     1GB     8GB |    20-32   128-256     0
-----------------------------------------------------------------------------------------------------------------------------------------------------
 hns          no     ||      0    143 |    432   5164   13190 |      5     12      9 ||      2h      7d |     6GB    26GB |   20-256  128-1536   0-4
 dpetrov      no     ||      0      4 |     42    128       0 |      0      0      0 ||      2h      7d |     8GB     8GB |       32       256     0
-----------------------------------------------------------------------------------------------------------------------------------------------------
 owners       no     ||      7   1740 |   7297  63960    5512 |    826    884      1 ||      2h      2d |     4GB    48GB |   20-256  128-4096   0-8
-----------------------------------------------------------------------------------------------------------------------------------------------------
```




