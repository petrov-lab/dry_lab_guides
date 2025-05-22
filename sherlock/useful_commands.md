# Useful commands

Here's a list of commands that are useful when working on Sherlock. 

### Aliases and functions 

Personally, I alias some of these commands, meaning I have a specific keyword that will run the entire command. This is done by modifying your `.bashrc` (found at `/home/users/[stanford_id].bashrc`). Here is a simple example:
```bash
#in .bashrc
alias scratch="cd $SCRATCH"
```
With this alias, when I type `scratch` it executes the command `cd $SCRATCH` to take me to my scratch directory.

Functions are basically aliases for more complex operations and they can also be added to `.bashrc`. Here is a complex example function I got from some stackoverflow comment that will broadly uncompress different filetypes:
```
function extract {
 if [ -z "$1" ]; then
    # display usage if no parameters given
    echo "Usage: extract <path/file_name>.<zip|rar|bz2|gz|tar|tbz2|tgz|Z|7z|xz|ex|tar.bz2|tar.gz|tar.xz>"
    echo "       extract <path/file_name_1.ext> [path/file_name_2.ext] [path/file_name_3.ext]"
    return 1
 else
    for n in $@
    do
      if [ -f "$n" ] ; then
          case "${n%,}" in
            *.tar.bz2|*.tar.gz|*.tar.xz|*.tbz2|*.tgz|*.txz|*.tar) 
                         tar xvf "$n"       ;;
            *.lzma)      unlzma ./"$n"      ;;
            *.bz2)       bunzip2 ./"$n"     ;;
            *.rar)       unrar x -ad ./"$n" ;;
            *.gz)        gunzip ./"$n"      ;;
            *.zip)       unzip ./"$n"       ;;
            *.z)         uncompress ./"$n"  ;;
            *.7z|*.arj|*.cab|*.chm|*.deb|*.dmg|*.iso|*.lzh|*.msi|*.rpm|*.udf|*.wim|*.xar)
                         7z x ./"$n"        ;;
            *.xz)        unxz ./"$n"        ;;
            *.exe)       cabextract ./"$n"  ;;
            *)
                         echo "extract: '$n' - unknown archive method"
                         return 1
                         ;;
          esac
      else
          echo "'$n' - file does not exist"
          return 1
      fi
    done
fi
}
```

### Commands for SLURM info

`sh_part` -- Tells you resource information on all of the partitions you have access to.

`squeue -u $USER` -- Shows your queue of jobs submitted to the cluster.

`squeue -o \"%.18i %.9P %.8j %.8u %.8g %.2c %.2t %.10M %.10l %.6D %R\" | grep -E \"JOBID|dpetrov\"` -- Shows the entire queue of jobs for the lab. Not super useful other than just being generally nosy.

### Other useful commands

Request an interactive compute node which will allow you to work on the command line with more resources. 
```
srun -c $cpus -p dpetrov,hns,normal --time $time --mem $gigs_ram --pty bash
```
