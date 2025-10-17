Working on Sherlock means that you will be using bash commands. You will have to use them to move around the filesystems, to look inside your directories, and to move files around. The bash basics can be pretty quickly learned from any sort of bash tutorial, and if you don't feel super comfortable moving through Sherlock, I would recommend looking at some.

The purpose of this directory is less to teach you how to do anything and more to show the kinds of things you can do with bash code. Bash code is incredibly powerful for file manipulation, whether its converted files from one type to another, reading data from files, or writing data to files. As I have become more comfortable with bash, I've found that I write a lot less Python or R code with the intent of manipulating files. With bash you can generally write much simpler code that runs much faster.

Powerhouse commands for bash include `grep`, `cut`, `cat`, `echo`,`awk`, `sed`, `sort`. Maybe the most useful command of all is `man`, which will bring up the documentation (manual) for any bash command (`man cat`). 

`grep` -- pattern matching in files. Super useful for finding specific text in files. Let's say you want to know how many sequences are in a fasta file. Every unique ID in a fasta file starts with `>`. We can simply count every line that has a `>`.

```
# -c will return the number of lines that have a match.
# if we omit -c, then it will return the lines themselves.
grep -c ">" test.fasta
```

`cut` -- get columns in files based on a delimiter. Lets say we have a bed file, which records genomic locations of annotations, and this bed file is storing gene locations. If we want to get the names of all the genes in our file, we can use `cut`.

Our bed file has the format `[chrom]\t[start]\t[end]\t[ID]`:
```
2L	2653045	2653926	FBgn0031442
2L	2977121	2978500	FBgn0031484
2L	3784946	3785622	FBgn0051960
2L	16799024	16801584	FBgn0040260
2R	11341713	11342760	FBgn0011555
```
We can just cut out the fourth column (based on tabs as the delimiter):
```
> cut -f4 genes.bed
FBgn0031442
FBgn0031484
FBgn0051960
FBgn0040260
FBgn0011555
```
If we had a .csv file, we can split on commas with `cut -f4 -d ',' test.csv`

`cat` -- concatenate files and print. This is most commonly used to print out a file's contents to the terminal: `cat test.txt`. If you list multiple files, it will print them all out one after the other. You can then redirect the output to a new file. `cat chr1.fasta chr2.fasta > genome.fasta`

`echo` -- display text. Can be useful if you have bash variables and need to figure out what is in them: `echo $SCRATCH`

