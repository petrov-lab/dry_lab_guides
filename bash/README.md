Working on Sherlock means that you will be using bash commands. You will have to use them to move around the filesystems, to look inside your directories, and to move files around. The bash basics can be pretty quickly learned from any sort of bash tutorial, and if you don't feel super comfortable moving through Sherlock, I would recommend looking at some.

The purpose of this directory is less to teach you how to do anything and more to show the kinds of things you can do with bash code. Bash code is incredibly powerful for file manipulation, whether its converted files from one type to another, reading data from files, or writing data to files. As I have become more comfortable with bash, I've found that I write a lot less Python or R code with the intent of manipulating files. With bash you can generally write much simpler code that runs much faster.

Powerhouse commands for bash include `grep`, `cut`, `cat`, `echo`, `sort`,`awk`, `sed`. Maybe the most useful command of all is `man`, which will bring up the documentation (manual) for any bash command (`man cat`). Most bash commands have a whole host of flags that will modify the function of the command in a potentially useful way. The man page will list all of these flags. It's also very helpful if you're reading someone else's code and you are trying to figure out what a command does when there are all these extra flags tacked on (`sort` vs `sort -t'_' -k1,1n -k2,2hr`).

`grep` -- pattern matching in files. Super useful for finding specific text in files. Let's say you want to know how many sequences are in a fasta file. Every unique ID in a fasta file starts with `>`. We can simply count every line that has a `>`.

`cut` -- get columns in files based on a delimiter. Lets say we have a bed file, which records genomic locations of annotations, and this bed file is storing gene locations. If we want to get the names of all the genes in our file, we can use `cut`.

`cat` -- concatenate files and print. This is most commonly used to print out a file's contents to the terminal: `cat test.txt`. If you list multiple files, it will print them all out one after the other. You can then redirect the output to a new file. `cat chr1.fasta chr2.fasta > genome.fasta`

`echo` -- display text. Can be useful if you have bash variables and need to figure out what is in them: `echo $SCRATCH`

`sort` -- sort files by a column. Very useful command for ordering values. Can sort numbers and text. Can also be used to get unique values with `sort -u`. 

`awk` -- This probably needs its own page. Awk is a scripting language used for text processing. If you need to do complex rearrangements of columns and/or had to change delimiters between files, awk is a good place to start.

`sed` -- This probably could also have its own page. sed is a utility in bash that edits text. It is exceptional at pattern matching (think find and replace). `sed` allows you to delete or replace specific patterns in files (among a bunch of other things).

There are many more bash commands than these few that I listed here. I find that I use a combination of these most commonly, however when doing more specific tasks it can be worthwhile googling similar problems to see if a command already exists to solve your issue.

Pipes (`|`) and input/output redirects (`<` or `>` or `>>`) are also powerful tools in bash scripting. Pipes allow you to chain commands together. Redirects allow you to save to or read from files, as well as chain commands together.
```
# Find the unique values in the fourth column of a comma-separated file, then save to a new file
cut -f4 -d ',' test.csv | sort -u > col4.txt
# First, we read in test.csv and we define columns by using the comma "," as a delimiter. We keep the fourth column and print it to stdout (standard output; the terminal)
# The pipe "|" captures the output and instead of printing it, it feeds the fourth column as an input to our next command,
# which first sorts the input and removes all duplicate values, and then prints to stdout
# The output redirect ">" captures the output and instead of printing it, writes it to a file named col4.txt
# Note that if col4.txt already exists, it will be overwritten. If we used ">>" instead, then the output would
# be appended to an existing col4.txt.
```
