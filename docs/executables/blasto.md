## Intro

This module will load a fasta formatted file and query each fasta sequence for blast.
The user may add blast parameters as space separated list after the sequence name. All queries are
listed into a log table. The user can either let the program running while waiting for the results
using the -C option, or quit and check if the results are ready later using -W -t <queryTable.tsv>

## Examples

```bash
$cat input.fa

>sequence1
GCGAAGCCCAAGAGGATGAAGCCAGAGATGGTGTTGGAGTTGCTGGGGCTGCTGAGGGTATTGATCTGTCTGTGACCTGCGATAGCATCAGAAGTTGTTTCACATTCTAGTTATAGCTGAGGGAGGTTATGTTTTGAGCAAGCAGGAAAC
>Sequence2
AGCTCCTGAGAAACTTGGGGGGCGCGACACAGATAGGGTGAAAGCAGAGTGATAGACCTGGGATGGTTACGGGACCAAGGGAAGACCAGGCTGGTTGGCATACACCGGTGAACGGATGGGAGTCCTAGGGAAAGATGATGCGCCTAACAG
>sequence2_filtered database='nt' filter="T" nucl_penalty=-5 gapcosts='1,11'
AGCTCCTGAGAAACTTGGGGGGCGCGACACAGATAGGGTGAAAGCAGAGTGATAGACCTGGGATGGTTACGGGACCAAGGGAAGACCAGGCTGGTTGGCATACACCGGTGAACGGATGGGAGTCCTAGGGAAAGATGATGCGCCTAACAG
>sequence3
TCGTTTGATTCTGCAAGCAGCACCTACTGTGGGGTATTGATAAGATCTCTGATGGCGTCTGAAATTCTTCTGAGATTAGAGGAAGATCAGGTGTGTTTTAATGTCGAGCAGGTGTTTCCCCAAGATTAGTGGGGGGATTCGGTTTTTCCT

$blasto -S -f /usr/home/JDoe/project1/input.fa -o /usr/home/JDoe/project1/run1
$blasto -W -t /usr/home/JDoe/project1/run1.queryTable.tsv -o /usr/home/JDoe/project1/run1
```

## Help
```bash

usage: blasto [-h] [-S] [-C] [-W] [-f INPUTFASTA] [-t INPUTTSV]
              [-o OUTPUTPREFIX] [--format_type FORMAT_TYPE]
              [--sleepTime SLEEPTIME] [--description]

This module will load a fasta formatted file and query each fasta sequence for
blast The user may add blast parameters as space separated list after the
sequence name. All queries are listed into a log table. The user can either
let the program running while waiting for the results using the -C option, or
quit and check if the results are ready later using -W -t <queryTable.tsv>

optional arguments:
  -h, --help            show this help message and exit
  -S, --submitFromFasta
                        Read in fasta file and submit blast queries. Write out
                        submitted query IDs. (default: False)
  -C, --continueThrough
                        Read from fasta file, submit and continue checking.
                        Write results when they are ready and exit after all
                        results are finished. (default: False)
  -W, --checkAndWriteResults
                        Read query IDs from tsv and check status. If results
                        are ready, collect and safe. (default: False)
  -f INPUTFASTA, --inputFasta INPUTFASTA
                        Fasta formatted input file containing one or more
                        input sequences. The sequence name may contain
                        additional blast paramers, (default: )
  -t INPUTTSV, --inputTsv INPUTTSV
                        Tab separated input file containing sequence IDs,
                        output prefix, query IDs, query arguments. (default: )
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Output prefix. All files will start with this prefix,
                        blast output files will be written two
                        <prefix>_<sequenceID>.<format_type> (default: )
  --format_type FORMAT_TYPE
                        format of the blast output (default: Tabular)
  --sleepTime SLEEPTIME
                        time to wait before checking again if your jobs are
                        done, only active if -C is on (default: 60)
  --description         Get a description of what this script does. (default:
                        False)
```
