# mapgtf

## About
A tool to annotate genome maps by gtf feature attributes.

## Usage

Install the `AGEpy` python package, which provides this tool.

Run `mapgtf --help` from the command line to show the usage:

```
usage: mapgtf [-h] [-f FEATURE] [-s SOURCE] [-a ATTRIBUTE [ATTRIBUTE ...]]
              [-A OTHER [OTHER ...]] [-k] [-o OUTPUT] [-c MAP_CHR_COL]
              [-p MAP_POS_COL] [-e MAP_END_COL] [-n NCPUS] [-v]
              GTFFILE MAPFILE

positional arguments:
  GTFFILE               annotation gtf file (gff version 2); see
                        http://www.ensembl.org/info/website/upload/gff.html
  MAPFILE               genome map file to convert; VCF files need uncommented
                        header line!

optional arguments:
  -h, --help            show this help message and exit
  -f FEATURE            optionally select feature from gtf file, e.g. exon
                        (default: None)
  -s SOURCE             optionally select source to subset gtf file, e.g.
                        ensembl (default: None)
  -a ATTRIBUTE [ATTRIBUTE ...]
                        select attribute tags to include as annotation columns
                        (default: ['gene_id', 'gene_name', 'gene_biotype',
                        'transcript_id', 'transcript_name',
                        'transcript_biotype'])
  -A OTHER [OTHER ...]  select additional columns to be added (default:
                        ['feature'])
  -k                    keep non-intersecting lines from file (default: False)
  -o OUTPUT             output file name (default: annotated.map)
  -c MAP_CHR_COL        map column index [1-based] for chromosome ids
                        (default: 1)
  -p MAP_POS_COL        map column index [1-based] for positions (default: 2)
  -e MAP_END_COL        map column index [1-based] for end positions, if
                        MAPFILE contains ranges (default: None)
  -n NCPUS              select number of processes for parallel computing
                        (default: 1)
  -v                    be more verbose (default: False)
```

## Input format

`GTFFILE` is a gtf (gff version 2) genome annotation formatted file.
See the specification at:
http://www.ensembl.org/info/website/upload/gff.html

`MAPFILE` is a **tab separated** text file.
It must contain a **header line**,
and at least a chromosome ID and position column.
In addition, another position column for range ends can be defined.

See the examples:

### Map file with single positions

*File content*

```
start   chr stats   notes
1       I   0.12    NA
10      I   0.44    NA
5       II  0.12    NA
12      II  0.01    important
10      III 0.59    NA
240     III 0.81    NA
```

*Command*

```bash
mapgtf -c 2 -p 1 /path/to/gtf /path/to/map
```

### Map file with regions (similar to bed file formats)

*File content*

```
chr start   end stats
I   1       5   0.12
I   10      100 0.44
II  5       9   0.12
II  12      17  0.01
III 10      190 0.59
III 240     500 0.81
```

*Command*

```bash
mapgtf -e 3 /path/to/gtf /path/to/map
```

### VCF (variant calling format) file

*File content*

See the specification at:
https://samtools.github.io/hts-specs/VCFv4.2.pdf

Since the header line in a vcf file starts with a comment character, it needs to be removed.

*Command*

```bash
sed 's/^#CHROM/CHROM/' /path/to/vcf > /path/to/vcf.mod
mapgtf /path/to/gtf /path/to/vcf.mod
```

