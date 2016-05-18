# david

*david* - A command line tool to perform enrichment analysis of gene lists
via the [DAVID][david] web API. Query gene lists from tabular files 
(xls[x] sheets, tsv files, csv files, etc.) for any database 
(e.g. GO terms, KEGG pathways) available from DAVID.

It is part of the [AGEpy][agepy] python package and installed as 
standalone-script.

## Installation

See the installation instructions from the [AGEpy][agepy] python package.

The script is installed to a path defined by the shell variable `PYTHONUSERBASE`
for Unix systems. This might be e.g. `$HOME/.local` or `/usr/local`. 
Add the `bin/` folder within this directory to your shell variable `PATH` and
export it, e.g. 

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then you can access *david* from the command line.

## Usage

### Options summary

To invoke the usage description of the script, call it with the help option as
in the following:

```bash
david --help
```

This shows you the following description:

```
usage: david [-h] [-d CHAR] [-c N] [-e N] [-s NAME [NAME ...]] [-g GTF]
             [-t CATEGORY[,CATEGORY,...]] [-i ID] [-u EMAIL] [-p P] [-n N]
             [-o FOLDER] [-x] [-v]
             TABLE [TABLE ...]

david - perform DAVID enrichment analysis of gene sets

positional arguments:
  TABLE                 Input file name[s], format is text or xls[x], same for
                        each file and auto-detected.

optional arguments:
  -h, --help            show this help message and exit
  -d CHAR, --delimiter CHAR
                        For text input files, a character used as column
                        delimiter. E.g. '\t' for tsv, ',' for csv, ';' for
                        csv2, etc. (default: \t)
  -c N, --column N      Select a column index (0-based) with ids to query.
                        None means all. (default: None)
  -e N, --expression-column N
                        Select a column index (0-based) with gene expression
                        values (default: None)
  -s NAME [NAME ...], --sheet NAME [NAME ...]
                        If input is an xls[x] file, select sheets by name.
                        None means all. (default: None)
  -g GTF, --gtf GTF     A gtf file to convert the ids in the input tables from
                        gene_name to gene_id (as in ENSEMBL gtf files). If no
                        file is provided, no conversion is performed.
                        (default: None)
  -t CATEGORY[,CATEGORY,...], --DAVIDtypes CATEGORY[,CATEGORY,...]
                        DAVID category types to enrich for. Comma separated
                        string. (default: GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_M
                        F_FAT,KEGG_PATHWAY,BIOCARTA,PFAM,PROSITE)
  -i ID, --DAVIDid ID   DAVID data format id. (default: ENSEMBL_GENE_ID)
  -u EMAIL, --DAVIDuser EMAIL
                        DAVID user id (email address). (default: None)
  -p P, --max-pvalue P  Maximum p-value for enrichment (default: 0.05)
  -n N, --min-genes N   Minimum number of genes in term (default: 2)
  -o FOLDER, --output-folder FOLDER
                        Output folder. (default: david_output)
  -x, --output-xlsx     Also return xlsx output files, with all categories
                        merged into a single workbook. (default: False)
  -v, --verbose         Be verbose. (default: False)

This program is part of the python package AGEpy. Author: Sven E. Templer
<sven.templer@gmail.com>. Copyright (c) 2016 - Bioinformatics Core Facility at
the Max Planck Institute for Biology of Ageing, Cologne, Germany. Please file
issues at https://github.com/mpg-age-bioinformatics/AGEpy/issues and find
documentation at http://agepy.readthedocs.io/en/latest/david/
```

### Input options

`TABLE` is a file path to your input table, for example a differential gene
expression result. The input file format (text or excel) is auto detected, and
must be identical for all supplied files. Input files can either be

* one or more xlsx workbook(s) with one or multiple sheets, or
* one or more tab or other delimited file(s) (e.g. csv).

Each table (text file, sheet of excel workbook) must be of the same format,
containing a header row with column names and can be of one of the following
format:

* Containing a single column with gene ids (option `-c`) and optionally a column
  with gene expression values (e.g. log fold changes, option `-e`), e.g.:

```
Genes log2FC  Description
ABC   1.1     Nothing
XYZ   3.6     important
FOO   -0.2    for david
...   ...     ...
BAR   0.0     in this column
```

* Containg one or multiple columns, each with a list of gene ids, e.g.:

```
Wildtype Mutant_1 Mutant_2
ABC      ABC      FOO
XYZ      EFG      HIJ
FOO      BAR      KLM
BAR               NOP
                  QRS
                  XYZ
```

`-d CHAR` is used to specify a character `CHAR` used as column 
delimiter, if the input format of `TABLE` is text.

`-c N` specifies a column index `N` (0-based), if only a single column of the
input tables contains gene identifiers to be queried. If this is not provided,
then all columns are used for the query.

`-e N` defines a column index `N` (0-based), if there exists a column with gene
expression values that should be attached as `, `-separated list to each term.

`-s NAME [NAME ...]` allows to select one (or more) sheets from excel input. By
default every sheet is used.

`-g GTF` can be used to specify a gtf (gff version 2, as provided by Ensembl)
file, which is used to convert gene symbols (attribute `gene_name`) to gene ids
(attribute `gene_id`). If supplied, both types are attached per term as 
`, `-separated list.

### Query options

`-t CATEGORY[,CATEGORY,...]` is used to select the categories on which to 
perform enrichment tests for the gene sets. It is a comma separated list.
Find all available categories at the ['annot' tag table][david_api] of the DAVID
web API description.

`-i ID` specifies the input format of gene identifiers. Available ID types can
be checked at the ['type' tag table][david_api] of the DAVID web API 
description.

`-u EMAIL` must be set to specify the user ID (email address), which needs to be
registered at the [DAVID webservice][david_register].

`-p P` defines the maximum p-value to select and return the terms by their
significance in the enrichment test.

`-n N` defnines the minimum number of genes that needs to be in a term to be
returned in the results.

### Output options

The output of the *david* tool is a single tab separated value (tsv) file per
input table, column and enrichment category. Optionally, a single excel workbook
can be returned containing each category. The general file naming scheme for tsv
result files is specified as follows:

```
<OUTPUT_FOLDER>/<INPUT_TABLE_BASE_NAME>-<SHEET_NAME>-<COLUMN_NAME>-<CATEGORY_NAME>.tsv
```

For excel result format, the name will be:

```
<OUTPUT_FOLDER>/<INPUT_TABLE_BASE_NAME>-<SHEET_NAME>-<COLUMN_NAME>.xlsx
```

`-o FOLDER` defines the output directory, and will be created if it does not
already exist.

`-x` enables the additional output in an excel workbook as described above.

`-v` can be used to create a more verbose printing of the tool to the standard
output on what is currently done.

## Legal

```
This program is part of the python package AGEpy.

Author: Sven E. Templer <sven.templer@gmail.com>

Copyright (c) 2016 - Bioinformatics Core Facility at the
                     Max Planck Institute for Biology of Ageing,
                     Cologne, Germany
```


[agepy]: https://github.com/mpg-age-bioinformatics/AGEpy
[david]: https://david.ncifcrf.gov/
[david_register]: https://david.ncifcrf.gov/webservice/register.htm
[david_api]: https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html
