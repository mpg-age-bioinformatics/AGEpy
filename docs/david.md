# david

A command line tool to perform enrichment analysis of gene lists via the [DAVID][david] web API. Query gene lists from tabular files (xls[x] sheets, tsv files, etc.) for any database (e.g. GO terms, KEGG pathways) available from DAVID.
Returns an xlsx file per gene list, with one sheet for each database.

It is part of the [AGEpy][agepy] python package and installed as standalone-script.

## Installation

See the installation instructions from the [AGEpy][agepy] python package.

The script is installed to a path defined by the shell variable `PYTHONUSERBASE` for Unix systems. This might be e.g. `$HOME/.local` or `/usr/local`. Add the `bin/` folder within this directory to your shell variable `PATH` and export it, e.g. 

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then you can access `david` from the command line.

## Usage

The basic syntax is:

```
david table.xlsx
```

Where `table.xlsx` is a file name to a table, in xls[x], or separated value (e.g. tsv, csv, csv2) format.
If you provide multiple files, they all must be of the same type and the input options (e.g. delimiter,
sheet name, column index, see below) are similar for each file.

David query options:

* `-t LIST`: a comma separated list with data bases to perform enrichment on, e.g. `GOTERM_BP_FAT`
or `KEGG_PATHWAY`
* `-i ID`: the id type of the gene list, e.g. `ENSEMBL_GENE_ID`. Check the tables of 
[available ID types](https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html).
* `-u USER`: the DAVID web API user name (email address), which needs to be registered at
the [DAVID webservice][david_register] beforehand.
* `-p P`: maximum p-value to select significant enriched terms
* `-n N`: minimum number of genes within a term

Input options:

* `-d DELIMITER`: define a delimiter for text separated input files
* `-c N`: define a single column (0-based) index, from which the gene list is to be taken.
Default is to iterate over all columns.
* `-s SHEET`: define a single sheet from an xls[x] file to extract gene lists from.
Default is to iterate over all sheets.
* `-g GTF_TABLE`: Provide a gtf file (optimised for ENSEMBL gtf files) if your gene list
contains gene symbols (gtf attribute `gene_name`) and they should be converted to ENSEMBL
gene ids (attribute `gene_id`) before querying.

Output options:

* `-o OUTPUT_FOLDER`: determine the output folder. Output files are named as base name
of the input tables + sheet name + column name + `.xlsx`.

Other:

* Run `david --help` to see the full list and explanation of options.

Input data format:

* tables need a column header
* all columns are used by default
* column indexing is 0-based, so `david -c 0` takes the 1st, `david -c 1` the 2nd column, and so on.
* all sheets are used by default (if input is in xls[x] format)
* sheet indexing is based sheet names (spaces in names need a quote)

Single column example from a text file.  
e.g. `david -c 0 -g my.gtf table.tsv`

```
Genes Values Description
ABC   1.1    Nothing
XYZ   3.6    important
FOO  -0.2    for david
...   ...    ...
BAR   0.0    in this column
```

Multiple column example from the sheet named 'Gene List' in an excel file.  
e.g. `david -s "Gene List" -g my.gtf table.xlsx`

```
Wildtype Mutant_1 Mutant_2
ABC      ABC      FOO
XYZ      EFG      HIJ
FOO      BAR      KLM
BAR               NOP
                  QRS
                  XYZ
```





[agepy]: https://github.com/mpg-age-bioinformatics/AGEpy
[david]: https://david.ncifcrf.gov/
[david_register]: https://david.ncifcrf.gov/webservice/register.htm
