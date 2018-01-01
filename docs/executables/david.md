## Intro

Queries the DAVID database for an enrichment analysis and plots CellPlots as well as SymPlots (see plots).

## Examples

```bash
$ cat input.tsv

ensembl_gene_id  log2(fold_change)
ENSG00000272449           1.859500
ENSG00000130762           0.601051
ENSG00000083444          -0.881957
ENSG00000162493          -0.638433
ENSG00000253368           0.654517

$ david -i input.tsv -o /usr/home/JDoe/project1/datasetA -d ENSEMBL_GENE_ID -u 'email.registered@david.com'
```

## Help

```bash
$ david --help

usage: david [-h] [-i INPUT] [-o OUTPUT] [-d DATABASE] [-c CATEGORIES]
             [-u USER] [-v] [-p PVALUE] [-n NGENES] [-b BACKGROUND]

Queries the DAVID database for an enrichment analysis and plots CellPlots as
well as SymPlots (see plots). Check
https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html for database ==
'type' tag and categories == 'annot' tag.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        A file with tab separated values where the first
                        column contains the identifiers to be queried and the
                        second column the respective log2fc for each
                        identifier. (default: None)
  -o OUTPUT, --output OUTPUT
                        /path/to/output/prefix (default: None)
  -d DATABASE, --database DATABASE
                        a string for the database to query, e.g.
                        'WORMBASE_GENE_ID'. (default: None)
  -c CATEGORIES, --categories CATEGORIES
                        a comma separated list of categories. (default: GOTERM
                        _BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,KEGG_PATHWAY,BIOCA
                        RTA,PFAM,PROSITE)
  -u USER, --user USER  a user ID registered at DAVID for querying (default:
                        None)
  -v, --verbose         Print more. (default: None)
  -p PVALUE, --pvalue PVALUE
                        Maximum p value for enrichment of a term. (default:
                        0.1)
  -n NGENES, --ngenes NGENES
                        Minimum number of genes within a term. (default: 2)
  -b BACKGROUND, --background BACKGROUND
                        A file with tab separated values where the first
                        column contains the identifiers to used as a
                        background. None for whole DAVID database as
                        background. (default: None)
```
