## Intro

`obo2tsv` parses a gene ontology obo file to tsv. It will include for each term a column for parent terms as well as child terms.

## Examples

```
$ obo2tsv -u http://geneontology.org/ontology/go-basic.obo -o go-basic.tsv -c 4 --organism http://geneontology.org/gene-associations/gene_association.fb.gz
```

Links to other `--organism` can be found on [http://geneontology.org/page/download-annotations](http://geneontology.org/page/download-annotations).

## Help

```
$ obo2tsv --help

usage: obo2tsv [-h] [-i INPUT] [-u URL] [-o OUTPUT] [-c CPUS]
               [--organism ORGANISM]

obo to tsv parser

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        go-basic.obo file. Files can be downloaded from
                        http://geneontology.org/page/download-ontology.
                        (default: None)
  -u URL, --url URL     If no go-basic.obo input file is specified, a url to a
                        target obo file can be specified instead. (default:
                        http://geneontology.org/ontology/go-basic.obo)
  -o OUTPUT, --output OUTPUT
                        Name of output tab separated file. (default: go-
                        basic.tsv)
  -c CPUS, --cpus CPUS  Number of cpus. (default: 36)
  --organism ORGANISM   Optional, merge GO obo.tsv with a GO annotation for an
                        organism: either a link to a file on geneontology.org
                        eg. http://geneontology.org/gene-
                        associations/gene_association.fb.gz or the path for
                        the respective downloded .gz file. (default: None)
```
