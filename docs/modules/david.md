## ___DAVIDenrich___

Queries the DAVID database for an enrichment analysis.
Check https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html for database == "type" tag and categories ==  "annot" tag.

**`DAVIDenrich(database, categories, user, ids, ids_bg = None, name = '', name_bg = '', verbose = False, p = 0.1, n = 2)`**

* **`database`** A string for the database to query, e.g. 'WORMBASE_GENE_ID'
* **`categories`** A comma separated string with databases
* **`user`**  A user ID registered at DAVID for querying
* **`ids`**  A list with identifiers
* **`name`**  A string with the name for the query set
* **`ids_bg`**  A list with the background identifiers to enrich against, 'None' for whole set
* **`name_bg`**  A string with the name for the background set
* **`p`**  Maximum p value for enrichment of a term
* **`n`**  Minimum number of genes within a term
* **`returns`**  None if no ids match the queried database, or a pandas dataframe with results

```python
>>> import AGEpy as age
>>> print sigGenes[:10]

[u'WBGene00022275', u'WBGene00004418', u'WBGene00018774',
 u'WBGene00018772', u'WBGene00018958', u'WBGene00021662',
 u'WBGene00255594', u'WBGene00021658', u'WBGene00021026',
 u'WBGene00022042']

>>> categories=['GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', 'KEGG_PATHWAY','BIOCARTA', 'PFAM', 'PROSITE' ]
>>> DAVIDdf=age.DAVIDenrich('WORMBASE_GENE_ID', categories, 'email.registered@david.com', sigGenes)
>>> print DAVIDdf.head()

categoryName                                     termName listHits  \
0  GOTERM_BP_FAT                       GO:0006412~translation      177   
1  GOTERM_BP_FAT         GO:0006518~peptide metabolic process      198   
2  GOTERM_BP_FAT      GO:0043043~peptide biosynthetic process      177   
3  GOTERM_BP_FAT        GO:0043604~amide biosynthetic process      180   
4  GOTERM_BP_FAT  GO:0043603~cellular amide metabolic process      206   

     percent               ease  \
0  5.85704831238  4.32627669357e-43   
1  6.55195234944  1.36601477909e-42   
2  5.85704831238  4.04090150003e-42   
3  5.95632031767  1.05565138148e-40   
4  6.81667769689  3.74871147863e-40   

                                         geneIds listTotals popHits  \
0  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     379   
1  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     455   
2  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     384   
3  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     402   
4  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     499   

popTotals foldEnrichment         bonferroni          benjamini  \
0     11221  2.79042292227  1.28576943333e-39  1.28576943333e-39   
1     11221  2.60009830425  4.05979592345e-39  2.02989796172e-39   
2     11221  2.75408929047  1.20095592581e-38  4.00318641936e-39   
3     11221   2.6753612131  3.13739590576e-37  7.84348976441e-38   
4     11221  2.46662227543  1.11411705145e-36   2.2282341029e-37   

            afdr  
0  7.78207551683e-40  
1  2.45717759656e-39  
2  7.26874466353e-39  
3  1.89889814083e-37  
4   6.7431553467e-37  
```
___

## ___DAVIDgetGeneAttribute___

Returns a list of gene names for given gene ids.

**`DAVIDgetGeneAttribute(x, df, refCol="ensembl_gene_id", fieldTOretrieve="gene_name")`**

* **`x`** a string with the list of IDs separated by ', '
* **`df`**  a dataframe with the reference column and a the column to retrieve
* **`refCol`** the header of the column containing the identifiers
* **`fieldTOretrieve`** the field to retrieve from parsedGTF eg. 'gene_name'
* **`returns`** list of fieldTOretrieve separeted by ', ' in the same order as the given in x

```python
>>> import AGEpy as age
>>> print df.head()

ensembl_gene_id            gene            locus sample_1 sample_2 status  \
0  WBGene00022275        Y74C9A.1    I:43732-44677       N2     daf2     OK   
1  WBGene00004418  F53G12.9,rpl-7  I:111037-113672       N2     daf2     OK   
2  WBGene00018774  F53G12.9,rpl-7  I:111037-113672       N2     daf2     OK   
3  WBGene00018772        F53G12.4  I:134336-137282       N2     daf2     OK   
4  WBGene00018958        F56C11.6  I:171339-175991       N2     daf2     OK   

     value_1      value_2  log2(fold_change)  test_stat  p_value   q_value  \
0     0.195901     0.986634           2.332390    2.32959  0.00570  0.031216   
1  3354.820000  2463.480000          -0.445539   -2.71381  0.00005  0.000556   
2  3354.820000  2463.480000          -0.445539   -2.71381  0.00005  0.000556   
3     1.235670     2.992460           1.276040    3.16508  0.00005  0.000556   
4     2.651180     3.795600           0.517696    1.73994  0.00410  0.024157   

significant                                              GO_id  \
0         yes                                                NaN   
1         yes  GO:0003735; GO:0000463; GO:0044822; GO:0002181...   
2         yes                                                NaN   
3         yes                                                NaN   
4         yes                 GO:0016787; GO:0005615; GO:0004104   

                                           GO_term    gene_biotype  \
0                                                NaN  protein_coding   
1  structural constituent of ribosome; maturation...  protein_coding   
2                                                NaN  protein_coding   
3                                                NaN  protein_coding   
4  hydrolase activity; extracellular space; choli...  protein_coding   

  NormInt evidence  
0 -0.356904       no  
1  3.458609       no  
2  3.458609       no  
3  0.283965       no  
4  0.501360       no  

>>> print DAVIDdf.head()

categoryName                                     termName listHits  \
0  GOTERM_BP_FAT                       GO:0006412~translation      177   
1  GOTERM_BP_FAT         GO:0006518~peptide metabolic process      198   
2  GOTERM_BP_FAT      GO:0043043~peptide biosynthetic process      177   
3  GOTERM_BP_FAT        GO:0043604~amide biosynthetic process      180   
4  GOTERM_BP_FAT  GO:0043603~cellular amide metabolic process      206   

     percent               ease  \
0  5.85704831238  4.32627669357e-43   
1  6.55195234944  1.36601477909e-42   
2  5.85704831238  4.04090150003e-42   
3  5.95632031767  1.05565138148e-40   
4  6.81667769689  3.74871147863e-40   

                                         geneIds listTotals popHits  \
0  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     379   
1  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     455   
2  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     384   
3  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     402   
4  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     499   

popTotals foldEnrichment         bonferroni          benjamini  \
0     11221  2.79042292227  1.28576943333e-39  1.28576943333e-39   
1     11221  2.60009830425  4.05979592345e-39  2.02989796172e-39   
2     11221  2.75408929047  1.20095592581e-38  4.00318641936e-39   
3     11221   2.6753612131  3.13739590576e-37  7.84348976441e-38   
4     11221  2.46662227543  1.11411705145e-36   2.2282341029e-37   

            afdr  
0  7.78207551683e-40  
1  2.45717759656e-39  
2  7.26874466353e-39  
3  1.89889814083e-37  
4   6.7431553467e-37  

>>> gene_names=df[["ensembl_gene_id","gene"]].drop_duplicates()
>>> DAVIDdf["gene_names"]=DAVIDdf["geneIds"].apply(lambda x: \
                              age.DAVIDgetGeneAttribute(x,\
                              gene_names,\
                              refCol="ensembl_gene_id",\
                              fieldTOretrieve="gene"))
>>> print DAVIDdf.head()

categoryName                                     termName listHits  \
0  GOTERM_BP_FAT                       GO:0006412~translation      177   
1  GOTERM_BP_FAT         GO:0006518~peptide metabolic process      198   
2  GOTERM_BP_FAT      GO:0043043~peptide biosynthetic process      177   
3  GOTERM_BP_FAT        GO:0043604~amide biosynthetic process      180   
4  GOTERM_BP_FAT  GO:0043603~cellular amide metabolic process      206   

     percent               ease  \
0  5.85704831238  4.32627669357e-43   
1  6.55195234944  1.36601477909e-42   
2  5.85704831238  4.04090150003e-42   
3  5.95632031767  1.05565138148e-40   
4  6.81667769689  3.74871147863e-40   

                                         geneIds listTotals popHits  \
0  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     379   
1  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     455   
2  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     384   
3  WBGENE00002063, WBGENE00013678, WBGENE00006626...       1878     402   
4  WBGENE00002063, WBGENE00006626, WBGENE00007584...       1878     499   

popTotals foldEnrichment         bonferroni          benjamini  \
0     11221  2.79042292227  1.28576943333e-39  1.28576943333e-39   
1     11221  2.60009830425  4.05979592345e-39  2.02989796172e-39   
2     11221  2.75408929047  1.20095592581e-38  4.00318641936e-39   
3     11221   2.6753612131  3.13739590576e-37  7.84348976441e-38   
4     11221  2.46662227543  1.11411705145e-36   2.2282341029e-37   

            afdr                                         gene_names  
0  7.78207551683e-40  ife-5, Y105E8A.20, tsn-1, yars-1, ife-3, C14C1...  
1  2.45717759656e-39  ife-5, tsn-1, C14C10.1, ife-3, rps-30, iff-2, ...  
2  7.26874466353e-39  ife-5, Y105E8A.20, tsn-1, yars-1, ife-3, C14C1...  
3  1.89889814083e-37  ife-5, Y105E8A.20, tsn-1, yars-1, ife-3, C14C1...  
4   6.7431553467e-37  ife-5, tsn-1, C14C10.1, ife-3, rps-30, Y51H4A....  
```
___

## ***DAVIDplot***

Queries the DAVID database for an enrichment analysis and plots CellPlots as
well as SymPlots (see plots) using the 20 most significant terms.
Check https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html for database == "type" tag and categories ==  "annot" tag.

**`DAVIDplot(database, categories, user, df_ids, output, df_ids_bg = None, name = '', name_bg = '', verbose = False, p = 0.1, n = 2)`**

* **`database`** a string for the database to query, e.g. 'WORMBASE_GENE_ID'
* **`categories`** a comma separated string with databases
* **`user`** a user ID registered at DAVID for querying
* **`df_ids`** a dataframe where the first column contains the identifiers
    to be queried and the second column the respective log2fc for each identifier.
* **`output`** /path/to/output/prefix
* **`df_ids_bg`** a dataframe where the first column contains the identifiers to be used as background. 'None' for whole set
* **`name`** a string with the name for the query set
* **`name_bg`** a string with the name for the background set
* **`p`** Maximum p value for enrichment of a term
* **`n`** Minimum number of genes within a term

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> print df.head()

ensembl_gene_id  log2(fold_change)
0  ENSG00000272449           1.859500
1  ENSG00000130762           0.601051
2  ENSG00000083444          -0.881957
3  ENSG00000162493          -0.638433
4  ENSG00000253368           0.654517

>>> categories=['GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', 'KEGG_PATHWAY','BIOCARTA', 'PFAM', 'PROSITE' ]
>>> DAVIDdf=DAVIDplot('ENSEMBL_GENE_ID', categories, 'email.registered@david.com', df, "/usr/home/JDoe/mydataset")
```
___
