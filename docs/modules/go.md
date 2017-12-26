## ___getGeneAssociation___

This function collects GO annotation from http://geneontology.org/page/download-annotations.

* **`URL_or_file`** either a link to a file on geneontology.org eg. http://geneontology.org/gene-associations/gene_association.fb.gz or the path for the respective  downloded .gz file.

### *`getGeneAssociation(URL_or_file)`*

* **`returns`** a Pandas dataframe with the parsed table.

```python
>>> import pandas as pd
>>> gA=age.getGeneAssociation("http://geneontology.org/gene-associations/gene_association.wb.gz")
>>> print gA.head()

DB    DB_Object_ID DB_Object_Symbol Qualifier       GO ID  \
0  WB  WBGene00000001            aap-1            GO:0005942   
1  WB  WBGene00000001            aap-1            GO:0005942   
2  WB  WBGene00000001            aap-1            GO:0008286   
3  WB  WBGene00000001            aap-1            GO:0008286   
4  WB  WBGene00000001            aap-1            GO:0008286   

                        DB:Reference Evidence      With (or) From Aspect  \
0                        GO_REF:0000002      IEA  InterPro:IPR001720      C   
1  WB_REF:WBPaper00005614|PMID:12393910      IDA                          C   
2  WB_REF:WBPaper00005614|PMID:12393910      IGI   WB:WBGene00000090      P   
3  WB_REF:WBPaper00005614|PMID:12393910      IGI   WB:WBGene00000898      P   
4  WB_REF:WBPaper00005614|PMID:12393910      IMP                          P   

DB_Object_Name DB_Object_Synonym DB_Object_Type       Taxon      Date  \
0                       Y110A7A.10           gene  taxon:6239  20170321   
1                       Y110A7A.10           gene  taxon:6239  20151214   
2                       Y110A7A.10           gene  taxon:6239  20151214   
3                       Y110A7A.10           gene  taxon:6239  20151214   
4                       Y110A7A.10           gene  taxon:6239  20060302   

Assigned_by Annotation Extension Gene Product Form ID  
0          WB                                            
1          WB                                            
2          WB                                            
3          WB                                            
4          WB                                            
```
___
