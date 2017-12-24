## ___organismsKEGG___

Lists all organisms present in the KEGG database.

* **`returns`** a dataframe containing one organism per row.

```python
>>> import AGEpy as age
>>> KEGGo=age.organismsKEGG()
>>> print KEGGo.head()

0    1                                                  2  \
0  T01001  hsa                               Homo sapiens (human)   
1  T01005  ptr                       Pan troglodytes (chimpanzee)   
2  T02283  pps                              Pan paniscus (bonobo)   
3  T02442  ggo  Gorilla gorilla gorilla (western lowland gorilla)   
4  T01416  pon                  Pongo abelii (Sumatran orangutan)   

                                 3  
0  Eukaryotes;Animals;Vertebrates;Mammals  
1  Eukaryotes;Animals;Vertebrates;Mammals  
2  Eukaryotes;Animals;Vertebrates;Mammals  
3  Eukaryotes;Animals;Vertebrates;Mammals  
4  Eukaryotes;Animals;Vertebrates;Mammals
```
___

## ___databasesKEGG___

Finds KEGG database identifiers for a respective organism given example ensembl ids.

* **`organism`** an organism as listed in organismsKEGG()
* **`ens_ids`** a list of ensenbl ids of the respective organism
* **`returns`** nothing if no database was found, or a string if a database was found

```python
>>> import AGEpy as age
>>> print sigGenes[:10]

['ENSG00000272449', 'ENSG00000130762', 'ENSG00000083444', 'ENSG00000162493', 'ENSG00000253368', 'ENSG00000235912', 'ENSG00000169174', 'ENSG00000240563', 'ENSG00000200174', 'ENSG00000162607']

>>> KEGGd=age.databasesKEGG("hsa",sigGenes)

For hsa the following db was found: Ensembl

>>> print KEGGd.head()

Ensembl
```
___

## ***ensembl_to_kegg***

Looks up KEGG mappings of KEGG ids to ensembl ids.

* **`organism`** an organisms as listed in organismsKEGG()
* **`kegg_db`** a matching KEGG db as reported in databasesKEGG
* **`returns`** a Pandas dataframe of with 'KEGGid' and 'ENSid'.

```python
>>> import AGEpy as age
>>> KEGGi=age.ensembl_to_kegg("hsa","Ensembl")

KEGG API: http://rest.genome.jp/link/Ensembl/hsa

>>> print KEGGi.head()

KEGGid            ENSid
0      hsa:1  ENSG00000121410
1     hsa:10  ENSG00000156006
2    hsa:100  ENSG00000196839
3   hsa:1000  ENSG00000170558
4  hsa:10000  ENSG00000117020
```
___

## ***ecs_idsKEGG***

Uses KEGG to retrieve all ids and respective ecs for a given KEGG organism.

* **`organism`** an organisms as listed in organismsKEGG()
* **`returns`** a Pandas dataframe of with 'ec' and 'KEGGid'.

```python
>>> import AGEpy as age
>>> KEGGie=age.ecs_idsKEGG("hsa")
>>> print KEGGie.head()

ec    KEGGid
0   ec:2.7.11.1  hsa:9344
1   ec:2.7.11.1  hsa:5894
2   ec:2.7.11.1   hsa:673
3   ec:2.7.12.2  hsa:5607
4  ec:2.7.11.24  hsa:5598
```
___

## ***idsKEGG***

Uses KEGG to retrieve all ids for a given KEGG organism.

* **`organism`** an organism as listed in organismsKEGG()
* **`returns`** a Pandas dataframe of with 'gene_name' and 'KEGGid'.

```python
>>> import AGEpy as age
>>> KEGGids = age.idsKEGG("hsa")
>>> print KEGGids.head()

gene_name         KEGGid
0      uncharacterized LOC100287010  hsa:100287010
1      uncharacterized LOC100288846  hsa:100288846
2                      DKFZp434L192     hsa:222029
3  uncharacterized protein FLJ30679     hsa:146512
4      uncharacterized LOC100128288  hsa:100128288
```
___

## ***pathwaysKEGG***

Retrieves all pathways for a given organism.

* **`organism`** an organism as listed in organismsKEGG()
* **`returns df`** a Pandas dataframe with the columns 'KEGGid','pathIDs', and 'pathName'.
* **`df_`** a Pandas dataframe with a columns for 'KEGGid', and one column for each pathway with the corresponding gene ids below

```python
>>> import AGEpy as age
>>> KEGGp, KEGGp_ =age.pathwaysKEGG("hsa")
>>> print KEGGp.head()

>>> print KEGGp_.head()


```
___
