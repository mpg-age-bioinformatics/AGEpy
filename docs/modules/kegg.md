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
* **`returns df_`** a Pandas dataframe with a columns for 'KEGGid', and one column for each pathway with the corresponding gene ids below

```python
>>> import AGEpy as age
>>> KEGGp, KEGGp_ =age.pathwaysKEGG("hsa")
>>> print KEGGp.head()

KEGGid                                             pathID  \
0     hsa:10  path:hsa00232, path:hsa01100, path:hsa00983, p...   
1    hsa:100        path:hsa05340, path:hsa01100, path:hsa00230   
2   hsa:1000                       path:hsa04514, path:hsa05412   
3  hsa:10000  path:hsa04211, path:hsa04630, path:hsa04152, p...   
4  hsa:10005        path:hsa01100, path:hsa00120, path:hsa04146   

                                      pathName  
0  Caffeine metabolism - Homo sapiens (human), Dr...  
1  Primary immunodeficiency - Homo sapiens (human...  
2  Arrhythmogenic right ventricular cardiomyopath...  
3  Longevity regulating pathway - multiple specie...  
4  Metabolic pathways - Homo sapiens (human), Pri...  

>>> print KEGGp_.head()

KEGGid path:hsa04151 path:hsa00250 path:hsa05230 path:hsa00790  \
0  hsa:10327           NaN           NaN           NaN           NaN   
1    hsa:124           NaN           NaN           NaN           NaN   
2    hsa:125           NaN           NaN           NaN           NaN   
3    hsa:126           NaN           NaN           NaN           NaN   
4    hsa:127           NaN           NaN           NaN           NaN   

path:hsa04610 path:hsa00020 path:hsa04612 path:hsa01524 path:hsa01521  \
0           NaN           NaN           NaN           NaN           NaN   
1           NaN           NaN           NaN           NaN           NaN   
2           NaN           NaN           NaN           NaN           NaN   
3           NaN           NaN           NaN           NaN           NaN   
4           NaN           NaN           NaN           NaN           NaN   

 ...      path:hsa04672 path:hsa00730 path:hsa04670 path:hsa05168  \
0      ...                NaN           NaN           NaN           NaN   
1      ...                NaN           NaN           NaN           NaN   
2      ...                NaN           NaN           NaN           NaN   
3      ...                NaN           NaN           NaN           NaN   
4      ...                NaN           NaN           NaN           NaN   

path:hsa00531 path:hsa00532 path:hsa00533 path:hsa00534 path:hsa04914  \
0           NaN           NaN           NaN           NaN           NaN   
1           NaN           NaN           NaN           NaN           NaN   
2           NaN           NaN           NaN           NaN           NaN   
3           NaN           NaN           NaN           NaN           NaN   
4           NaN           NaN           NaN           NaN           NaN   

path:hsa04340  
0           NaN  
1           NaN  
2           NaN  
3           NaN  
4           NaN  

>>> print KEGGp_[["path:hsa04151"]].dropna().head()

path:hsa04151
27	hsa:2538
41	hsa:5105
42	hsa:5106
58	hsa:57818
66	hsa:92579

```
___

## ***KEGGmatrix***

Looks for all KEGG annotatios of an organism in biomart and the respective pathways in KEGG. It can also retrieve links to pathways figures with red labeled genes provided in a dataframe.

* **`organism`** a KEGG organism identifier
* **`dataset`** a biomaRt dataset
* **`query_attributes`** biomaRt query attributes, the name can change but the output should stay in the same order ie. 'ensembl_gene_id','kegg_enzyme'
* **`host`** biomart_host
* **`links`** if True, returns df_links
* **`dfexp`** a Pandas dataframe with the following columns: 'ensembl_gene_id', 'log2FC'
* **`kegg_db`** a KEGG database as recovered by the databasesKEGG function
* **`database`** a biomaRt database, depecrated, default=None.
* **`returns df`** a Pandas dataframe with the 'KEGGid','pathsIDs','pathName','ensembl_gene_id','kegg_enzyme'
* **`returns df_`** a matrix with a column for each KEGG pathway for a given organism and the expression values in the respective dfexp in parameter
* **`returns fullmatrix`** a matrix with a column for each KEGG pathway for a given organism
* **`returns df_links`** a dataframe with links for each pathway and the links in the dfexp highlighted red (if df_links.

```python
>>> import AGEpy as age
>>> print dge.head()

ensembl_gene_id    log2FC
0  ENSG00000272449  1.859500
1  ENSG00000130762  0.601051
2  ENSG00000083444 -0.881957
3  ENSG00000162493 -0.638433
4  ENSG00000253368  0.654517

>>> results=age.KEGGmatrix("hsa",'hsapiens_gene_ensembl', dfexp=dge, kegg_db="Ensembl" )

>>> print results[0]

index    KEGGid                                            pathIDs  \
0   58519  hsa:8813                       path:hsa00510, path:hsa01100   
1  120375  hsa:8711                                                NaN   
2  120515  hsa:6725                                                NaN   
3   48229  hsa:6714  path:hsa04360, path:hsa04611, path:hsa05100, p...   
4   20996  hsa:2534  path:hsa04725, path:hsa04360, path:hsa05416, p...   

                                        pathName  ensembl_gene_id  \
0  Metabolic pathways - Homo sapiens (human), N-G...  ENSG00000000419   
1                                                NaN  ENSG00000000938   
2                                                NaN  ENSG00000000938   
3  Mitophagy - animal - Homo sapiens (human), EGF...  ENSG00000000938   
4  Axon guidance - Homo sapiens (human), Adherens...  ENSG00000000938   

kegg_enzyme  
0  ec:2.4.1.83  
1  ec:2.7.10.2  
2  ec:2.7.10.2  
3  ec:2.7.10.2  
4  ec:2.7.10.2  

>>> print results[1]

ensembl_gene_id  kegg_enzyme     KEGGid    log2FC  path:hsa04151  \
26   ENSG00000149930  ec:2.7.11.1  hsa:10000  0.456776       0.456776   
251  ENSG00000183421  ec:2.7.11.1  hsa:10000  0.407559       0.407559   
476  ENSG00000112079  ec:2.7.11.1  hsa:10000 -0.498607      -0.498607   
701  ENSG00000196730  ec:2.7.11.1  hsa:10000  0.920654       0.920654   
926  ENSG00000123572  ec:2.7.11.1  hsa:10000  2.108620       2.108620   

path:hsa00250  path:hsa05230  path:hsa00790  path:hsa04610  \
26             NaN       0.456776            NaN            NaN   
251            NaN       0.407559            NaN            NaN   
476            NaN      -0.498607            NaN            NaN   
701            NaN       0.920654            NaN            NaN   
926            NaN       2.108620            NaN            NaN   

path:hsa00020      ...        path:hsa04672  path:hsa00730  \
26             NaN      ...                  NaN            NaN   
251            NaN      ...                  NaN            NaN   
476            NaN      ...                  NaN            NaN   
701            NaN      ...                  NaN            NaN   
926            NaN      ...                  NaN            NaN   

path:hsa04670  path:hsa05168  path:hsa00531  path:hsa00532  \
26             NaN            NaN            NaN            NaN   
251            NaN            NaN            NaN            NaN   
476            NaN            NaN            NaN            NaN   
701            NaN            NaN            NaN            NaN   
926            NaN            NaN            NaN            NaN   

path:hsa00533  path:hsa00534  path:hsa04914  path:hsa04340  
26             NaN            NaN       0.456776            NaN  
251            NaN            NaN       0.407559            NaN  
476            NaN            NaN      -0.498607            NaN  
701            NaN            NaN       0.920654            NaN  
926            NaN            NaN       2.108620            NaN  

>>> print results[2]

KEGGid path:hsa04151 path:hsa00250 path:hsa05230 path:hsa00790  \
0  hsa:10327           NaN           NaN           NaN           NaN   
1    hsa:124           NaN           NaN           NaN           NaN   
2    hsa:125           NaN           NaN           NaN           NaN   
3    hsa:126           NaN           NaN           NaN           NaN   
4    hsa:127           NaN           NaN           NaN           NaN   

path:hsa04610 path:hsa00020 path:hsa04612 path:hsa01524 path:hsa01521  \
0           NaN           NaN           NaN           NaN           NaN   
1           NaN           NaN           NaN           NaN           NaN   
2           NaN           NaN           NaN           NaN           NaN   
3           NaN           NaN           NaN           NaN           NaN   
4           NaN           NaN           NaN           NaN           NaN   

 ...      path:hsa04672 path:hsa00730 path:hsa04670 path:hsa05168  \
0      ...                NaN           NaN           NaN           NaN   
1      ...                NaN           NaN           NaN           NaN   
2      ...                NaN           NaN           NaN           NaN   
3      ...                NaN           NaN           NaN           NaN   
4      ...                NaN           NaN           NaN           NaN   

path:hsa00531 path:hsa00532 path:hsa00533 path:hsa00534 path:hsa04914  \
0           NaN           NaN           NaN           NaN           NaN   
1           NaN           NaN           NaN           NaN           NaN   
2           NaN           NaN           NaN           NaN           NaN   
3           NaN           NaN           NaN           NaN           NaN   
4           NaN           NaN           NaN           NaN           NaN   

path:hsa04340  
0           NaN  
1           NaN  
2           NaN  
3           NaN  
4           NaN  

>>> print results[3]

URL   pathway
0  http://www.kegg.jp/kegg-bin/show_pathway?hsa04...  hsa04151
1  http://www.kegg.jp/kegg-bin/show_pathway?hsa00...  hsa00250
2  http://www.kegg.jp/kegg-bin/show_pathway?hsa05...  hsa05230
3  http://www.kegg.jp/kegg-bin/show_pathway?hsa00...  hsa00790
4  http://www.kegg.jp/kegg-bin/show_pathway?hsa04...  hsa04610

>>> print results[3].loc[0,"URL"]

http://www.kegg.jp/kegg-bin/show_pathway?hsa04151/ec:2.7.11.1%09red/ec:2.7.11.22%09red/ec:2.7.10.1%09red/ec:3.1.3.16%09red/ec:2.7.10.2%09red/ec:2.3.2.27%09red/ec:1.14.13.39%09red/ec:2.7.12.2%09red/ec:3.1.3.48%09red
```
___
