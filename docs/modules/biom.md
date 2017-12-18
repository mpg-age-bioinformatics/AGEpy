## ___datasetsBM___

Lists BioMart datasets.

* **`host`** address of the host server, default='http://www.ensembl.org/biomart'
* **`returns`** nothing

```python
>>> import AGEpy as age
>>> age.datasetsBM()

u'acarolinensis_gene_ensembl'	Anole lizard genes (AnoCar2.0),
u'acarolinensis_genomic_sequence'	Anole lizard sequences (AnoCar2.0),
u'amelanoleuca_gene_ensembl'	Panda genes (ailMel1),
u'amelanoleuca_genomic_sequence'	Panda sequences (ailMel1),
u'amexicanus_gene_ensembl'	Cave fish genes (AstMex102),
u'amexicanus_genomic_sequence'	Cave fish sequences (AstMex102),
u'anancymaae_gene_ensembl'	Ma's night monkey genes (Anan_2.0),
u'anancymaae_genomic_sequence'	Ma's night monkey sequences (Anan_2.0),
u'aplatyrhynchos_gene_ensembl'	Duck genes (BGI_duck_1.0),
u'aplatyrhynchos_genomic_sequence'	Duck sequences (BGI_duck_1.0),
u'btaurus_gene_ensembl'	Cow genes (UMD3.1),
u'btaurus_genomic_sequence'	Cow sequences (UMD3.1),
u'btaurus_marker_end'	marker_feature_end,
u'btaurus_marker_start'	marker_feature,
u'btaurus_qtl_feature'	qtl_feature,
.
.
.
```
___

## ___filtersBM___

Lists BioMart filters for a specific dataset.

* **`dataset`** dataset to list filters of
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> age.filtersBM('hsapiens_gene_ensembl')

u'affy_hc_g110'	'AFFY HC G110 probe ID(s) [e.g. 266_s_at]' (type	id_list, values	[]),
u'affy_hg_focus'	'AFFY HG Focus probe ID(s) [e.g. 212481_s_at]' (type	id_list, values	[]),
u'affy_hg_u133_plus_2'	'AFFY HG U133 Plus 2 probe ID(s) [e.g. 1553551_s_at]' (type	id_list, values	[]),
u'affy_hg_u133a'	'AFFY HG U133A probe ID(s) [e.g. 211600_at]' (type	id_list, values	[]),
u'affy_hg_u133a_2'	'AFFY HG U133A 2 probe ID(s) [e.g. 211600_at]' (type	id_list, values	[]),
u'affy_hg_u133b'	'AFFY HG U133B probe ID(s) [e.g. 224321_at]' (type	id_list, values	[]),
u'affy_hg_u95a'	'AFFY HG U95A probe ID(s) [e.g. 33866_at]' (type	id_list, values	[]),
u'affy_hg_u95av2'	'AFFY HG U95Av2 probe ID(s) [e.g. 33866_at]' (type	id_list, values	[]),
u'affy_hg_u95b'	'AFFY HG U95B probe ID(s) [e.g. 48794_s_at]' (type	id_list, values	[]),
u'affy_hg_u95c'	'AFFY HG U95C probe ID(s) [e.g. 66888_at]' (type	id_list, values	[]),
u'affy_hg_u95d'	'AFFY HG U95D probe ID(s) [e.g. 70806_at]' (type	id_list, values	[]),
u'affy_hg_u95e'	'AFFY HG U95E probe ID(s) [e.g. 88289_at]' (type	id_list, values	[]),
u'affy_hta_2_0'	'AFFY HTA 2 0 probe ID(s) [e.g. TC04001102.hg]' (type	id_list, values	[]),
u'affy_huex_1_0_st_v2'	'AFFY HuEx 1 0 st v2 probe ID(s) [e.g. 4037584]' (type	id_list, values	[]),
u'affy_hugene_1_0_st_v1'	'AFFY HuGene 1 0 st v1 probe ID(s) [e.g. 8165644]' (type	id_list, values	[]),
u'affy_hugene_2_0_st_v1'	'AFFY HuGene 2 0 st v1 probe ID(s) [e.g. 17100641]' (type	id_list, values	[]),
u'affy_hugenefl'	'AFFY HuGeneFL probe ID(s) [e.g. Z70759_at]' (type	id_list, values	[]),
u'affy_primeview'	'AFFY PrimeView probe ID(s) [e.g. 11761516_x_at]' (type	id_list, values	[]),
.
.
.

```
___

## ___attributesBM___

Lists BioMart attributes for a specific dataset.

* **`dataset`** dataset to list attributes of
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> age.attributesBM('hsapiens_gene_ensembl')

u'3_utr_end'	'3' UTR end' (default	False),
 u'3_utr_start'	'3' UTR start' (default	False),
 u'3utr'	'3' UTR' (default	False),
 u'5_utr_end'	'5' UTR end' (default	False),
 u'5_utr_start'	'5' UTR start' (default	False),
 u'5utr'	'5' UTR' (default	False),
 u'acarolinensis_homolog_associated_gene_name'	'Anole lizard gene name' (default	False),
 u'acarolinensis_homolog_canonical_transcript_protein'	'Query protein or transcript ID' (default	False),
 u'acarolinensis_homolog_chrom_end'	'Anole lizard chromosome/scaffold end (bp)' (default	False),
 u'acarolinensis_homolog_chrom_start'	'Anole lizard chromosome/scaffold start (bp)' (default	False),
 u'acarolinensis_homolog_chromosome'	'Anole lizard chromosome/scaffold name' (default	False),
 u'acarolinensis_homolog_dn'	'dN with Anole lizard' (default	False),
 u'acarolinensis_homolog_ds'	'dS with Anole lizard' (default	False),
 u'acarolinensis_homolog_ensembl_gene'	'Anole lizard gene stable ID' (default	False),
 .
 .
 .

```
___

## ___queryBM___

Queries BioMart.

* **`query_attributes`** list of attributes to recover from BioMart
* **`query_dataset`** dataset to query
* **`query_filter`** one BioMart filter associated with the items being queried
* **`query_items`** list of items to be queried (must assoiate with given filter)
* **`query_querydic`** for complex queries this option should be used instead of 'filters' and 'items' and a dictionary of filters provided here eg. querydic={"filter1":["item1","item2"],"filter2":["item3","item4"]}. If using querydic, don't query more than 350 items at once.
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** a Pandas dataframe of the queried attributes

```python
>>> import AGEpy as age
>>> queryDf=queryBM(query_attributes=["ensembl_gene_id","external_gene_name", \
                                  "go_id","name_1006","definition_1006"],\
                query_dataset='hsapiens_gene_ensembl')
>>> print queryDf.head()

ensembl_gene_id external_gene_name       go_id            name_1006  \
0  ENSG00000283891             MIR628  GO:0005615  extracellular space   
1  ENSG00000251931          RNU6-871P                                    
2  ENSG00000207766             MIR626                                    
3  ENSG00000275323         AC012314.7  GO:0003723          RNA binding   
4  ENSG00000275323         AC012314.7  GO:0005634              nucleus   

                                     definition_1006  
0  "That part of a multicellular organism outside..."  
1                                                     
2                                                     
3  "Interacting selectively and non-covalently wi..."  
4  "A membrane-bounded organelle of eukaryotic ce..."  
```
___

## ___FilterGOstring___

Filters GO terms based on given strings using ENSEMBL's biomart homology mapping.

* **`names_filter`** list of substrings to filter GO names on. Default=["age-", "aging", "aged", 'aging', 'aging.', 'aging,']
* **`exclude_names`** list of substrings to be used for exclusion of GO names. Default=["packaging","voltage","cleavage-",
                       "stage-1","cage-like","message-specific",
                       "damage-associated","stage-specific","foraging",
                       "DNA-damaging","engaging","damaged","packaged"]
* **`defs_filter`** list of substrings to filter GO defenitions on. Default=[" age-", " aging", " aged", ' aging', ' aging.', ' aging,']
* **`exclude_defs`** list of substrings to be used for exclustion of GO defenitions. Default=["packaging","voltage","cleavage-",
                         "stage-1","cage-like","message-specific",
                         "damage-associated","stage-specific","foraging",
                         "DNA-damaging","engaging","damaged","packaged"]
* **`host`** biomart host server, default="http://www.ensembl.org/biomart"
* **`HSA`** retrieved hsa dataframe
* **`MUS`** retrieved mus dataframe
* **`CEL`** retrieved cel dataframe
* **`DMEL`** retrieved dmel dataframe

* **`returns`**  homology_df, HSA, MUS, CEL, DMEL

```python
>>> import AGEpy as age
>>> homology_df, HSA, MUS, CEL, DMEL=age.FilterGOstring()
>>> print homology_df.head()

HSA_ensembl_gene_id HSA_external_gene_name  \
0     ENSG00000000003                 TSPAN6   
1     ENSG00000000005                   TNMD   
2     ENSG00000000460               C1orf112   
3     ENSG00000000971                    CFH   
4     ENSG00000002079                  MYH16   

                                         HSA_go_id  \
0  GO:0039532, , GO:0070062, GO:0016021, GO:00160...   
1  GO:0005737, , GO:0016020, GO:0035990, GO:00717...   
2                                                NaN   
3  , GO:0030449, GO:0070062, GO:0045087, GO:00725...   
4                                                NaN   

                                     HSA_name_1006  \
0  , negative regulation of NIK/NF-kappaB signali...   
1  , nuclear envelope, cytoplasm, negative regula...   
2                                                NaN   
3  , innate immune response, heparan sulfate prot...   
4                                                NaN   

                               HSA_definition_1006 MUS_ensembl_gene_id  \
0  "The component of a membrane consisting of the..."  ENSMUSG00000067377   
1  "The component of a membrane consisting of the..."  ENSMUSG00000031250   
2                                                NaN   ENSMUSG00000041406   
3  "Interacting selectively and non-covalently wi..."                 NaN   
4                                                NaN                  NaN   

CEL_ensembl_gene_id DMEL_ensembl_gene_id MUS_external_gene_name  \
0                 NaN                  NaN                 Tspan6   
1                 NaN                  NaN                   Tnmd   
2                 NaN                  NaN               BC055324   
3                 NaN                  NaN                   None   
4                 NaN                  NaN                   None   

                                         MUS_go_id   ...     \
0  GO:0039532, , GO:0070062, GO:0016021, GO:00160...   ...      
1  GO:0016020, GO:0035990, GO:0071773, GO:0016021...   ...      
2               GO:0005575, GO:0008150, GO:0003674,    ...      
3                                               None   ...      
4                                               None   ...      

                               MUS_definition_1006 CEL_external_gene_name  \
0  "The component of a membrane consisting of the..."                   None   
1  "The component of a membrane consisting of the..."                   None   
2  "Elemental activities, such as catalysis or bi..."                   None   
3                                               None                    None   
4                                               None                    None   

CEL_go_id CEL_name_1006 CEL_definition_1006 DMEL_external_gene_name  \
0      None          None                None                    None   
1      None          None                None                    None   
2      None          None                None                    None   
3      None          None                None                    None   
4      None          None                None                    None   

DMEL_go_id DMEL_name_1006 DMEL_definition_1006 evidence  
0       None           None                 None      NaN  
1       None           None                 None      NaN  
2       None           None                 None      NaN  
3       None           None                 None      NaN  
4       None           None                 None      NaN  

```

**evidence** indicates from which organisms there is evidence of the intended string
