### ___datasetsBM___

Lists BioMart datasets.

* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> datasetsBM()
```
___

### ___datasetsBM___

Lists BioMart filters for a specific dataset.

* **`dataset`** dataset to list filters of
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> filtersBM()

```
___

### ___attributesBM___

Lists BioMart attributes for a specific dataset.

* **`dataset`** dataset to list attributes of
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** nothing

```python
>>> import AGEpy as age
>>> attributesBM()

```
___

### ___queryBM___

Queries BioMart.

* **`query_filter`** one BioMart filter associated with the items being queried
* **`query_items`** list of items to be queried (must assoiate with given filter)
* **`query_querydic`** for complex queries this option should be used instead of 'filters' and 'items' and a dictionary of filters provided here eg. querydic={"filter1":["item1","item2"],"filter2":["item3","item4"]}. If using querydic, don't query more than 350 items at once.
* **`query_attributes`** list of attributes to recover from BioMart
* **`query_dataset`** dataset to query
* **`host`** address of the host server, default='http://www.ensembl.org/biomart'

* **`returns`** a Pandas dataframe of the queried attributes

```python
>>> import AGEpy as age
>>> queryDf=queryBM()


```
___

### ___FilterGOstring___

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
>>> homology_df=FilterGOstring()


```
