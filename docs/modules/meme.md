## ___filterMotifs___

Selectes motifs from a meme file based on the number of sites.
* **`memeFile`** MEME file to be read
* **`outFile`** MEME file to be written
* **`minSites`** minimum number of sites each motif needs to have to be valid
* **`returns`** nothing

```python
>>> import AGEpy as age
>>> age.filterMotifs("/path/to/input.meme","/path/to/output.meme", 15)
```
___
