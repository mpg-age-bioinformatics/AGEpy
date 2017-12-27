## ___getFasta___

Retrieves a sequence from an opened multifasta file.

**`getFasta(opened_file, sequence_name)`**

* **`opened_file`** an opened multifasta file eg. opened_file=open("/path/to/file.fa",'r+')
* **`sequence_name`** the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
* **`returns`** a string with the sequence of interest

```python
>>> import AGEpy as age
>>> fafile="/path/to/GRCm38.dna.primary_assembly.fa"
>>> with open(fafile, "r") as fastafile:
...     chr2=age.getFasta(fastafile, "2")
>>> print len(chr2)

182113224

>>> print chr2[82113224:82113284]

AGGGTGAATGATGTTTCTGGTACAGTGTACCAGTAAACCTAGCAGTAGGAGCATCAGTAT
```
___

## ___writeFasta___

Writes a fasta sequence into a file.

**`writeFasta(sequence, sequence_name, output_file)`**

* **`sequence`** a string with the sequence to be written
* **`sequence_name`** name of the the fasta sequence
* **`output_file`** /path/to/file.fa to be written
* **`returns`** nothing

```python
>>> import AGEpy as age
>>> print len(chr2)

182113224

>>> print chr2[82113224:82113284]

AGGGTGAATGATGTTTCTGGTACAGTGTACCAGTAAACCTAGCAGTAGGAGCATCAGTAT

>>> age.writeFasta(chr2,"2 my version of this sequence","/path/to/out/file.fa")
```
___

## ___rewriteFasta___

Rewrites a specific sequence in a multifasta file while keeping the sequence header.

**`rewriteFasta(sequence, sequence_name, fasta_in, fasta_out)`**

* **`sequence`** a string with the sequence to be written
* **`sequence_name`** the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
* **`fasta_in`** /path/to/original.fa
* **`fasta_out`** /path/to/destination.fa
* **`returns`** nothing

```python
>>> import AGEpy as age
>>> fafile="/path/to/GRCm38.dna.primary_assembly.fa"
>>> with open(fafile, "r") as fastafile:
...     chr2=age.getFasta(fastafile, "2")
>>> chr2=chr2.strip("N")
>>> age.rewriteFasta(chr2, "2", fafile, "/path/to/modified/file.fa")
```
___
