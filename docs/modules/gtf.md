## ___readGTF___

Reads a GTF file and labels the respective columns in agreement with GTF file standards:
'seqname','source','feature','start','end','score','strand','frame','attribute'.

**`readGTF(infile)`**

* **`infile`** /path/to/file.gtf
* **`returns`** a Pandas dataframe of the respective GTF

```python
>>> import AGEpy as age
>>> GTF=age.readGTF("gencode.v24.primary_assembly.annotation.gtf")
>>> print GTF.head()

seqname  source     feature  start    end score strand frame  \
0    chr1  HAVANA        gene  11869  14409     .      +     .   
1    chr1  HAVANA  transcript  11869  14409     .      +     .   
2    chr1  HAVANA        exon  11869  12227     .      +     .   
3    chr1  HAVANA        exon  12613  12721     .      +     .   
4    chr1  HAVANA        exon  13221  14409     .      +     .   

                                         attribute  
0  gene_id "ENSG00000223972.5"; gene_type "transc..."  
1  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."
```
___

## ***retrieve_GTF_field***

Returns a field of choice from the attribute column of the GTF.

**`retrieve_GTF_field(field,gtf)`**

* **`field`** field to be retrieved
* **`returns`** a Pandas dataframe with one column containing the field of choice

```python
>>> import AGEpy as age
>>> GTF=age.readGTF("/gencode.v24.primary_assembly.annotation.gtf")
>>> print GTF.head()

seqname  source     feature  start    end score strand frame  \
0    chr1  HAVANA        gene  11869  14409     .      +     .   
1    chr1  HAVANA  transcript  11869  14409     .      +     .   
2    chr1  HAVANA        exon  11869  12227     .      +     .   
3    chr1  HAVANA        exon  12613  12721     .      +     .   
4    chr1  HAVANA        exon  13221  14409     .      +     .   

                                         attribute  
0  gene_id "ENSG00000223972.5"; gene_type "transc..."  
1  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."

>>> GTF["gene_id"]=age.retrieve_GTF_field("gene_id",GTF)
>>> print GTF.head()

seqname  source     feature  start    end score strand frame  \
0    chr1  HAVANA        gene  11869  14409     .      +     .   
1    chr1  HAVANA  transcript  11869  14409     .      +     .   
2    chr1  HAVANA        exon  11869  12227     .      +     .   
3    chr1  HAVANA        exon  12613  12721     .      +     .   
4    chr1  HAVANA        exon  13221  14409     .      +     .   

                                         attribute            gene_id  
0  gene_id "ENSG00000223972.5"; gene_type "transc..."  ENSG00000223972.5  
1  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5  
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5  
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5  
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5
```
___

## ___attributesGTF___

List the type of attributes in a the attribute section of a GTF file

**`attributesGTF(inGTF)`**

* **`inGTF`** GTF dataframe to be analysed
* **`returns`** a list of attributes present in the attribute section

```python
>>> import AGEpy as age
>>> attributes=age.attributesGTF(GTF)
>>> print attributes

['gene_status', 'havana_gene', 'transcript_support_level', 'level', 'transcript_type', 'tag', 'protein_id', 'gene_id', 'exon_id', 'transcript_id', 'exon_number', 'ont', 'havana_transcript', 'ccdsid', 'transcript_name', 'gene_type', 'transcript_status', 'gene_name']
```
___
## ___parseGTF___

Reads an extracts all attributes in the attributes section of a GTF and constructs a new dataframe wiht one collumn per attribute instead of the attributes column.

**`parseGTF(inGTF)`**

* **`inGTF`** GTF dataframe to be parsed
* **`returns`** a dataframe of the orignal input GTF with attributes parsed

```python
>>> GTF=age.readGTF("gencode.v24.primary_assembly.annotation.gtf")
>>> print GTF.head()

seqname  source     feature  start    end score strand frame  \
0    chr1  HAVANA        gene  11869  14409     .      +     .   
1    chr1  HAVANA  transcript  11869  14409     .      +     .   
2    chr1  HAVANA        exon  11869  12227     .      +     .   
3    chr1  HAVANA        exon  12613  12721     .      +     .   
4    chr1  HAVANA        exon  13221  14409     .      +     .   

                                         attribute  
0  gene_id "ENSG00000223972.5"; gene_type "transc..."  
1  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."  

>>> GTFpa=age.parseGTF(GTF)
>>> print GTFpa.head()

seqname  source     feature  start    end score strand frame gene_status  \
0    chr1  HAVANA        gene  11869  14409     .      +     .       KNOWN   
1    chr1  HAVANA  transcript  11869  14409     .      +     .       KNOWN   
2    chr1  HAVANA        exon  11869  12227     .      +     .       KNOWN   
3    chr1  HAVANA        exon  12613  12721     .      +     .       KNOWN   
4    chr1  HAVANA        exon  13221  14409     .      +     .       KNOWN   

            havana_gene    ...               exon_id      transcript_id  \
0  OTTHUMG00000000961.2    ...                   NaN                NaN   
1  OTTHUMG00000000961.2    ...                   NaN  ENST00000456328.2   
2  OTTHUMG00000000961.2    ...     ENSE00002234944.1  ENST00000456328.2   
3  OTTHUMG00000000961.2    ...     ENSE00003582793.1  ENST00000456328.2   
4  OTTHUMG00000000961.2    ...     ENSE00002312635.1  ENST00000456328.2   

  exon_number  ont     havana_transcript ccdsid transcript_name  \
0         NaN  NaN                   NaN    NaN             NaN   
1         NaN  NaN  OTTHUMT00000362751.1    NaN     DDX11L1-002   
2           1  NaN  OTTHUMT00000362751.1    NaN     DDX11L1-002   
3           2  NaN  OTTHUMT00000362751.1    NaN     DDX11L1-002   
4           3  NaN  OTTHUMT00000362751.1    NaN     DDX11L1-002   

                            gene_type transcript_status gene_name  
0  transcribed_unprocessed_pseudogene               NaN   DDX11L1  
1  transcribed_unprocessed_pseudogene             KNOWN   DDX11L1  
2  transcribed_unprocessed_pseudogene             KNOWN   DDX11L1  
3  transcribed_unprocessed_pseudogene             KNOWN   DDX11L1  
4  transcribed_unprocessed_pseudogene             KNOWN   DDX11L1
```
___

## ___writeGTF___

Write a GTF dataframe into a file.

**`writeGTF(inGTF,file_path)`**

* **`inGTF`** GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
* **`file_path`** /path/to/the/file.gtf
* **`returns`** nothing

```python
>>> import AGEpy as age
>>> writeGTF(GTFpa,"/path/to/new/file.gtf")
```
___

## ___MAPGenoToTrans___

Gets all positions of all bases in an exon.

**`MAPGenoToTrans(parsedGTF,feature)`**

* **`df`** a Pandas dataframe with 'start','end', and 'strand' information for each entry. df must contain ['seqname','feature','start','end','strand','frame','gene_id',  'transcript_id','exon_id','exon_number']
* **`feature`** feature upon wich to generate the map, eg. 'exon' or 'transcript'
* **`returns`** a dictionary with a string with the comma separated positions of all bases in the exon

```python
>>> import AGEpy as age
>>> print GTF.head()

seqname  source     feature  start    end score strand frame  \
0    chr1  HAVANA        gene  11869  14409     .      +     .   
1    chr1  HAVANA  transcript  11869  14409     .      +     .   
2    chr1  HAVANA        exon  11869  12227     .      +     .   
3    chr1  HAVANA        exon  12613  12721     .      +     .   
4    chr1  HAVANA        exon  13221  14409     .      +     .   

                                         attribute            gene_id  \
0  gene_id "ENSG00000223972.5"; gene_type "transc..."  ENSG00000223972.5   
1  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   

     transcript_id            exon_id exon_number  
0                NaN                NaN         NaN  
1  ENST00000456328.2                NaN         NaN  
2  ENST00000456328.2  ENSE00002234944.1           1  
3  ENST00000456328.2  ENSE00003582793.1           2  
4  ENST00000456328.2  ENSE00002312635.1           3

>>> GtoT=age.MAPGenoToTrans(GTF,"exon")
>>> print GtoT

{ENST23923910:'234,235,236,1021,..'}
```
___

## ___GTFtoBED___

Transform a GTF dataframe into a bed dataframe

**`GTFtoBED(inGTF,name)`**

* **`inGTF`** GTF dataframe for transformation
* **`name`** field of the GTF data frame to be use for the bed 'name' positon
* **`returns`** a bed dataframe with the corresponding bed fiels: 'chrom','chromStart','chromEnd','name','score','strand'

```python
>>> import AGEpy as age
>>> bed = age.GTFtoBED(GTF, "gene_id")
```
___
## ___GetTransPosition___

Maps a genome position to transcript positon.

**`GetTransPosition(df, field, dic, refCol="transcript_id")`**

* **`df`** a Pandas dataframe
* **`field`** the head of the column containing the genomic position
* **`dic`** a dictionary containing for each transcript the respective bases eg. {ENST23923910:'234,235,236,1021,..'}. See *MAPGenoToTrans*.
* **`refCol`** header of the reference column with IDs, eg. 'transcript_id'

```python
>>> import AGEpy as age
>>> print GTF_.head()

seqname  source feature  start    end score strand frame  \
2    chr1  HAVANA    exon  11869  12227     .      +     .   
3    chr1  HAVANA    exon  12613  12721     .      +     .   
4    chr1  HAVANA    exon  13221  14409     .      +     .   
6    chr1  HAVANA    exon  12010  12057     .      +     .   
7    chr1  HAVANA    exon  12179  12227     .      +     .   

                                         attribute            gene_id  \
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
6  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
7  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   

     transcript_id            exon_id exon_number  target  
2  ENST00000456328.2  ENSE00002234944.1           1   12000  
3  ENST00000456328.2  ENSE00003582793.1           2   12617  
4  ENST00000456328.2  ENSE00002312635.1           3   14000  
6  ENST00000450305.2  ENSE00001948541.1           1   12040  
7  ENST00000450305.2  ENSE00001671638.2           2   12210  

>>> GTF_["transcript target"]=GTF_.apply(age.GetTransPosition, \
                                     args=("target",GtoT),axis=1)
>>> print GTF_.head()

seqname  source feature  start    end score strand frame  \
2    chr1  HAVANA    exon  11869  12227     .      +     .   
3    chr1  HAVANA    exon  12613  12721     .      +     .   
4    chr1  HAVANA    exon  13221  14409     .      +     .   
6    chr1  HAVANA    exon  12010  12057     .      +     .   
7    chr1  HAVANA    exon  12179  12227     .      +     .   

                                         attribute            gene_id  \
2  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
3  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
4  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
6  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   
7  gene_id "ENSG00000223972.5"; transcript_id "EN..."  ENSG00000223972.5   

     transcript_id            exon_id exon_number  target  transcript target  
2  ENST00000456328.2  ENSE00002234944.1           1   12000                132  
3  ENST00000456328.2  ENSE00003582793.1           2   12617                364  
4  ENST00000456328.2  ENSE00002312635.1           3   14000               1248  
6  ENST00000450305.2  ENSE00001948541.1           1   12040                 31  
7  ENST00000450305.2  ENSE00001671638.2           2   12210                 80  
```
___
