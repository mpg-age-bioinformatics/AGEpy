### ___GetBEDnarrowPeakgz___

Reads a gz compressed BED narrow peak file from a web address or local file and returns a pandas dataframe.

**URL_or_PATH_TO_file:** Source of input bed. Either a web link or a path to a local file.

```python
>>> import AGEpy as age

>>> eCLIP_1_bednarrowPeak="https://www.encodeproject.org/files/ENCFF066PCT/@@download/ENCFF066PCT.bed.gz"
>>> bed=age.GetBEDnarrowPeakgz(eCLIP_1_bednarrowPeak)
>>> print bed.head()

chrom chromStart   chromEnd    name score strand       signalValue  \
0  chr7  139371278  139371296  Peak_0  1000      +  5.09062636514014   
1  chr7  139371257  139371278  Peak_1  1000      +   5.0840236303159   
2  chr7  155781335  155781431  Peak_2  1000      +  3.70481328524336   
3  chr7   87156569   87156676  Peak_3  1000      +  3.95023151551588   
4  chr7  105073472  105073521  Peak_4  1000      +  4.14556204062503   

     -log10(pValue) -log10(qvalue) peak  
0  48.9834262537309             -1   -1  
1  48.7463712698062             -1   -1  
2  42.6519289009201             -1   -1  
3  37.7848384917051             -1   -1  
4  34.0756845242392             -1   -1
```
___

### ___writeBED___

Writes a bed dataframe into a bed file.

**inBED:** A pandas dataframe with the contents of the bed file to be written.

**file_path:** Path to target file.

```python
>>> import AGEpy as age

>>> print bed.head()

chrom chromStart   chromEnd    name score strand       signalValue  \
0  chr7  139371278  139371296  Peak_0  1000      +  5.09062636514014   
1  chr7  139371257  139371278  Peak_1  1000      +   5.0840236303159   
2  chr7  155781335  155781431  Peak_2  1000      +  3.70481328524336   
3  chr7   87156569   87156676  Peak_3  1000      +  3.95023151551588   
4  chr7  105073472  105073521  Peak_4  1000      +  4.14556204062503   

     -log10(pValue) -log10(qvalue) peak  
0  48.9834262537309             -1   -1  
1  48.7463712698062             -1   -1  
2  42.6519289009201             -1   -1  
3  37.7848384917051             -1   -1  
4  34.0756845242392             -1   -1

>>> age.writeBED(bed,"/path/to/file.bed")
```
___

### ___dfTObedtool___

Transforms a pandas dataframe into a bedtool.
