## ___readSAM___

Reads and parses a sam file.

**`readSAM(SAMfile,header=False)`**

* **`SAMfile`** /path/to/file.sam
* **`header`** logical, if True, reads the header information
* **`returns`** a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' and a list of the headers if header=True

```python
>>> import AGEpy as age
>>> SAMdf=age.readSAM("sample1.sam")
>>> print SAMdf.head()

CIGAR  \
0  J00137:91:HJG75BBXX:6:1101:27458:1244    4     *         0    0       *   
1   J00137:91:HJG75BBXX:6:1101:2483:1226    4     *         0    0       *   
2   J00137:91:HJG75BBXX:6:1101:6593:1244   16    II  11210427  255  2S146M   
3   J00137:91:HJG75BBXX:6:1101:9293:1244    0     I  10433525  255    150M   
4  J00137:91:HJG75BBXX:6:1101:13271:1244   16   III   5277278  255    150M   

RNEXT PNEXT TLEN                                                SEQ  \
0     *     0    0  CCAAAATCAGTTACAAAAAAATTAAATATCGAGTTCCTCCCCCAGA...   
1     *     0    0  ACGTGACCGATGGTTGGCATGGCACGCATACCACGGAAGCGTCTGC...   
2     *     0    0  AACAACAGCAGCAGCAGATTTACCAAAGGTTCCCAGCAAGACTAAT...   
3     *     0    0  CTTGATTGTACTGCTGTGGTGGACCGCGTGGTCCTCCTTGTTGGTT...   
4     *     0    0  GGACATGATGATCATGGCCACGACTCTCATGGACATAGTCATGATC...   

                                              QUAL  
0  AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ...  
1  A-<AF7FF-AF-A<--F-7<FAA7FFF7-<-<F7JF---7--77J7...  
2  7-AA7JF)JJJF<--7FF--7--AAA--<A7<--<<A-7-A-<<FJ...  
3  AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ...  
4  JFFAJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJ...
```
___

## ___writeSAM___
Writes a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' into a sam file.

**`writeSAM(sam, SAMfile, header=None)`**

* **`sam`** pandas dataframe to be writen
* **`SAMfile`** /path/to/file.sam
* **`returns`** nothing

```
>>> import AGEpy as age
>>> age.writeSAM(SAMdf,"modified.sam")
```
___

## ___SAMflags___
Explains a SAM flag.

**`SAMflags(x)`**

* **`x`** flag
* **`returns`** complete SAM flag explanation

```
>>> import AGEpy as age
>>> print age.SAMflags(64)
```
["0: Read unpaired",
 "0:  Read not mapped in proper pair",
 "0:  Read mapped",
 "0:  Mate mapped",
 "0:  Read direct strand",
 "0:  Mate direct strand",
 "1:  First in pair",
 "0:  First in pair",
 "0:  Primary alignment",
 "0:  Read passes platform/vendor quality checks",
 "0:  Read is not PCR or optical duplicate",
 "0:  Not supplementary alignment"]

 ___
