## ___GetBEDnarrowPeakgz___

Reads a gz compressed BED narrow peak file from a web address or local file and returns a pandas dataframe.

`writeBED(inBED, file_path)`

* **`URL_or_PATH_TO_file`** source of input bed. Either a web link or a path to a local file.

* **`returns`** a pandas dataframe of the inpud bed.

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

## ___writeBED___

Writes a bed dataframe into a bed file.

* **`inBED`** a pandas dataframe with the contents of the bed file to be written.
* **`file_path`** path to target file.

* **`returns`** nothing.

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

## ___dfTObedtool___

Transforms a pandas dataframe into a bedtool. Requires `bedtools` to be in your `path`.

* **`df`** a pandas dataframe.
* **`returns`** a bedtool.

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

>>> bedtool=age.dfTObedtool(bed)
>>> print bedtool.head()

chr7	139371278	139371296	Peak_0	1000	+	5.09062636514014	48.9834262537309	-1	-1
chr7	139371257	139371278	Peak_1	1000	+	5.0840236303159	48.7463712698062	-1	-1
chr7	155781335	155781431	Peak_2	1000	+	3.70481328524336	42.6519289009201	-1	-1
chr7	87156569	87156676	Peak_3	1000	+	3.95023151551588	37.7848384917051	-1	-1
chr7	105073472	105073521	Peak_4	1000	+	4.14556204062503	34.0756845242392	-1	-1
chr7	128761857	128761952	Peak_5	1000	+	4.02131461357736	33.9350181783027	-1	-1
chr7	121296414	121296454	Peak_6	1000	+	3.50632247892067	30.2512926812531	-1	-1
chr7	139368342	139368352	Peak_7	1000	+	4.41912711395099	29.6666535015756	-1	-1
chr7	87155583	87155635	Peak_8	1000	+	4.08769554637519	29.3752024210392	-1	-1
chr7	105540000	105540028	Peak_9	1000	+	4.2212263105571	29.0451450847765	-1	-1

>>> print type(bed)

<class 'pandas.core.frame.DataFrame'>

>>> print type(bedtool)

<class 'pybedtools.bedtool.BedTool'>
```

___

## ___GetPeaksExons___

Annotates a bedtool, BED narrow peak.

* **`bed`** a pandas dataframe in bed format
* **`parsedGTF`** a parsed GTF file as outputed by parseGTF()

* **`returns`** a Pandas dataframe

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

>>> GTF=age.readGTF("/beegfs/group_bit/data/projects/departments/Bioinformatics/bit_RNAseq_eCLIP/downloads/gencode.v24.primary_assembly.annotation.gtf")
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

>>> bedAn=age.GetPeaksExons(bed,GTFpa)
>>> print bedAn.head()

chrom  chromStart   chromEnd     name  score strand  signalValue  \
0  chr7   155781335  155781431   Peak_2   1000      +     3.704813   
1  chr7   155781335  155781431   Peak_2   1000      +     3.704813   
2  chr7   121296414  121296454   Peak_6   1000      +     3.506322   
3  chr7    87155538   87155583  Peak_16   1000      +     4.077391   
4  chr7   107904733  107904812  Peak_17   1000      +     3.674368   

   -log10(pValue)  -log10(qvalue)  peak         ...           \
0       42.651929              -1    -1         ...            
1       42.651929              -1    -1         ...            
2       30.251293              -1    -1         ...            
3       22.798739              -1    -1         ...            
4       21.118496              -1    -1         ...            

              gene_id  exon_id_count  exon_id norm. mean -log10(pValue)  \
0  ENSG00000184863.10              1                          42.651929   
1  ENSG00000184863.10              1                          42.651929   
2  ENSG00000106034.17              1                          30.251293   
3  ENSG00000135164.18              3                        2951.868281   
4  ENSG00000091140.12              1                          21.118496   

  exon_id signalValue transcript_id_count  \
0            3.704813                   1   
1            3.704813                   1   
2            3.506322                   1   
3           42.703999                   3   
4            3.674368                   1   

  transcript_id norm. mean -log10(pValue)  transcript_id signalValue  \
0                               42.651929                   3.704813   
1                               42.651929                   3.704813   
2                               30.251293                   3.506322   
3                             2951.868281                  42.703999   
4                               21.118496                   3.674368   

  gene_id_count gene_id norm. mean -log10(pValue)  gene_id signalValue  
0             4                        116.619012            17.830941  
1             4                        116.619012            17.830941  
2             2                         30.251293             2.144090  
3             8                       3300.707425            73.902289  
4             5                        135.139064            22.210269  
```
**gene_id_count**: number of intervals overlapping this gene

**transcript_id_count**: number of intervals overlapping this transcript

**exon_id_count**: number of intervals overlapping this exon
___
