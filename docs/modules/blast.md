## ___BLASTquery___

Performs a blast query online. As in https://ncbi.github.io/blast-cloud/

**`BLASTquery(query,database,program,filter=None, format_type=None, expect=None, nucl_reward=None, nucl_penalty=None, gapcosts=None, matrix=None, hitlist_size=None, descriptions=None, alignments=None, ncbi_gi=None, threshold=None, word_size=None, composition_based_statistics=None, organism=None, others=None, num_threads=None, baseURL="http://blast.ncbi.nlm.nih.gov", verbose=False)`**

* **`query`** Search query. Allowed values: Accession, GI, or FASTA.
* **`database`** BLAST database. Allowed values: nt, nr, refseq_rna, refseq_protein, swissprot, pdbaa, pdbnt
* **`program`** BLAST program. Allowed values:  blastn, megablast, blastp, blastx, tblastn, tblastx
* **`filter`** Low complexity filtering. Allowed values: F to disable. T or L to enable. Prepend “m” for mask at lookup (e.g., mL)
* **`format_type`** Report type. Allowed values: HTML, Text, XML, XML2, JSON2, or Tabular. HTML is the default.
* **`expect`** Expect value. Allowed values: Number greater than zero.
* **`nucl_reward`** Reward for matching bases (BLASTN and megaBLAST). Allowed values: Integer greater than zero.
* **`nucl_penalty`** Cost for mismatched bases (BLASTN and megaBLAST). Allowed values: Integer less than zero.
* **`gapcosts`** Gap existence and extension costs. Allowed values: Pair of positive integers separated by a space such as “11 1”.
* **`matrix`** Scoring matrix name. Allowed values: One of BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM250, PAM30 or PAM70. Default: BLOSUM62 for all applicable programs.
* **`hitlist_size`** Number of databases sequences to keep. Allowed values: Integer greater than zero.
* **`descriptions`** Number of descriptions to print (applies to HTML and Text). Allowed values: Integer greater than zero.
* **`alignments`** Number of alignments to print (applies to HTML and Text). Allowed values: Integer greater than zero.
* **`ncbi_gi`** Show NCBI GIs in report. Allowed values: T or F.
* **`threshold`** Neighboring score for initial words. Allowed values: Positive integer (BLASTP default is 11). Does not apply to BLASTN or MegaBLAST).
* **`word_size`** Size of word for initial matches. Allowed values: Positive integer.
* **`composition_based_statistics`** Composition based statistics algorithm to use. Allowed values: One of 0, 1, 2, or 3. See comp_based_stats command line option in the BLAST+ user manual for details.
* **`organism`** an organism as in https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
* **`others`** here you can add other parameters as seen in a blast bookmarked page. Define you query in https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
        Once your query is defined click on "Bookmark" on right upper side of the page. You can copy fragments of the URL
        which define the query. Eg. For organism "Homo sapiens (taxid:9606)" you will see the string "EQ_MENU=Homo%20sapiens%20%28taxid%3A9606%29" - this is
        the string you can use here in others.
* **`num_threads`** Number of virtual CPUs to use. Allowed values: Integer greater than zero (default is 1). Supported only on the cloud.
* **`verbose`** print more

* **`returns`** BLAST search request identifier

```python
>>> import AGEpy as age
>>> seq="CTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTAC"
>>> RID=age.BLASTquery(seq,"nt","blastn")
>>> print RID

4MS2JV8T014
```
___

## ___BLASTcheck___

Checks the status of a query.

**`BLASTcheck(rid,baseURL="http://blast.ncbi.nlm.nih.gov")`**

* **`rid`**  BLAST search request identifier. Allowed values: The Request ID (RID) returned when the search was submitted
* **`baseURL`** server url. Default=http://blast.ncbi.nlm.nih.gov

* **`returns status`** status for the query.
* **`returns therearehist`** yes or no for existing hits on a finished query.

```python
>>> import AGEpy as age
>>> status, therearehits=age.BLASTcheck(RID)

RID: 4MRYDZSC014; status:READY; hits: yes

>>> print status, therearehits

READY yes
```
___

## ___BLASTresults___

Retrieves results for an RID.

**`BLASTresults(rid, format_type="Tabular", hitlist_size= None, alignments=None, ncbi_gi = None, format_object=None, baseURL="http://blast.ncbi.nlm.nih.gov")`**

* **`rid`** BLAST search request identifier. Allowed values: The Request ID (RID) returned when the search was submitted
* **`format_type`** Report type. Allowed values: HTML, Text, XML, XML2, JSON2, or Tabular.
* **`hitlist_size`** Number of databases sequences to keep. Allowed values: Integer greater than zero.
* **`alignments`** Number of alignments to print (applies to HTML and Text). Allowed values: Integer greater than zero.
* **`ncbi_gi`** Show NCBI GIs in report. Allowed values: T or F.
* **`format_object`** Object type. Allowed values: SearchInfo (status check) or Alignment (report formatting).
* **`baseURL`** server url. Default=http://blast.ncbi.nlm.nih.gov

* **`returns`** the result of a BLAST query. If format_type="Tabular" it will parse the content into a Pandas dataframe.

```python
>>> import AGEpy as age
>>> r=age.BLASTresults(RID)
>>> print r.head()

query id                                        subject ids  \
0  Query_17381                       gi|1012955506|gb|JN214348.1|   
1  Query_17381                       gi|631786534|tpe|HG975427.1|   
2  Query_17381                        gi|369762889|gb|JN900492.1|   
3  Query_17381                   gi|371502118|ref|NM_001126118.1|   
4  Query_17381  gi|371502115|ref|NM_001126112.2|;gi|454521556|...   

query acc.ver subject acc.ver % identity alignment length mismatches  \
0   Query_17381      JN214348.1    100.000             1190          0   
1   Query_17381      HG975427.1    100.000             1190          0   
2   Query_17381      JN900492.1    100.000             1190          0   
3   Query_17381  NM_001126118.1    100.000             1190          0   
4   Query_17381  NM_001126112.2    100.000             1190          0   

gap opens q. start q. end s. start s. end evalue bit scor  
0         0        1   1190      614   1803    0.0     2147  
1         0        1   1190      766   1955    0.0     2147  
2         0        1   1190      877   2066    0.0     2147  
3         0        1   1190      888   2077    0.0     2147  
4         0        1   1190      768   1957    0.0     2147
```
