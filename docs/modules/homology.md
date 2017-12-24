## ___getHomoloGene___

Returns NBCI's Homolog Gene tables.

* **`taxfile`** path to local file or to baseURL/taxfile, default="build_inputs/taxid_taxname",
* **`genefile`** path to local file or to baseURL/genefile, defult="homologene.data"
* **`proteinsfile`** path to local file or to baseURL/proteinsfile, default="build_inputs/all_proteins.data"
* **`proteinsclusterfile`** path to local file or to baseURL/proteinsclusterfile, default="build_inputs/proteins_for_clustering.data"
* **`baseURL`** baseURL for downloading files, default="http://ftp.ncbi.nih.gov/pub/HomoloGene/current/"
* **`returns genedf`** Homolog gene Pandas dataframe
* **`returns protclusdf`** Pandas dataframe. Lists one protein per gene that were used for homologene clustering.
                    If a gene has multiple protein accessions derived from alternative splicing,
                    only one protein isoform that give most protein alignment to proteins in other species
                    was selected for clustering and it is listed in this file.
* **`returns proteinsdf`** Pandas dataframe. Lists all proteins and their gene information.
                    If a gene has multple protein accessions derived from alternative splicing event,
                    each protein accession is list in a separate line.

```python
>>> import AGEpy as age
>>> genedf, protclusdf, proteinsdf = age.getHomoloGene()
>>> print genedf.head()

HID Taxonomy ID Gene ID Gene Symbol Protein gi Protein accession  \
0   3        9606      34       ACADM    4557231       NP_000007.1   
1   3        9598  469356       ACADM  160961497    NP_001104286.1   
2   3        9544  705168       ACADM  109008502    XP_001101274.1   
3   3        9615  490207       ACADM  545503811    XP_005622188.1   
4   3        9913  505968       ACADM  115497690    NP_001068703.1   

                organism  
0            Homo sapiens  
1         Pan troglodytes  
2          Macaca mulatta  
3  Canis lupus familiaris  
4              Bos taurus  

>>> print protclusdf.head()

taxid entrez GeneID gene symbol gene description protein accession.ver  \
0  3702      10723019   AT1G27045        AT1G27045        NP_001185103.1   
1  3702      10723020   AT2G41231        AT2G41231        NP_001189726.1   
2  3702      10723023   AT1G24095        AT1G24095        NP_001185076.1   
3  3702      10723026   AT1G12855        AT1G12855        NP_001184976.1   
4  3702      10723027   AT4G22758        AT4G22758        NP_001190802.1   

  mrna accession.ver length of protein  listed in column 5  \
0     NM_001198174.1                                   227   
1     NM_001202797.1                                    99   
2     NM_001198147.1                                   213   
3     NM_001198047.1                                   462   
4     NM_001203873.1                                   255   

  -11) contains data about gene location on the genome  \
0                                          240254421     
1                                          240254678     
2                                          240254421     
3                                          240254421     
4                                          240256243     

  starting position of gene in 0-based coordinate  \
0                                         9391608   
1                                        17195291   
2                                         8523246   
3                                         4382159   
4                                        11958309   

  end position of the gene in 0-based coordinate strand  \
0                                        9393018      +   
1                                       17195914      +   
2                                        8524928      +   
3                                        4383610      +   
4                                       11960035      +   

  nucleotide gi of genomic sequence where this gene is annotated  \
0                                          AT1G27045               
1                                          AT2G41231               
2                                          AT1G24095               
3                                          AT1G12855               
4                                          AT4G22758               

               organism  
0  Arabidopsis thaliana  
1  Arabidopsis thaliana  
2  Arabidopsis thaliana  
3  Arabidopsis thaliana  
4  Arabidopsis thaliana

>>> print proteinsdf.head()

taxid entrez GeneID gene symbol gene description protein accession.ver  \
0  3702      10723019   AT1G27045        AT1G27045        NP_001185103.1   
1  3702      10723020   AT2G41231        AT2G41231        NP_001189725.1   
2  3702      10723020   AT2G41231        AT2G41231        NP_001189726.1   
3  3702      10723023   AT1G24095        AT1G24095        NP_001185076.1   
4  3702      10723026   AT1G12855        AT1G12855        NP_001184976.1   

 mrna accession.ver length of protein  listed in column 5  \
0     NM_001198174.1                                   227   
1     NM_001202796.1                                   104   
2     NM_001202797.1                                    99   
3     NM_001198147.1                                   213   
4     NM_001198047.1                                   462   

 -11) contains data about gene location on the genome  \
0                                          240254421     
1                                          240254678     
2                                          240254678     
3                                          240254421     
4                                          240254421     

 starting position of gene in 0-based coordinate  \
0                                         9391608   
1                                        17195291   
2                                        17195291   
3                                         8523246   
4                                         4382159   

 end position of the gene in 0-based coordinate strand  \
0                                        9393018      +   
1                                       17195914      +   
2                                       17195914      +   
3                                        8524928      +   
4                                        4383610      +   

 nucleotide gi of genomic sequence where this gene is annotated  \
0                                          AT1G27045               
1                                          AT2G41231               
2                                          AT2G41231               
3                                          AT1G24095               
4                                          AT1G12855               

              organism  
0  Arabidopsis thaliana  
1  Arabidopsis thaliana  
2  Arabidopsis thaliana  
3  Arabidopsis thaliana  
4  Arabidopsis thaliana
```
___
