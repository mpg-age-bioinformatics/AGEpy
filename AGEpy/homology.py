import pandas as pd
import os

def getHomoloGene(taxfile="build_inputs/taxid_taxname",\
                  genefile="homologene.data",\
                  proteinsfile="build_inputs/all_proteins.data",\
                  proteinsclusterfile="build_inputs/proteins_for_clustering.data",\
                  baseURL="http://ftp.ncbi.nih.gov/pub/HomoloGene/current/"):
    """
    Returns NBCI's Homolog Gene tables.

    :param taxfile: path to local file or to baseURL/taxfile
    :param genefile: path to local file or to baseURL/genefile
    :param proteinsfile: path to local file or to baseURL/proteinsfile
    :param proteinsclusterfile: path to local file or to baseURL/proteinsclusterfile
    :param baseURL: baseURL for downloading files

    :returns genedf: Homolog gene Pandas dataframe
    :returns protclusdf: Pandas dataframe. Lists one protein per gene that were used for homologene clustering.
                        If a gene has multiple protein accessions derived from alternative splicing,
                        only one protein isoform that give most protein alignment to proteins in other species
                        was selected for clustering and it is listed in this file.
    :returns proteinsdf: Pandas dataframe. Lists all proteins and their gene information.
                        If a gene has multple protein accessions derived from alternative splicing event,
                        each protein accession is list in a separate line.
    """

    def getDf(inputfile):
        if os.path.isfile(inputfile):
            df=pd.read_table(inputfile,header=None)
        else:
            df = urllib2.urlopen(baseURL+inputfile)
            df=df.read().split("\n")
            df=[ s for s in df if len(s) > 0 ]
            df=[s.split("\t") for s in df]
            df=pd.DataFrame(df)
        return df

    taxdf=getDf(taxfile)
    taxdf.set_index([0],inplace=True)
    taxdi=taxdf.to_dict().get(1)

    genedf=getDf(genefile)
    genecols=["HID","Taxonomy ID","Gene ID","Gene Symbol","Protein gi","Protein accession"]
    genedf.columns=genecols
    genedf["organism"]=genedf["Taxonomy ID"].apply(lambda(x):taxdi.get(x))

    proteinsdf=getDf(proteinsfile)
    proteinscols=["taxid","entrez GeneID","gene symbol","gene description","protein accession.ver","mrna accession.ver",\
                 "length of protein  listed in column 5","-11) contains data about gene location on the genome",\
                  "starting position of gene in 0-based coordinate",\
                  "end position of the gene in 0-based coordinate","strand","nucleotide gi of genomic sequence where this gene is annotated"]
    proteinsdf.columns=proteinscols
    proteinsdf["organism"]=proteinsdf["taxid"].apply(lambda(x):taxdi.get(x))

    protclusdf=getDf(proteinsclusterfile)
    protclustercols=["taxid","entrez GeneID","gene symbol","gene description","protein accession.ver","mrna accession.ver",\
                 "length of protein  listed in column 5","-11) contains data about gene location on the genome",\
                  "starting position of gene in 0-based coordinate",\
                  "end position of the gene in 0-based coordinate","strand","nucleotide gi of genomic sequence where this gene is annotated"]
    protclusdf.columns=proteinscols
    protclusdf["organism"]=protclusdf["taxid"].apply(lambda(x):taxdi.get(x))

    return genedf, protclusdf, proteinsdf
