import pandas as pd
import numpy as np
from collections import OrderedDict
import pybedtools
import csv

def readGTF(infile):
    """
    Reads a GTF file and labels the respective columns in agreement with GTF file standards:
    'seqname','source','feature','start','end','score','strand','frame','attribute'.

    :param infile: path/to/file.gtf
    :returns: a Pandas dataframe of the respective GTF

    """
    df=pd.read_table(infile, sep='\t', comment="#", header=None, dtype=str)
    df.columns=['seqname','source','feature','start','end','score','strand','frame','attribute']
    #df = df.astype(str) # from DTC
    return df

def retrieve_GTF_field(field,gtf):
    """
    Returns a field of choice from the attribute column of the GTF

    :param field: field to be retrieved
    :returns: a Pandas dataframe with one columns containing the field of choice

    """
    inGTF=gtf.copy()
    def splits(x):
        l=x.split(";")
        l=[ s.split(" ") for s in l]
        res=np.nan
        for s in l:
            if field in s:
                if '"' in s[-1]:
                    res=s[-1][1:-1]
                else:
                    res=s[-1]
        return res

    inGTF[field]=inGTF['attribute'].apply(lambda x: splits(x))
    return inGTF[[field]]

def attributesGTF(inGTF):
    """
    List the type of attributes in a the attribute section of a GTF file

    :param inGTF: GTF dataframe to be analysed
    :returns: a list of attributes present in the attribute section

    """
    df=pd.DataFrame(inGTF['attribute'].str.split(";").tolist())
    desc=[]
    for i in df.columns.tolist():
        val=df[[i]].dropna()
        val=pd.DataFrame(val[i].str.split(' "').tolist())[0]
        val=list(set(val))
        for v in val:
            if len(v) > 0:
                l=v.split(" ")
                if len(l)>1:
                    l=l[1]
                else:
                    l=l[0]
                desc.append(l)
    desc=list(set(desc))
    finaldesc=[]
    for d in desc:
        if len(d) > 0:
            finaldesc.append(d)
    return finaldesc

def parseGTF(inGTF):
    """
    Reads an extracts all attributes in the attributes section of a GTF and constructs a new dataframe wiht one collumn per attribute instead of the attributes column

    :param inGTF: GTF dataframe to be parsed
    :returns: a dataframe of the orignal input GTF with attributes parsed.

    """

    desc=attributesGTF(inGTF)
    ref=inGTF.copy()
    ref.reset_index(inplace=True, drop=True)
    df=ref.drop(['attribute'],axis=1).copy()
    for d in desc:
        field=retrieve_GTF_field(d,ref)
        df=pd.concat([df,field],axis=1)
    return df

def writeGTF(inGTF,file_path):
    """
    Write a GTF dataframe into a file

    :param inGTF: GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
    :param file_path: path/to/the/file.gtf
    :returns: nothing
    """
    cols=inGTF.columns.tolist()
    if len(cols) == 9:
        if 'attribute' in cols:
            df=inGTF
    else:
        df=inGTF[cols[:8]]
        df['attribute']=""
        for c in cols[8:]:
            if c == cols[len(cols)-1]:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'";'
            else:
                df['attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'"; '
    df.to_csv(file_path, sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)

def GTFtoBED(inGTF,name):
    """
    Transform a GTF dataframe into a bed dataframe

    :param inGTF: GTF dataframe for transformation
    :param name: field of the GTF data frame to be use for the bed 'name' positon

    returns: a bed dataframe with the corresponding bed fiels: 'chrom','chromStart','chromEnd','name','score','strand'
    """

    bed=inGTF.copy()
    bed.reset_index(inplace=True, drop=True)
    if name not in bed.columns.tolist():
        field=retrieve_GTF_field(name, bed)
        bed=pd.concat([bed,field],axis=1)
    bed=bed[['seqname','start','end',name,'score','strand']]
    bed.columns=['chrom','chromStart','chromEnd','name','score','strand']
    bed.drop_duplicates(inplace=True)
    bed.reset_index(inplace=True,drop=True)
    return bed

def MAPGenoToTrans(parsedGTF,feature):
    """
    Gets all positions of all bases in an exon

    :param df: a Pandas dataframe with 'start','end', and 'strand' information for each entry.
                df must contain 'seqname','feature','start','end','strand','frame','gene_id',
                'transcript_id','exon_id','exon_number']
    :param feature: feature upon wich to generate the map, eg. 'exon' or 'transcript'

    :returns: a string with the comma separated positions of all bases in the exon
    """
    GenTransMap=parsedGTF[parsedGTF["feature"]==feature]
    def getExonsPositions(df):
        start=int(df["start"])
        stop=int(df["end"])
        strand=df["strand"]
        r=range(start,stop+1)
        if strand=="-":
            r.sort(reverse=True)
        r=[ str(s) for s in r]
        return ",".join(r)

    GenTransMap["feature_bases"]=GenTransMap.apply(getExonsPositions, axis=1)
    GenTransMap=GenTransMap.sort_values(by=["transcript_id","exon_number"],ascending=True)
    def CombineExons(df):
        return pd.Series(dict( feature_bases = ','.join(df['feature_bases']) ) )
    GenTransMap=GenTransMap.groupby("transcript_id").apply(CombineExons)
    GenTransMap=GenTransMap.to_dict().get("feature_bases")

    return GenTransMap

def GetTransPosition(df,field,dic,refCol="transcript_id"):
    """
    Maps a genome position to transcript positon"

    :param df: a Pandas dataframe
    :param field: the head of the column containing the genomic position
    :param dic: a dictionary containing for each transcript the respective bases eg. {ENST23923910:'234,235,236,1021,..'}
    :param refCol: header of the reference column with IDs, eg. 'transcript_id'

    :returns: position on transcript
    """
    try:
        gen=str(int(df[field]))
        transid=df[refCol]
        bases=dic.get(transid).split(",")
        bases=bases.index(str(gen))+1
    except:
        bases=np.nan
    return bases

def getPromotersBed(gtf,fa,upstream=2000,downstream=200):
    """
    Reads a gtf file and returns a bed file for the promoter coordinates.
    
    :param gtf: path/to/file.gtf. Must be an ensembl gtf.
    :param fa: path/to/fasta.fa. Must be an ensembl fasta file.
    :param upstream: number of bases upstream of transcript start sites the promoter should start
    :param downstream: number of bases downstream of transcript start sites the promoter should end
    :returns: a pandas dataframe in bed format

    """
    chrsizes={}
    with open(fa, "r") as f:
        for line in f.readlines():
            if line[0] == ">":
                l=line.split(" ")
                seqname=l[0][1:]
                size=int(l[2].split(":")[-2])
                chrsizes[seqname]=size
    gtf=readGTF(gtf)
    gtf=gtf[gtf["feature"]=="transcript"]
    gtf.reset_index(inplace=True, drop=True)

    gtf["gene_id"]=retrieve_GTF_field(field="gene_id",gtf=gtf)
    gtf["gene_name"]=retrieve_GTF_field(field="gene_name",gtf=gtf)

    def getcoord(df):
        seqname=df["seqname"]
        strand=df["strand"]
        if strand == "+":
            tss=int(df["start"])
            promoter_start=tss-upstream
            promoter_end=tss+downstream
        else:
            tss=int(df["end"])
            promoter_start=tss-downstream
            promoter_end=tss+upstream

        if promoter_start < 0:
            promoter_start=0
        if promoter_end > chrsizes[seqname]:
            promoter_end=chrsizes[seqname]

        return str(promoter_start)+","+str(promoter_end)

    gtf["promoter"]=gtf.apply(getcoord, axis=1)
    gtf["start"]=gtf["promoter"].apply(lambda x: int(x.split(",")[0]) )
    gtf["end"]=gtf["promoter"].apply(lambda x: int(x.split(",")[1]) )
    
    gtf["id, name"]=gtf["gene_id"]+", "+gtf["gene_name"]
    gtf_=gtf.drop(["source","feature","attribute","promoter","gene_id","gene_name"],axis=1)
    gtf_=gtf_.drop_duplicates()
    gtf_counts=gtf_[["id, name"]]
    gtf_counts["#"]=1
    gtf_counts=gtf_counts.groupby(["id, name"]).sum()
    beds=gtf_[["seqname","start","end","id, name","score","strand"]]
    beds.columns=['chrom', 'start', 'stop', 'name', 'score', 'strand']
    beds=beds[beds["name"].isin( gtf_counts[gtf_counts["#"]==1].index.tolist() )]
    genes=[ s for s in list(set(gtf_counts[gtf_counts["#"]>1].index.tolist())) if str(s).lower() != "nan" ]

    for gene_id in genes:
        tmp=gtf[gtf["id, name"]==gene_id]
        strand=tmp["strand"].tolist()[0]
        bed=GTFtoBED(inGTF=tmp,name="id, name")
        bed = pybedtools.BedTool.from_dataframe(bed)
        bed=bed.sort()
        bed=bed.merge()
        bed = pd.read_table(bed.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
        bed["name"]=gene_id
        bed["score"]="."
        bed["strand"]=strand
        beds=pd.concat([beds,bed])
    
    beds = pybedtools.BedTool.from_dataframe(beds)
    beds = beds.sort()
    beds = pd.read_table(beds.fn, names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    
    beds.reset_index(inplace=True, drop=True)
    beds["i"]=beds.index.tolist()
    beds["i"]=beds["i"].astype(str)
    beds["name"]=beds["i"]+": "+beds["name"]
    beds=beds.drop(["i"],axis=1)

    return beds
