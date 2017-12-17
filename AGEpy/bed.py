import pandas as pd
import os
from urllib import urlopen
import urllib2
import StringIO
import gzip
import pybedtools
from pybedtools import BedTool
from gtf import GTFtoBED

def writeBED(inBED, file_path):
    """
    Writes a bed dataframe into a bed file.
    Bed format: 'chrom','chromStart','chromEnd','name','score','strand'

    :param inBED: bed dataframe to be written.
    :param file_path: /path/to/file.bed

    :returns: nothing

    """
    inBED.to_csv(file_path,index=None,sep="\t",header=None)

def GetBEDnarrowPeakgz(URL_or_PATH_TO_file):
    """
    Reads a gz compressed BED narrow peak file from a web address or local file

    :param URL_or_PATH_TO_file: web address of path to local file

    :returns: a Pandas dataframe
    """

    if os.path.isfile(URL_or_PATH_TO_file):
        response=open(URL_or_PATH_TO_file, "r")
        compressedFile = StringIO.StringIO(response.read())
    else:
        response = urllib2.urlopen(URL_or_PATH_TO_file)
        compressedFile = StringIO.StringIO(response.read())
    decompressedFile = gzip.GzipFile(fileobj=compressedFile)
    out=decompressedFile.read().split("\n")
    out=[ s.split("\t") for s in out]
    out=pd.DataFrame(out)
    out.columns=["chrom","chromStart","chromEnd","name","score","strand","signalValue","-log10(pValue)","-log10(qvalue)","peak"]
    out["name"]=out.index.tolist()
    out["name"]="Peak_"+out["name"].astype(str)
    out=out[:-1]
    return out

def dfTObedtool(df):
    """
    Transforms a pandas dataframe into a bedtool

    :param df: Pandas dataframe

    :returns: a bedtool
    """

    df=df.astype(str)
    df=df.drop_duplicates()
    df=df.values.tolist()
    df=["\t".join(s) for s in df ]
    df="\n".join(df)
    df=BedTool(df, from_string=True)
    return df

def GetPeaksExons(bed,parsedGTF):
    """
    Annotates a bedtool, BED narrow peak

    :param bed: a pandas dataframe in bed format
    :param parsedGTF: a parsed GTF file as outputed by parseGTF() with the following columns

    :returns: a Pandas dataframe
    """

    bedtool_AB=dfTObedtool(bed)
        
    exonsGTF=parsedGTF[parsedGTF["feature"]=="exon"]
    exonsGTF.reset_index(inplace=True, drop=True)

    exonsBED=GTFtoBED(exonsGTF, "exon_id")
    exonsBED.columns=['chrom', 'chromStart', 'chromEnd', 'exon_id', 'score', 'strand']
    exonsBEDcols=exonsBED.columns.tolist()

    bedcols=bed.columns.tolist()
    exonsBEDcols_=[]
    for c in exonsBEDcols:
        if c in bedcols:
            exonsBEDcols_.append(c+"_exon")
        else:
            exonsBEDcols_.append(c)

    cols=[bedcols,exonsBEDcols_,["overlap"] ]
    cols=[item for sublist in cols for item in sublist]
   
    bedtool_exons=dfTObedtool(exonsBED)

    bedtool_target_exons=bedtool_AB.intersect(bedtool_exons, wo=True, s=True)
    dfTargetE=pd.read_table(bedtool_target_exons.fn, names=cols)
    ExonsTransGenes=parsedGTF[["exon_id","transcript_id","gene_id"]].drop_duplicates()
    dfTargets=pd.merge(dfTargetE,ExonsTransGenes,on=["exon_id"],how="left")
    dfTargets["count"]=1

    def getCounts(df,field):
        """
        For each field in a bed narrow peak returns the number or times that field is present,\
        the normalized mean of the '-log10(pValue)' and normalized mean of the signal value.

        :param df: a Pandas dataframe of a bed narrow peak
        :param field: field to analyse, ie. exons or transcripts

        :returns: a Pandas dataframe
        """

        tmp=df[[field,'name',"count"]].drop_duplicates()
        tmp=tmp.drop(["name"],axis=1)
        tmp["count"]=tmp["count"].astype(int)
        tmp.columns=[field,"%s_count" %str(field)]
        tmp=tmp.groupby(field, as_index=False).sum()
        df=pd.merge(df,tmp,on=field,how="left")

        tmp=df[[field,'name',"-log10(pValue)"]].drop_duplicates()
        tmp=tmp.drop(["name"],axis=1)
        tmp["-log10(pValue)"]=tmp["-log10(pValue)"].astype(float)
        tmp=tmp.groupby(field).apply(lambda l: reduce(lambda x, y: x*y, l["-log10(pValue)"]) )
        tmp=pd.DataFrame(tmp)
        tmp.reset_index(inplace=True,drop=False)
        tmp.columns=[field,"%s norm. mean -log10(pValue)" %str(field)]
        df=pd.merge(df,tmp,on=field,how="left")

        tmp=df[[field,'name',"signalValue"]].drop_duplicates()
        tmp=tmp.drop(["name"],axis=1)
        tmp["signalValue"]=tmp["signalValue"].astype(float)
        tmp=tmp.groupby(field).apply(lambda l: reduce(lambda x, y: x*y, l["signalValue"]) )
        tmp=pd.DataFrame(tmp)
        tmp.reset_index(inplace=True,drop=False)
        tmp.columns=[field,"%s signalValue" %str(field)]
        df=pd.merge(df,tmp,on=field,how="left")

        return df

    for f in ["exon_id","transcript_id"]:
        dfTargets=getCounts(dfTargets,f)

    def getCounts_GeneIDs(df):
        """
        For each gene id in a bed narrow peak returns the number or times that field is present,\
        the normalized mean of the '-log10(pValue)' and normalized mean of the signal value.

        :param df: a Pandas dataframe of a bed narrow peak

        :returns: a Pandas dataframe
        """

        field="gene_id"

        tmp=df[[field,"transcript_id","transcript_id_count"]].drop_duplicates()
        tmp=tmp.drop(["transcript_id"],axis=1)
        tmp["transcript_id_count"]=tmp["transcript_id_count"].astype(int)
        tmp.columns=[field,"%s_count" %str(field)]
        tmp=tmp.groupby(field, as_index=False).sum()
        df=pd.merge(df,tmp,on=field,how="left")

        tmp=df[[field,'transcript_id',"transcript_id norm. mean -log10(pValue)"]].drop_duplicates()
        tmp=tmp.drop(["transcript_id"],axis=1)
        tmp["transcript_id norm. mean -log10(pValue)"]=tmp["transcript_id norm. mean -log10(pValue)"].astype(float)
        tmp.columns=[field,"%s norm. mean -log10(pValue)" %str(field)]
        tmp=tmp.groupby(field, as_index=False).sum()
        df=pd.merge(df,tmp,on=field,how="left")



        tmp=df[[field,'transcript_id',"transcript_id signalValue"]].drop_duplicates()
        tmp=tmp.drop(["transcript_id"],axis=1)
        tmp["transcript_id signalValue"]=tmp["transcript_id signalValue"].astype(float)
        tmp.columns=[field,"%s signalValue" %str(field)]
        tmp=tmp.groupby(field, as_index=False).sum()
        df=pd.merge(df,tmp,on=field,how="left")

        return df

    dfTargets=getCounts_GeneIDs(dfTargets)


    dfTargets=dfTargets.drop(["count"],axis=1)
    return dfTargets
