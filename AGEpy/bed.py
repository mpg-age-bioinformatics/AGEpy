import pandas as pd
import numpy as np
import os
#from urllib import urlopen # python2 
#import urllib2 # python2
import urllib.request as urllib2
#import StringIO python2
from io import StringIO
import gzip
import pybedtools
from pybedtools import BedTool
from .gtf import GTFtoBED
from .gtf import readGTF
from .gtf import retrieve_GTF_field
import sys

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

def AnnotateBED(bed, GTF, genome_file, bedcols=None, promoter=[1000,200]):
    """
    Annotates a bed file.

    :param bed: either a /path/to/file.bed or a Pandas dataframe in bed format. /path/to/file.bed implies bedcols.
    :param GTF: /path/to/file.gtf
    :param genome_file: /path/to/file.genome - a tab separated values of chr name and size information
    :param bedcols: a comma separated string of column headers to use when reading in a bed file. eg: "chr,start,end,name"
    :param promoter: a list containing the upstream start of the promoter region from the TSS and the downstream end of the promoter region from the TSS.

    :returns: a Pandas dataframe with the annotated bed file. exons and promoters will be reported as well in the annotated_gene_features column.
    """
    if type(bed) == type("string"):
        bed=pd.read_table(bed,header=None)
        bed.columns=bedcols.split(",")

    print("Reading GTF file.")
    sys.stdout.flush()

    GTF=readGTF(GTF)
    GTF["gene_name"]=retrieve_GTF_field("gene_name", GTF)
    GTF["gene_id"]=retrieve_GTF_field("gene_id", GTF)
    GTF["gene_name"]=GTF["gene_name"]+"/"+GTF["gene_id"]
    GTF=GTF.drop(["gene_id"],axis=1)

    print("Generating promoters annotation.")
    sys.stdout.flush()

    promoters=GTF[GTF["feature"]=="transcript"]
    promoters_plus=promoters[promoters["strand"]=="+"]
    promoters_minus=promoters[promoters["strand"]=="-"]

    upstream=promoter[0]
    downstream=promoter[1]

    promoters_plus.loc[:,"promoter_start"]=promoters_plus.loc[:,"start"].astype(int)-upstream
    promoters_plus.loc[:,"promoter_end"]=promoters_plus.loc[:,"start"].astype(int)+downstream

    promoters_minus.loc[:,"promoter_start"]=promoters_minus["end"].astype(int)-downstream
    promoters_minus.loc[:,"promoter_end"]=promoters_minus["end"].astype(int)+upstream

    promoters=pd.concat([promoters_plus,promoters_minus])

    promoters=promoters[["seqname","feature","promoter_start","promoter_end","gene_name"]]
    promoters.columns=["seqname","feature","start","end","gene_name"]

    promoters.loc[:,"feature"]="promoter"
    promoters.drop_duplicates(inplace=True)
    promoters.reset_index(inplace=True, drop=True)

    chr_sizes=pd.read_table(genome_file,header=None)
    chr_sizes.columns=["seqname","size"]
    chr_sizes.loc[:,"seqname"]=chr_sizes["seqname"].astype(str)
    promoters.loc[:,"seqname"]=promoters["seqname"].astype(str)

    promoters=pd.merge(promoters,chr_sizes,how="left",on=["seqname"])
    def CorrectStart(df):
        s=df["start"]
        if s < 0:
            s=0
        return s

    def CorrectEnd(df):
        s=df["end"]
        e=df["size"]
        if s > e:
            s=e
        return s

    promoters.loc[:,"start"]=promoters.apply(CorrectStart,axis=1)
    promoters.loc[:,"end"]=promoters.apply(CorrectEnd,axis=1)

    promoters.drop(["size"],axis=1, inplace=True)

    GTFs=GTF[["seqname","feature","start","end","gene_name"]]
    GTFs=GTFs[ GTFs["feature"]!= "gene"]

    GTFs.drop_duplicates(inplace=True)
    GTFs.reset_index(inplace=True, drop=True)

    GTFs=pd.concat([GTFs,promoters])

    def NewName(df):
        name=df["gene_name"]
        feature=df["feature"]
        if feature == "transcript":
            res=name
        else:
            res=name+":"+feature
        return res

    GTFs.loc[:,"gene_name"]=GTFs.apply(NewName, axis=1)
    GTFs=GTFs[["seqname","start","end","gene_name"]]

    print( "Intersecting annotation tables and bed." )
    sys.stdout.flush()

    refGTF=dfTObedtool(GTFs)
    pos=dfTObedtool(bed)

    colsGTF=GTFs.columns.tolist()
    newCols=bed.columns.tolist()

    for f in colsGTF:
        newCols.append(f+"_")
    newCols_=[ s for s in newCols if s not in ["seqname_","start_", "end_"]]

    pos=pos.intersect(refGTF, loj=True)
    pos=pd.read_table(pos.fn , names=newCols)
    pos=pos[newCols_]

    print("Merging features.")
    sys.stdout.flush()

    def GetFeature(x):
        if ":" in x:
            res=x.split(":")[1]
        else:
            res=np.nan
        return res

    def GetName(x):
        if ":" in x:
            res=x.split(":")[0]
        elif type(x) == type("string"):
            if x != ".":
                res=x
            else:
                res=np.nan
        else:
            res=np.nan
        return res

    pos["gene_feature_"]=pos["gene_name_"].apply( lambda x: GetFeature(x) )
    pos["gene_name_"]=pos["gene_name_"].apply( lambda x: GetName(x) )

    refcol=pos.columns.tolist()
    refcol=[ s for s in refcol if s != "gene_feature_" ]

    def CombineAnn(df):
        def JOIN(x):
            return ', '.join([ str(s) for s in list(set(df[x]))  if str(s) != "nan" ] )
        return pd.Series(dict( gene_feature_ = JOIN("gene_feature_") ) )

    pos_=pos.groupby(refcol).apply(CombineAnn)
    pos_.reset_index(inplace=True, drop=False)

    def MergeNameFeatures(df):
        name=df["gene_name_"]
        feature=df["gene_feature_"]
        if (type(name) == type("string")) & (name != ".") :
            if type(feature) == type("string"):
                if len(feature) > 0:
                    res=name+": "+feature
                else:
                    res=name
            else:
                res=name
        else:
            res=np.nan
        return res

    pos_["annotated_gene_features"]=pos_.apply(MergeNameFeatures,axis=1)

    pos_=pos_.drop(["gene_name_","gene_feature_"],axis=1)

    def CombineAnn(df):
        def JOIN(x):
            return '; '.join([ str(s) for s in list(set(df[x]))  if str(s) != "nan" ] )
        return pd.Series(dict( annotated_gene_features = JOIN("annotated_gene_features") ) )

    refcol=[ s for s in refcol if s != "gene_name_" ]
    pos_=pos_.groupby(refcol).apply(CombineAnn)
    pos_.reset_index(inplace=True, drop=False)

    return pos_
