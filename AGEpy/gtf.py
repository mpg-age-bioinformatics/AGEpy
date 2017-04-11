import pandas as pd

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
