import pandas as pd
import numpy as np

def readSAM(SAMfile,header=False):
    """
    Reads and parses a sam file.

    :param SAMfile: /path/to/file.sam
    :param header: logical, if True, reads the header information

    :returns: a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' and a list of the headers if header=True

    """
    if header==True:
        f=open(SAMfile,"r+")
        head=[]
        for line in f.readlines():
            if line[0]=="@":
                head.append(line)
            else:
                continue
        f.close()

    sam=pd.read_table(SAMfile,sep="this_gives_one_column",comment="@",header=None)
    sam=pd.DataFrame(sam[0].str.split("\t").tolist())
    acols=[0,1,2,3,4,5,6,7,8,9]
    sam_=sam[acols]
    samcols=sam.columns.tolist()
    bcols=[ s for s in samcols if s not in acols ]
    sam_[10]=sam[bcols[0]]
    if len(bcols) > 1:
        for c in bcols[1:]:
            sam_[10]=sam_[10].astype(str)
            sam[c]=sam[c].astype(str)
            sam_[10]=sam_[10]+"\t"+sam[c]

    sam_.columns=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL']

    if header==True:
        return sam_, head
    else:
        return sam_

def writeSAM(sam,SAMfile,header=None):
    """
    Writes a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' into a sam file

    :param sam: pandas dataframe to be writen
    :param SAMfile: /path/to/file.sam

    :returns: nothing
    """
    def toNone(x):
        if x=="None":
            x=np.nan
        return x

    sam.reset_index(inplace=True,drop=True)
    QUAL=pd.DataFrame(sam['QUAL'].str.split("\t").tolist())
    cols=QUAL.columns.tolist()

    for c in cols:
        QUAL[c]=QUAL[c].apply(lambda x: toNone(x))

    sam=sam.drop(['QUAL'],axis=1)
    sam=pd.concat([sam,QUAL],axis=1)
    sam=sam.astype(str)
    sam=sam.as_matrix()

    tfile=open(SAMfile, "w+")

    if header != None:
        for l in header:
            tfile.write(l)

    for l in sam:
        l=[ s for s in l if s not in ['nan'] ]
        l="\t".join(l)
        tfile.write(l+"\n")

    tfile.close()

def SAMflags(x):
    """
    Explains a SAM flag.

    :param x: flag

    :returns: complete SAM flag explanaition
    """
    flags=[]

    if x & 1:
        l="1: Read paired"
    else:
        l="0: Read unpaired"
    flags.append(l)

    if x & 2 :
        l="1: Read mapped in proper pair"
    else:
        l="0: Read not mapped in proper pair"
    flags.append(l)

    if x & 4 :
        l="1: Read unmapped"
    else:
        l="0: Read mapped"
    flags.append(l)

    if x & 8 :
        l="1: Mate unmapped"
    else:
        l="0: Mate mapped"
    flags.append(l)

    if x & 16 :
        l="1: Read reverse strand"
    else:
        l="0: Read direct strand"
    flags.append(l)

    if x & 32 :
        l="1: Mate reverse strand"
    else:
        l="0: Mate direct strand"
    flags.append(l)

    if x & 64 :
        l="1: First in pair"
    else:
        l="0: Second in pair"
    flags.append(l)

    if x & 128 :
        l="1: Second in pair"
    else:
        l="0: First in pair"
    flags.append(l)

    if x & 256 :
        l="1: Not primary alignment"
    else:
        l="0: Primary alignment"
    flags.append(l)

    if x & 512 :
        l="1: Read fails platform/vendor quality checks"
    else:
        l="0: Read passes platform/vendor quality checks"
    flags.append(l)

    if x & 1024 :
        l="1: Read is PCR or optical duplicate"
    else:
        l="0: Read is not PCR or optical duplicate"
    flags.append(l)

    if x & 2048 :
        l="1: Supplementary alignment"
    else:
        l="0: Not supplementary alignment"
    flags.append(l)

    return flags
