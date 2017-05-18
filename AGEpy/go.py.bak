import pandas as pd
import urllib2
import StringIO
import gzip


def getGeneAssociation(URL_or_file):
    """
    This function collects GO annotation from http://geneontology.org/page/download-annotations.

    :param URL_or_file: either a link to a file on geneontology.org eg. http://geneontology.org/gene-associations/gene_association.fb.gz or the path for the respective  downloded .gz file.
    :returns: a Pandas dataframe with the parsed table.
    """
    if URL_or_file[:4] == "http":
        response = urllib2.urlopen(URL_or_file)
        compressedFile = StringIO.StringIO(response.read())
        decompressedFile = gzip.GzipFile(fileobj=compressedFile)
    else:
        decompressedFile = gzip.GzipFile(URL_or_file)
    out=decompressedFile.read().split("\n")

    version=[ s for s in out if len(s) > 0 ]
    version=[ s for s in version if s[0] == '!' ]
    version=[ s for s in version if "!gaf-version:" in s ]
    version=version[0]


    if version=="!gaf-version: 2.0":
        reco=version
    else:
        reco=None

    out=[ s for s in out if len(s) > 0 ]
    out=[ s for s in out if s[0] != "!" ]
    out=[s.split("\t") for s in out]
    out=pd.DataFrame(out)
    mgi_cols=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier (this field is optional)","GO ID","DB:Reference","Evidence Code","Evidence Code Qualifier (optional)",\
     "Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon","Date","Assigned_by"]
    fb_cols=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO ID","DB:Reference","Evidence",\
     "With (or) From","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon","Date","Assigned_by","Annotation Extension",\
     "Gene Product Form ID"]
    gaf_20=["DB","DB Object ID","DB Object Symbol","Qualifier","GO ID","DB:Reference (|DB:Reference)","Evidence Code","With (or) From","Aspect","DB Object Name",\
     "DB Object Synonym (|Synonym)","DB Object Type","Taxon(|taxon)","Date","Assigned By","Annotation Extension","Gene Product Form ID"]
    cols={"fb":fb_cols,"wb":fb_cols,"mgi":fb_cols,"!gaf-version: 2.0":gaf_20}
    colsType=URL_or_file.split(".")
    colsType=colsType[len(colsType)-2]
    if colsType=="gaf":
        colsType=reco
    if colsType in cols.keys():
        try:
            cols=cols.get(colsType)
            out.columns=cols
        except ValueError as err:
            print "Could not fit headers."
            print err
            sys.stdout.flush()
    else:
        print "Could not find headers for %s." %colsType
        sys.stdout.flush()
    return out
