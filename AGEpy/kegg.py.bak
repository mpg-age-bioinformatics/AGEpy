import pandas as pd
import sys
from urllib import urlopen
try:
    from rpy2.robjects.packages import importr
    try:
        biomaRt = importr("biomaRt")
    except:
        print "rpy2 could be loaded but 'biomaRt' could not be found.\nIf you want to use 'biomaRt' related functions please install 'biomaRt' in R.\n\n$ R\n> source('http://bioconductor.org/biocLite.R')\n> biocLite()\n> biocLite('biomaRt')\n> quit()"
        sys.stdout.flush()
except:
    print "Failed to import rpy2 module.\nPlease make sure you are using the same version of R you had when AGEpy was installed."
    sys.stdout.flush()




def organismsKEGG():
    """
    Lists all organisms present in the KEGG database.

    :returns: a list of lists containing one organism per list.

    """
    organisms=urlopen("http://rest.kegg.jp/list/organism").read()
    organisms=organisms.split("\n")
    #for o in organisms:
    #    print o
    #    sys.stdout.flush()
    return organisms


def databasesKEGG(organism,ens_ids):
    """
    Finds KEGG database identifiers for a respective organism given example ensembl ids.


    :param organism: an organism as listed in organismsKEGG()
    :param ens_ids: a list of ensenbl ids of the respective organism

    :returns: nothing if no database was found, or a string if a database was found

    """
    all_genes=urlopen("http://rest.kegg.jp/list/"+organism).read()
    all_genes=all_genes.split("\n")
    dbs=[]
    while len(dbs) == 0:
        for g in all_genes:
            if len(dbs) == 0:
                kid = g.split("\t")[0]
                gene=urlopen("http://rest.kegg.jp/get/"+kid).read()
                DBLINKS=gene.split("\n")
                DBLINKS=[ s for s in DBLINKS if ":" in s ]
                for d in DBLINKS:
                    test=d.split(" ")
                    test=test[len(test)-1]
                    if test in ens_ids:
                        DBLINK=[ s for s in DBLINKS if test in s ]
                        DBLINK=DBLINK[0].split(":")
                        DBLINK=DBLINK[len(DBLINK)-2]
                        dbs.append(DBLINK)
            else:
                break
    ens_db=dbs[0].split(" ")
    ens_db=ens_db[len(ens_db)-1]
    test_db=urlopen("http://rest.genome.jp/link/"+ens_db+"/"+organism).read()
    test_db=test_db.split("\n")
    if len(test_db) == 1:
        print "For "+organism+" the following db was found: "+ens_db
        print "This database does not seem to be valid KEGG-linked database identifier"
        print "For \n'hsa' use 'ensembl-hsa'\n'mmu' use 'ensembl-mmu'\n'cel' use 'EnsemblGenomes-Gn'\n'dme' use 'FlyBase'"
        sys.stdout.flush()
        ens_db = None
    else:
        print "For "+organism+" the following db was found: "+ens_db
        sys.stdout.flush()
    return ens_db


def ensembl_to_kegg(organism,kegg_db):
    """
    Looks up KEGG mappings of KEGG ids to ensembl ids

    :param organism: an organisms as listed in organismsKEGG()
    :param kegg_db: a matching KEGG db as reported in databasesKEGG

    :returns: a Pandas dataframe of with 'KEGGid' and 'ENSid'.

    """
    print "KEGG API: http://rest.genome.jp/link/"+kegg_db+"/"+organism
    sys.stdout.flush()
    kegg_ens=urlopen("http://rest.genome.jp/link/"+kegg_db+"/"+organism).read()
    kegg_ens=kegg_ens.split("\n")
    final=[]
    for i in kegg_ens:
        final.append(i.split("\t"))
    df=pd.DataFrame(final[0:len(final)-1])[[0,1]]
    ens_id=pd.DataFrame(df[1].str.split(":").tolist())[1]
    df=pd.concat([df,ens_id],axis=1)
    df.columns=['KEGGid','ensDB','ENSid']
    df=df[['KEGGid','ENSid']]
    return df

def ecs_idsKEGG(organism):
    """
    Uses KEGG to retrieve all ids and respective ecs for a given KEGG organism

    :param organism: an organisms as listed in organismsKEGG()

    :returns: a Pandas dataframe of with 'ec' and 'KEGGid'.

    """
    kegg_ec=urlopen("http://rest.kegg.jp/link/"+organism+"/enzyme").read()
    kegg_ec=kegg_ec.split("\n")
    final=[]
    for k in kegg_ec:
        final.append(k.split("\t"))
    df=pd.DataFrame(final[0:len(final)-1])[[0,1]]
    df.columns=['ec','KEGGid']
    return df


def idsKEGG(organism):
    """
    Uses KEGG to retrieve all ids for a given KEGG organism

    :param organism: an organism as listed in organismsKEGG()

    :returns: a Pandas dataframe of with 'gene_name' and 'KEGGid'.

    """
    ORG=urlopen("http://rest.kegg.jp/list/"+organism).read()
    ORG=ORG.split("\n")
    final=[]
    for k in ORG:
        final.append(k.split("\t"))
    df=pd.DataFrame(final[0:len(final)-1])[[0,1]]
    df.columns=['KEGGid','description']
    field = pd.DataFrame(df['description'].str.split(';',1).tolist())[0]
    field = pd.DataFrame(field)
    df = pd.concat([df[['KEGGid']],field],axis=1)
    df.columns=['KEGGid','gene_name']
    df=df[['gene_name','KEGGid']]
    return df

def pathwaysKEGG(organism):
    """
    Retrieves all pathways for a given organism.

    :param organism: an organism as listed in organismsKEGG()

    :returns df: a Pandas dataframe with the columns 'KEGGid','pathIDs', and 'pathName'.
    :returns df_: a Pandas dataframe with a columns for 'KEGGid', and one column for each pathway with the corresponding gene ids below
    """

    #print "KEGG API: http://rest.kegg.jp/list/pathway/"+organism
    #sys.stdout.flush()
    kegg_paths=urlopen("http://rest.kegg.jp/list/pathway/"+organism).read()
    kegg_paths=kegg_paths.split("\n")
    final=[]
    for k in kegg_paths:
        final.append(k.split("\t"))
    df=pd.DataFrame(final[0:len(final)-1])[[0,1]]
    df.columns=['pathID','pathName']
    print "Collecting genes for " + str(len(df)) + " pathways"
    sys.stdout.flush()
    df_pg=pd.DataFrame()
    for i in df['pathID'].tolist():
        print i+" ",
        sys.stdout.flush()
        path_genes=urlopen("http://rest.kegg.jp/link/genes/"+i).read()
        path_genes=path_genes.split("\n")
        final=[]
        for k in path_genes:
            final.append(k.split("\t"))
        final.pop()
        if len(final[0]) is not 2:
          print " skipped (no genes found in data)"
          continue
        print str(len(final)) + " gene(s) assigned"
        df_tmp=pd.DataFrame(final[0:len(final)])[[0,1]]
        df_tmp.columns=['pathID','KEGGid']
        df_pg=pd.concat([df_pg,df_tmp])
    df=pd.merge(df,df_pg,on=["pathID"], how="outer")
    KEGGids=list(set(df['KEGGid'].tolist()))
    print "Merging final table"
    sys.stdout.flush()
    df_final=pd.DataFrame()
    for k in KEGGids:
        df_tmp=df[df['KEGGid']==k]
        pathIDs=", ".join(df_tmp['pathID'].tolist())
        pathName=", ".join(df_tmp['pathName'].tolist())
        d={'KEGGid':k,'pathIDs':pathIDs,'pathName':pathName}
        df_tmp=pd.DataFrame(d,index=[0])
        df_final=pd.concat([df_final,df_tmp])
    df_final=df_final.dropna()

    df_=df[['KEGGid']].drop_duplicates()
    for i in list(set(df['pathID'].tolist())):
        tmp=df[df['pathID']==i]
        N=": ".join([i,tmp["pathName"].tolist()[0]])
        tmp=tmp[['KEGGid']]
        tmp.columns=[N]
        df_=pd.merge(df_,tmp,left_on=['KEGGid'],right_on=[N],how="outer")

    return df_final, df_


def biomaRtTOkegg(df):
    """
    Transforms a pandas dataframe with the columns 'ensembl_gene_id','kegg_enzyme'
    to dataframe ready for use in ...

    :param df: a pandas dataframe with the following columns: 'ensembl_gene_id','kegg_enzyme'

    :returns: a pandas dataframe with the following columns: 'ensembl_gene_id','kegg_enzyme'
    """
    df=df.dropna()
    ECcols=df.columns.tolist()
    df.reset_index(inplace=True,drop=True)
    # field = ECsb[['kegg_enzyme']]
    field = pd.DataFrame(df['kegg_enzyme'].str.split('+',1).tolist())[1]
    field = pd.DataFrame(field)
    df=pd.concat([df[['ensembl_gene_id']],field],axis=1)
    df.columns=ECcols
    df.drop_duplicates(inplace=True)
    df.reset_index(inplace=True,drop=True)
    plus=df['kegg_enzyme'].tolist()
    plus=[ s for s in plus if "+" in s ]
    noPlus=df[~df['kegg_enzyme'].isin(plus)]
    plus=df[df['kegg_enzyme'].isin(plus)]
    noPlus.reset_index(inplace=True, drop=True)
    plus.reset_index(inplace=True, drop=True)
    for p in range(0,len(plus)):
        enz=plus.ix[p]['kegg_enzyme']
        enz=enz.split("+")
        enz=pd.DataFrame(enz)
        enz.colums=['kegg_enzyme']
        enz['ensembl_gene_id']=plus.ix[p]['kegg_enzyme']
        noPlus=pd.concat([noPlus,enz])
    noPlus=noPlus.drop_duplicates()
    noPlus=noPlus[['ensembl_gene_id','kegg_enzyme']]
    noPlus['fake']='ec:'
    noPlus['kegg_enzyme']=noPlus['fake']+noPlus['kegg_enzyme']
    noPlus=noPlus[['ensembl_gene_id','kegg_enzyme']]

    return noPlus



def expKEGG(organism,names_KEGGids):
    """
    Gets all KEGG pathways for an organism

    :param organism: an organism as listed in organismsKEGG()
    :param names_KEGGids: a Pandas dataframe with the columns 'gene_name': and  'KEGGid' as reported from idsKEGG(organism) (or a subset of it).

    :returns df: a Pandas dataframe with 'KEGGid','pathID(1):pathNAME(1)', 'pathID(n):pathNAME(n)'
    :returns paths: a list of retrieved KEGG pathways
    """
    #print "KEGG API: http://rest.kegg.jp/list/pathway/"+organism
    #sys.stdout.flush()
    kegg_paths=urlopen("http://rest.kegg.jp/list/pathway/"+organism).read()
    kegg_paths=kegg_paths.split("\n")
    final=[]
    for k in kegg_paths:
        final.append(k.split("\t"))
    df=pd.DataFrame(final[0:len(final)-1])[[0,1]]

    df.columns=['pathID','pathName']
    print "Collecting genes for pathways"
    sys.stdout.flush()
    df_pg=pd.DataFrame()
    for i in df['pathID'].tolist():
        print i
        sys.stdout.flush()
        path_genes=urlopen("http://rest.kegg.jp/link/genes/"+i).read()
        path_genes=path_genes.split("\n")
        final=[]
        for k in path_genes:
            final.append(k.split("\t"))
        if len(final[0]) > 1:
            df_tmp=pd.DataFrame(final[0:len(final)-1])[[0,1]]
            df_tmp.columns=['pathID','KEGGid']
            df_pg=pd.concat([df_pg,df_tmp])
    df=pd.merge(df,df_pg,on=["pathID"], how="outer")
    df=df[df['KEGGid'].isin(names_KEGGids['KEGGid'].tolist())]
    df=pd.merge(df,names_KEGGids,how='left',on=['KEGGid'])
    df_fA=pd.DataFrame(columns=['KEGGid'])
    paths=[]
    for k in df[['pathID']].drop_duplicates()['pathID'].tolist():
        df_tmp=df[df['pathID']==k]
        pathName=df_tmp['pathName'].tolist()[0]
        pathName=" : ".join([k,pathName])
        keggIDs_in_path=df_tmp[['KEGGid']].drop_duplicates()['KEGGid'].tolist()
        a={pathName:keggIDs_in_path}
        a=pd.DataFrame(a,index=range(len(keggIDs_in_path)))
        a['KEGGid']=a[pathName].copy()
        df_fA=pd.merge(df_fA,a,how='outer',on=['KEGGid'])
        paths.append(pathName)

    return df_fA, paths

rbiomart_host="www.ensembl.org"

def KEGGmatrix(organism, dataset, database, query_attributes=['ensembl_gene_id','kegg_enzyme'], host=rbiomart_host,links=True,dfexp=None,kegg_db=None ):
    """
    This looks for all KEGG annotatios of an organism in biomaRt and the respective pathways in KEGG.

    :param dfexp: a Pandas dataframe with the following columns: 'ensembl_gene_id', 'log2FC'
    :param organism: a KEGG organism identifier
    :param dataset: a biomaRt dataset
    :param database: a biomaRt database
    :param query_attributes: biomaRt query attributes, the name can change but the output should stay in the same order ie. 'ensembl_gene_id','kegg_enzyme'
    :param host: biomaRt_host
    :param links: if True, returns df_links
    :param dfexp: a Pandas dataframe with the folowing columns 'KEGGid' and 'log2FC'
    :param kegg_db: a KEGG database as recovered by the databasesKEGG function


    :returns df: a Pandas dataframe with the 'KEGGid','pathsIDs','pathName','ensembl_gene_id','kegg_enzyme'
    :returns df_: a matrix with a column for each KEGG pathway for a given organism and the expression values in the respective dfexp in parameter
    :returns fullmatrix: a matrix with a column for each KEGG pathway for a given organism
    :returns df_links: a dataframe with links for each pathway and the links in the dfexp highlighted red (if df_links.

    """
    try:
        # Get all ensembl gene ids and keeg enzyme labels from biomaRt
        biomaRt = importr("biomaRt")
        ensemblMart=biomaRt.useMart(database, host=host)
        ensembl=biomaRt.useDataset(dataset, mart=ensemblMart)
        biomaRt_output=biomaRt.getBM(attributes=query_attributes,mart=ensembl)
        biomaRt_output = [tuple([biomaRt_output[j][i] for j in range(biomaRt_output.ncol)]) for i in range(biomaRt_output.nrow)]
        biomaRt_output = pd.DataFrame(biomaRt_output)
        biomaRt_output.columns = ['ensembl_gene_id','kegg_enzyme']
        biomaRt_output = biomaRt_output[biomaRt_output['kegg_enzyme']!='']
        biomaRt_output.reset_index(inplace=True,drop=True)
        biomaRt_output=biomaRtTOkegg(biomaRt_output)
    except:
        # Do it wiht KEGG
        ec_KEGGid=ecs_idsKEGG(organism)
        KEGGid_ENSid=ensembl_to_kegg(organism,kegg_db)
        biomaRt_output=pd.merge(ec_KEGGid,KEGGid_ENSid,on=['KEGGid'],how="outer")
        biomaRt_output=biomaRt_output.drop(['KEGGid'],axis=1)
        biomaRt_output=biomaRt_output[['ENSid','ec']]
        biomaRt_output.columns=['ensembl_gene_id','kegg_enzyme']


    # Gett all pathways
    df, df_=pathwaysKEGG(organism)
    fullmatrix=df_.copy()

    # Get all KEGG ecs from KEGG
    ecs=ecs_idsKEGG(organism)

    biomaRt_output=pd.merge(biomaRt_output,ecs,left_on=['kegg_enzyme'],right_on=['ec'],how="outer")
    biomaRt_output=biomaRt_output.drop(['ec'],axis=1)

    df=pd.merge(df, biomaRt_output, how="outer",on=['KEGGid']).drop_duplicates()
    df=df[df['KEGGid'].astype(str)!="nan"]
    df=df.sort('ensembl_gene_id')
    df=df.drop_duplicates(subset=['KEGGid','pathIDs','pathName','kegg_enzyme' ])
    df=df.reset_index()
    if not isinstance(dfexp, pd.DataFrame):
        print "Returning df and fullmatrix"
        sys.stdout.flush()
        return df, fullmatrix
    else:
        expDic=dfexp[['ensembl_gene_id','log2FC']].dropna()
        expDic=expDic.set_index(['ensembl_gene_id'])
        expDic=expDic.to_dict().get("log2FC")

        dfexp=pd.merge(biomaRt_output, dfexp, how="right",on=['ensembl_gene_id'])

        #expDic=dfexp[['KEGGid','log2FC']].dropna()
        #expDic=expDic.set_index(['KEGGid'])
        #expDic=expDic.to_dict().get("log2FC")

        cols=df_.columns.tolist()
        cols=[ s for s in cols if "KEGGid" not in s ]
        #for c in cols:
        #    df_[c]=df_[c].apply(lambda x: expDic.get(x))

        df_=pd.merge(dfexp,df_,on=['KEGGid'],how="left")
        df_=df_.dropna(subset=['KEGGid','ensembl_gene_id'])
        df_=df_.sort(columns=cols)

        def get_expression(df_,col,expDic=expDic):
            v=df_[col]
            if str(v) != "nan":
                eID=df_['ensembl_gene_id']
                ex=expDic.get(eID)
            else:
                ex=v
            return ex

        for c in cols:
            df_[c]=df_.apply(get_expression, args=(c,),axis=1)

        if links==True:
            df_links=pd.DataFrame()
            for p in cols:
                dfT=df_.dropna(subset=[p])
                dfT=dfT.dropna(subset=['kegg_enzyme'])
                print dfT.head()
                dfT=dfT.drop_duplicates(subset=['kegg_enzyme'])
                if len(dfT) > 0:
                    pathway=p.split(":")[1]
                    URL="http://www.kegg.jp/kegg-bin/show_pathway?"+pathway
                    for i in dfT['kegg_enzyme'].tolist():
                        gKEGG=i
                        color="red"
                        text="/"+gKEGG+"%09"+color
                        URL=URL+text
                    print URL
                    d={"pathway":pathway, "URL":URL}
                    d=pd.DataFrame(d,index=[0])
                    df_links=pd.concat([df_links,d])
            df_links.reset_index(inplace=True, drop=True)

            return df, df_, fullmatrix, df_links

        else:
            return df, df_, fullmatrix
