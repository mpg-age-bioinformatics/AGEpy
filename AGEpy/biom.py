import biomart
import sys
import pandas as pd
import numpy as np
from biomart import BiomartServer
#from cStringIO import StringIO # python2
from io import BytesIO as cStringIO
from io import StringIO

biomart_host="http://www.ensembl.org/biomart"

def datasetsBM(host=biomart_host):
    """
    Lists BioMart datasets.

    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    stdout_ = sys.stdout #Keep track of the previous value.
    stream = StringIO()
    sys.stdout = stream   
    server = BiomartServer(biomart_host)
    server.show_datasets()
    sys.stdout = stdout_ # restore the previous stdout.
    variable = stream.getvalue() 
    v=variable.replace("{"," ") 
    v=v.replace("}"," ") 
    v=v.replace(": ","\t")
    print(v)

def filtersBM(dataset,host=biomart_host):
    """
    Lists BioMart filters for a specific dataset.

    :param dataset: dataset to list filters of.
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    stdout_ = sys.stdout #Keep track of the previous value.
    stream = StringIO()
    sys.stdout = stream   
    server = BiomartServer(host)
    d=server.datasets[dataset]
    d.show_filters()
    sys.stdout = stdout_ # restore the previous stdout.
    variable = stream.getvalue() 
    v=variable.replace("{"," ") 
    v=v.replace("}"," ") 
    v=v.replace(": ","\t")
    print(v)

def attributesBM(dataset,host=biomart_host):
    """
    Lists BioMart attributes for a specific dataset.

    :param dataset: dataset to list attributes of.
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    stdout_ = sys.stdout #Keep track of the previous value.
    stream = StringIO()
    sys.stdout = stream   
    server = BiomartServer(host)
    d=server.datasets[dataset]
    d.show_attributes()
    sys.stdout = stdout_ # restore the previous stdout.
    variable = stream.getvalue() 
    v=variable.replace("{"," ") 
    v=v.replace("}"," ") 
    v=v.replace(": ","\t")
    print(v)

def queryBM(query_attributes,query_dataset,query_filter=None,query_items=None,query_dic=None,host=biomart_host):
    """
    Queries BioMart.

    :param query_attributes: list of attributes to recover from BioMart
    :param query_dataset: dataset to query
    :param query_filter: one BioMart filter associated with the items being queried
    :param query_items: list of items to be queried (must assoiate with given filter)
    :param query_dic: for complex queries this option should be used instead of 'filters' and 'items' and a dictionary of filters provided here eg. querydic={"filter1":["item1","item2"],"filter2":["item3","item4"]}. If using querydic, don't query more than 350 items at once.
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: a Pandas dataframe of the queried attributes

    """
    server = BiomartServer(host)
    d=server.datasets[query_dataset]
    res=[]

    if not query_dic:
        if query_items:
            chunks=[query_items[x:x+350] for x in xrange(0, len(query_items), 350)]
            for c in chunks:
                response=d.search({'filters':{query_filter:c},'attributes':query_attributes})
                for line in response.iter_lines():
                    line = line.decode('utf-8')
                    res.append(line.split("\t"))
        else:
            response=d.search({'attributes':query_attributes})
            for line in response.iter_lines():
                line = line.decode('utf-8')
                res.append(line.split("\t"))
    
    elif query_dic:
        response=d.search({'filters':query_dic,'attributes':query_attributes})
        for line in response.iter_lines():
            line = line.decode('utf-8')
            res.append(line.split("\t"))
    res=pd.DataFrame(res)
    res.columns=query_attributes
    return(res)

def FilterGOstring(names_filter=["age-", "aging", "aged", 'aging', 'aging.', 'aging,'],\
                   exclude_names=["packaging","voltage","cleavage-",\
                       "stage-1","cage-like","message-specific",\
                       "damage-associated","stage-specific","foraging",\
                       "DNA-damaging","engaging","damaged","packaged"],\
                   defs_filter=[" age-", " aging", " aged", ' aging', ' aging.', ' aging,'],\
                   exclude_defs=["packaging","voltage","cleavage-",\
                         "stage-1","cage-like","message-specific",\
                         "damage-associated","stage-specific","foraging",\
                         "DNA-damaging","engaging","damaged","packaged"],\
                   host=biomart_host,\
                   HSA=None,MUS=None,CEL=None,DMEL=None):
    
    """
    Filters GO terms based on given strings using ENSEMBL's biomart homology mapping.
    
    :param names_filter: list of substrings to filter GO names on
    :param exclude_names: list of substrings to be used for exclusion of GO names
    :param defs_filter: list of substrings to filter GO defenitions on
    :param exclude_defs: list of substrings to be used for exclustion of GO defenitions
    :param host: biomart host server, default="http://www.ensembl.org/biomart"
    :param HSA: retrieved hsa dataframe
    :param MUS: retrieved mus dataframe
    :param CEL: retrieved cel dataframe
    :param DMEL: retrieved dmel dataframe
    
    :returns homology_df, HSA, MUS, CEL, DMEL   
 
    
    """
    
    if type(HSA) == type(None):
        queries={'hsapiens_gene_ensembl':["ensembl_gene_id","external_gene_name", \
                                          "go_id","name_1006","definition_1006"],\
                 "mmusculus_gene_ensembl":["ensembl_gene_id","external_gene_name", \
                                           "go_id","name_1006","definition_1006"],\
                 "celegans_gene_ensembl":["ensembl_gene_id","external_gene_name", \
                                          "go_id","name_1006","definition_1006"],\
                 "dmelanogaster_gene_ensembl":["ensembl_gene_id","external_gene_name", \
                                               "go_id","name_1006","definition_1006"]}

        def QueryBioMart(dataset,attributes,host=host):
            #print dataset
            #sys.stdout.flush()
            server = BiomartServer( host )
            organism=server.datasets[dataset]
            response=organism.search({'attributes':attributes})
            response=response.content
            response=response.decode()
            response=response.split("\n")
            response=[s.split("\t") for s in response ]
            response=pd.DataFrame(response,columns=attributes)
            return response

        homology=[ "ensembl_gene_id","celegans_homolog_ensembl_gene","dmelanogaster_homolog_ensembl_gene","mmusculus_homolog_ensembl_gene"]
        hsa_homology=QueryBioMart('hsapiens_gene_ensembl',homology)

        HSA=QueryBioMart('hsapiens_gene_ensembl',queries['hsapiens_gene_ensembl'])
        MUS=QueryBioMart('mmusculus_gene_ensembl',queries['mmusculus_gene_ensembl'])
        CEL=QueryBioMart('celegans_gene_ensembl',queries['celegans_gene_ensembl'])
        DMEL=QueryBioMart('dmelanogaster_gene_ensembl',queries['dmelanogaster_gene_ensembl'])    

        HSA=pd.merge(HSA,hsa_homology,on=["ensembl_gene_id"],how="outer")
        HSA.columns=['HSA_ensembl_gene_id', 'HSA_external_gene_name','HSA_go_id', 'HSA_name_1006', 'HSA_definition_1006',\
                     'CEL_ensembl_gene_id', 'DMEL_ensembl_gene_id', 'MUS_ensembl_gene_id']
        MUS.columns=['MUS_ensembl_gene_id','MUS_external_gene_name',"MUS_go_id",'MUS_name_1006','MUS_definition_1006']
        CEL.columns=['CEL_ensembl_gene_id','CEL_external_gene_name',"CEL_go_id",'CEL_name_1006','CEL_definition_1006']
        DMEL.columns=['DMEL_ensembl_gene_id','DMEL_external_gene_name',"DMEL_go_id",'DMEL_name_1006','DMEL_definition_1006']
        
    HSA_gos=HSA[['HSA_go_id','HSA_name_1006','HSA_definition_1006']]
    MUS_gos=MUS[['MUS_go_id','MUS_name_1006','MUS_definition_1006']]
    CEL_gos=CEL[['CEL_go_id','CEL_name_1006','CEL_definition_1006']]
    DMEL_gos=DMEL[['DMEL_go_id','DMEL_name_1006','DMEL_definition_1006']]

    GOS=pd.DataFrame()
    for tmp in [ HSA_gos, MUS_gos, CEL_gos, DMEL_gos]:
        tmp.columns=['go_id','name_1006','definition_1006']
        GOS=pd.concat([GOS,tmp])
    GOS=GOS.drop_duplicates()
    GOS=GOS.dropna()
    GOS.reset_index(inplace=True, drop=True)

    names=GOS["name_1006"].tolist()
    defs=GOS["definition_1006"].tolist()

    filtered_names=[]
    for age in names_filter:
        tmp=list(set( [ s for s in names if age in s ] ))
        exclude=exclude_names
        for e in exclude:
            tmp=[ s for s in tmp if e not in s ]
        filtered_names.append(tmp)   

    filtered_defs=[]
    for age in defs_filter:
        tmp=list(set( [ s for s in defs if age in s ] ))
        exclude=exclude_defs
        for e in exclude:
            tmp=[ s for s in tmp if e not in s ]
        filtered_defs.append(tmp)

    # def checkStrings(filtered_names,names_filter):
    #     print_names=" ".join(filtered_names)
    #     print_names=print_names.split(" ")
    #     print_names_=[]
    #     for age in names_filter:
    #         tmp=[s for s in print_names if age in s ]
    #         print_names_.append(tmp)
    #     print_names_=[item for sublist in print_names_ for item in sublist]
    #     print_names_=list(set(print_names_))
    #     return print_names_

    filtered_names = [item for sublist in filtered_names for item in sublist]
    filtered_defs = [item for sublist in filtered_defs for item in sublist]

    # names_=checkStrings(filtered_names,names_filter)
    # defs_=checkStrings(filtered_defs,defs_filter)

    print("\nStrings being used for filtering in names section:")
    for f in filtered_names:
        print(f)

    print("\nStrings being used for filtering in defenitions section:")
    for f in filtered_defs:
        print("\n"+f)

    def CHECK_AGE(x,l):
        if x in l:
            res=x
        else:
            res=np.nan
        return res

    GOS["names_filter"]=GOS["name_1006"].apply(lambda x: CHECK_AGE(x,filtered_names) )
    GOS["defs_filter"]=GOS["definition_1006"].apply(lambda x: CHECK_AGE(x,filtered_defs) )
    AGEGOS=GOS[['go_id',"names_filter","defs_filter"]].dropna(thresh=2)
    AGEGOS=list(set(AGEGOS["go_id"].tolist()))

    def CombineAnnHSA(df):
         return pd.Series(dict(HSA_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['HSA_ensembl_gene_id']))  if str(s) != "nan" ] ) ,\
                               HSA_go_id = ', '.join([ str(s) for s in list(set(df['HSA_go_id'])) if str(s) != "nan" ]), \
                               HSA_external_gene_name = ', '.join([ str(s) for s in list(set(df['HSA_external_gene_name'])) if str(s) != "nan" ] ) ,\
                               HSA_name_1006 = ', '.join([ str(s) for s in list(set(df['HSA_name_1006'])) if str(s) != "nan" ] ) ,\
                               HSA_definition_1006 = ', '.join([ str(s) for s in list(set(df['HSA_definition_1006'])) if str(s) != "nan" ] ) ,\
                               CEL_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['CEL_ensembl_gene_id'])) if str(s) != "nan" ] ) ,\
                               DMEL_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['DMEL_ensembl_gene_id'])) if str(s) != "nan"] ) ,\
                               MUS_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['MUS_ensembl_gene_id'])) if str(s) != "nan" ] ) ,\
                              ) ) 

    def CombineAnnMUS(df):
         return pd.Series(dict(MUS_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['MUS_ensembl_gene_id']))  if str(s) != "nan" ] ) ,\
                               MUS_external_gene_name = ', '.join([ str(s) for s in list(set(df['MUS_external_gene_name'])) if str(s) != "nan" ]), \
                               MUS_go_id = ', '.join([ str(s) for s in list(set(df['MUS_go_id'])) if str(s) != "nan" ] ) ,\
                               MUS_name_1006 = ', '.join([ str(s) for s in list(set(df['MUS_name_1006'])) if str(s) != "nan" ] ) ,\
                               MUS_definition_1006 = ', '.join([ str(s) for s in list(set(df['MUS_definition_1006'])) if str(s) != "nan" ] ) ,\
                              ) ) 

    def CombineAnnCEL(df):
         return pd.Series(dict(CEL_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['CEL_ensembl_gene_id']))  if str(s) != "nan" ] ) ,\
                               CEL_external_gene_name = ', '.join([ str(s) for s in list(set(df['CEL_external_gene_name'])) if str(s) != "nan" ]), \
                               CEL_go_id = ', '.join([ str(s) for s in list(set(df['CEL_go_id'])) if str(s) != "nan" ] ) ,\
                               CEL_name_1006 = ', '.join([ str(s) for s in list(set(df['CEL_name_1006'])) if str(s) != "nan" ] ) ,\
                               CEL_definition_1006 = ', '.join([ str(s) for s in list(set(df['CEL_definition_1006'])) if str(s) != "nan" ] ) ,\
                              ) ) 
    def CombineAnnDMEL(df):
         return pd.Series(dict(DMEL_ensembl_gene_id = ', '.join([ str(s) for s in list(set(df['DMEL_ensembl_gene_id']))  if str(s) != "nan" ] ) ,\
                               DMEL_external_gene_name = ', '.join([ str(s) for s in list(set(df['DMEL_external_gene_name'])) if str(s) != "nan" ]), \
                               DMEL_go_id = ', '.join([ str(s) for s in list(set(df['DMEL_go_id'])) if str(s) != "nan" ] ) ,\
                               DMEL_name_1006 = ', '.join([ str(s) for s in list(set(df['DMEL_name_1006'])) if str(s) != "nan" ] ) ,\
                               DMEL_definition_1006 = ', '.join([ str(s) for s in list(set(df['DMEL_definition_1006'])) if str(s) != "nan" ] ) ,\
                              ) )  

    HSA=HSA.groupby(by=["HSA_ensembl_gene_id","MUS_ensembl_gene_id",\
                        "DMEL_ensembl_gene_id","CEL_ensembl_gene_id"], as_index=False).apply(CombineAnnHSA)
    MUS=MUS.groupby(by="MUS_ensembl_gene_id", as_index=False).apply(CombineAnnMUS)
    DMEL=DMEL.groupby(by="DMEL_ensembl_gene_id", as_index=False).apply(CombineAnnDMEL)
    CEL=CEL.groupby(by="CEL_ensembl_gene_id", as_index=False).apply(CombineAnnCEL)

    MUS.reset_index(inplace=True,drop=True)
    HSA.reset_index(inplace=True,drop=True)
    DMEL.reset_index(inplace=True,drop=True)
    CEL.reset_index(inplace=True,drop=True)
    
    HSA=HSA[['HSA_ensembl_gene_id', 'HSA_external_gene_name', 'HSA_go_id', 'HSA_name_1006','HSA_definition_1006',\
             'MUS_ensembl_gene_id', 'CEL_ensembl_gene_id', 'DMEL_ensembl_gene_id']]
    MUS=MUS[['MUS_ensembl_gene_id', 'MUS_external_gene_name', 'MUS_go_id', 'MUS_name_1006','MUS_definition_1006']]
    CEL=CEL[['CEL_ensembl_gene_id', 'CEL_external_gene_name', 'CEL_go_id', 'CEL_name_1006','CEL_definition_1006']]
    DMEL=DMEL[['DMEL_ensembl_gene_id', 'DMEL_external_gene_name', 'DMEL_go_id', 'DMEL_name_1006','DMEL_definition_1006']]

    homDF=pd.merge(HSA,MUS,on=["MUS_ensembl_gene_id"], how="outer")
    homDF=pd.merge(homDF,CEL,on=["CEL_ensembl_gene_id"], how="outer")
    homDF=pd.merge(homDF,DMEL,on=["DMEL_ensembl_gene_id"], how="outer")

    def AGE_ORG(df,AGEGOS=AGEGOS):
        HSA_go_id=str(df["HSA_go_id"]).split(", ")
        MUS_go_id=str(df["MUS_go_id"]).split(", ")
        CEL_go_id=str(df["CEL_go_id"]).split(", ")
        DMEL_go_id=str(df["DMEL_go_id"]).split(", ")

        res=[]

        for n,l in zip(["HSA","MUS","CEL","DMEL"],\
                       [HSA_go_id,MUS_go_id,CEL_go_id,DMEL_go_id]):
            r=[ s for s in l if s in AGEGOS ]
            if len(r) > 0:
                res.append(n)
        if len(res) > 0:
            res.sort()
            res=", ".join(res)
        else:
            res=np.nan
        return res

    homDF.reset_index(inplace=True, drop=True)
    homDF["evidence"]=homDF.apply(AGE_ORG, axis=1)

    return homDF,HSA,MUS,CEL,DMEL
