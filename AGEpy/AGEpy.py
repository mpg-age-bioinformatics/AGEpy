#!/usr/bin/env python

"""Bioinformatics tools developed at the Max Planck Institute for Biology of Ageing"""

import pandas as pd
import numpy as np
import os
import csv
import sys
from suds.client import Client as sudsclient
import os
import ssl
from biomart import BiomartServer
from urllib import urlopen
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec


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
    label=field
    field = pd.DataFrame(gtf['attribute'].str.split(field).tolist())[1]
    field = field.astype(str)
    field = pd.DataFrame(field.str.split(';',1).tolist())
    field = pd.DataFrame(field[0].str.split('"').tolist())[1]
    field = pd.DataFrame(field)
    field.columns=[label]
    return field

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
    df=inGTF.drop(['attribute'],axis=1)
    for d in desc:
        field=retrieve_GTF_field(d,inGTF)
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
    if name not in inGTF.columns.tolist():
        field=retrieve_GTF_field(name, inGTF)
        inGTF=pd.concat([inGTF,field],axis=1)
    bed=inGTF[['seqname','start','end',name,'score','strand']]
    bed.columns=['chrom','chromStart','chromEnd','name','score','strand']
    bed.drop_duplicates(inplace=True)
    bed.reset_index(inplace=True,drop=True)
    return bed

def writeBED(inBED, file_path):
    """
    Writes a bed dataframe into a bed file.
    Bed format: 'chrom','chromStart','chromEnd','name','score','strand'
 
    :param inBED: bed dataframe to be written.
    :param file_path: /path/to/file.bed
    
    :returns: nothing    

    """
    inBED.to_csv(file_path,index=None,sep="\t",header=None)

biomart_host="http://www.ensembl.org/biomart"
rbiomart_host="www.ensembl.org"

def databasesBM(host=biomart_host):
    """
    Lists BioMart databases.
    
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host) 
    server.show_databases()

def RdatabasesBM(host=rbiomart_host):
    """
    Lists BioMart databases through a RPY2 connection.
    
    :param host: address of the host server, default='www.ensembl.org'

    :returns: nothing

    """
    biomaRt = importr("biomaRt")
    print(biomaRt.listMarts(host=host))


def datasetsBM(host=biomart_host):
    """
    Lists BioMart datasets.
    
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host)
    server.show_datasets()

def RdatasetsBM(database,host=rbiomart_host):
    """
    Lists BioMart datasets through a RPY2 connection.
    
    :param database: a database listed in RdatabasesBM()
    :param host: address of the host server, default='www.ensembl.org'

    :returns: nothing

    """
    biomaRt = importr("biomaRt")
    ensemblMart=biomaRt.useMart(database, host=host)
    print(biomaRt.listDatasets(ensemblMart))

def filtersBM(dataset,host=biomart_host):
    """
    Lists BioMart filters for a specific dataset.
    
    :param dataset: dataset to list filters of.
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host)
    d=server.datasets[dataset]
    d.show_filters()

def RfiltersBM(dataset,database,host=rbiomart_host):
    """
    Lists BioMart filters through a RPY2 connection.
    
    :param dataset: a dataset listed in RdatasetsBM()
    :param database: a database listed in RdatabasesBM()
    :param host: address of the host server, default='www.ensembl.org'

    :returns: nothing

    """
    biomaRt = importr("biomaRt")
    ensemblMart=biomaRt.useMart(database, host=host)
    ensembl=biomaRt.useDataset(dataset, mart=ensemblMart)
    print(biomaRt.listFilters(ensembl))

def attributesBM(dataset,host=biomart_host):
    """
    Lists BioMart attributes for a specific dataset.
    
    :param dataset: dataset to list attributes of.
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host)
    d=server.datasets[dataset]
    d.show_attributes()

def RattributesBM(dataset,database,host=rbiomart_host):
    """
    Lists BioMart attributes through a RPY2 connection.
    
    :param dataset: a dataset listed in RdatasetsBM()
    :param database: a database listed in RdatabasesBM()
    :param host: address of the host server, default='www.ensembl.org'
    
    :returns: nothing

    """
    biomaRt = importr("biomaRt")
    ensemblMart=biomaRt.useMart(database, host=rbiomart_host)
    ensembl=biomaRt.useDataset(dataset, mart=ensemblMart)
    print(biomaRt.listAttributes(ensembl))

def queryBM(query_filter,query_items,query_attributes,query_dataset,query_dic=None,host=biomart_host):
    """
    Queries BioMart.

    :param query_filtery: one BioMart filter associated with the items being queried
    :param query_items: list of items to be queried (must assoiate with given filter)
    :param query_querydic: for complex queries this option should be used instead of 'filters' and 'items' and a dictionary of filters provided here eg. querydic={"filter1":["item1","item2"],"filter2":["item3","item4"]}. If using querydic, don't query more than 350 items at once. 
    :param query_attributes: list of attributes to recover from BioMart  
    :param query_dataset: dataset to query
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: a Pandas dataframe of the queried attributes

    """
    server = BiomartServer(host)
    d=server.datasets[query_dataset]
    res=[]
    
    if query_dic is None:
        chunks=[query_items[x:x+350] for x in xrange(0, len(query_items), 350)]
        for c in chunks:
            response=d.search({'filters':{query_filter:c},'attributes':query_attributes})
            for line in response.iter_lines():
                line = line.decode('utf-8')
                res.append(line.split("\t")) 
    else:
        response=d.search({'filters':query_dic,'attributes':query_attributes})
        for line in response.iter_lines():
            line = line.decode('utf-8')
            res.append(line.split("\t"))
    res=pd.DataFrame(res)
    res.columns=query_attributes
    return res
 
def RqueryBM(query_filter,query_items,query_attributes,dataset,database,host=rbiomart_host):
    """
    Queries BioMart.

    :param query_filtery: one BioMart filter associated with the items being queried
    :param query_items: list of items to be queried (must assoiate with given filter)
    :param query_attributes: list of attributes to recover from BioMart  
    :param dataset: dataset to query
    :param database: database to query
    :param host: address of the host server, default='www.ensembl.org'

    return: a Pandas dataframe of the queried attributes

    """
 
    biomaRt = importr("biomaRt")
    ensemblMart=biomaRt.useMart(database, host=rbiomart_host)
    ensembl=biomaRt.useDataset(dataset, mart=ensemblMart)
    df=biomaRt.getBM(attributes=query_attributes, filters=query_filter, values=query_items, mart=ensembl)
    output = [tuple([df[j][i] for j in range(df.ncol)]) for i in range(df.nrow)]
    output = pd.DataFrame(output)
    output.columns = query_attributes
    return output


david_categories = [
  'GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', 'KEGG_PATHWAY',
  'BIOCARTA', 'PFAM', 'PROSITE' ]

david_fields = [
  'categoryName', 'termName', 'listHits', 'percent',
  'ease', 'geneIds', 'listTotals', 'popHits', 'popTotals',
  'foldEnrichment', 'bonferroni', 'benjamini', 'afdr']
# include:
# 'fisher'
# 'termName' to 'term' and 'term_name'

def DAVIDenrich(database, categories, user, ids, ids_bg = None, name = '', name_bg = '', verbose = False, p = 0.1, n = 2):

    """
    Queries the DAVID database for an enrichment analysis
    Check https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html for database == "type" tag and categories ==  "annot" tag.

    :param database: A string for the database to query, e.g. 'WORMBASE_GENE_ID'
    :param categories: A comma separated string with databases
    :param user: A user ID registered at DAVID for querying
    :param ids: A list with identifiers
    :param name: A string with the name for the query set
    :param ids_bg: A list with the background identifiers to enrich against,
      'None' for whole set
    :param name_bg: A string with the name for the background set
    :param p: Maximum p value for enrichment of a term
    :param n: Minimum number of genes within a term
    :param ct: Maybe another threshold

    :returns: None if no ids match the queried database, or a pandas data frame with results
    """

    ids = ','.join(ids)
    use_bg = 0
    if ids_bg is not None:
      ids_bg = ','.join(ids_bg)
    ssl._create_default_https_context = ssl._create_unverified_context
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    client = sudsclient(url)
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
    client_auth = client.service.authenticate(user)
    if verbose:
      print 'User Authentication:', client_auth
      sys.stdout.flush()
    size = client.service.addList(ids, database, name, 0) #| inputListIds,idType,listName,listType)
    if verbose:
      print 'Mapping rate of ids: ', str(size)
      sys.stdout.flush()
    if not float(size) > float(0):
      return None
    if ids_bg is not None:
      size_bg = client.service.addList(ids_bg, database, name_bg, 1)
      if verbose:
        print 'Mapping rate of background ids: ', str(size_bg)
        sys.stdout.flush()
    client_categories = client.service.setCategories(categories)
    if verbose:
      print 'Categories used: ', client_categories
      sys.stdout.flush()
    client_report = client.service.getChartReport(p, n)
    size_report = len(client_report)
    if verbose:
      print 'Records reported: ', str(size_report)
      sys.stdout.flush()

    if size_report > 0:
        df = []
        for r in client_report:
            d = dict(r)
            line = []
            for f in david_fields:
                line.append(str(d[f]))
            df.append(line)
        df = pd.DataFrame(df)
        df.columns=david_fields
    else:
        df=None
    
    return df


def id_nameDAVID(df,GTF=None,name_id=None):
    """
    Given a DAVIDenrich output it converts ensembl gene ids to genes names and adds this column to the output

    :param df: a dataframe output from DAVIDenrich
    :param GTF: a GTF dataframe from readGTF()
    :param name_id: instead of a gtf dataframe a dataframe with the columns 'gene_name' and 'gene_id' can be given as input
    
    :returns: a pandas dataframe with a gene name column added to it.
    """
    if name_id is None:
        gene_name=retrieve_GTF_field('gene_name',GTF)
        gene_id=retrieve_GTF_field('gene_id', GTF)
        GTF=pd.concat([gene_name,gene_id],axis=1)
    else:
        GTF=name_id.copy()
    df['Gene_names']="genes"
    terms=df['termName'].tolist()
    enrichN=pd.DataFrame()
    for term in terms:
        tmp=df[df['termName']==term]
        tmp=tmp.reset_index(drop=True)
        ids=tmp.xs(0)['geneIds']
        ids=pd.DataFrame(data=ids.split(", "))
        ids.columns=['geneIds']
        ids['geneIds']=ids['geneIds'].map(str.lower)
        GTF['gene_id']=GTF['gene_id'].map(str.lower)
        ids=pd.merge(ids, GTF, how='left', left_on='geneIds', right_on='gene_id')
        names=ids['gene_name'].tolist()
        names= ', '.join(names)
        tmp["Gene_names"]=names
        #tmp=tmp.replace(to_replace=tmp.xs(0)['Gene_names'], value=names)
        enrichN=pd.concat([enrichN, tmp])
    enrichN=enrichN.reset_index(drop=True)
    
    gene_names=enrichN[['Gene_names']]
    gpos=enrichN.columns.get_loc("geneIds")
    enrichN=enrichN.drop(['Gene_names'],axis=1)
    cols=enrichN.columns.tolist()
    enrichN=pd.concat([enrichN[cols[:gpos+1]],gene_names,enrichN[cols[gpos+1:]]],axis=1)   

    return enrichN


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
    #print "KEGG API: http://rest.genome.jp/link/"+ens_db+"/"+organism
    #sys.stdout.flush()
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


def getFileFormat (path):
  
  """
  Return the file format

  :param path: The path to the file
  
  :returns: None, if file is missing, else one of the strings 'xlsx', 'xls', 'txt'
  """

  xlsx_sign = b'\x50\x4B\x05\06'
  xls_sign = b'\x09\x08\x10\x00\x00\x06\x05\x00'
  ret = None
  signs = [(0, 512, 8), (2, -22, 4)] # whence, offset, size
  for w, o, s in signs:
    with open(path, 'rb') as f:
      f.seek(o, w)
      b = f.read(s)
      if b == xls_sign:
        ret = 'xls'
      elif b == xlsx_sign:
        ret = 'xlsx'
      else:
        ret = 'txt'
  return ret


def readDataFrame (path, sheet = None, sep = '\t'):

  """
  Returns a pandas data frame

  :param path: The path to the file
  :param sheet: Sheet name or integer 0-based index for xls[x] files, None for all
  :param sep: A separator for text format

  :returns: A pandas data frame
  """

  ff = getFileFormat(path)
  if ff is None:
    print 'error: file not matching format xls[x]/txt'
    sys.exit()
  if ff in ['xls', 'xlsx']:
    df = pd.read_excel(path, sheetname = sheet)
  else:
    df = pd.read_table(path, sep = sep)
  return df

def getFasta(opened_file, sequence_name):
    """
    Retrieves a sequence from an opened multifasta file

    :param opened_file: an opened multifasta file eg. opened_file=open("/path/to/file.fa",'r+')
    :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)

    returns: a string with the sequence of interest
    """

    lines = opened_file.readlines()
    seq=str("")
    for i in range(0, len(lines)):
        line = lines[i]
        if line[0] == ">":
            fChr=line.split(" ")[0].split("\n")[0]
            fChr=fChr[1:]
            if fChr == sequence_name:
                s=i
                code=['N','A','C','T','G']
                firstbase=lines[s+1][0]
                while firstbase in code:
                    s=s + 1
                    seq=seq+lines[s]
                    firstbase=lines[s+1][0]
    
    if len(seq)==0:
        seq=None
    else:        
        seq=seq.split("\n")
        seq="".join(seq)
    
    return seq  


def writeFasta(sequence, sequence_name, output_file):
    """
    Writes a fasta sequence into a file.

    :param sequence: a string with the sequence to be written
    :param sequence_name: name of the the fasta sequence
    :param output_file: /path/to/file.fa to be written

    :returns: nothing
    """
    i=0
    f=open(output_file,'w')
    f.write(">"+str(sequence_name)+"\n")
    while i <= len(sequence):
        f.write(sequence[i:i+60]+"\n")
        i=i+60
    f.close()  

def rewriteFasta(sequence, sequence_name, fasta_in, fasta_out):
    """
    Rewrites a specific sequence in a multifasta file while keeping the sequence header.

    :param sequence: a string with the sequence to be written  
    :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
    :param fasta_in: /path/to/original.fa
    :param fasta_out: /path/to/destination.fa
    
    :returns: nothing
    """
    f=open(fasta_in, 'r+')
    f2=open(fasta_out,'w')
    lines = f.readlines()
    i=0
    while i < len(lines):
        line = lines[i]
        if line[0] == ">":
            f2.write(line)
            fChr=line.split(" ")[0]
            fChr=fChr[1:]
            if fChr == sequence_name:
                code=['N','A','C','T','G']
                firstbase=lines[i+1][0]
                while firstbase in code:
                    i=i+1
                    firstbase=lines[i][0]
                s=0
                while s <= len(sequence):
                    f2.write(sequence[s:s+60]+"\n")
                    s=s+60                
            else:
                i=i+1 
        else:
            f2.write(line)
            i=i+1

    f2.close
    f.close

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

def CellPlot(df, output_file=None, gene_expression="log2FC", figure_title="CellPlot", pvalCol="elimFisher", lowerLimit=None, upperLimit=None, colorBarType='Spectral'):
    """
    Python implementation of the CellPlot from the CellPlot package for R.
    -inf or inf enrichments will come out as min found float or max found float, respectively.

    :param df: pandas dataframe with the following columns - 'Enrichment', 'Term', and 'log2fc'.
               For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term. 
               eg. '-inf,-1,2,3.4,3.66,inf'
    :param output_file: prefix for an output file. If given it will create output_file.CellPlot.svg and output_file.CellPlot.png 
    :param gene_expression: label for the color gradiant bar.
    :param figure_title: Figure title.
    :param pvalCol: name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking      
    :param lowerLimit: lower limit for the heatmap bar (default is the 0.1 percentile)
    :param upperLimit: upper limit for the heatmap bar (default is the 0.9 percentile)      
    :param colorBarType: type of heatmap, 'Spectral' is dafault, alternative eg. 'seismic'
    :returns: a matplotlib figure
    """    
    limits=pd.DataFrame(df['log2fc'].str.split(",").tolist())
    limits=limits.as_matrix().flatten()
    limits=pd.DataFrame(limits)
    limits[0]=limits[0].astype(str)
    limits=limits[limits[0]!="None"][limits[0]!=""]
    
    try:
        limits=[float(x) for x in limits[0].tolist()]
        
    except ValueError,e:
        print "error",e,"on line"
 
    if upperLimit:
        maxFC=upperLimit
    else:
        maxFC=np.percentile(limits,90)
    
    if lowerLimit:
        minFC=lowerLimit
    else:
        minFC=np.percentile(limits,10)

    #maxFC=np.percentile(limits,90)
    #minFC=np.percentile(limits,10)

    cmap = matplotlib.cm.get_cmap(colorBarType)
    norm = matplotlib.colors.Normalize(vmin=minFC, vmax=maxFC)

    if len(df)>10:
        siz=len(df)*3/10
    elif len(df)==1:
        siz=1
    elif len(df)==2:
        siz=2
    else:
        siz=3

    
    fig = plt.figure(figsize=(8, siz))
    #fig.suptitle(figure_title, fontsize=24, fontweight='bold')

    ax1 = fig.add_axes([0.05, 3.5/( float(siz)*float(10)/float(3) ), 0.9, 2])
    ax2 = fig.add_axes([0.05, 1.5/( float(siz)*float(10)/float(3) ), 0.9, 1.5/( float(siz)*float(10)/float(3) )])
    arrangment=np.arange(len(df))+.5
    enr=df['Enrichment'].tolist()
    enr=[x for x in enr if str(x) != str(float("inf"))]
    enr=[x for x in enr if str(x) != str(float("-inf"))]

    m=max(enr)
    
    maxE=max(enr)
    minE=min(enr)

    def getINFs(x,maxFC=maxFC,minFC=minFC):
        if x == str(float("inf")):
            return maxFC
        elif x == str(float("-inf")):
            return minFC
        else:
            return x

    def fix_enrichment(x,minE=minE,maxE=maxE):
        if str(x) == str(float("inf")):
            return maxE
        elif str(x) == str(float("-inf")):
            return minE
        else:
            return x

    df['Enrichment']=df['Enrichment'].apply(lambda x: fix_enrichment(x))


    ax1.barh(arrangment, df['Enrichment'].tolist(), color='white', edgecolor='black')#range(0,len(test))
    for i,pos in zip(df.index.tolist(),arrangment):
        fcs=df.ix[i,'log2fc'].split(",")
        fcs=pd.DataFrame(fcs)
        fcs[0]=fcs[0].astype(str)
        fcs[0]=fcs[0].apply(lambda x: getINFs(x))
        
        #fcs=fcs[fcs[0]!=""].astype(float)[0].tolist()
        fcs=fcs.astype(float)[0].tolist()        

        try:
            w=float(df.ix[i,'Enrichment'])/float(len(fcs))
        except:
            print df.ix[i,]
        p=0
        fcs.sort(key=float)
        for f in fcs:
            #if float(f) > maxFC:
            #    f=maxFC
            #if float(f) < minFC:
            #    f=minFC
            ax1.barh(pos, w, left=p, color=cmap(norm(float(f))), edgecolor='black')
            p=p+w
        if pvalCol:
            if df.ix[i,pvalCol] < 0.05:
                barAn=len(fcs)
            else:
                barAn=str(len(fcs))+" (NS)"
        else:
            barAn=len(fcs)
        ax1.text(df.ix[i,'Enrichment']+m*.02, pos+0.25, barAn, ha='left', va='bottom')

    ax1.set_yticks(arrangment+0.4)
    ax1.set_yticklabels(df['Term'].tolist())

    ax1.tick_params(
        axis='y',          
        which='both',      
        left='off',      
        right='off',         
        labelleft='on')

    ax1.tick_params(
        axis='x',          
        which='both',     
        bottom='off',      
        top='on',         
        labelbottom='off',
        labeltop='on') 
    
    ax1.set_ylim(ymax = max(arrangment) + 1.5 ) #1.5
    ax1.set_xlabel("GO Term Enrichment")
    ax1.xaxis.set_label_position('top') 

    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

    cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm, orientation='horizontal')
    if not lowerLimit:
        if not upperLimit:
            cb1.set_label(gene_expression+'\n(0.1-0.9 percentiles)\n\n\n'+figure_title)
        else:
            cb1.set_label(gene_expression+'\n\n\n'+figure_title)
    else:
        cb1.set_label(gene_expression+'\n\n\n'+figure_title)

    #plt.subplots_adjust(top=8.5)

    #ax1.annotate('a fractional title', xy=(.025, .975), xycoords='figure fraction', horizontalalignment='center', verticalalignment='top', fontsize=20)
    #plt.subplots_adjust(top=500)
    #ax1.text(2, 2, 'right bottom',
    #    horizontalalignment='center',
    #    verticalalignment='center',
    #    transform=ax1.transAxes)

    if output_file:
        plt.savefig(output_file+".CellPlot.png",dpi=300,bbox_inches='tight', pad_inches=0.1,format='png')
        plt.savefig(output_file+".CellPlot.svg",dpi=300,bbox_inches='tight', pad_inches=0.1,format='svg')

    return fig

def SymPlot(df,output_file=None,figure_title="SymPlot",pvalCol="elimFisher"):
    """
    Python implementation of the SymPlot from the CellPlot package for R.
    -inf or inf enrichments will come out as min found float or max found float, respectively.    

    :param df: pandas dataframe with the following columns - 'Enrichment', 'Significant', 'Annotated', 'Term', and 'log2fc'.
               For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term. 
               eg. '-inf,-1,2,3.4,3.66,inf'
    :param output_file: prefix for an output file. If given it witll create output_file.SymPlot.svg and output_file.SymPlot.png 
    :param figure_title: Figure title.
    :param pvalCol: name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking 
    :returns: a matplotlib figure
    """
    maxAn=df['Annotated'].max()

    arrangment=np.arange(len(df))+.5

    def getINFs(x):
        if x == str(float("inf")):
            return 1
        elif x == str(float("-inf")):
            return -1
        else:
            return x

    enr=df['Enrichment'].tolist()
    enr=[x for x in enr if str(x) != str(float("inf"))]
    enr=[x for x in enr if str(x) != str(float("-inf"))]

    maxE=max(enr)
    minE=min(enr)

    def fix_enrichment(x,minE=minE,maxE=maxE):
        if str(x) == str(float("inf")):
            return maxE
        elif str(x) == str(float("-inf")):
            return minE
        else:
            return x

    df['Enrichment']=df['Enrichment'].apply(lambda x: fix_enrichment(x))

    limits=df['Enrichment'].tolist()
    maxFC=np.percentile(limits,90)
    minFC=np.percentile(limits,10)

    cmap = matplotlib.cm.get_cmap('Spectral')
    norm = matplotlib.colors.Normalize(vmin=minFC, vmax=maxFC)

    if len(df) >= 5:
        siz=len(df)*4/10
    else:
        size=5*4/10


    fig = plt.figure(figsize=(8, siz))
    #fig.suptitle(figure_title, fontsize=24, fontweight='bold')
    gs = gridspec.GridSpec(1, 3, width_ratios=[2,0.75,2])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = fig.add_axes([0.125, 0.11/100*len(df), 0.775, 0.075*10/len(df)])#/100.0075

    alldown=[]
    allup=[]

    for i,pos in zip(df.index.tolist(),arrangment):
        f=df.ix[i,'Enrichment']
        sigN=df.ix[i,'Significant']
        ann=float(df.ix[i,'Annotated'])

        if ann!=maxAn:
            p=float(maxAn-ann)/2
        else:
            p=0
        ax2.barh(pos, ann, left=p, color=cmap(norm(float(f))),edgecolor='black')#

        fcs=df.ix[i,'log2fc'].split(",")
        fcs=pd.DataFrame(fcs)
        fcs[0]=fcs[0].astype(str)
        fcs[0]=fcs[0].apply(lambda x: getINFs(x))
        #fcs=fcs[fcs[0]!=""].astype(float)     
        fcs=fcs.astype(float)  
        down=len(fcs[fcs[0]<0])/ann*100
        up=len(fcs[fcs[0]>0])/ann*100
        alldown.append(down)
        allup.append(up)

        ax1.barh(pos, down, color="blue",edgecolor='blue')
        ax3.barh(pos, up, color="red",edgecolor='red')

    ax1.spines['top'].set_visible(True)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax1.tick_params(axis='x',which='both',bottom='off', top='on',labelbottom='off',labeltop='on')
    ax1.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')

    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off',labeltop='off')
    ax2.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

    ax3.spines['top'].set_visible(True)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    ax3.tick_params(axis='x',which='both',bottom='off',top='on',labelbottom='off',labeltop='on')
    ax3.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')

    fa=10*0.1/len(df)+1
    fb=10*0.08/len(df)+1
    
    ax1.set_title('Downregulated (%)',y=fa)#
    ax2.set_title('Annotated\n(max=%s)' %str(maxAn),y=fb)#
    ax3.set_title('Upregulated (%)',y=fa)

    ax1.set_xlim(max(max(alldown),max(allup)), 0)
    ax2.set_xlim(0, maxAn)
    ax3.set_xlim(0, max(max(alldown),max(allup)))

    ax1.set_ylim(ymax = max(arrangment)+1.5)
    ax2.set_ylim(ymax = max(arrangment)+1.5)
    ax3.set_ylim(ymax = max(arrangment)+1.5)
    
    
    ax1.set_yticks(arrangment+0.4)
    def get_label_with_sig (df):
        termLabel=df['Term']
        if pvalCol:
            pvalue=df[pvalCol]
            if pvalue > 0.05:
                return "(NS) "+termLabel
            else:
                return termLabel
        else:
            return termLabel

    df['newLabels']=df.apply(get_label_with_sig, axis=1)  

    ax1.set_yticklabels(df['newLabels'].tolist())

    cb1 = matplotlib.colorbar.ColorbarBase(ax4, cmap=cmap,norm=norm, orientation='horizontal')
    cb1.set_label('GO Term Enrichment (0.1-0.9 percentiles)\n\n\n'+figure_title)
    
    fig.subplots_adjust(wspace=0)

    if output_file:
        plt.savefig(output_file+".SymPlot.png",dpi=300,bbox_inches='tight', pad_inches=0.1,format='png')
        plt.savefig(output_file+".SymPlot.svg",dpi=300,bbox_inches='tight', pad_inches=0.1,format='svg')
    
    return fig

DNAcode={'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TAG': 'STOP', 'GGA': 'G', 'TAA': 'STOP', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'TGA': 'STOP', 'GAC': 'D', 'CGT': 'R', 'TGG': 'W', 'GAA': 'E', 'CGC': 'R'}




if __name__ == '__main__':
    print "AGEpy"
    
