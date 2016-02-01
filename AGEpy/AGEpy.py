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

def databasesBM(host=biomart_host):
    """
    Lists BioMart datasets.
    
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host) 
    server.show_databases()


def datasetsBM(host=biomart_host):
    """
    Lists BioMart datasets.
    
    :param host: address of the host server, default='http://www.ensembl.org/biomart'

    :returns: nothing

    """
    server = BiomartServer(host)
    server.show_datasets()


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
    client = sudsclient('https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl')
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
        df.append(david_fields)
        for r in client_report:
            d = dict(r)
            line = []
            for f in david_fields:
                line.append(str(d[f]))
                df.append(line)
        df = pd.DataFrame(df)
    else:
        df=None
    
    return df


def organismsKEGG():
    """
    Lists all organisms present in the KEGG database.

    :returns: nothing

    """
    organisms=urlopen("http://rest.kegg.jp/list/organism").read()
    organisms=organisms.split("\n")
    for o in organisms:
        o=o.split("\t")
        print str(o[1]), str(o[2])
        sys.stdout.flush()


def databasesKEGG(organism,ens_ids):
    """
    Finds KEGG database identifiers for a respective organism given example ensembl ids.
    

    :param organism: an organism as listed in organismsKEGG()
    :param ens_ids: a list of ensenbl ids of the respective organism

    :returns: nothing

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

    :returns: a Pandas dataframe with the columns 'KEGGid','pathIDs', and 'pathName'.
    
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
    return df_final




def expKEGG(organism,names_KEGGids):
    """
    Gets all KEGG pathways for an organism

    :param organism: an organism as listed in organismsKEGG()
    :param names_KEGGids: a Pandas dataframe with the columns 'gene_name' and  'KEGGid' as reported from idsKEGG(organism) (or a subset of it).

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


def getFileFormat (path):
  
  """
  Return the file format

  :params path: The path to the file
  
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

  :params path: The path to the file
  :params sheet: Sheet name or integer 0-based index for xls[x] files, None for all
  :params sep: A separator for text format

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



if __name__ == '__main__':
    print "AGEpy"
    
