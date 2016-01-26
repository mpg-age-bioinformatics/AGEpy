#!/usr/bin/env python

""" Bioinformatics tools developed at the Max Planck Institute for Biology of Ageing"""

import pandas as pd
import numpy as np
import os
import csv
import sys
from suds.client import Client as sudsclient
import os
import ssl

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
   
    :param inBED: bed dataframe to be written.
    :param file_path: /path/to/file.bed
    
    :returns: nothing    

    """
    inBED.to_csv(file_path,index=None,sep="\t",header=None)


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

    :returns: None if no ids match the queried database, or a list with results
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
      size_bg = cliet.service.addList(ids, database, name_bg, 1)
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
    return client_report


def write_DAVID(report, path, sep = '\t'):

  """
  Write the DAVIDenrich list to a file

  :param report: A report list as returned from DAVIDenrich
  :param path: A file path to write to
  :param sep: A character used as separator

  :returns: None
  """

  n = len(report)
  if n < 1:
    print 'warning: report of length 0'
    #sys.exit()
  with open(path, 'w') as out:
    out.write(sep.join(david_fields) + '\n')
    for r in report:
      d = dict(r)
      line = []
      for f in david_fields:
        line.append(str(d[f]))
      out.write(sep.join(line) + '\n')


def pandas_DAVID(report):

  """
  Return a pandas data frame from a DAVIDenrich

  :param report: A report list as returned from DAVIDenrich

  :returns: A pandas data frame
  """

  df = []
  df.append(david_fields)
  for r in report:
    d = dict(r)
    line = []
    for f in david_fields:
      line.append(str(d[f]))
    df.append(line)
  pdf = pd.DataFrame(df)
  return pdf


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
    
