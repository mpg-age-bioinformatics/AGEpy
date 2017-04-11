import pandas as pd
import sys
from suds.client import Client as sudsclient
import ssl

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

    ids = ','.join([str(i) for i in ids])
    use_bg = 0
    if ids_bg is not None:
      ids_bg = ','.join([str(i) for i in ids_bg])
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
                line.append(unicode(d[f]).encode('ascii','ignore'))
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

def DAVIDgetGeneAttribute(x,parsedGTF,refCol="ensembl_gene_id",fieldTOretrieve="gene_name"):
    """
    Returns a list of gene names for given gene ids.

    :param x: a string with the list of IDs separated by ', '
    :param parsedGTF: a parsed GTF with at least refcCol and fieldTOretrieve as retrieved from parsedGTF()
    :param refCol: the header of the column containing the identifiers
    :param fieldTOretrieve: the field to retrieve from parsedGTF eg. 'gene_name'

    :returns: list of fieldTOretrieve separeted by ', ' in the same order as the given in x
    """

    l=x.split(", ")
    l=[ s.upper() for s in l ]
    tmpdf=pd.DataFrame({refCol:l},index=range(len(l)))
    df_fix=parsedGTF[[refCol,fieldTOretrieve]].drop_duplicates()
    df_fix[refCol]=df_fix[refCol].apply(lambda x: x.upper())
    ids=pd.merge(tmpdf,df_fix,how="left",on=[refCol])
    ids=ids[fieldTOretrieve].tolist()
    ids=[ str(s) for s in ids ]
    ids=", ".join(ids)
    return ids
