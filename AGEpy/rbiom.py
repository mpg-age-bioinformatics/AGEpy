import pandas as pd
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

rbiomart_host="www.ensembl.org"

def RdatabasesBM(host=rbiomart_host):
    """
    Lists BioMart databases through a RPY2 connection.

    :param host: address of the host server, default='www.ensembl.org'

    :returns: nothing

    """
    biomaRt = importr("biomaRt")
    print(biomaRt.listMarts(host=host))

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
