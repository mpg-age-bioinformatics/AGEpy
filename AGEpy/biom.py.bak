import pandas as pd
from biomart import BiomartServer

biomart_host="http://www.ensembl.org/biomart"

def databasesBM(host=biomart_host):
    """
    Lists BioMart databases.

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
