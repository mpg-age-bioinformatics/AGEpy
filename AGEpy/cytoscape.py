import pandas as pd
import numpy as np
import os
import sys
import requests
import os
import sys
from time import sleep
import json
import urllib2
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from wand.image import Image as WImage
import paramiko
import tempfile
from shutil import copyfile

cytoscape_host="192.168.50.110"
cytoscape_port="1234"

def CheckResponse(r):
    status=str(r.status_code)
    if status == "200":
        #print "Finished"
        sys.stdout.flush()
    elif status == "201":
        sys.stdout.flush()
    else:
        print r, r.status_code
        sys.stdout.flush()
        
def checkCytoscapeVersion(host=cytoscape_host,port=cytoscape_port):
    """
    Checks cytoscape version
    
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    
    :returns: cytoscape and api version
    """
    
    URL="http://"+str(host)+":"+str(port)+"/v1/version/"
    r = requests.get(url = URL)
    r=json.loads(r.content)
    for k in r.keys():
        print k, r[k]
        
def cytoscape(namespace,command="",PARAMS={},host=cytoscape_host,port=cytoscape_port,method="POST",verbose=False):
    """
    General function for interacting with Cytoscape API.
    
    :param namespace: namespace where the request should be executed. eg. "string"
    :param commnand: command to execute. eg. "protein query"
    :param PARAMs: a dictionary with the parameters. Check your swagger normaly running on
    http://localhost:1234/v1/swaggerUI/swagger-ui/index.html?url=http://localhost:1234/v1/commands/swagger.json
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    :param method: type of http call, ie. "POST" or "GET" or "HELP".
    :param verbose: print more information
    
    :returns: For "POST" the data in the content's response. For "GET" None. 
    
    eg.

    cytoscape("G","string","pubmed query",{"pubmed":"p53 p21","limit":"50"})


    """     

    if (method == "GET") or (method == "G"):
        P=[]
        for p in PARAMS.keys():
            v=str(PARAMS[p])
            v=v.replace(" ","%20")  
            P.append(str(p)+"="+str(PARAMS[p]))
        P="&".join(P)
        URL="http://"+str(host)+":"+str(port)+"/v1/commands/"+str(namespace)+"/"+str(command)+"?"+P
        if verbose:
            print "'"+URL+"'"
            sys.stdout.flush()
        r = requests.get(url = URL)
        CheckResponse(r)
        res=None
    
    elif (method == "POST") or (method == "P"):
        URL="http://"+str(host)+":"+str(port)+"/v1/commands/"+str(namespace)+"/"+str(command)
        r = requests.post(url = URL, json = PARAMS)
        CheckResponse(r)
        res=r.content
        if verbose:
            print res
        res=json.loads(res)
        res=res["data"]
        
    elif (method=="HTML") or (method == "H") or (method=="HELP"):
        P=[]
        for p in PARAMS.keys():
            v=str(PARAMS[p])
            v=v.replace(" ","%20")  
            P.append(str(p)+"="+str(PARAMS[p]))
        P="&".join(P)
        URL="http://"+str(host)+":"+str(port)+"/v1/commands/"+str(namespace)+"/"+str(command)+"?"+P
        if verbose:
            print "'"+URL+"'"
            sys.stdout.flush()
        response = urllib2.urlopen(URL)
        #print response
        res = response.read()
        print res
        sys.stdout.flush()
    
    return res

def getTableColumns(table, columns, namespace = "default", network = "current", host=cytoscape_host,port=cytoscape_port,verbose=False):
    """
    Gets tables from cytoscape
    
    :param table: table to retrieve eg. node
    :param columns: columns to retrieve
    :param namespace: namepsace, default="default"
    :param network: a network name or id, default="current"
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    :param verbose: print more information

    :returns: a pandas dataframe
    """
    
    
    if type(network) != int:
        network=cytoscape("network", "get attribute",\
                      {"network":network,\
                       "namespace":namespace,\
                       "columnList":"SUID"},host=host,port=port)
        network=network[0]["SUID"]

    
    df=pd.DataFrame()
    def target(column):
        URL="http://"+str(host)+":"+str(port)+"/v1/networks/"+str(network)+"/tables/"+namespace+table+"/columns/"+column  
        if verbose:
            print "'"+URL+"'"
            sys.stdout.flush()
        response = urllib2.urlopen(URL)
        response = response.read()
        colA=json.loads(response)

        col=pd.DataFrame()    
        colHeader=colA["name"]
        colValues=colA["values"]
        col[colHeader]=colValues
        return col
    
    ncols=["name"]
    for c in columns:
        ncols.append(c.replace(" ","%20") )
    for c in ncols:
        col=target(c)
        df=pd.concat([df,col],axis=1)
        
    df.index=df["name"].tolist()
    df=df.drop(["name"],axis=1)
    
    return df

def loadTableData(df, df_key='index',table="node", \
                   table_key_column = "name",
                   network="current",\
                   namespace="default",\
                   host=cytoscape_host,port=cytoscape_port,verbose=False):

    """
    Loads tables into cytoscape
    
    :param df: a pandas dataframe to load
    :param df_key: key column in df
    :param table: target table, default="node"
    :param table_key_column: table key column, default="name"
    :param network: a network name or id, default="current"
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    :param verbose: print more information

    :returns: output of put request
    """
    
    
    
    if type(network) != int:
        networkID=cytoscape("network", "get attribute",\
                      {"network":network,\
                       "namespace":namespace,\
                       "columnList":"SUID"},host=host,port=port)
        networkname=cytoscape("network", "get attribute",\
                      {"network":network,\
                       "namespace":namespace,\
                       "columnList":"name"},host=host,port=port)
        
        network=networkID[0]["SUID"]
        networkname=networkname[0]["name"]
        
    tmp=df.copy()
    if df_key!="index":
        tmp.index=tmp[df_key].tolist()
        tmp=tmp.drop([df_key],axis=1)
            
    tablen=networkname+" default node"
    
    data=[]
    
    for c in tmp.columns.tolist():
        tmpcol=tmp[[c]].dropna()
        for r in tmpcol.index.tolist():
            cell={}
            cell[str(table_key_column)]=str(r) # {"name":"p53"}
            val=tmpcol.loc[r,c]
            cell[str(c)]=val
            data.append(cell)
    
    
    upload={"key":table_key_column,"dataKey":table_key_column,\
            "data":data}
    
    
    URL="http://"+str(host)+":"+str(port)+"/v1/networks/"+str(network)+"/tables/"+namespace+table  
    if verbose:
        print "'"+URL+"'"
        sys.stdout.flush()
    r = requests.put(url = URL, json = upload)
    if verbose:
        print r
    CheckResponse(r)
    res=r.content
    return res

def simple_defaults(defaults_dic):
    """
    Simplifies defaults.
    
    :param defaults_dic: a dictionary of the form { visualProperty_A:value_A, visualProperty_B:value_B, ..}
    
    :returns: a list of dictionaries with each item corresponing to a given key in defaults_dic
    """
     
    defaults=[]
    for d in defaults_dic.keys():
        dic={}
        dic["visualProperty"]=d
        dic["value"]=defaults_dic[d]
        defaults.append(dic)
    return defaults
    
def update_style(title,defaults=None,mappings=None,host=cytoscape_host,port=cytoscape_port,verbose=False):
    """
    Updates a visual style
    
    :param title: title of the visual style
    :param defaults: a list of dictionaries for each visualProperty
    :param mappings: a list of dictionaries for each visualProperty
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    
    :retunrs: nothing
    """
    
    URL="http://"+str(host)+":"+str(port)+"/v1/styles/"+str(title)
    if verbose:
        print URL
        sys.stdout.flush()

    response = urllib2.urlopen(URL)    
    response = response.read()
    response = json.loads(response)
    
    olddefaults=response["defaults"]
    oldmappings=response["mappings"]

    if mappings:
        mappings_visual_properties=[ m["visualProperty"] for m in mappings ]
        newmappings=[ m for m in oldmappings if m["visualProperty"] not in mappings_visual_properties ]
        for m in mappings:
            newmappings.append(m)
    else:
        newmappings=oldmappings

    if defaults:
        defaults_visual_properties=[ m["visualProperty"] for m in defaults ]
        newdefaults=[ m for m in olddefaults if m["visualProperty"] not in defaults_visual_properties ]
        for m in defaults:
            newdefaults.append(m)
    else:
        newdefaults=olddefaults
        
    r=requests.delete(URL)
    CheckResponse(r)    

    URL="http://"+str(host)+":"+str(port)+"/v1/styles"
    PARAMS={"title":title,\
           "defaults":newdefaults,\
           "mappings":newmappings}
    r = requests.post(url = URL, json = PARAMS)
    CheckResponse(r)
    
def create_styles(title,defaults=None,mappings=None,host=cytoscape_host,port=cytoscape_port):
    """
    Creates a new visual style
    
    :param title: title of the visual style
    :param defaults: a list of dictionaries for each visualProperty
    :param mappings: a list of dictionaries for each visualProperty
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    
    :retunrs: nothing
    """   
    
    try:
        update_style(title,defaults=defaults,mappings=mappings,host=host,port=port)
        print "Existing style was updated."
        sys.stdout.flush()
    except:
        print "Creating new style."
        sys.stdout.flush()
        URL="http://"+str(host)+":"+str(port)+"/v1/styles"
        PARAMS={"title":title,\
               "defaults":defaults,\
               "mappings":mappings}
        r = requests.post(url = URL, json = PARAMS)
        CheckResponse(r)

    
def mapVisualProperty(visualProperty,mappingType,mappingColumn,\
                      lower=None,center=None,upper=None,\
                      discrete=None,\
                      network="current",table="node",\
                      namespace="default",host=cytoscape_host,\
                      port=cytoscape_port, verbose=False):
    """"
    Generates a dictionary for a given visual property
    
    :param visualProperty: visualProperty
    :param mappingType: mappingType
    :param mappingColumn: mappingColumn
    :param lower: for "continuous" mappings a list of the form [value,rgb_string]
    :param center: for "continuous" mappings a list of the form [value,rgb_string]
    :param upper: for "continuous" mappings a list of the form [value,rgb_string]
    :param discrete: for discrete mappings, a list of lists of the form [ list_of_keys, list_of_values ]
    :param network: a network name or id, default="current"
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234
    
    :retunrs: a dictionary for the respective visual property

    """
    if type(network) != int:
        networkID=cytoscape("network", "get attribute",\
                      {"network":network,\
                       "namespace":namespace,\
                       "columnList":"SUID"},host=host,port=port)
        networkname=cytoscape("network", "get attribute",\
                      {"network":network,\
                       "namespace":namespace,\
                       "columnList":"name"},host=host,port=port)
        
        network=networkID[0]["SUID"]
        networkname=networkname[0]["name"]

        URL="http://"+str(host)+":"+str(port)+"/v1/networks/"+str(network)+"/tables/"+namespace+table+"/columns/"
        if verbose:
            print URL 
            sys.stdout.flush()
        response = urllib2.urlopen(URL)
        response = response.read()
        response = json.loads(response)

    for r in response:
        if r["name"]==mappingColumn:
            mappingColumnType=r["type"]
            break
            
    PARAMS={"mappingType" : mappingType,\
            "mappingColumn" : mappingColumn,
            "mappingColumnType" : mappingColumnType,
            "visualProperty" : visualProperty}
    if mappingType == "continuous":
        PARAMS["points"]=[{"value" : lower[0],\
                           "lesser" : lower[1],\
                           "equal" : lower[1],\
                           "greater" : lower[1]},\
                          {"value" : center[0],
                           "lesser" : center[1],
                           "equal" : center[1],
                           "greater" : center[1] },\
                          {"value" : upper[0],\
                           "lesser" : upper[1],\
                           "equal" : upper[1],\
                           "greater" : upper[1]}]
        
    if discrete:
        PARAMS["map"]=[]
        for k,v in zip(discrete[0],discrete[1]):
            PARAMS["map"].append({ "key":k,"value":v})        
    return PARAMS

def result(filetype="PNG",saveas=None, host=cytoscape_host,port=cytoscape_port):
    """
    Checks the current network. 
    
    Note: works only on localhost
    
    :param outfile: name of the temporary outfile
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234

    :returns: nothing
    """
    sleep(1)

    def MAKETMP():
        (fd, tmp_file) = tempfile.mkstemp()
        tmp_file="/tmp/"+tmp_file.split("/")[-1]
        return tmp_file
    
    outfile=MAKETMP()
    
    extensions={"PNG":".png","PDF":".pdf","CYS":".cys","CYJS":".cyjs"}
    ext=extensions[filetype]
    
    response=cytoscape("view","fit content",host=host,port=port)
    response=cytoscape("view", "export" , \
                       {"options":filetype,\
                        "OutputFile":outfile},\
                      host=host,port=port)
    if host!='localhost':
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(host)
        ftp_client=ssh.open_sftp()
        ftp_client.get(outfile+ext,outfile+ext)
        ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+outfile+ext )
        
    img = WImage(filename=outfile+ext)
    if saveas:
        copyfile(outfile+ext,saveas)
    os.remove(outfile+ext)
    return img
