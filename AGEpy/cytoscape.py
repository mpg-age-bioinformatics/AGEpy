import pandas as pd
import numpy as np
import os
import sys
import requests
import os
import sys
from time import sleep
import json
#import urllib2 # python2
import urllib.request as urllib2
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
        print(r, r.status_code)
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
        print(k, r[k])
        
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

    cytoscape("string","pubmed query",{"pubmed":"p53 p21","limit":"50"})


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
            print("'"+URL+"'")
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
            print(res)
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
            print("'"+URL+"'")
            sys.stdout.flush()
        response = urllib2.urlopen(URL)
        #print response
        res = response.read()
        print(res)
        sys.stdout.flush()
    
    return res

def getTableColumns(table, columns, namespace = "default", network = "current", host=cytoscape_host,port=cytoscape_port,verbose=False):
    """
    Gets tables from cytoscape
    
    :param table: table to retrieve eg. node
    :param columns: columns to retrieve in list format
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
            print("'"+URL+"'")
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
        try:
            col=target(c)
            df=pd.concat([df,col],axis=1)
        except:
            print("Could not find "+c)
            sys.stdout.flush()        

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
    :param df_key: key column in df, default="index"
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
            check=tmpcol[tmpcol.index.isin([r])]
            
            cell={}
            cell[str(table_key_column)]=str(r) # {"name":"p53"}
            if len(check) == 1:
                val=tmpcol.loc[r,c]
            
                if type(val) != str:
                    val=float(val)
            else:
                print(check)
                val=""
            cell[str(c)]=val
            data.append(cell)
    
    
    upload={"key":table_key_column,"dataKey":table_key_column,\
            "data":data}
    
    
    URL="http://"+str(host)+":"+str(port)+"/v1/networks/"+str(network)+"/tables/"+namespace+table  
    if verbose:
        print("'"+URL+"'")
        sys.stdout.flush()
    r = requests.put(url = URL, json = upload)
    if verbose:
        print(r)
    CheckResponse(r)
    res=r.content
    return res

def simple_defaults(defaults_dic):
    """
    Simplifies defaults.
    
    :param defaults_dic: a dictionary of the form { visualProperty_A:value_A, visualProperty_B:value_B, ..}
    
    :returns: a list of dictionaries with each item corresponding to a given key in defaults_dic
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

    if defaults: 
        defaults_=[]    
        for d in defaults:
            if d:
                defaults_.append(d)
        defaults=defaults_

    if mappings:
        mappings_=[]
        for m in mappings:
            if m:
                mappings_.append(m)
        mappings=mappings_


    URL="http://"+str(host)+":"+str(port)+"/v1/styles/"+str(title)
    if verbose:
        print(URL)
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
    
    if defaults:
        defaults_=[]
        for d in defaults:
            if d:
                defaults_.append(d)
        defaults=defaults_

    if mappings:
        mappings_=[]
        for m in mappings:
            if m:
                mappings_.append(m)
        mappings=mappings_

    try:
        update_style(title,defaults=defaults,mappings=mappings,host=host,port=port)
        print("Existing style was updated.")
        sys.stdout.flush()
    except:
        print("Creating new style.")
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
            print(URL)
            sys.stdout.flush()
        response = urllib2.urlopen(URL)
        response = response.read()
        response = json.loads(response)

    mappingColumnType=None
    for r in response:
        if r["name"]==mappingColumn:
            mappingColumnType=r["type"]
            break
    if not mappingColumnType:
        print("For mappingType: "+mappingType+" it was not possible to find a  mappingColumnType.")
        sys.stdout.flush()
    
        
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
    
    if not mappingColumnType:
        res=None
    else:
        res=PARAMS    

    return res

def result(filetype="PNG",saveas=None, host=cytoscape_host,port=cytoscape_port):
    """
    Checks the current network. 
    
    Note: works only on localhost
   
    :param filetype: file type, default="PNG" 
    :param saveas: /path/to/non/tmp/file.prefix
    :param host: cytoscape host address, default=cytoscape_host
    :param port: cytoscape port, default=1234

    :returns: an image
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

def aDiffCytoscape(df,aging_genes,target,species="caenorhabditis elegans",limit=None, cutoff=0.4,\
                  taxon=None,host=cytoscape_host,port=cytoscape_port):
    """
    Plots tables from aDiff/cuffdiff into cytoscape using String protein queries. 
    Uses top changed genes as well as first neighbours and difusion fo generate subnetworks. 
    
    :param df: df as outputed by aDiff for differential gene expression
    :param aging_genes: ENS gene ids to be labeled with a diagonal
    :param species: species for string app query. eg. "caenorhabditis elegans", "drosophila melanogaster", "mus musculus", "homo sapiens"
    :param limit: limit for string app query. Number of extra genes to recover. If None, limit=N(query_genes)*.25
    :param cuttoff: confidence cuttoff for sting app query. Default=0.4
    :param taxon: taxon id for string app query. For the species shown above, taxon id will be automatically identified
    :param cytoscape_host: host address for cytoscape
    :param cytoscape_port: cytoscape port
    :param target: target destination for saving files without prefix. eg. "/beegfs/group_bit/home/JBoucas/test/N2_vs_daf2"

    :returns: nothing
    
    """

    ##### TEMPORARY FIX - STRING APP NOT ACCEPTING QUERIES ABOVE 2000 GENES ####
    df=df.sort_values(by=["q_value"],ascending=True)
    df.reset_index(inplace=True, drop=True)
    tmp=df[:1999]
    df=tmp.copy()
    ##### END OF TEMPORARY FIX #####

    query_genes=df["ensembl_gene_id"].tolist()
    df["NormInt"]=df["value_1"]*df["value_2"]
    df["NormInt"]=df["NormInt"].apply(lambda x: np.log10(np.sqrt(x)) )

    if not limit:
        limit=int(len(query_genes)*.25)

    # Annotate aging evindence
    def CheckEvidence(x,aging_genes=aging_genes):
        if x in aging_genes:
            res="aging_gene"
        else:
            res="no"
        return res
    df["evidence"]=df["ensembl_gene_id"].apply(lambda x:CheckEvidence(x) )

    # fix infinit values
    def FixInfs(x):
        if str(x) in ["-inf","inf"]:
            res=np.nan
        else:
            res=x
            return res
    df["NormInt"]=df["NormInt"].apply( lambda x:  FixInfs(x) )
    df["log2(fold_change)"]=df["log2(fold_change)"].apply( lambda x:  FixInfs(x) )
        
    taxons={"caenorhabditis elegans":"6239","drosophila melanogaster":"7227",\
       "mus musculus":"10090","homo sapiens":"9606"}
    if not taxon: 
        taxon=taxons[species]


    # destroy any existing network still present in cytoscape
    response=cytoscape("network", "list",\
                       host=host, port=port)
    if "networks" in response.keys():
        response=response["networks"]
        #print response
        if len(response) > 0:
            for r in response:
                rr=cytoscape("network", "destroy",{"network":"SUID:"+str(r)},\
                             host=host, port=port)
    
    # String protein query
    query_genes=[ str(s) for s in query_genes ]

    response=cytoscape("string", "protein query",\
                      {"query":",".join(query_genes),\
                       "cutoff":str(cutoff),\
                       "species":species,\
                       "limit":str(limit),\
                       "taxonID":taxon},\
                       host=host, port=port)

    print("giving some time to cytoscape..")
    sys.stdout.flush()
    sleep(10)

    # apply new layout
    response=cytoscape("layout", "force-directed",\
                      {"defaultSpringCoefficient":".000004",\
                       "defaultSpringLength":"5"},\
                       host=host, port=port)
    
    # redefine defaults for node visualization
    response=loadTableData(df[["ensembl_gene_id","log2(fold_change)","NormInt","evidence"]].dropna(),\
                       df_key="ensembl_gene_id",table_key_column="query term",\
                      host=host, port=port)
    defaults_dic={"NODE_SHAPE":"ellipse",\
               "NODE_SIZE":60,\
               "NODE_FILL_COLOR":"#AAAAAA",\
               "EDGE_TRANSPARENCY":120}
    defaults_list=simple_defaults(defaults_dic)

    # apply mappings - blue / white / red - from -4 to +4 log2FC
    NODE_LABEL=mapVisualProperty("NODE_LABEL","passthrough","display name",host=host, port=port)

    create_styles("dataStyle",defaults_list,[NODE_LABEL],host=host, port=port)
    response=cytoscape("vizmap", "apply", {"styles":"dataStyle"},host=host, port=port)
    cmap = matplotlib.cm.get_cmap("bwr")
    norm = matplotlib.colors.Normalize(vmin=-4, vmax=4)
    min_color=matplotlib.colors.rgb2hex(cmap(norm(-4)))
    center_color=matplotlib.colors.rgb2hex(cmap(norm(0)))
    max_color=matplotlib.colors.rgb2hex(cmap(norm(4)))  

    NODE_FILL_COLOR=mapVisualProperty('NODE_FILL_COLOR','continuous','log2(fold_change)',\
                         lower=[-4,min_color],center=[0.0,center_color],upper=[4,max_color],\
                                     host=host, port=port)
    
    # apply diamond shape and increase node size to nodes with aging evidence
    NODE_SHAPE=mapVisualProperty('NODE_SHAPE','discrete','evidence',discrete=[ ["aging_gene","no"], ["DIAMOND", "ellipse"] ],\
                                host=host, port=port)
    NODE_SIZE=mapVisualProperty('NODE_SIZE','discrete','evidence',discrete=[ ["aging_gene","no"], ["100.0","60.0"] ],\
                               host=host, port=port)
    update_style("dataStyle",mappings=[NODE_SIZE,NODE_SHAPE,NODE_FILL_COLOR],\
                host=host, port=port)
    response=cytoscape("vizmap", "apply", {"styles":"dataStyle"},\
                      host=host, port=port)
    
    # apply mappings - reds - to Normalized expression (as in MA plots) to border color and border size
    NormIntDf = getTableColumns('node',['NormInt'],host=host, port=port)
    if 'NormInt' in NormIntDf.columns.tolist():
        min_NormInt = min(NormIntDf.dropna()['NormInt'].tolist())
        max_NormInt = max(NormIntDf.dropna()['NormInt'].tolist())
        cent_NormInt = np.mean([min_NormInt,max_NormInt])

        cmap = matplotlib.cm.get_cmap("Reds")
        norm = matplotlib.colors.Normalize(vmin=min_NormInt, vmax=max_NormInt)
        min_color=matplotlib.colors.rgb2hex(cmap(norm(np.mean([min_NormInt,max_NormInt]))))
        center_color=matplotlib.colors.rgb2hex(cmap(norm(cent_NormInt)))
        max_color=matplotlib.colors.rgb2hex(cmap(norm(max_NormInt)))  

        NODE_BORDER_PAINT=mapVisualProperty('NODE_BORDER_PAINT','continuous','NormInt',\
                             lower=[min_NormInt,min_color],center=[np.mean([min_NormInt,max_NormInt]),center_color],upper=[max_NormInt,max_color],\
                                           host=host, port=port)
        update_style("dataStyle",mappings=[NODE_BORDER_PAINT],\
                    host=host, port=port)
        response=cytoscape("vizmap", "apply", {"styles":"dataStyle"},\
                          host=host, port=port)

        NODE_BORDER_WIDTH=mapVisualProperty('NODE_BORDER_WIDTH','continuous','NormInt',\
                          lower=[min_NormInt,2],center=[np.mean([min_NormInt,max_NormInt]),4],upper=[max_NormInt,8],\
                                           host=host, port=port)
        update_style("dataStyle",mappings=[NODE_BORDER_WIDTH],\
                    host=host, port=port)
        response=cytoscape("vizmap", "apply", {"styles":"dataStyle"},\
                          host=host, port=port)

        response=cytoscape("network","rename",\
                           {"name":'main String network'},\
                          host=host, port=port)
    
    # create network with edges only
    response=cytoscape("network","select",\
                   {"edgeList":"all",\
                   "extendEdges":"true"},\
                  host=host, port=port)
    
    response=cytoscape("network","create",\
                   {"source":"current",\
                    "nodeList":"selected"},\
                  host=host, port=port)
    
    response=cytoscape("network","rename",\
                       {"name":'main String network (edges only)'},\
                      host=host, port=port)
    
    # top 10 changed genes > first neighbours

    response=cytoscape("network","set current",
                       {"network":"main String network (edges only)"},\
                      host=host, port=port)
    log2fcDf = getTableColumns('node',['log2(fold_change)'],host=host, port=port)
    if 'log2(fold_change)' in log2fcDf.columns.tolist():
        log2fcDf['log2(fold_change)']=log2fcDf['log2(fold_change)'].apply(lambda x: abs(x))
        log2fcDf=log2fcDf.sort_values(by=['log2(fold_change)'],ascending=False)
        top_nodes=log2fcDf.index.tolist()[:int(len(log2fcDf)*.10)]
        response=cytoscape("network","select",
                           {"nodeList":"name:"+",".join(top_nodes)},\
                          host=host, port=port)
        response=cytoscape("network","select",
                           {"firstNeighbors":"",\
                            "direction":"any",\
                            "network":"current"},\
                           host=host, port=port)
        response=cytoscape("network","create",
                           {"source":"current",\
                            "nodeList":"selected"},\
                          host=host, port=port)
        response=cytoscape("network","select",
                           {"edgeList":"all",\
                           "extendEdges":"true"},\
                          host=host, port=port)
        response=cytoscape("network","delete",
                           {"nodeList":"unselected"},\
                          host=host, port=port)
        response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                          host=host,port=port)
        response=cytoscape("layout", "force-directed",\
                          host=host, port=port)
        response=cytoscape("network","rename",\
                           {"name":'top '+str(int(len(log2fcDf)*.10))+' changed firstNeighbors'},\
                          host=host, port=port)
    
        #top 10 changed genes difusion
        response=cytoscape("network","set current",
                           {"network":"main String network (edges only)"},\
                          host=host, port=port)
        response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                          host=host, port=port)
        response=cytoscape("network","select",
                           {"nodeList":"name:"+",".join(top_nodes)},\
                          host=host, port=port)
        response=cytoscape("diffusion","diffuse",host=host, port=port)
        response=cytoscape("network","create",
                           {"source":"current",\
                          "nodeList":"selected"},\
                          host=host, port=port)
        response=cytoscape("network","select",
                           {"edgeList":"all",\
                           "extendEdges":"true"},\
                          host=host, port=port)
        response=cytoscape("network","delete",
                           {"nodeList":"unselected"},\
                          host=host, port=port)
        response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                          host=host, port=port)
        response=cytoscape("layout", "force-directed",host=host, port=port)
        response=cytoscape("network","rename",\
                           {"name":'top '+str(int(len(log2fcDf)*.10))+' changed diffusion'},\
                          host=host, port=port)
    
    def MAKETMP():
        (fd, f) = tempfile.mkstemp()
        f="/tmp/"+f.split("/")[-1]
        return f


    cys=MAKETMP()
    cyjs=MAKETMP()
    main_png=MAKETMP()
    main_pdf=MAKETMP()
    edg_png=MAKETMP()
    edg_pdf=MAKETMP()
    neig_png=MAKETMP()
    neig_pdf=MAKETMP()
    dif_png=MAKETMP()
    dif_pdf=MAKETMP()
    
    response=cytoscape("session", "save as" , \
                   {"file":cys},\
                 host=host, port=port)

    response=cytoscape("network", "export" , \
                       {"options":'CYJS',\
                        "OutputFile":cyjs},\
                      host=host, port=port)

    response=cytoscape("network","set current",
                       {"network":"main String network"},\
                      host=host, port=port)
    response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                      host=host, port=port)
    sleep(5)
    response=cytoscape("view", "export" , \
                       {"options":"PNG",\
                        "OutputFile":main_png},\
                      host=host, port=port)

    response=cytoscape("view", "export" , \
                       {"options":"PDF",\
                        "OutputFile":main_pdf},\
                      host=host, port=port)

    
    response=cytoscape("network","set current",
                   {"network":"main String network (edges only)"},\
                  host=host, port=port)
    response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                      host=host, port=port)
    sleep(5)
    response=cytoscape("view", "export" , \
                       {"options":"PNG",\
                        "OutputFile":edg_png},\
                      host=host, port=port)

    response=cytoscape("view", "export" , \
                       {"options":"PDF",\
                        "OutputFile":edg_pdf},\
                      host=host, port=port)
    

    try:
        response=cytoscape("network","set current",
                           {"network":'top '+str(int(len(log2fcDf)*.10))+' changed firstNeighbors'},\
                          host=host, port=port)
        response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                          host=host, port=port)
        sleep(5)

        response=cytoscape("view", "export" , \
                           {"options":"PNG",\
                            "OutputFile":neig_png},\
                          host=host, port=port)

        response=cytoscape("view", "export" , \
                           {"options":"PDF",\
                            "OutputFile":neig_pdf},\
                          host=host, port=port)
    except:
        print("No "+"changed firstNeighbors")
        sys.stdout.flush()

    try:
        response=cytoscape("network","set current",
                           {"network":'top '+str(int(len(log2fcDf)*.10))+' changed diffusion'},\
                          host=host, port=port)
        response=cytoscape("network","deselect",{"edgeList":"all", "nodeList":"all"},\
                          host=host, port=port)
        sleep(5)

        response=cytoscape("view", "export" , \
                           {"options":"PNG",\
                            "OutputFile":dif_png},\
                          host=host, port=port)

        response=cytoscape("view", "export" , \
                           {"options":"PDF",\
                            "OutputFile":dif_pdf},\
                          host=host, port=port)
    except:
        print("No "+"changed diffusion")
        sys.stdout.flush()    

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host)

    ftp_client=ssh.open_sftp()

    for f, extension, local in zip([cys,cyjs,main_png,main_pdf,edg_png,edg_pdf,neig_png,neig_pdf,dif_png,dif_pdf],\
                                    [".cys",".cyjs",".png",".pdf",".png",".pdf",".png",".pdf",".png",".pdf" ],\
                                    [target+".cys",target+".cyjs",target+".main.png",target+".main.pdf",\
                                    target+".main.edges.png",target+".main.edges.pdf",\
                                    target+".topFirstNeighbors.png",target+".topFirstNeighbors.pdf",\
                                    target+".topDiffusion.png",target+".topDiffusion.pdf"]):
        try:
            ftp_client.get(f+extension,local)
            ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+f+extension )
        except:
            print("No "+local)
            sys.stdout.flush()
    
