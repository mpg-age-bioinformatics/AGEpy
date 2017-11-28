import requests
import sys
import json

def CheckResponse(r):
    if str(r.status_code) == "200":
        print "Finished"
        sys.stdout.flush()
    else:
        print r
        sys.stdout.flush()
        
def cytoscape(t,namespace,command,PARAMS,host="localhost",port=1234):
    """
    General function for interacting with Cytoscape API.
    
    :param t: type of http call, ie. "POST" or "GET".
    :param namespace: namespace where the request should be executed. eg. "string"
    :param commnand: command to execute. eg. "protein query"
    :param PARAM: a dictionary with the parameters. Check your swagger normaly running on
    http://localhost:1234/v1/swaggerUI/swagger-ui/index.html?url=http://localhost:1234/v1/commands/swagger.json
    
    :returns: For "POST" the data in the content's response. For "GET" None. 
    
    eg.

    cytoscape("G","string","pubmed query",{"pubmed":"p53 p21","limit":"50"})


    """     
    
    if (t == "GET") or (t == "G"):
        P=[]
        for p in PARAMS.keys():
            v=str(PARAMS[p])
            v=v.replace(" ","%0A")
            P.append(str(p)+"="+v)
        P="&".join(P)
        URL="http://"+str(host)+":"+str(port)+"/v1/commands/"+str(namespace)+"/"+str(command)+"?"+P
        print "'"+URL+"'"
        sys.stdout.flush()
        r = requests.get(url = URL)
        CheckResponse(r)
        res=None
    
    elif (t == "POST") or (t == "P"):
        URL="http://"+str(host)+":"+str(port)+"/v1/commands/"+str(namespace)+"/"+str(command)
        r = requests.post(url = URL, json = PARAMS)
        CheckResponse(r)
        res=r.content
        res=json.loads(res)
        res=res["data"]
    
    return res
