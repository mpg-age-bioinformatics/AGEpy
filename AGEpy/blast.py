import requests
import itertools
import pandas as pd
import sys

def variablename(var):
    """
    Returns the string of a variable name.
    """
    s=[tpl[0] for tpl in itertools.ifilter(lambda x: var is x[1], globals().items())]
    s=s[0].upper()
    return s

def BLASTquery(query,database,program,filter=None,\
               format_type=None, expect=None,\
               nucl_reward=None, nucl_penalty=None,\
               gapcosts=None, matrix=None,\
               hitlist_size=None, descriptions=None,\
               alignments=None,\
               ncbi_gi=None, threshold=None,\
               word_size=None, composition_based_statistics=None,\
               organism=None, others=None,\
               num_threads=None, baseURL="http://blast.ncbi.nlm.nih.gov",\
              verbose=False):
    """
    Performs a blast query online.

    As in https://ncbi.github.io/blast-cloud/

    :param query: Search query. Allowed values: Accession, GI, or FASTA.
    :param database: BLAST database. Allowed values: nt, nr, refseq_rna, refseq_protein, swissprot, pdbaa, pdbnt
    :param program: BLAST program. Allowed values:  blastn, megablast, blastp, blastx, tblastn, tblastx
    :param filter: Low complexity filtering. Allowed values: F to disable. T or L to enable. Prepend "m" for mask at lookup (e.g., mL)
    :param format_type: Report type. Allowed values: HTML, Text, XML, XML2, JSON2, or Tabular. HTML is the default.
    :param expect: Expect value. Allowed values: Number greater than zero.
    :param nucl_reward: Reward for matching bases (BLASTN and megaBLAST). Allowed values: Integer greater than zero.
    :param nucl_penalty: Cost for mismatched bases (BLASTN and megaBLAST). Allowed values: Integer less than zero.
    :param gapcosts: Gap existence and extension costs. Allowed values: Pair of positive integers separated by a space such as "11 1".
    :param matrix: Scoring matrix name. Allowed values: One of BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM250, PAM30 or PAM70. Default: BLOSUM62 for all applicable programs.
    :param hitlist_size: Number of databases sequences to keep. Allowed values: Integer greater than zero.
    :param descriptions: Number of descriptions to print (applies to HTML and Text). Allowed values: Integer greater than zero.
    :param alignments: Number of alignments to print (applies to HTML and Text). Allowed values: Integer greater than zero.
    :param ncbi_gi: Show NCBI GIs in report. Allowed values: T or F.
    :param threshold: Neighboring score for initial words. Allowed values: Positive integer (BLASTP default is 11). Does not apply to BLASTN or MegaBLAST).
    :param word_size: Size of word for initial matches. Allowed values: Positive integer.
    :param composition_based_statistics: Composition based statistics algorithm to use. Allowed values: One of 0, 1, 2, or 3. See comp_based_stats command line option in the BLAST+ user manual for details.
    :param organism: an organism as in https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
    :param others: here you can add other parameters as seen in a blast bookmarked page. Define you query in https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
            Once your query is defined click on "Bookmark" on right upper side of the page. You can copy fragments of the URL
            which define the query. Eg. For organism "Homo sapiens (taxid:9606)" you will see the string "EQ_MENU=Homo%20sapiens%20%28taxid%3A9606%29" - this is
            the string you can use here in others.
    :param num_threads: Number of virtual CPUs to use. Allowed values: Integer greater than zero (default is 1). Supported only on the cloud.
    :param verbose: print more

    :returns: BLAST search request identifier
    """

    if organism:
        organism=organism.replace(" ", "%20").replace("(", "%28").replace(")", "%29").replace(":", "%3A")
        EQ_MENU=organism
    else:
        EQ_MENU=None

    URL=baseURL+"/Blast.cgi?"
    URL=URL+"QUERY="+str(query)+"&DATABASE="+str(database)+"&PROGRAM="+str(program)
    for o,varname in zip([filter, format_type, expect, nucl_reward, nucl_penalty,\
              gapcosts, matrix, hitlist_size, descriptions, alignments,\
              ncbi_gi, threshold, word_size, composition_based_statistics,\
              EQ_MENU, num_threads],\
              ['FILTER' , 'FORMAT_TYPE', 'EXPECT', 'NUCL_REWARD', 'NUCL_PENALTY',\
              'GAPCOSTS', 'MATRIX', 'HITLIST_SIZE', 'DESCRIPTIONS', 'ALIGNMENTS',\
              'NCBI_GI', 'THRESHOLD', 'WORD_SIZE', 'COMPOSITION_BASED_STATISTICS',\
              'EQ_MENU', 'NUM_THREADS']):
        if o:
            URL=URL+"&"+ varname +"="+str(o)

    if others:
        URL=URL+"&"+others

    URL=URL+"&CMD=Put"

    if verbose:
        print(URL)
        sys.stdout.flush()

    response=requests.get(url = URL)
    r=response.content.split("\n")
    RID=[ s for s in r if "RID = " in s ]
    if len(RID) > 0:
        RID=RID[0].split(" ")[-1]
    else:
        print("Could not return an RID for this query.")
        RID=None
    return RID

def BLASTcheck(rid,baseURL="http://blast.ncbi.nlm.nih.gov"):
    """
    Checks the status of a query.

    :param rid: BLAST search request identifier. Allowed values: The Request ID (RID) returned when the search was submitted
    :param baseURL: server url. Default=http://blast.ncbi.nlm.nih.gov

    :returns status: status for the query.
    :returns therearehist: yes or no for existing hits on a finished query.
    """

    URL=baseURL+"/Blast.cgi?"
    URL=URL+"FORMAT_OBJECT=SearchInfo&RID="+rid+"&CMD=Get"
    response=requests.get(url = URL)
    r=response.content.split("\n")
    try:
        status=[ s for s in r if "Status=" in s ][0].split("=")[-1]
        ThereAreHits=[ s for s in r if "ThereAreHits=" in s ][0].split("=")[-1]
    except:
        status=None
        ThereAreHits=None

    print(rid, status, ThereAreHits)
    sys.stdout.flush()

    return status, ThereAreHits

def BLASTresults(rid, format_type="Tabular", \
                 hitlist_size= None, alignments=None, \
                 ncbi_gi = None, format_object=None,\
                 baseURL="http://blast.ncbi.nlm.nih.gov"):
    """
    Retrieves results for an RID.

    :param rid: BLAST search request identifier. Allowed values: The Request ID (RID) returned when the search was submitted
    :param format_type: Report type. Allowed values: HTML, Text, XML, XML2, JSON2, or Tabular. Tabular is the default.
    :param hitlist_size: Number of databases sequences to keep. Allowed values: Integer greater than zero.
    :param alignments: Number of alignments to print (applies to HTML and Text). Allowed values: Integer greater than zero.
    :param ncbi_gi: Show NCBI GIs in report. Allowed values: T or F.
    :param format_object: Object type. Allowed values: SearchInfo (status check) or Alignment (report formatting).
    :param baseURL: server url. Default=http://blast.ncbi.nlm.nih.gov

    :returns: the result of a BLAST query. If format_type="Tabular" it will parse the content into a Pandas dataframe.
    """

    URL=baseURL+"/Blast.cgi?"
    URL=URL+"RID="+str(rid)+"&FORMAT_TYPE="+str(format_type)
    for o in [ hitlist_size, alignments,\
              ncbi_gi, format_object]:
        if o:
            URL=URL+"&"+ variablename(var) +"="+str(o)
    URL=URL+"&CMD=Get"
    response=requests.get(url = URL)
    response=response.content

    if format_type=="Tabular":
        result=response.split("\n")
        result=[ s.split("\t") for s in result][6:]
        header=result[:7]
        content=result[7:]
        fields=header[5][0].strip("# Fields: ").split(", ")
        result=pd.DataFrame(content,columns=fields)
        response=result[:int(header[-1][0].split(" ")[1])]

    return response
