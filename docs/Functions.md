Help on module AGEpy.AGEpy in AGEpy:

NAME
    AGEpy.AGEpy - Bioinformatics tools developed at the Max Planck Institute for Biology of Ageing

FUNCTIONS

    CellPlot(df, output_file=None, gene_expression='log2FC', figure_title='CellPlot', pvalCol='elimFisher', lowerLimit=None, upperLimit=None, colorBarType='Spectral')
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
    
    DAVIDenrich(database, categories, user, ids, ids_bg=None, name='', name_bg='', verbose=False, p=0.1, n=2)
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
    
    GTFtoBED(inGTF, name)
        Transform a GTF dataframe into a bed dataframe
        
        :param inGTF: GTF dataframe for transformation
        :param name: field of the GTF data frame to be use for the bed 'name' positon
        
        returns: a bed dataframe with the corresponding bed fiels: 'chrom','chromStart','chromEnd','name','score','strand'
    
    KEGGmatrix(organism, dataset, database, query_attributes=['ensembl_gene_id', 'kegg_enzyme'], host='www.ensembl.org', links=True, dfexp=None, kegg_db=None)
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
    
    RattributesBM(dataset, database, host='www.ensembl.org')
        Lists BioMart attributes through a RPY2 connection.
        
        :param dataset: a dataset listed in RdatasetsBM()
        :param database: a database listed in RdatabasesBM()
        :param host: address of the host server, default='www.ensembl.org'
        
        :returns: nothing
    
    RdatabasesBM(host='www.ensembl.org')
        Lists BioMart databases through a RPY2 connection.
        
        :param host: address of the host server, default='www.ensembl.org'
        
        :returns: nothing
    
    RdatasetsBM(database, host='www.ensembl.org')
        Lists BioMart datasets through a RPY2 connection.
        
        :param database: a database listed in RdatabasesBM()
        :param host: address of the host server, default='www.ensembl.org'
        
        :returns: nothing
    
    RfiltersBM(dataset, database, host='www.ensembl.org')
        Lists BioMart filters through a RPY2 connection.
        
        :param dataset: a dataset listed in RdatasetsBM()
        :param database: a database listed in RdatabasesBM()
        :param host: address of the host server, default='www.ensembl.org'
        
        :returns: nothing
    
    RqueryBM(query_filter, query_items, query_attributes, dataset, database, host='www.ensembl.org')
        Queries BioMart.
        
        :param query_filtery: one BioMart filter associated with the items being queried
        :param query_items: list of items to be queried (must assoiate with given filter)
        :param query_attributes: list of attributes to recover from BioMart  
        :param dataset: dataset to query
        :param database: database to query
        :param host: address of the host server, default='www.ensembl.org'
        
        return: a Pandas dataframe of the queried attributes
    
    SAMflags(x)
        Explains a SAM flag.
        
        :param x: flag
        
        :returns: complete SAM flag explanaition
    
    SymPlot(df, output_file=None, figure_title='SymPlot', pvalCol='elimFisher')
        Python implementation of the SymPlot from the CellPlot package for R.
        -inf or inf enrichments will come out as min found float or max found float, respectively.    
        
        :param df: pandas dataframe with the following columns - 'Enrichment', 'Significant', 'Annotated', 'Term', and 'log2fc'.
                   For log2fc each cell must contain a comma separated string with the log2fc for the genes enriched in the respective term. 
                   eg. '-inf,-1,2,3.4,3.66,inf'
        :param output_file: prefix for an output file. If given it witll create output_file.SymPlot.svg and output_file.SymPlot.png 
        :param figure_title: Figure title.
        :param pvalCol: name of the column containing the p values to determine if the terms should be marked as NS - not significant, use None for no marking 
        :returns: a matplotlib figure
    
    attributesBM(dataset, host='http://www.ensembl.org/biomart')
        Lists BioMart attributes for a specific dataset.
        
        :param dataset: dataset to list attributes of.
        :param host: address of the host server, default='http://www.ensembl.org/biomart'
        
        :returns: nothing
    
    attributesGTF(inGTF)
        List the type of attributes in a the attribute section of a GTF file
        
        :param inGTF: GTF dataframe to be analysed
        :returns: a list of attributes present in the attribute section
    
    biomaRtTOkegg(df)
        Transforms a pandas dataframe with the columns 'ensembl_gene_id','kegg_enzyme' 
        to dataframe ready for use in ...
        
        :param df: a pandas dataframe with the following columns: 'ensembl_gene_id','kegg_enzyme' 
        
        :returns: a pandas dataframe with the following columns: 'ensembl_gene_id','kegg_enzyme'
    
    databasesBM(host='http://www.ensembl.org/biomart')
        Lists BioMart databases.
        
        :param host: address of the host server, default='http://www.ensembl.org/biomart'
        
        :returns: nothing
    
    databasesKEGG(organism, ens_ids)
        Finds KEGG database identifiers for a respective organism given example ensembl ids.
        
        
        :param organism: an organism as listed in organismsKEGG()
        :param ens_ids: a list of ensenbl ids of the respective organism
        
        :returns: nothing if no database was found, or a string if a database was found
    
    datasetsBM(host='http://www.ensembl.org/biomart')
        Lists BioMart datasets.
        
        :param host: address of the host server, default='http://www.ensembl.org/biomart'
        
        :returns: nothing
    
    ecs_idsKEGG(organism)
        Uses KEGG to retrieve all ids and respective ecs for a given KEGG organism
        
        :param organism: an organisms as listed in organismsKEGG()
        
        :returns: a Pandas dataframe of with 'ec' and 'KEGGid'.
    
    ensembl_to_kegg(organism, kegg_db)
        Looks up KEGG mappings of KEGG ids to ensembl ids
        
        :param organism: an organisms as listed in organismsKEGG()
        :param kegg_db: a matching KEGG db as reported in databasesKEGG
        
        :returns: a Pandas dataframe of with 'KEGGid' and 'ENSid'.
    
    expKEGG(organism, names_KEGGids)
        Gets all KEGG pathways for an organism
        
        :param organism: an organism as listed in organismsKEGG()
        :param names_KEGGids: a Pandas dataframe with the columns 'gene_name': and  'KEGGid' as reported from idsKEGG(organism) (or a subset of it).
        
        :returns df: a Pandas dataframe with 'KEGGid','pathID(1):pathNAME(1)', 'pathID(n):pathNAME(n)'
        :returns paths: a list of retrieved KEGG pathways
    
    filtersBM(dataset, host='http://www.ensembl.org/biomart')
        Lists BioMart filters for a specific dataset.
        
        :param dataset: dataset to list filters of.
        :param host: address of the host server, default='http://www.ensembl.org/biomart'
        
        :returns: nothing
    
    getFasta(opened_file, sequence_name)
        Retrieves a sequence from an opened multifasta file
        
        :param opened_file: an opened multifasta file eg. opened_file=open("/path/to/file.fa",'r+')
        :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
        
        returns: a string with the sequence of interest
    
    getFileFormat(path)
        Return the file format
        
        :param path: The path to the file
        
        :returns: None, if file is missing, else one of the strings 'xlsx', 'xls', 'txt'
    
    id_nameDAVID(df, GTF=None, name_id=None)
        Given a DAVIDenrich output it converts ensembl gene ids to genes names and adds this column to the output
        
        :param df: a dataframe output from DAVIDenrich
        :param GTF: a GTF dataframe from readGTF()
        :param name_id: instead of a gtf dataframe a dataframe with the columns 'gene_name' and 'gene_id' can be given as input
        
        :returns: a pandas dataframe with a gene name column added to it.
    
    idsKEGG(organism)
        Uses KEGG to retrieve all ids for a given KEGG organism
        
        :param organism: an organism as listed in organismsKEGG()
        
        :returns: a Pandas dataframe of with 'gene_name' and 'KEGGid'.
    
    organismsKEGG()
        Lists all organisms present in the KEGG database.
        
        :returns: a list of lists containing one organism per list.
    
    parseGTF(inGTF)
        Reads an extracts all attributes in the attributes section of a GTF and constructs a new dataframe wiht one collumn per attribute instead of the attributes column
        
        :param inGTF: GTF dataframe to be parsed
        :returns: a dataframe of the orignal input GTF with attributes parsed.
    
    pathwaysKEGG(organism)
        Retrieves all pathways for a given organism.
        
        :param organism: an organism as listed in organismsKEGG()
        
        :returns df: a Pandas dataframe with the columns 'KEGGid','pathIDs', and 'pathName'.
        :returns df_: a Pandas dataframe with a columns for 'KEGGid', and one column for each pathway with the corresponding gene ids below
    
    queryBM(query_filter, query_items, query_attributes, query_dataset, query_dic=None, host='http://www.ensembl.org/biomart')
        Queries BioMart.
        
        :param query_filtery: one BioMart filter associated with the items being queried
        :param query_items: list of items to be queried (must assoiate with given filter)
        :param query_querydic: for complex queries this option should be used instead of 'filters' and 'items' and a dictionary of filters provided here eg. querydic={"filter1":["item1","item2"],"filter2":["item3","item4"]}. If using querydic, don't query more than 350 items at once. 
        :param query_attributes: list of attributes to recover from BioMart  
        :param query_dataset: dataset to query
        :param host: address of the host server, default='http://www.ensembl.org/biomart'
        
        :returns: a Pandas dataframe of the queried attributes
    
    readDataFrame(path, sheet=None, sep='\t')
        Returns a pandas data frame
        
        :param path: The path to the file
        :param sheet: Sheet name or integer 0-based index for xls[x] files, None for all
        :param sep: A separator for text format
        
        :returns: A pandas data frame
    
    readGTF(infile)
        Reads a GTF file and labels the respective columns in agreement with GTF file standards:
        'seqname','source','feature','start','end','score','strand','frame','attribute'.
        
        :param infile: path/to/file.gtf
        :returns: a Pandas dataframe of the respective GTF
    
    readSAM(SAMfile, header=False)
        Reads and parses a sam file.
        
        :param SAMfile: /path/to/file.sam
        :param header: logical, if True, reads the header information
        
        :returns: a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' and a list of the headers if header=True
    
    retrieve_GTF_field(field, gtf)
        Returns a field of choice from the attribute column of the GTF
        
        :param field: field to be retrieved
        :returns: a Pandas dataframe with one columns containing the field of choice
    
    rewriteFasta(sequence, sequence_name, fasta_in, fasta_out)
        Rewrites a specific sequence in a multifasta file while keeping the sequence header.
        
        :param sequence: a string with the sequence to be written  
        :param sequence_name: the name of the sequence to be retrieved eg. for '>2 dna:chromosome chromosome:GRCm38:2:1:182113224:1 REF' use: sequence_name=str(2)
        :param fasta_in: /path/to/original.fa
        :param fasta_out: /path/to/destination.fa
        
        :returns: nothing
    
    writeBED(inBED, file_path)
        Writes a bed dataframe into a bed file.
        Bed format: 'chrom','chromStart','chromEnd','name','score','strand'
        
        :param inBED: bed dataframe to be written.
        :param file_path: /path/to/file.bed
        
        :returns: nothing
    
    writeFasta(sequence, sequence_name, output_file)
        Writes a fasta sequence into a file.
        
        :param sequence: a string with the sequence to be written
        :param sequence_name: name of the the fasta sequence
        :param output_file: /path/to/file.fa to be written
        
        :returns: nothing
    
    writeGTF(inGTF, file_path)
        Write a GTF dataframe into a file
        
        :param inGTF: GTF dataframe to be written. It should either have 9 columns with the last one being the "attributes" section or more than 9 columns where all columns after the 8th will be colapsed into one.
        :param file_path: path/to/the/file.gtf
        :returns: nothing
    
    writeSAM(sam, SAMfile, header=None)
        Writes a pandas dataframe with the respective SAM columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL' into a sam file
        
        :param sam: pandas dataframe to be writen
        :param SAMfile: /path/to/file.sam
        
        :returns: nothing

DATA
    DNAcode = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T',...
    __warningregistry__ = {("Not importing directory 'pandas': missing __i...
    biomart_host = 'http://www.ensembl.org/biomart'
    david_categories = ['GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT',...
    david_fields = ['categoryName', 'termName', 'listHits', 'percent', 'ea...
    rbiomart_host = 'www.ensembl.org'


