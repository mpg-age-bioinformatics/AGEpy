import pandas as pd
import os
import sys
import AGEpy as age
import numpy as np
import requests
from io import StringIO
from itertools import islice
from urllib.parse import urlencode
from biomart import BiomartServer
import time
import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objects as go
import networkx as nx
import matplotlib
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch
import plotly.graph_objects as go

StringDB="/nexus/posix0/MAGE-flaski/service/databases/StringDB/"
StringDB_proteins_file="protein.info.v12.0.txt"
biomarthost="http://dec2021.archive.ensembl.org/biomart/"
taxons={"caenorhabditis elegans":"6239","drosophila melanogaster":"7227",\
        "mus musculus":"10090","homo sapiens":"9606", "saccharomyces cerevisiae": "4932", "danio rerio":"7955"}
tags={"caenorhabditis elegans":"CEL","drosophila melanogaster":"DMEL",\
      "mus musculus":"MUS","homo sapiens":"HSA"}


def get_evidence (species, tags=tags) :
    """
    Retrieve aging-related gene evidence for a given species.

    This function loads (or creates, if missing) a species-specific aging
    evidence table derived from GO annotations via AGEpy/biomart, then
    extracts the list of genes associated with aging evidence.

    Parameters
    ----------
    species : str
        Species name (e.g. "caenorhabditis elegans", "mus musculus").
        Must be present as a key in the 'tags' dictionary.

    tags : dict, optional
        Mapping from species names to organism-specific STRING/AGEpy tags
        (default: global 'tags').

    Returns
    -------
    aging_genes : list of str
        List of gene names for the given species that have aging-related
        evidence.
    """
    
    if species in tags.keys():
        organismtag=tags[species]
        taxon_id=taxons[species]
        species_name=species.replace(' ', '_')
        
        if not os.path.isfile(f"{StringDB}{species_name}_{taxon_id}_homdf.txt"):
            print("Could not find ageing evidence table. Using biomart to create one.")
            sys.stdout.flush()
            homdf,HSA,MUS,CEL,DMEL=age.FilterGOstring(host=biomarthost)
            homdf.to_csv(f"{StringDB}{species_name}_{taxon_id}_homdf.txt", sep="\t", index=None)
        else:
            print("Found existing ageing evidence table.")
            sys.stdout.flush()
            homdf=pd.read_csv(f"{StringDB}{species_name}_{taxon_id}_homdf.txt", sep="\t", low_memory=False)
        aging_genes=homdf[[organismtag+"_ensembl_gene_id","evidence"]].dropna()
        aging_genes=aging_genes[aging_genes[organismtag+"_ensembl_gene_id"]!=""]
        aging_genes=aging_genes[organismtag+"_ensembl_gene_id"].tolist()
    
        aging_genes=homdf[[organismtag+"_external_gene_name","evidence"]].dropna()
        aging_genes=aging_genes[aging_genes[organismtag+"_external_gene_name"]!="None"]
        aging_genes=aging_genes[organismtag+"_external_gene_name"].tolist()

    return aging_genes


def prepare_query (fin, taxon_id, aging_genes=None, ntop=1000):
    """
    Prepare query genes for STRING network analysis from a results file.

    Parameters
    ----------
    fin : str
        Path to a tab-separated differential expression results file
        containing 'ensembl_gene_id', 'gene_name', 'padj', and 'baseMean'.

    taxon_id : str
        NCBI taxonomy identifier used to select the STRING protein
        annotation file.

    aging_genes : iterable of str, optional
        Set of Ensembl gene IDs considered as aging-related. These genes
        are labeled as "aging_gene" in the output.

    ntop : int, optional (default=1000)
        Maximum number of significant genes (sorted by 'padj') to return
        for STRING querying.

    Returns
    -------
    dfin : pandas.DataFrame
        Input DataFrame with added 'evidence' and log10-transformed
        'baseMean' columns.

    query_genes : list of str
        Gene names passing 'padj < 0.05', present in STRING, limited to
        the top 'ntop' entries.
    """
        
    def CheckEvidence(x):
        if aging_genes is not None and x in aging_genes:
            return "aging_gene"
        return "no"

    dfin=pd.read_csv(fin, sep="\t")
    dfin["evidence"]=dfin["ensembl_gene_id"].apply(lambda x:CheckEvidence(x) )
    dfin["baseMean"]=dfin["baseMean"].apply(lambda x: np.log10(x))
    qdf=dfin[dfin["padj"]<0.05]
    if qdf.shape[0] == 0:
            sys.exit()
    qdf=qdf.sort_values(by=["padj"],ascending=True)
    query_genes=qdf["gene_name"].dropna().tolist()

    StringDB_proteins=pd.read_csv(f"{StringDB}{taxon_id}.{StringDB_proteins_file}", sep="\t", low_memory=False).rename(columns={'#string_protein_id':'string_protein_id'})
    string_gene_set = StringDB_proteins['preferred_name'].tolist()
    query_genes = [g for g in query_genes if g in string_gene_set][:ntop]
    
    return dfin, query_genes



# API query
# Query STRING DB interactions

def API_string_network(query_genes, taxon_id, batch_size=50, limit=10, min_score=400):
    """
    Retrieve a protein–protein interaction network from the STRING database
    using the STRING REST API (`network` endpoint), querying identifiers in
    batches to avoid request size limits.

    Parameters
    ----------
    query_genes : iterable of str
        List or iterable of gene_names

    taxon_id : int
        NCBI taxonomy identifier specifying the species to query
        (e.g. 9606 for Homo sapiens).

    batch_size : int, optional (default=50)
        Number of identifiers to include per API request. Identifiers are
        split into successive batches of this size and queried separately.

    limit : int
        Maximum number of interaction partners returned per protein, as
        defined by the STRING API `limit` parameter.

    min_score : int, optional (default=400)
        Minimum required interaction score (STRING combined score),
        corresponding to the API parameter `required_score`.

    Returns
    -------
    pandas.DataFrame
        A concatenated DataFrame containing all interaction edges returned
        by STRING across all batches, in the original STRING TSV format
        (one row per interaction).
    nodes_api : set[str]
        All nodes admitted to the network.
    

    Notes
    -----
    - This function queries the STRING endpoint:
      https://string-db.org/api/tsv/network
    - Identifiers within each batch are joined using the STRING-required
      `%0d` delimiter.
    - Each batch response is appended and concatenated into a single
      DataFrame.
    - HTTP errors are raised immediately if a request fails.
    """
    base_url = "https://string-db.org/api/tsv/network"
    all_dfs = []
    
    def batch(iterable, n):
        """Yield successive n-sized batches from iterable."""
        it = iter(iterable)
        while True:
            chunk = list(islice(it, n))
            if not chunk:
                return
            yield chunk

    for i, batch_ids in enumerate(batch(query_genes, batch_size), 1):
        identifiers = "%0d".join(batch_ids)
        params = {
            "identifiers": identifiers,
            "species": taxon_id,
            "limit": limit,
            "required_score": min_score,
            "caller_identity": "python_script"
        }
        r = requests.get(base_url, params=params)
        r.raise_for_status()
        df = pd.read_csv(StringIO(r.text), sep="\t", low_memory=False)
        all_dfs.append(df)

    all_dfs = [df for df in all_dfs if not df.empty]
    edges_api = pd.concat(all_dfs, ignore_index=True)
    nodes_api = (set(edges_api["preferredName_A"]) | set(edges_api["preferredName_B"]))
    
    return edges_api, nodes_api 

## Example usage: edges_api, nodes_api = API_string_network(query_genes, taxon_id, batch_size=2, limit=limit)



# Read local detailed table (create if does not exist)
def read_local_links(db_file, taxon_id):
    """
    Load a local STRING links file, creating an annotated version if it does
    not already exist.

    If the specified file does not exist, this function reads the corresponding
    raw STRING links file (expected to be the same path without the
    '_annotated' suffix), joins it with the global `StringDB_proteins`
    annotation table to recover preferred gene names, normalizes all STRING
    scores by dividing them by 1000, selects a standard set of STRING evidence
    columns, and writes the annotated result to disk.

    If the file already exists, it is loaded directly.
    
    Parameters
    ----------
    db_file : file path for links database.

    taxon_id : int
        NCBI taxonomy identifier (e.g. 9606 for Homo sapiens).

    """
    if not os.path.exists(db_file):
        print("Create and load annotated link database...")
        StringDB_proteins=pd.read_csv(f"{StringDB}{taxon_id}.{StringDB_proteins_file}", sep="\t", low_memory=False).rename(columns={'#string_protein_id':'string_protein_id'})
        StringDB_links = pd.read_csv(db_file.replace("_annotated.tsv", ".txt.gz"), sep=" ", low_memory=False)

        df_links_db = (
            StringDB_links
            .merge(
                StringDB_proteins,
                left_on="protein1",
                right_on="string_protein_id",
                how="left"
            )
            .rename(columns={"preferred_name": "preferredName_A"})
            .drop(columns="string_protein_id")
            .merge(
                StringDB_proteins,
                left_on="protein2",
                right_on="string_protein_id",
                how="left"
            )
            .rename(columns={"preferred_name": "preferredName_B"})
            .drop(columns="string_protein_id")
        )

        df_links_db["ncbiTaxonId"] = taxon_id

        # Rename columns to STRING-style names
        df_links_db = df_links_db.rename(
            columns={
                "protein1": "stringId_A", "protein2": "stringId_B", "neighborhood": "nscore",
                "fusion": "fscore", "cooccurence": "pscore", "coexpression": "ascore",
                "experimental": "escore", "database": "dscore", "textmining": "tscore",
                "combined_score": "score",
            }
        )

        # Normalize all score columns
        score_cols = [
            "nscore", "fscore", "pscore", "ascore",
            "escore", "dscore", "tscore", "score"
        ]
        df_links_db[score_cols] = df_links_db[score_cols] / 1000.0

        df_links_db = df_links_db[
            [
                "stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "ncbiTaxonId", "nscore",
                "fscore", "pscore", "ascore", "escore", "dscore", "tscore", "score",
            ]
        ]

        df_links_db.to_csv(db_file, sep="\t", index=None)

    else:
        print("Load annotated database:", db_file)
        df_links_db = pd.read_csv(db_file, sep="\t", low_memory=False)

    return df_links_db

## Example usage: db_links = read_local_links(f"{StringDB}{taxon_id}.protein.links.detailed.v12.0_annotated.tsv", taxon_id)



def local_string_network(df_links_db, query_genes, min_score=0.4, limit=10):
    """
    Approximate STRING /network API output from local flat files.

    Parameters
    ----------
    df_links_db : pd.DataFrame
        STRING links table, must contain 'preferredName_A', 'preferredName_B', 'score'.
    query_genes : list[str]
        List of query gene names.
    min_score : int, optional
        Minimum score threshold (0-1000), default 0.4.
    limit : int, optional
        Maximum number of neighbors per query gene (like API), default 50.

    Returns
    -------
    edges_local : pd.DataFrame
        All edges between admitted nodes with score >= min_score.
    nodes_local : set[str]
        All nodes admitted to the network.
    """

    # 1️⃣ Global score filter
    df = df_links_db[df_links_db["score"] >= min_score]

    # 2️⃣ Phase 1: Node admission per query gene
    nodes_local = set(query_genes)
    edges_seed_blocks = []

    for g in query_genes:
    # All candidate edges involving g
        df_g = df[
            (df["preferredName_A"] == g) |
            (df["preferredName_B"] == g)
        ].sort_values("score", ascending=False)
    
        new_nodes_added = 0
        retained_edges = []
    
        for _, row in df_g.iterrows():
            u = row["preferredName_A"]
            v = row["preferredName_B"]
    
            # Check whether this edge introduces a new node
            introduces_new = (
                (u not in nodes_local) or
                (v not in nodes_local)
            )
    
            if introduces_new:
                retained_edges.append(row)
                nodes_local.update([u, v])
                new_nodes_added += 1
    
                if new_nodes_added >= limit:
                    break
    
        if retained_edges:
            edges_seed_blocks.append(pd.DataFrame(retained_edges))

    # 3️⃣ Edge set: Complete edges among all admitted nodes
    edges_local = df[
        df["preferredName_A"].isin(nodes_local) &
        df["preferredName_B"].isin(nodes_local)
    ].drop_duplicates(subset=["preferredName_A", "preferredName_B"]).reset_index(drop=True)

    return edges_local, nodes_local

## Example usage: edges_local, nodes_local = local_string_network(df_links_db, query_genes, min_score=0.4, limit=limit)



def plot_network_RNAseq(
    edges_df,
    dfin,
    pdf_file,
    png_file,
    graphml_file,
    directions=None  # Optional DataFrame with "from" and "to"
):
    """
    Plot a gene interaction network with Cytoscape-compatible visual attributes
    and export it as both a static PDF and a GraphML file.

    This function builds an undirected NetworkX graph from an edge list,
    assigns node attributes based on differential expression data, computes
    a spring layout, and renders a publication-quality static network plot.
    Visual properties (node colors, sizes, shapes, positions, and edge styles)
    are stored as native node and edge attributes so that the exported GraphML
    file can be loaded into Cytoscape with matching aesthetics.

    Optionally, a separate set of directed relationships can be provided to
    draw arrows on top of the undirected network without altering the base
    graph structure.

    Parameters
    ----------
    edges_df : pandas.DataFrame
        DataFrame describing undirected network edges. Expected columns include:
        - `preferredName_A` : source node name
        - `preferredName_B` : target node name
        - `score` (optional) : edge weight
        - `dashed` (optional) : "yes"/"no" flag to control edge linestyle

    dfin : pandas.DataFrame
        DataFrame containing node-level attributes indexed by gene name.
        Required columns:
        - `gene_name`
        - `log2FoldChange`
        - `baseMean`
        - `evidence` ("yes"/"no")

    pdf_file : str
        Path to the output PDF file for the static network visualization.

    png_file : str
        Path to the output PNG file for the static network visualization.

    graphml_file : str
        Path to the output GraphML file for use in Cytoscape or other
        network visualization tools.

    directions : pandas.DataFrame, optional
        Optional DataFrame specifying directed relationships to be visualized
        as arrows. Must contain two columns:
        - `from` : source node name
        - `to`   : target node name
        These arrows are visual only and do not change the underlying
        undirected graph topology.

    Returns
    -------
    None
        This function does not return a value. It produces side effects by
        writing a PDF, a PNG and a GraphML file to disk.
    """
    # -----------------------------
    # 1️⃣ Build graph - create edges
    # -----------------------------
    G = nx.Graph()

    for _, r in edges_df.iterrows():
        u = str(r["preferredName_A"])
        v = str(r["preferredName_B"])
        val = r.get("score", 1.0)
        score = float(val.iloc[0] if hasattr(val, "iloc") else val)
        dashed = str(r["dashed"]) if "dashed" in r else "no" # only if 'dashed' column present in the dfin dataframe
        G.add_edge(u, v, score=score, dashed=dashed, arrow=False)

    # -----------------------------
    # 2️⃣ Node attributes
    # -----------------------------
    expr = (
        dfin
        .set_index("gene_name")[["log2FoldChange", "baseMean", "evidence"]]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0)
    )

    # attach information to the nodes (log2FC, baseMean and evidence) indexed by gene_name
    for n in G.nodes():
        n_str = str(n)
        if n_str in expr.index:
            lfc = expr.loc[n_str, "log2FoldChange"]
            bm  = expr.loc[n_str, "baseMean"]
            ev  = expr.loc[n_str, "evidence"]

            G.nodes[n]["log2FoldChange"] = float(lfc.iloc[0] if hasattr(lfc, "iloc") else lfc)
            G.nodes[n]["baseMean"] = float(bm.iloc[0] if hasattr(bm, "iloc") else bm)
            G.nodes[n]["evidence"] = str(ev.iloc[0] if hasattr(ev, "iloc") else ev)
        else:
            G.nodes[n]["log2FoldChange"] = 0.0
            G.nodes[n]["baseMean"] = 0.0
            G.nodes[n]["evidence"] = "no"

        G.nodes[n]["linewidth"] = 2.0

    # -----------------------------
    # 3️⃣ Layout
    # -----------------------------
    pos = nx.spring_layout(G, seed=42, iterations=50)
    pos_arr = np.array([pos[n] for n in G.nodes()])

    # -----------------------------
    # 4️⃣ Colormaps
    # -----------------------------
    # set scale for filling the nodes
    fc_norm = mcolors.Normalize(vmin=-4, vmax=4)
    fc_cmap = matplotlib.colormaps.get_cmap("bwr")

    # set scale for node borders based on baseMean value range
    bm_vals = np.array([G.nodes[n]["baseMean"] for n in G.nodes()])
    bm_norm = mcolors.Normalize(vmin=bm_vals.min(), vmax=bm_vals.max())
    bm_cmap = matplotlib.colormaps.get_cmap("Reds")

    # Actual node fill and border color
    fc_colors = np.array([
        fc_cmap(fc_norm(G.nodes[n]["log2FoldChange"]))
        for n in G.nodes()
    ])

    ec_colors = np.array([
        bm_cmap(bm_norm(G.nodes[n]["baseMean"]))
        for n in G.nodes()
    ])

    # Node size and shape based on ageing evidence, ellipse 60 if no evidence and diamond 120 if evidence
    sizes = np.array([
        120 if G.nodes[n]["evidence"] == "yes" else 60
        for n in G.nodes()
    ])

    markers = np.array([
        'diamond' if G.nodes[n]["evidence"] == "yes" else 'ellipse'
        for n in G.nodes()
    ])

    # -----------------------------
    # 5️⃣ Store Cytoscape attributes
    # -----------------------------
    for i, n in enumerate(G.nodes()):
        G.nodes[n]["fill_color"] = matplotlib.colors.to_hex(fc_colors[i])
        G.nodes[n]["stroke_color"] = matplotlib.colors.to_hex(ec_colors[i])
        G.nodes[n]["size"] = float(sizes[i])
        G.nodes[n]["shape"] = str(markers[i])
        G.nodes[n]["linewidth"] = 2.0
        G.nodes[n]["x"] = float(pos_arr[i, 0])
        G.nodes[n]["y"] = float(pos_arr[i, 1])

    # -----------------------------
    # 6️⃣ Draw edges
    # -----------------------------
    plt.figure(figsize=(12, 12))

    # Draw dashed or plain edges if 'dashed' column exists
    for u, v, d in G.edges(data=True):
        style = 'dashed' if d.get("dashed", "no") == "yes" else 'solid'
        plt.plot(
            [pos[u][0], pos[v][0]],
            [pos[u][1], pos[v][1]],
            color="black",
            linewidth=0.5,
            linestyle=style,
            zorder=1
        )

    # -----------------------------
    # 7️⃣ Draw arrows - only if directions df is provided
    # -----------------------------
    if directions is not None:
        for _, r in directions.iterrows():
            u, v = r["from"], r["to"]
            if u in pos and v in pos:
                plt.gca().add_patch(
                    FancyArrowPatch(
                        posA=pos[u],
                        posB=pos[v],
                        arrowstyle='-|>',
                        color='black',
                        linewidth=0.5,
                        mutation_scale=20,
                        zorder=2
                    )
                )

    # -----------------------------
    # 8️⃣ Draw nodes
    # -----------------------------
    for marker_type in np.unique(markers):
        idx = np.where(markers == marker_type)[0]
        plt.scatter(
            pos_arr[idx, 0],
            pos_arr[idx, 1],
            s=sizes[idx],
            marker='D' if marker_type == 'diamond' else 'o',
            c=fc_colors[idx],
            edgecolors=ec_colors[idx],
            linewidths=2,
            zorder=2
        )

    # -----------------------------
    # 9️⃣ Labels & export
    # -----------------------------
    nx.draw_networkx_labels(G, pos, font_size=6)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(pdf_file, dpi=300)
    plt.savefig(png_file, dpi=300)
    plt.close()

    nx.write_graphml(G, graphml_file)

    print("Network saved:\n", pdf_file, "\n", graphml_file)


def plot_network_interactive(edges_df):
    """
    Create a lightweight interactive network visualization from an edge list.

    This function builds an undirected graph from a simple edge DataFrame,
    computes a force-directed layout, and displays an interactive Plotly
    network. All nodes share the same visual properties and are labeled
    by their names.

    Parameters
    ----------
    edges_df : pandas.DataFrame
        DataFrame describing undirected network edges.
        Must contain the columns:
        - `preferredName_A`
        - `preferredName_B`

    Returns
    -------
    None
        This function does not return a value.
        It displays an interactive Plotly network.
    """
    import networkx as nx
    import plotly.graph_objects as go

    # -----------------------------
    # Build graph
    # -----------------------------
    G = nx.Graph()

    for _, r in edges_df.iterrows():
        u = str(r["preferredName_A"])
        v = str(r["preferredName_B"])
        G.add_edge(u, v)

    pos = nx.spring_layout(G, seed=42, k=0.7, iterations=100)

    # -----------------------------
    # Edges
    # -----------------------------
    edge_x, edge_y = [], []
    for u, v in G.edges():
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(color="black", width=0.5),
        hoverinfo="none"
    )

    # -----------------------------
    # Nodes
    # -----------------------------
    node_trace = go.Scatter(
        x=[pos[n][0] for n in G.nodes()],
        y=[pos[n][1] for n in G.nodes()],
        mode="markers+text",
        text=list(G.nodes()),
        textposition="middle center",
        hoverinfo="text",
        marker=dict(
            size=12,
            color="#cccccc",
            line=dict(width=1, color="black"),
            symbol="circle"
        )
    )

    fig = go.Figure([edge_trace, node_trace])
    fig.update_layout(
        width=1000,
        height=1000,
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=20, b=20)
    )

    fig.show()


def plot_network_RNAseq_interactive(edges_local, dfin, directions=None):
    """
    Create an interactive gene interaction network integrating RNA-seq data.

    This function builds an undirected NetworkX graph from a local edge list,
    attaches RNA-seq–derived attributes (log2FoldChange, baseMean expression,
    and aging evidence) to each node, computes a force-directed layout, and
    renders an interactive Plotly visualization.

    Node fill color encodes log2 fold change (blue–white–red scale),
    node border color encodes base mean expression (Reds scale),
    node size and shape encode aging evidence.
    Gene names are displayed as centered labels, while quantitative values
    are shown only on mouse hover.

    Parameters
    ----------
    edges_local : pandas.DataFrame
        DataFrame describing undirected network edges.
        Must contain the columns:
        - `preferredName_A`
        - `preferredName_B`

    dfin : pandas.DataFrame
        DataFrame containing RNA-seq results and node-level annotations.
        Required columns:
        - `gene_name`
        - `log2FoldChange`
        - `baseMean`
        - `evidence` ("yes"/"no")

    directions : pandas.DataFrame, optional
        Reserved for future use. Intended to specify directed relationships
        between nodes. Currently not used in the interactive visualization.

    Returns
    -------
    None
        This function does not return any object.
        It displays an interactive Plotly network in the current environment.
    """
    # -----------------------------
    # 1️⃣ Build graph
    # -----------------------------
    G = nx.Graph()

    for _, r in edges_local.iterrows():
        u = str(r["preferredName_A"])
        v = str(r["preferredName_B"])
        G.add_edge(u, v)

    expr = (
        dfin
        .set_index("gene_name")[["log2FoldChange", "baseMean", "evidence"]]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0)
    )

    # attach information to the nodes (log2FC, baseMean and evidence) indexed by gene_name
    for n in G.nodes():
        if n in expr.index:
            lfc = expr.loc[n, "log2FoldChange"]
            bm  = expr.loc[n, "baseMean"]
            ev  = expr.loc[n, "evidence"]

            G.nodes[n]["log2FoldChange"] = float(lfc.iloc[0] if hasattr(lfc, "iloc") else lfc)
            G.nodes[n]["baseMean"] = float(bm.iloc[0] if hasattr(bm, "iloc") else bm)
            G.nodes[n]["evidence"] = str(ev.iloc[0] if hasattr(ev, "iloc") else ev)
        else:
            G.nodes[n]["log2FoldChange"] = 0.0
            G.nodes[n]["baseMean"] = 0.0
            G.nodes[n]["evidence"] = "no"

    # -----------------------------
    # 2️⃣ Layout (spread out)
    # -----------------------------
    pos = nx.spring_layout(G, seed=42, k=0.7, iterations=100)

    # -----------------------------
    # 3️⃣ fill nodes and node borders according to log2FC and baseMean, respectively
    # -----------------------------
    # set scale for filling the nodes
    fc_norm = mcolors.Normalize(vmin=-4, vmax=4)
    fc_cmap = matplotlib.colormaps.get_cmap("bwr")

    # set scale for node borders based on baseMean value range
    bm_vals = np.array([G.nodes[n]["baseMean"] for n in G.nodes()])
    bm_norm = mcolors.Normalize(vmin=bm_vals.min(), vmax=bm_vals.max())
    bm_cmap = matplotlib.colormaps.get_cmap("Reds")

    # Actual node fill and border color
    node_fill = [
        matplotlib.colors.to_hex(fc_cmap(fc_norm(G.nodes[n]["log2FoldChange"])))
        for n in G.nodes()
    ]

    node_edge = [
        matplotlib.colors.to_hex(bm_cmap(bm_norm(G.nodes[n]["baseMean"])))
        for n in G.nodes()
    ]

    # Node size and shape based on ageing evidence, circle 12 if no evidence and diamond 24 if evidence
    node_size = [12 if G.nodes[n]["evidence"] == "no" else 24 for n in G.nodes()]
    node_symbol = ["circle" if G.nodes[n]["evidence"] == "no" else "diamond" for n in G.nodes()]

    # add the log2FC and baseMean values as labels when overing the nodes with the mouse
    hover_text = [
        f"{n}<br>log2FC={G.nodes[n]['log2FoldChange']:.2f}"
        f"<br>baseMean={G.nodes[n]['baseMean']:.2f}"
        for n in G.nodes()
    ]

    # -----------------------------
    # 4️⃣ Edges
    # -----------------------------
    edge_x, edge_y = [], []
    for u, v in G.edges():
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(color="black", width=0.5),
        hoverinfo="none"
    )

    # -----------------------------
    # 5️⃣ Nodes
    # -----------------------------
    node_trace = go.Scatter(
        x=[pos[n][0] for n in G.nodes()],
        y=[pos[n][1] for n in G.nodes()],
        mode="markers+text",
        text=list(G.nodes()),
        textposition="middle center",
        hoverinfo="text",
        hovertext=hover_text,
        marker=dict(
            size=node_size,
            color=node_fill,
            line=dict(width=2, color=node_edge),
            symbol=node_symbol
        )
    )

    fig = go.Figure([edge_trace, node_trace])
    fig.update_layout(
        width=1000,
        height=1000,
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=20, b=20)
    )

    fig.show()


# Main function
def build_ppi_RNAseq (fin, organism, min_score=0.4, limit=5, ntop=1000, directions=None, API=True) :
    """
    Build ppi network from dataframe with ensembl gene id, gene names, log2FoldChange values, baseMean and adj. p-values.
    Gene names are queried to StringDB using the API or local flat files.
    Nodes are colored by log2FoldChange values, their stoke colored according to baseMean values.
    Edges can be directional, information to provide in directions.
    Write two files: .graphml and .pdf containing the network plot.
    
    Parameters
    ----------
    fin : path
        Experiment result table, must contain 'ensembl_gene_id', 'gene_name', 'log2FoldChange' 'baseMean', 'padj'.
    organism : species name among 'caenorhabditis elegans', 'danio rerio', 'drosophila melanogaster',
        'homo sapiens', 'mus musculus', 'saccharomyces cerevisiae'
    min_score : int, optional
        Minimum score threshold (0-1000), default 0.4.
    limit : int, optional
        Maximum number of neighbors per query gene, default 5.
    ntop : int, optional
        Maximum number of query genes, take the ntop genes sorted by adj. p-value.
    directions : pd.DataFrame, optional
        Dataframe with columns 'from' and 'to', for oriented edges between pairs of specific genes.
    API : bolean, optional
        Whether the network is built from the StringDB API (True) or the local flat files (False). Default True.

    Returns
    -------
    df_network : pd.DataFrame
        All edges between admitted nodes with score >= min_score.
    nodes_local : set[str], optional
        All nodes admitted to the network for local-based network.
    """
    taxon_id=taxons[organism]
    dirname=os.path.dirname(fin)
    name=os.path.basename(fin).replace('results.tsv','')

    # build gene_list and add evidence
    aging_genes = get_evidence(organism)
    print("Prepare query gene list...")
    dfin, gene_list = prepare_query(fin, taxon_id, aging_genes, ntop=ntop)
    print("Number of query genes:", len(gene_list))
    
    if API==True :
        print("Query StringDB API...")
        min_score = int(min_score * 1000)
        df_network, nodes = API_string_network(gene_list, taxon_id, batch_size=2, limit=limit, min_score=min_score)
        
    else :
        print("Query local link file...")
        db_links = read_local_links(f"{StringDB}{taxon_id}.protein.links.detailed.v12.0_annotated.tsv", taxon_id)
        df_network, nodes = local_string_network(db_links, gene_list, min_score=min_score, limit=limit)

    df_network = df_network.loc[:, ~df_network.columns.duplicated()]
    print("Plot network...")
    plot_network_RNAseq(df_network, dfin,
                pdf_file=f"{dirname}{name}network.pdf",
                png_file=f"{dirname}{name}network.png",
                graphml_file=f"{dirname}{name}network.graphml",
                directions=directions)


    return df_network, nodes
    
## Example usage: edges_local, nodes_local = bp.build_ppi_RNAseq(deseq2_file_path, 'caenorhabditis elegans', limit=5, API=False)


# Example on how to query StringDB to get the edges from a list of genes of interest
'''
taxon_id=taxons[organism] ## you need to search for the ensembl taxon_id of your organism, also available at the top of network.py

# query StringDB API
edges_df, nodes = API_string_network(gene_list)
    
# query the local flat files
db_links = read_local_links(f"{StringDB}{taxon_id}.protein.links.detailed.v12.0_annotated.tsv", taxon_id)
edges_df, nodes = local_string_network(db_links, gene_list)
'''

# Example on how to plot the network from edges_df, adjust to you needs (node colors...)
'''
import networkx as nx
import matplotlib.pyplot as plt

# -----------------------------
# 1️⃣ Build graph
# -----------------------------
G = nx.Graph()

for _, r in edges_df.iterrows():
    u = str(r["preferredName_A"])
    v = str(r["preferredName_B"])

    score = r.get("score", 1.0)
    if hasattr(score, "iloc"):
        score = score.iloc[0]

    dashed = r.get("dashed", "no")

    G.add_edge(
        u,
        v,
        score=float(score),
        dashed=str(dashed)
    )

# -----------------------------
# 2️⃣ Layout
# -----------------------------
pos = nx.spring_layout(G, seed=42, iterations=50)

# -----------------------------
# 3️⃣ Draw
# -----------------------------
plt.figure(figsize=(10, 10))

# Edges
for u, v, d in G.edges(data=True):
    x = [pos[u][0], pos[v][0]]
    y = [pos[u][1], pos[v][1]]
    style = "dashed" if d.get("dashed") == "yes" else "solid"
    plt.plot(x, y, color="black", linewidth=0.6, linestyle=style)

# Nodes
nx.draw_networkx_nodes(
    G,
    pos,
    node_size=300,
    node_color="#cccccc",
    edgecolors="black",
    linewidths=1.0,
)

# Labels
nx.draw_networkx_labels(G, pos, font_size=8)

plt.axis("off")
plt.tight_layout()

# -----------------------------
# 4️⃣ Save outputs
# -----------------------------
plt.savefig(pdf_file, dpi=300)
plt.savefig(png_file, dpi=300)
plt.close()

nx.write_graphml(G, graphml_file)

print(
    "Network saved:\n"
    f"  PDF     : {pdf_file}\n"
    f"  PNG     : {png_file}\n"
    f"  GraphML : {graphml_file}"
)
'''

