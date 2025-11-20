import urllib.request
import gzip
import io
import pandas as pd
from io import StringIO
import requests
import re

def geo_to_sra_samples(gse_id: str, email=None, api_key=None):
    """
    Given a GEO series ID (e.g. 'GSE129642'), download the series matrix,
    extract sample metadata, infer SRA Experiments, detect BioProject,
    fetch SRA RunInfo automatically, and return a merged GEO+SRA table.

    Parameters
    ----------
    gse_id : str
        GEO Series ID (e.g. "GSE129642")
    email : str, optional
        optional email passed to NCBI E-utilities
    api_key : str, optional
        NCBI API key for higher request limits

    Returns
    -------
    merged_samples_df : DataFrame
        Final GEO+SRA combined sample table
    groups_df : DataFrame
        GEO sample-label groups
    runinfo_df : DataFrame
        The SRA RunInfo table retrieved for the BioProject
    """

    # -------------------------------------------------------------
    # 1. Download GEO Series Matrix
    # -------------------------------------------------------------
    prefix_folder = gse_id[:6] + "nnn"
    url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{prefix_folder}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    )

    print(f"Downloading GEO matrix: {url}")
    try:
        with urllib.request.urlopen(url) as response:
            compressed_data = response.read()
    except Exception as e:
        raise RuntimeError(f"Failed to download {url}: {e}")

    with gzip.GzipFile(fileobj=io.BytesIO(compressed_data)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    lines = text.splitlines()

    # -------------------------------------------------------------
    # 2. Extract Series BioProject + SRA Study
    # -------------------------------------------------------------
    series_bioproject = None
    series_sra = None

    for line in lines:
        if not line.startswith("!Series_relation"):
            continue

        val = line.split("\t", 1)[-1].strip().strip('"')

        if val.startswith("BioProject:") and "bioproject/" in val:
            series_bioproject = val.split("bioproject/")[-1]

        elif val.startswith("SRA:") and "sra?term=" in val:
            series_sra = val.split("sra?term=")[-1]

    if series_bioproject is None:
        raise RuntimeError("Unable to extract BioProject from GEO Series.")

    print(f"Detected BioProject: {series_bioproject}")

    # -------------------------------------------------------------
    # 3. Parse GEO sample metadata
    # -------------------------------------------------------------
    start_idx = next(i for i, l in enumerate(lines) if l.startswith("!Sample_title"))
    try:
        end_idx = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    except StopIteration:
        end_idx = len(lines)

    sample_lines = lines[start_idx:end_idx]
    meta_df = pd.read_csv(StringIO("\n".join(sample_lines)), sep="\t", header=None, dtype=str)
    meta_df[0] = meta_df[0].str.lstrip("!")

    if (meta_df[0] == "Sample_geo_accession").any():
        samples = meta_df.loc[meta_df[0] == "Sample_geo_accession"].iloc[0, 1:]
        meta_df.columns = ["field"] + samples.tolist()
    else:
        meta_df.columns = ["field"] + [f"s{i}" for i in range(1, meta_df.shape[1])]

    # -------------------------------------------------------------
    # 4. Extract sample-level SRA Experiments
    # -------------------------------------------------------------
    sra = meta_df[meta_df["field"] == "Sample_relation"]
    samples_dict = {}

    for sample in sra.columns[1:]:
        entries = sra[sample].dropna().tolist()
        found = [
            v.split("sra?term=")[-1]
            for v in entries
            if isinstance(v, str) and v.startswith("SRA:")
        ]
        if found:
            samples_dict[sample] = found[0]

    # -------------------------------------------------------------
    # 5. Characteristics → labels
    # -------------------------------------------------------------
    ch1_rows = [f for f in meta_df["field"] if isinstance(f, str) and f.endswith("_ch1")]
    ch1_df = meta_df[meta_df["field"].isin(ch1_rows)].transpose()[1:]

    # remove invariant cols
    for col in list(ch1_df.columns):
        if len(set(ch1_df[col])) == 1:
            ch1_df = ch1_df.drop(columns=[col])

    ch1_df["label"] = ch1_df.apply(lambda r: "; ".join(r.astype(str)), axis=1)
    ch1_df = ch1_df[["label"]]

    # groups
    unique = list(dict.fromkeys(ch1_df["label"]))
    group_map = {lab: f"group_{i+1}" for i, lab in enumerate(unique)}
    ch1_df["group"] = ch1_df["label"].map(group_map)
    ch1_df["Experiment"] = ch1_df.index.map(samples_dict)
    ch1_df["sample"] = ch1_df.index

    # titles
    titles = meta_df[meta_df["field"] == "Sample_title"].drop(columns=["field"]).transpose()
    titles.columns = ["sample_title"]
    ch1_df = ch1_df.join(titles)

    # replicate file-name
    def add_reps(groups):
        ct = {}
        out = []
        for g in groups:
            ct[g] = ct.get(g, 0) + 1
            out.append(f"{g}.rep_{ct[g]}")
        return out

    ch1_df["file_name"] = add_reps(ch1_df["group"])

    samples_df = ch1_df[["sample", "Experiment", "sample_title", "group", "file_name", "label"]]
    groups_df = samples_df[["group", "label"]].drop_duplicates()

    samples_df["series_bioproject"] = series_bioproject
    samples_df["series_sra"] = series_sra
    groups_df["series_bioproject"] = series_bioproject
    groups_df["series_sra"] = series_sra

    # -------------------------------------------------------------
    # 6. Fetch SRA RunInfo automatically for the detected BioProject
    # -------------------------------------------------------------
    print("Fetching SRA RunInfo...")

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    esearch = {
        "db": "sra",
        "term": f"{series_bioproject}[bioproject]",
        "usehistory": "y",
        "retmax": 0,
        "retmode": "json",
    }
    if email:
        esearch["email"] = email
    if api_key:
        esearch["api_key"] = api_key

    r = requests.get(base + "esearch.fcgi", params=esearch)
    r.raise_for_status()

    js = r.json()["esearchresult"]
    if int(js["count"]) == 0:
        raise ValueError(f"No SRA results found for {series_bioproject}")

    runinfo_query = {
        "db": "sra",
        "query_key": js["querykey"],
        "WebEnv": js["webenv"],
        "rettype": "runinfo",
        "retmode": "text",
    }
    if api_key:
        runinfo_query["api_key"] = api_key

    r = requests.get(base + "efetch.fcgi", params=runinfo_query)
    r.raise_for_status()

    runinfo_text = r.text.strip()
    runinfo_df = pd.read_csv(StringIO(runinfo_text))

    # -------------------------------------------------------------
    # 7. Merge GEO sample table with SRA RunInfo
    # -------------------------------------------------------------
    runinfo_sub = runinfo_df[["Experiment", "Run", "BioSample", "LibraryLayout"]]

    merged = samples_df.merge(runinfo_sub, on="Experiment", how="outer")

    merged = merged[
        [
            "series_bioproject",
            "series_sra",
            "sample",
            "Experiment",
            "Run",
            "BioSample",
            "LibraryLayout",
            "sample_title",
            "group",
            "label",
        ]
    ]

    return merged, groups_df, runinfo_df


def gsm_to_gse(gsm_id: str):
    """
    Given a GEO Sample ID (e.g. 'GSM3717993'), return the GEO Series ID(s)
    (e.g. ['GSE123456', 'GSE789012']) that this sample belongs to.

    Strategy:
      1. Ask GEO for the 'series' target of this GSM in text/brief form.
      2. Extract all GSE accessions via regex (GSE\\d+).
      3. If that fails, fall back to generic HTML/text and again grep for GSE\\d+.
    """

    GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

    gse_pattern = re.compile(r"\bGSE\d+\b", re.IGNORECASE)
    series_ids: set[str] = set()

    # --- 1) Preferred: series view in text form ---
    params_series = {
        "acc": gsm_id,
        "targ": "series",   # show series that include this GSM
        "form": "text",
        "view": "brief",
    }
    r = requests.get(GEO_BASE, params=params_series)
    r.raise_for_status()
    text = r.text

    for match in gse_pattern.findall(text):
        series_ids.add(match.upper())

    if series_ids:
        return sorted(series_ids)

    # --- 2) Fallback: default (HTML/text) view of GSM page and scan for GSEs ---
    params_default = {
        "acc": gsm_id,
    }
    r = requests.get(GEO_BASE, params=params_default)
    r.raise_for_status()
    text = r.text

    for match in gse_pattern.findall(text):
        series_ids.add(match.upper())

    return sorted(series_ids)



def gsms_to_sra_samples(gsms):
    """
    Retrieve SRA sample information for a list of GSM accessions.

    This function maps GSM IDs to their corresponding GSE accession,
    fetches SRA sample metadata for the GSE, caches previously queried samples,
    and returns a combined DataFrame of SRA samples corresponding to the input GSMs.

    Parameters
    ----------
    gsms : list or iterable
        A collection of GSM accession identifiers (e.g., ["GSM12345", "GSM67890"]).

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing SRA sample records corresponding to the provided GSM accessions.
        Columns and structure are determined by `geo_to_sra_samples()`.
    """
    queries = pd.DataFrame(columns=["sample"])
    all_sra_samples = pd.DataFrame()

    for gsm in gsms:
        # Check if sample already cached
        if gsm in queries["sample"].tolist():
            sra_subset = queries[queries["sample"] == gsm]
        else:
            # Retrieve GSE accession and query SRA samples
            gse = gsm_to_gse(gsm)[0]
            sra_results , groups_df, runinfo_df= geo_to_sra_samples(gse)

            # Cache all results from this GSE
            queries = pd.concat([queries, sra_results], ignore_index=True)

            # Extract only the row(s) for the current GSM
            sra_subset = sra_results[sra_results["sample"] == gsm]

        # Accumulate all matching results
        all_sra_samples = pd.concat([all_sra_samples, sra_subset], ignore_index=True)

    return all_sra_samples