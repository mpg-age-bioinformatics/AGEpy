import requests
import pandas as pd
from io import StringIO
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional


NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _make_headers(email: Optional[str] = None) -> Dict[str, str]:
    ua = "python-requests (SRA metadata fetcher; contact: {})".format(
        email or "no-email-provided"
    )
    return {"User-Agent": ua}


def fetch_runinfo_for_bioproject(
    bioproject: str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Get NCBI SRA RunInfo table (CSV) for a given BioProject (PRJNA...).

    This is equivalent to: esearch -db sra -query <bioproject> | efetch -format runinfo
    """

    # 1) esearch to get WebEnv + QueryKey
    esearch_params = {
        "db": "sra",
        "term": bioproject,
        "usehistory": "y",
        "retmax": "0",
    }
    if api_key:
        esearch_params["api_key"] = api_key

    resp = requests.get(
        NCBI_EUTILS_BASE + "esearch.fcgi",
        params=esearch_params,
        headers=_make_headers(email),
    )
    resp.raise_for_status()

    # parse minimal XML for WebEnv & QueryKey
    root = ET.fromstring(resp.text)
    webenv = root.findtext("WebEnv")
    query_key = root.findtext("QueryKey")
    if not webenv or not query_key:
        raise RuntimeError(
            f"Could not get WebEnv/QueryKey for BioProject {bioproject}. "
            f"Raw response (first 500 chars):\n{resp.text[:500]}"
        )

    # 2) efetch in runinfo (CSV) mode
    efetch_params = {
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "retmode": "text",
        "rettype": "runinfo",
    }
    if api_key:
        efetch_params["api_key"] = api_key

    runinfo_resp = requests.get(
        NCBI_EUTILS_BASE + "efetch.fcgi",
        params=efetch_params,
        headers=_make_headers(email),
    )
    runinfo_resp.raise_for_status()
    csv_text = runinfo_resp.text.strip()

    if not csv_text or not csv_text.startswith("Run,"):
        raise RuntimeError(
            f"RunInfo CSV for {bioproject} looks wrong.\n"
            f"First 500 chars:\n{csv_text[:500]}"
        )

    runinfo_df = pd.read_csv(StringIO(csv_text))
    return runinfo_df


def fetch_biosample_attributes(
    biosample_ids: List[str],
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    For a list of BioSample accessions (e.g. SAMN...), fetch their Attributes via
    NCBI BioSample EFetch and return a wide DataFrame: one row per BioSample,
    columns = unique attribute names (tissue, cell_line, treatment, etc.).
    """
    biosample_ids = sorted({b for b in biosample_ids if isinstance(b, str) and b})
    records: Dict[str, Dict[str, str]] = {}

    for bs in biosample_ids:
        params = {
            "db": "biosample",
            "id": bs,
            "retmode": "xml",
        }
        if api_key:
            params["api_key"] = api_key

        resp = requests.get(
            NCBI_EUTILS_BASE + "efetch.fcgi",
            params=params,
            headers=_make_headers(email),
        )
        resp.raise_for_status()
        xml_text = resp.text

        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            # skip if something weird came back
            continue

        # BioSample XML: Attributes/Attribute elements
        # We try attribute_name > harmonized_name > display_name as column key.
        attr_dict: Dict[str, str] = {}
        for attr in root.findall(".//Attributes/Attribute"):
            name = (
                attr.get("attribute_name")
                or attr.get("harmonized_name")
                or attr.get("display_name")
            )
            if not name:
                continue
            col = name.strip().lower().replace(" ", "_")
            val = (attr.text or "").strip()
            if not val:
                continue
            # if multiple values with same attribute, join them
            if col in attr_dict:
                if val not in attr_dict[col].split(";"):
                    attr_dict[col] = attr_dict[col] + ";" + val
            else:
                attr_dict[col] = val

        if attr_dict:
            attr_dict["biosample"] = bs
            records[bs] = attr_dict

    if not records:
        return pd.DataFrame(columns=["biosample"])

    df = pd.DataFrame.from_dict(records, orient="index")
    df.reset_index(drop=True, inplace=True)
    return df


def fetch_sra_metadata_table_for_bioproject(
    bioproject: str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    High-level function:
    - Get SRA RunInfo table for BioProject
    - Get BioSample attributes for all involved BioSamples
    - Merge → SraRunTable-like metadata DataFrame.
    """
    runinfo_df = fetch_runinfo_for_bioproject(bioproject, email=email, api_key=api_key)

    # If there's no BioSample column, just return runinfo
    if "BioSample" not in runinfo_df.columns:
        return runinfo_df

    biosamples = [
        b for b in runinfo_df["BioSample"].astype(str).unique() if b and b != "nan"
    ]
    biosample_attr_df = fetch_biosample_attributes(
        biosamples, email=email, api_key=api_key
    )

    if biosample_attr_df.empty:
        return runinfo_df

    # Make merge keys consistent
    biosample_attr_df = biosample_attr_df.rename(columns={"biosample": "BioSample"})
    merged = runinfo_df.merge(biosample_attr_df, on="BioSample", how="left")

    return merged





