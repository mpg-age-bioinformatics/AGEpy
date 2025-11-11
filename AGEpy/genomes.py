import os
import sys
import re
import io
from pathlib import Path
from urllib.request import Request, urlopen, urlretrieve
from urllib.parse import urlparse

import numpy as np
import pandas as pd


def _is_url_file(url: str) -> bool:
    """
    Determine if a URL points to a file (HTTP/HTTPS/FTP) without relying on file extensions.
    """

    parsed = urlparse(url)
    scheme = parsed.scheme.lower()

    # ---------------- HTTP / HTTPS ----------------
    if scheme in ("http", "https"):
        # Try HEAD first
        try:
            req = Request(url, method="HEAD")
            with urlopen(req) as r:
                ctype = r.headers.get("Content-Type", "").lower()
                # HTML directory indexes show as HTML → directory
                if "text/html" in ctype:
                    return False
                return True
        except:
            # Fallback: GET a tiny amount and inspect
            try:
                req = Request(url, method="GET")
                with urlopen(req) as r:
                    ctype = r.headers.get("Content-Type", "").lower()
                    chunk = r.read(2048).lower()
                    if "text/html" in ctype or b"<html" in chunk or b"index of" in chunk:
                        return False
                    return True
            except:
                return False

    # ---------------- FTP ----------------
    if scheme == "ftp":
        host = parsed.hostname
        ftp_path = parsed.path

        try:
            ftp = FTP(host, timeout=10)
            ftp.login()  # anonymous
            # If size() works → FILE
            size = ftp.size(ftp_path)
            ftp.quit()
            return size is not None
        except:
            # size() fails → not a file → directory
            try:
                ftp.quit()
            except:
                pass
            return False

    # Unknown scheme → assume directory
    return False


def _ensure_trailing_slash(url: str) -> str:
    return url if url.endswith("/") else url + "/"


# --- main functions --------------------------------------------------------


def ensembl_releases(ensembl_url="ftp://ftp.ensembl.org/pub/"):
    with urlopen(ensembl_url) as response:
        lines = response.read().decode().splitlines()
    names = [line.split()[-1] for line in lines]
    releases = {x.split("/")[-1] for x in names if re.fullmatch(r"release-\d+", x.split("/")[-1])}
    return sorted(releases, key=lambda s: int(s.split("-")[1]))

def ensembl_release_organisms(release="last", ensembl_url = "ftp://ftp.ensembl.org/pub/"):

    if release == 'last' :
        release=ensembl_releases(ensembl_url = ensembl_url)[-1].split("-")[-1]
    
    ensembl_url_=f'{ensembl_url}release-{release}/gtf/'
    
    with urlopen(ensembl_url_) as response:
        lines = response.read().decode().splitlines()

    # Extract last path component
    names = [line.split()[-1] for line in lines]

    return names

def ensembl_gtf(
    organism="homo_sapiens",
    release="last",
    source="ftp://ftp.ensembl.org/pub/",
    dest_dir=".",
    filename=None,
):
    """
    Download a GTF file. 'source' can be:
      - A direct file URL (any extension) -> downloaded as-is.
      - A directory URL (ending with '/') pointing to an Ensembl GTF dir.
      - An Ensembl base (ftp/http) from which we construct the GTF directory.
    Returns: (local_path, download_url)
    """

    # Case 1: Direct file URL -> download immediately
    if _is_url_file(source):
        download_url = source
        if filename is None:
            filename = source.rstrip("/").split("/")[-1]
        dest = Path(dest_dir).expanduser()
        dest.mkdir(parents=True, exist_ok=True)
        local_path = dest / filename
        urlretrieve(download_url, local_path.as_posix())
        return local_path.as_posix(), download_url

    # Case 2/3: Directory or base -> figure out the correct directory
    base = _ensure_trailing_slash(source)
    # Heuristic: if it already looks like .../gtf/<organism>/ then treat as direct dir
    is_gtf_dir = base.rstrip("/").endswith(f"/gtf/{organism}")

    if release == "last" and not is_gtf_dir:
        release = ensembl_releases(ensembl_url=base)[-1].split("-")[-1]

    if not is_gtf_dir:
        base = f"{base}release-{release}/gtf/{organism}/"

    # List directory and choose a GTF
    with urlopen(base) as response:
        names = [line.split()[-1] for line in response.read().decode().splitlines()]

    # Try strict match when we know the release; otherwise accept any .gtf
    candidates = [n for n in names if n.endswith(f".{release}.gtf.gz")] if release else []
    if not candidates:
        candidates = [n for n in names if n.endswith(".gtf.gz") or n.endswith(".gtf")]

    # Prefer canonical (exclude abinitio and .chr)
    canonical = [n for n in candidates if "abinitio" not in n and ".chr." not in n]
    if canonical:
        gtf_file = sorted(canonical)[0]
    elif candidates:
        gtf_file = sorted(candidates)[0]
    else:
        raise FileNotFoundError(f"No GTF found in {base}")

    download_url = f"{base}{gtf_file}"

    if filename is None:
        filename = gtf_file
    dest = Path(dest_dir).expanduser()
    dest.mkdir(parents=True, exist_ok=True)
    local_path = dest / filename
    urlretrieve(download_url, local_path.as_posix())
    return local_path.as_posix(), download_url


def ensembl_dna(
    organism="homo_sapiens",
    release="last",
    seq_type="toplevel",  # or "primary_assembly"
    source="ftp://ftp.ensembl.org/pub/",
    dest_dir=".",
    filename=None,
):
    """
    Download a DNA FASTA file. 'source' can be:
      - A direct file URL (any extension) -> downloaded as-is.
      - A directory URL (ending with '/') pointing to an Ensembl dna/ dir.
      - An Ensembl base (ftp/http) from which we construct the dna/ directory.
    Returns: (local_path, download_url)
    """

    # Case 1: Direct file URL -> download immediately
    if _is_url_file(source):
        download_url = source
        if filename is None:
            filename = source.rstrip("/").split("/")[-1]
        dest = Path(dest_dir).expanduser()
        dest.mkdir(parents=True, exist_ok=True)
        local_path = dest / filename
        urlretrieve(download_url, local_path.as_posix())
        return local_path.as_posix(), download_url

    # Case 2/3: Directory or base
    base = _ensure_trailing_slash(source)
    # Heuristic: already a dna directory?
    is_dna_dir = base.rstrip("/").endswith("/dna")

    if release == "last" and not is_dna_dir:
        release = ensembl_releases(ensembl_url=base)[-1].split("-")[-1]

    if not is_dna_dir:
        base = f"{base}release-{release}/fasta/{organism}/dna/"

    # List directory and choose a FASTA according to seq_type
    with urlopen(base) as response:
        names = [line.split()[-1] for line in response.read().decode().splitlines()]

    # Prefer exact dna.<type> files; fallback to any .fa/.fa.gz
    pattern = f".dna.{seq_type}."
    matches = [n for n in names if pattern in n and (n.endswith(".fa.gz") or n.endswith(".fa"))]
    if not matches:
        matches = [n for n in names if n.endswith(".fa.gz") or n.endswith(".fa")]

    if not matches:
        raise FileNotFoundError(f"No FASTA found in {base} (seq_type='{seq_type}')")

    dna_file = sorted(matches)[0]
    download_url = f"{base}{dna_file}"

    if filename is None:
        filename = dna_file
    dest = Path(dest_dir).expanduser()
    dest.mkdir(parents=True, exist_ok=True)
    local_path = dest / filename
    urlretrieve(download_url, local_path.as_posix())
    return local_path.as_posix(), download_url