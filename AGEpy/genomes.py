import os
import sys
import re
import io
import subprocess
from pathlib import Path
from urllib.request import Request, urlopen, urlretrieve
from urllib.parse import urlparse
from ftplib import FTP  # needed by _is_url_file

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
        except Exception:
            # Fallback: GET a tiny amount and inspect
            try:
                req = Request(url, method="GET")
                with urlopen(req) as r:
                    ctype = r.headers.get("Content-Type", "").lower()
                    chunk = r.read(2048).lower()
                    if "text/html" in ctype or b"<html" in chunk or b"index of" in chunk:
                        return False
                    return True
            except Exception:
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
        except Exception:
            try:
                ftp.quit()
            except Exception:
                pass
            # size() fails → not a file → directory
            return False

    # Unknown scheme → assume directory
    return False


def _ensure_trailing_slash(url: str) -> str:
    return url if url.endswith("/") else url + "/"


def _log_and_check_checksum(
    base_url: str,
    filename: str,
    local_path: Path,
    dest_dir: Path,
    checksums_filename: str = "checksums.txt",
    raise_on_mismatch: bool = True,
):
    """
    Download CHECKSUMS from base_url (directory containing the file),
    find the entry for 'filename', compute local Unix 'sum' and compare.
    Log everything into dest_dir/checksums_filename.

    If the 'sum' command is not available or CHECKSUMS is missing,
    we just log what we can and do not raise (unless we have a clear mismatch).
    """

    # 1) Fetch remote CHECKSUMS
    checksums_url = f"{base_url}CHECKSUMS"
    remote_line = None
    remote_sum = None
    remote_blocks = None

    try:
        with urlopen(checksums_url) as r:
            lines = r.read().decode().splitlines()
        # Find line that ends with the filename (works for both full path or just name)
        for ln in lines:
            if ln.strip().endswith(filename):
                remote_line = ln.strip()
                break
        if remote_line:
            parts = remote_line.split()
            # Typical Ensembl format: "<sum> <blocks> <filename>"
            if len(parts) >= 2:
                remote_sum = parts[0]
                remote_blocks = parts[1]
    except Exception:
        # No CHECKSUMS or network issue – just skip remote info
        pass

    # 2) Compute local sum using Unix 'sum' if available
    local_sum = None
    local_blocks = None
    try:
        proc = subprocess.run(
            ["sum", str(local_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        out_line = proc.stdout.strip().splitlines()[0]
        parts = out_line.split()
        if len(parts) >= 1:
            local_sum = parts[0]
        if len(parts) >= 2:
            local_blocks = parts[1]
    except Exception:
        # 'sum' not available or failed; we'll just log remote info
        pass

    # 3) Decide status
    status = "NOT_CHECKED"
    if remote_sum is not None and local_sum is not None:
        if remote_sum == local_sum and (
            remote_blocks is None or local_blocks == remote_blocks
        ):
            status = "OK"
        else:
            status = "MISMATCH"

    # 4) Log to checksums.txt
    try:
        log_path = dest_dir / checksums_filename
        with log_path.open("a", encoding="utf-8") as log:
            log.write(
                "\t".join(
                    [
                        str(local_path),
                        filename,
                        remote_sum or "",
                        remote_blocks or "",
                        local_sum or "",
                        local_blocks or "",
                        status,
                    ]
                )
                + "\n"
            )
    except Exception:
        pass

    # 5) Optionally raise on mismatch
    if raise_on_mismatch and status == "MISMATCH":
        raise ValueError(
            f"Checksum mismatch for {local_path}: "
            f"remote={remote_sum} blocks={remote_blocks}, "
            f"local={local_sum} blocks={local_blocks}"
        )


def _download_readme(base_url: str, dest_dir: Path, organism: str, release: str, kind: str):
    """
    Download README from the directory and store as:
      <organism>.<release>.<kind>.readme.txt
    where kind is 'gtf' or 'dna'.
    """
    readme_url = f"{base_url}README"
    readme_name = f"{organism}.{release}.{kind}.readme.txt"
    readme_path = dest_dir / readme_name
    try:
        urlretrieve(readme_url, readme_path.as_posix())
    except Exception:
        # README may not exist; ignore silently
        pass


# --- main functions --------------------------------------------------------


def ensembl_releases(ensembl_url="ftp://ftp.ensembl.org/pub/"):
    with urlopen(ensembl_url) as response:
        lines = response.read().decode().splitlines()
    names = [line.split()[-1] for line in lines]
    releases = {x.split("/")[-1] for x in names if re.fullmatch(r"release-\d+", x.split("/")[-1])}
    return sorted(releases, key=lambda s: int(s.split("-")[1]))


def ensembl_release_organisms(release="last", ensembl_url="ftp://ftp.ensembl.org/pub/"):

    if release == "last":
        release = ensembl_releases(ensembl_url=ensembl_url)[-1].split("-")[-1]

    ensembl_url_ = f"{ensembl_url}release-{release}/gtf/"

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
    verify_checksums=True,
):
    """
    Download a GTF file. 'source' can be:
      - A direct file URL (any extension) -> downloaded as-is (no CHECKSUMS/README).
      - A directory URL (ending with '/') pointing to an Ensembl GTF dir.
      - An Ensembl base (ftp/http) from which we construct the GTF directory.

    Returns: (local_path, download_url)
    """

    dest = Path(dest_dir).expanduser()
    dest.mkdir(parents=True, exist_ok=True)

    # Case 1: Direct file URL -> download immediately (user gave full path)
    if _is_url_file(source):
        download_url = source
        if filename is None:
            filename = source.rstrip("/").split("/")[-1]
        local_path = dest / filename
        urlretrieve(download_url, local_path.as_posix())
        # Intentionally skip CHECKSUMS/README here, since user gave a full path
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

    local_path = dest / filename
    urlretrieve(download_url, local_path.as_posix())

    # CHECKSUMS + README (only when we derived the path ourselves)
    if verify_checksums:
        _log_and_check_checksum(base, gtf_file, local_path, dest)
    _download_readme(base, dest, organism, str(release), "gtf")

    return local_path.as_posix(), download_url


def ensembl_dna(
    organism="homo_sapiens",
    release="last",
    seq_type="toplevel",  # or "primary_assembly"
    source="ftp://ftp.ensembl.org/pub/",
    dest_dir=".",
    filename=None,
    verify_checksums=True,
):
    """
    Download a DNA FASTA file. 'source' can be:
      - A direct file URL (any extension) -> downloaded as-is (no CHECKSUMS/README).
      - A directory URL (ending with '/') pointing to an Ensembl dna/ dir.
      - An Ensembl base (ftp/http) from which we construct the dna/ directory.

    Returns: (local_path, download_url)
    """

    dest = Path(dest_dir).expanduser()
    dest.mkdir(parents=True, exist_ok=True)

    # Case 1: Direct file URL -> download immediately (user gave full path)
    if _is_url_file(source):
        download_url = source
        if filename is None:
            filename = source.rstrip("/").split("/")[-1]
        local_path = dest / filename
        urlretrieve(download_url, local_path.as_posix())
        # Intentionally skip CHECKSUMS/README here, since user gave a full path
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

    local_path = dest / filename
    urlretrieve(download_url, local_path.as_posix())

    # CHECKSUMS + README (only when we derived the path ourselves)
    if verify_checksums:
        _log_and_check_checksum(base, dna_file, local_path, dest)
    _download_readme(base, dest, organism, str(release), "dna")

    return local_path.as_posix(), download_url