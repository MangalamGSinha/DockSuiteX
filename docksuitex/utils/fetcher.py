"""Structure fetching from online databases.

This module provides functions for downloading protein and ligand structures
from public databases including RCSB PDB and PubChem. It automates the process
of retrieving structures for docking studies.

The fetching functions handle:
    - Protein structure download from RCSB PDB
    - Ligand structure download from PubChem
    - Batch downloading of multiple structures
    - Automatic file naming and organization

Example:
    Fetching protein structures::

        from docksuitex.utils import fetch_pdb

        # Download single protein
        pdb_path = fetch_pdb("1UBQ", save_to="structures")
        
        # Download multiple proteins
        paths = fetch_pdb(["1UBQ", "2HHB", "3CL0"], save_to="proteins")

    Fetching ligand structures::

        from docksuitex.utils import fetch_sdf

        # Download ligand from PubChem
        sdf_path = fetch_sdf("2244", save_to="ligands")  # Aspirin
        
        # Download multiple ligands
        paths = fetch_sdf([2244, 5090, 6323], save_to="ligands")

Note:
    These functions require an active internet connection and access to
    the respective databases. Rate limiting may apply for batch downloads.
"""

import requests
from pathlib import Path
from typing import Union, List


def fetch_pdb(pdbid: Union[str, List[str]], save_to: Union[str, Path] = ".") -> Union[Path, List[Path]]:
    """Download PDB structure file(s) from the RCSB Protein Data Bank.

    This function downloads the `.pdb` file(s) corresponding to the given
    4-character PDB ID(s).

    Args:
        pdbid (str | List[str]): The 4-character alphanumeric PDB ID (e.g., "1UBQ")
            or a list of PDB IDs.
        save_to (str | Path, optional): Directory to save the file(s).
            Defaults to the current directory.

    Returns:
        Path | List[Path]: The absolute path(s) to the downloaded `.pdb` file(s).

    Raises:
        ValueError: If `pdbid` is not a valid 4-character alphanumeric string.
        requests.RequestException: If the network request fails.
        RuntimeError: If the PDB file cannot be retrieved (e.g., invalid ID).
    """
    if isinstance(pdbid, list):
        return [fetch_pdb(pid, save_to) for pid in pdbid]

    pdbid = pdbid.upper().strip()
    if len(pdbid) != 4 or not pdbid.isalnum():
        raise ValueError(
            "❌ Invalid PDB ID. It must be a 4-character alphanumeric string.")

    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    save_path = Path(save_to).expanduser().resolve() / f"{pdbid}.pdb"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"❌ Failed to download PDB file from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    print(f"✅ Downloaded {pdbid}.pdb → saved to {save_path}")
    return save_path


def fetch_sdf(cid: Union[str, int, List[Union[str, int]]], save_to: Union[str, Path] = ".") -> Union[Path, List[Path]]:
    """Download 3D SDF structure file(s) from PubChem using Compound ID(s).

    This function downloads the `.sdf` file(s) corresponding to the given
    PubChem CID(s).

    Args:
        cid (str | int | List[str | int]): The numeric Compound ID (e.g., 2244 for Aspirin)
            or a list of CIDs.
        save_to (str | Path, optional): Directory to save the file(s).
            Defaults to the current directory.

    Returns:
        Path | List[Path]: The absolute path(s) to the downloaded `.sdf` file(s).

    Raises:
        ValueError: If `cid` is not a valid integer identifier.
        requests.RequestException: If the network request fails.
        RuntimeError: If the SDF file cannot be retrieved or is empty.
    """
    if isinstance(cid, list):
        return [fetch_sdf(c, save_to) for c in cid]

    cid = str(cid).strip()
    if not cid.isdigit():
        raise ValueError("❌ Invalid CID. It must be a numeric ID.")

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    save_path = Path(save_to).expanduser().resolve() / f"{cid}.sdf"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200 or not response.text.strip():
        raise RuntimeError(f"❌ Failed to download SDF file from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    print(f"✅ Downloaded {cid}.sdf → saved to {save_path}")
    return save_path