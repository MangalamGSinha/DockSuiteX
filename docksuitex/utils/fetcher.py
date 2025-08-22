import requests
from pathlib import Path
from typing import Union


def fetch_pdb(pdbid: str, save_dir: Union[str, Path] = ".") -> Path:
    """
    Downloads a .pdb file from the RCSB PDB database using a 4-character PDB ID.

    Args:
        pdbid (str): The 4-character alphanumeric PDB ID (e.g., '1UBQ').
        save_dir (str | Path, optional): Directory to save the file. Defaults to the current directory.

    Returns:
        Path: The full path to the downloaded .pdb file.

    Raises:
        ValueError: If the PDB ID is not a valid 4-character alphanumeric string.
        RuntimeError: If the download fails due to an invalid ID or network issue.
    """
    pdbid = pdbid.upper().strip()
    if len(pdbid) != 4 or not pdbid.isalnum():
        raise ValueError(
            "❌ Invalid PDB ID. It must be a 4-character alphanumeric string.")

    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    save_path = Path(save_dir).resolve() / f"{pdbid}.pdb"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"❌ Failed to download PDB file from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    return save_path


def fetch_sdf(cid: Union[str, int], save_dir: Union[str, Path] = ".") -> Path:
    """
    Downloads a 3D .sdf file from PubChem using the given Compound ID (CID).

    Args:
        cid (str | int): The numeric Compound ID for the molecule (e.g., 2244 for Aspirin).
        save_dir (str | Path, optional): Directory to save the file. Defaults to the current directory.

    Returns:
        Path: The full path to the downloaded .sdf file.

    Raises:
        ValueError: If the CID is not a valid numeric identifier.
        RuntimeError: If the download fails or the file is empty.
    """
    cid = str(cid).strip()
    if not cid.isdigit():
        raise ValueError("❌ Invalid CID. It must be a numeric ID.")

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    save_path = Path(save_dir).resolve() / f"{cid}.sdf"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200 or not response.text.strip():
        raise RuntimeError(f"❌ Failed to download SDF file from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    return save_path
