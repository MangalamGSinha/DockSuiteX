import requests
from pathlib import Path
from typing import Union


def fetch_pdb(pdbid: str, save_to: Union[str, Path] = ".") -> Path:
    """Download a PDB structure file from the RCSB Protein Data Bank.

    This function downloads the `.pdb` file corresponding to the given
    4-character PDB ID.

    Args:
        pdbid (str): The 4-character alphanumeric PDB ID (e.g., "1UBQ").
        save_to (str | Path, optional): Directory to save the file.
            Defaults to the current directory.

    Returns:
        Path: The absolute path to the downloaded `.pdb` file.

    Raises:
        ValueError: If `pdbid` is not a valid 4-character alphanumeric string.
        requests.RequestException: If the network request fails.
        RuntimeError: If the PDB file cannot be retrieved (e.g., invalid ID).
    """
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

    return save_path


def fetch_sdf(cid: Union[str, int], save_to: Union[str, Path] = ".") -> Path:
    """Download a 3D SDF structure file from PubChem using a Compound ID (CID).

    This function downloads the `.sdf` file corresponding to the given
    PubChem CID

    Args:
        cid (str | int): The numeric Compound ID (e.g., 2244 for Aspirin).
        save_to (str | Path, optional): Directory to save the file.
            Defaults to the current directory.

    Returns:
        Path: The absolute path to the downloaded `.sdf` file.

    Raises:
        ValueError: If `cid` is not a valid integer identifier.
        requests.RequestException: If the network request fails.
        RuntimeError: If the SDF file cannot be retrieved or is empty.
    """
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

    return save_path