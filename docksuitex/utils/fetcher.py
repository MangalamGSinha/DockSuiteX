import requests
from pathlib import Path
from typing import Union, List
import concurrent.futures


def _download_pdb(pid: str, save_to: Union[str, Path]) -> Path:
    """
    Download a single PDB file from RCSB.

    Args:
        pid (str): The 4-character PDB ID.
        save_to (Union[str, Path]): Directory to save the file.

    Returns:
        Path: The path to the downloaded file.

    Raises:
        ValueError: If the PDB ID is invalid.
        RuntimeError: If the download fails.
    """
    pid = pid.upper().strip()
    if len(pid) != 4 or not pid.isalnum():
        raise ValueError("❌ Invalid PDB ID. Must be 4-character alphanumeric.")

    url = f"https://files.rcsb.org/download/{pid}.pdb"
    save_path = Path(save_to).expanduser().resolve() / f"{pid}.pdb"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"❌ Failed to download PDB file from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    print(f"✅ Downloaded {pid}.pdb → {save_path}")
    return save_path


def fetch_pdb(pdbid: Union[str, List[str]], save_to: Union[str, Path] = ".", parallel: int = 4) -> Union[Path, List[Path]]:
    """
    Download PDB structure file(s) from the RCSB Protein Data Bank.

    This function downloads `.pdb` files for the given PDB ID(s). It supports
    parallel downloading when a list of IDs is provided.

    Args:
        pdbid (Union[str, List[str]]): A single 4-character PDB ID (e.g., "1UBQ")
            or a list of PDB IDs.
        save_to (Union[str, Path], optional): Directory to save the file(s).
            Defaults to the current directory.
        parallel (int, optional): Number of threads to use for parallel
            downloading. Defaults to 4.

    Returns:
        Union[Path, List[Path]]: The absolute path(s) to the downloaded `.pdb` file(s).
    """
    if isinstance(pdbid, list):
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = executor.map(lambda pid: _download_pdb(pid, save_to), pdbid)
            return list(futures)
    else:
        return _download_pdb(pdbid, save_to)


def _download_sdf(cid: str, save_to: Union[str, Path]) -> Path:
    """
    Download a single SDF file from PubChem.

    Args:
        cid (str): The Compound ID (CID).
        save_to (Union[str, Path]): Directory to save the file.

    Returns:
        Path: The path to the downloaded file.

    Raises:
        ValueError: If the CID is invalid.
        RuntimeError: If the download fails.
    """
    cid = str(cid).strip()
    if not cid.isdigit():
        raise ValueError("❌ Invalid CID. Must be numeric.")

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    save_path = Path(save_to).expanduser().resolve() / f"{cid}.sdf"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200 or not response.text.strip():
        raise RuntimeError(f"❌ Failed to download SDF from: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    print(f"✅ Downloaded {cid}.sdf → {save_path}")
    return save_path


def fetch_sdf(cid: Union[str, int, List[Union[str, int]]], save_to: Union[str, Path] = ".", parallel: int = 4) -> Union[Path, List[Path]]:
    """
    Download 3D SDF structure file(s) from PubChem.

    This function downloads `.sdf` files for the given Compound ID(s). It supports
    parallel downloading when a list of CIDs is provided.

    Args:
        cid (Union[str, int, List[Union[str, int]]]): A single numeric Compound ID
            (e.g., 2244) or a list of CIDs.
        save_to (Union[str, Path], optional): Directory to save the file(s).
            Defaults to the current directory.
        parallel (int, optional): Number of threads to use for parallel
            downloading. Defaults to 4.

    Returns:
        Union[Path, List[Path]]: The absolute path(s) to the downloaded `.sdf` file(s).
    """
    if isinstance(cid, list):
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = executor.map(lambda c: _download_sdf(c, save_to), cid)
            return list(futures)
    else:
        return _download_sdf(cid, save_to)
