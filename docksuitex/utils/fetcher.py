"""Molecule fetching utilities for DockSuiteX."""

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
    """Download PDB structure file(s) from the RCSB Protein Data Bank.

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


def _download_sdf(molecule_id: str, save_to: Union[str, Path]) -> Path:
    """
    Download a single SDF file from PubChem or ChEMBL.

    The source is auto-detected from the ID format:
    - Numeric IDs (e.g. ``2244``) → **PubChem** (3D conformer)
    - IDs starting with ``CHEMBL`` (e.g. ``CHEMBL25``) → **ChEMBL**

    Args:
        molecule_id (str): A PubChem CID (numeric) or ChEMBL ID (e.g. "CHEMBL25").
        save_to (Union[str, Path]): Directory to save the file.

    Returns:
        Path: The path to the downloaded file.

    Raises:
        ValueError: If the ID format is unrecognised.
        RuntimeError: If the download fails.
    """
    molecule_id = str(molecule_id).strip()

    # ── ChEMBL ─────────────────────────────────────────────────
    if molecule_id.upper().startswith("CHEMBL"):
        chembl_id = molecule_id.upper()
        if not chembl_id[6:].isdigit():
            raise ValueError(
                f"❌ Invalid ChEMBL ID '{chembl_id}'. "
                "Must be in the format CHEMBLnnn (e.g. CHEMBL25)."
            )
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.sdf"
        filename = f"{chembl_id}.sdf"
        source = "ChEMBL"

    # ── PubChem ────────────────────────────────────────────────
    elif molecule_id.isdigit():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{molecule_id}/SDF?record_type=3d"
        filename = f"{molecule_id}.sdf"
        source = "PubChem"

    else:
        raise ValueError(
            f"❌ Unrecognised molecule ID '{molecule_id}'. "
            "Provide a numeric PubChem CID or a ChEMBL ID (e.g. CHEMBL25)."
        )

    save_path = Path(save_to).expanduser().resolve() / filename
    save_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url)
    if response.status_code != 200 or not response.text.strip():
        raise RuntimeError(f"❌ Failed to download SDF from {source}: {url}")

    with open(save_path, "w") as f:
        f.write(response.text)

    print(f"✅ Downloaded {filename} ({source}) → {save_path}")
    return save_path


def fetch_sdf(
    molecule_id: Union[str, int, List[Union[str, int]]],
    save_to: Union[str, Path] = ".",
    parallel: int = 4,
) -> Union[Path, List[Path]]:
    """Download SDF structure file(s) from **PubChem** or **ChEMBL**.

    The source is auto-detected from each ID:

    - Numeric IDs (e.g. ``2244``) are fetched as 3D conformers from **PubChem**.
    - IDs starting with ``CHEMBL`` (e.g. ``CHEMBL25``) are fetched from **ChEMBL**.

    You can freely mix both types in a single call.

    Args:
        molecule_id (Union[str, int, List[Union[str, int]]]): A single ID or a
            list of IDs.  Each ID can be a numeric PubChem CID or a ChEMBL ID.
        save_to (Union[str, Path], optional): Directory to save the file(s).
            Defaults to the current directory.
        parallel (int, optional): Number of threads to use for parallel
            downloading. Defaults to 4.

    Returns:
        Union[Path, List[Path]]: The absolute path(s) to the downloaded ``.sdf`` file(s).
    """
    if isinstance(molecule_id, list):
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = executor.map(lambda mid: _download_sdf(mid, save_to), molecule_id)
            return list(futures)
    else:
        return _download_sdf(molecule_id, save_to)
