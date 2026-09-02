"""Format conversion utilities for DockSuiteX."""

import subprocess
import shutil
from pathlib import Path
from typing import Optional, Union
import os

# Locate obabel from PATH (installed via: pip install openbabel-wheel)
OBABEL_EXE = shutil.which("obabel")

# Commonly used molecular file formats and their Open Babel identifiers
COMMON_FORMATS = {
    # 3D structure formats
    "mol2", "pdb", "sdf", "mol", "xyz", "cif", "pdbqt",
    # Text / line-notation formats
    "smi", "smiles", "can", "inchi", "inchikey",
    # Sequence formats
    "fasta",
    # Simulation / other
    "gro", "ent", "mcif", "mmcif",
    # Depiction
    "svg",
    # Proprietary / legacy
    "cdx", "mol2", "rxn",
}


def convert(
    input: Union[str, Path, list[Union[str, Path]]],
    output_format: str,
    gen3d: bool = False,
    split: bool = False,
    output_dir: Optional[Union[str, Path]] = None,
    extra_args: Optional[list[str]] = None,
) -> Union[Path, list[Path]]:
    """Convert molecular file(s) from one format to another using Open Babel.

    Supports single files, lists of file paths, or a directory path.
    The output filename is automatically derived by replacing the input
    extension with the specified output_format.

    Args:
        input (str | Path | list): Path to the input file, a list of paths,
            or a directory containing molecular files.
        output_format (str): The target format identifier (e.g. "sdf", "pdbqt").
        gen3d (bool): Generate 3D coordinates. Essential when converting from
            SMILES or other 2D/text formats. Defaults to False.
        split (bool): Split multi-molecule/model files into individual files.
            Appends suffixes (1, 2, 3...) to the output filename. Defaults to False.
        output_dir (str | Path, optional): Directory to save the converted files.
            If not provided, files are saved in the same directory as the input.
        extra_args (list[str], optional): Additional command-line arguments
            passed directly to Open Babel (e.g. ["-h", "--minimize"]).

    Returns:
        Path | list[Path]: Path to the created output file (single input)
            or a list of paths (batch input).

    Raises:
        FileNotFoundError: If input files or Open Babel executable are missing.
        RuntimeError: If the Open Babel conversion fails.
    """
    if not OBABEL_EXE:
        raise FileNotFoundError(
            "Error: Open Babel executable not found on PATH.\n"
            "   Install it with: pip install openbabel-wheel"
        )

    # 1. Identify all input files
    input_paths: list[Path] = []
    is_batch = False

    if isinstance(input, list):
        input_paths = [Path(p).resolve() for p in input]
        is_batch = True
    else:
        p = Path(input).resolve()
        if p.is_dir():
            # Glob for common molecular files
            input_paths = [
                f for f in p.glob("*")
                if f.is_file() and f.suffix.lower().lstrip(".") in COMMON_FORMATS
            ]
            input_paths.sort()
            is_batch = True
        elif p.is_file():
            input_paths = [p]
            is_batch = False
        else:
            raise FileNotFoundError(f"❌ Input path not found: {input}")

    if not input_paths:
        raise ValueError(f"❌ No valid molecular files found in: {input}")

    results: list[Path] = []

    # 2. Process each file
    output_format = output_format.lower().lstrip(".")
    
    if output_dir:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    for inp in input_paths:
        # Determine specific output path for this file
        if output_dir:
            target_out = output_dir / inp.with_suffix(f".{output_format}").name
        else:
            target_out = inp.with_suffix(f".{output_format}")

        # Infer input format from extension
        current_in_format = inp.suffix.lower().lstrip(".")
        
        if not current_in_format:
            continue

        # Build command
        cmd = [
            OBABEL_EXE,
            "-i", current_in_format, str(inp),
            "-o", output_format, "-O", str(target_out),
        ]

        if split:
            cmd.append("-m")

        if gen3d:
            cmd.append("--gen3d")

        if extra_args:
            cmd += extra_args

        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            continue

        # Handle splitting vs normal output
        if split:
            # Open Babel with -m creates files like base1.ext, base2.ext, etc.
            i = 1
            found_any = False
            while True:
                indexed_out = target_out.parent / f"{target_out.stem}{i}{target_out.suffix}"
                if indexed_out.exists():
                    results.append(indexed_out)
                    found_any = True
                    i += 1
                else:
                    break
            
            if not found_any:
                continue

        else:
            if not target_out.exists():
                continue

            results.append(target_out)

    if not results:
        err_msg = f"❌ Conversion failed for {input_paths[0].name}"
        if 'result' in locals() and result.stderr:
            err_msg += f"\nOpen Babel Error: {result.stderr.strip()}"
        raise RuntimeError(err_msg)

    if not is_batch and not split:
        return results[0]

    return results




def get_supported_formats() -> dict[str, str]:
    """Query Open Babel for all supported file formats.

    Returns:
        dict[str, str]: Mapping of format identifier to description.

    Raises:
        FileNotFoundError: If the Open Babel executable is not found.
        RuntimeError: If Open Babel fails to list formats.
    """
    if not OBABEL_EXE:
        raise FileNotFoundError(
            "Error: Open Babel executable not found on PATH.\n"
            "   Install it with: pip install openbabel-wheel"
        )

    result = subprocess.run(
        [OBABEL_EXE, "-L", "formats"],
        capture_output=True, text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(f"Error: Failed to list formats:\n{result.stderr}")

    formats = {}
    for line in result.stdout.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        # Format lines look like: "sdf -- MDL MOL format"
        parts = line.split(" -- ", 1)
        if len(parts) == 2:
            fmt_id = parts[0].strip()
            description = parts[1].strip()
            formats[fmt_id] = description

    return formats
