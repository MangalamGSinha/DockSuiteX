"""Batch grid calculation module for DockSuiteX."""

import os
from pathlib import Path
from typing import List, Union, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback
import pprint

from ..grid_calculator import GridCalculator


class BatchGridCalculator:
    """Batch grid box calculation with parallel processing.

    Runs :class:`GridCalculator` on multiple receptor files in parallel.
    Supports both ``"blind"`` and ``"p2rank"`` modes.
    """

    def __init__(
        self,
        inputs: Union[List[Union[str, Path]], str, Path],
        mode: str = "p2rank",
        max_pockets: Optional[int] = None,
        padding: float = GridCalculator.DEFAULT_PADDING,
    ):
        """Initialise BatchGridCalculator.

        Args:
            inputs: A list of receptor file paths, a directory containing
                ``.pdb`` / ``.pdbqt`` files, or a single file path.
            mode: ``"blind"`` or ``"p2rank"``. Defaults to ``"p2rank"``.
            max_pockets: Maximum pockets to return per receptor in P2Rank
                mode. ``None`` returns all. Ignored in blind mode.
            padding: Bounding-box padding in Å. Defaults to 10.0.
        """
        if mode not in ("blind", "p2rank"):
            raise ValueError(f"Invalid mode '{mode}'. Choose 'blind' or 'p2rank'.")

        self.mode        = mode
        self.max_pockets = max_pockets
        self.padding     = padding

        if isinstance(inputs, (str, Path)):
            path = Path(inputs).resolve()
            if path.is_dir():
                self.files = sorted(
                    f.resolve()
                    for f in path.glob("*")
                    if f.suffix.lower() in (".pdb", ".pdbqt")
                )
            elif path.is_file():
                if path.suffix.lower() in (".pdb", ".pdbqt"):
                    self.files = [path]
                else:
                    raise ValueError(
                        f"❌ Invalid file type: {path.suffix}. Supported: .pdb, .pdbqt"
                    )
            else:
                raise ValueError(f"❌ Input path does not exist: {inputs}")
        elif isinstance(inputs, list):
            self.files = sorted(Path(f).resolve() for f in inputs)
        else:
            raise ValueError(
                "❌ Invalid input. Provide a list of files, a directory, or a single file."
            )

        if not self.files:
            raise ValueError("❌ No valid receptor files found.")

        self.results: Dict[str, List[Dict]] = {}

    # ── Worker ────────────────────────────────────────────────────────────────

    @staticmethod
    def _process_one(
        file_path: Union[str, Path],
        save_to: Path,
        cpu_per_worker: int,
        mode: str,
        padding: float,
    ) -> Dict:
        """Process one receptor file.

        This is an internal worker used for parallel processing.
        """
        try:
            protein_name  = Path(file_path).stem
            protein_out   = save_to / f"{protein_name}_grid"

            calc   = GridCalculator(
                file_path, mode=mode, _cpu=cpu_per_worker, padding=padding
            )
            result = calc.run(save_to=protein_out)

            return {
                "file":       str(file_path),
                "status":     "success",
                "result":     result,
                "output_dir": str(protein_out),
            }
        except Exception as e:
            return {
                "file":      str(file_path),
                "status":    "error",
                "error":     str(e),
                "traceback": traceback.format_exc(),
            }

    # ── Public API ────────────────────────────────────────────────────────────

    def run_all(
        self,
        save_to: Union[str, Path] = "p2rank_outputs",
        cpu: int = (os.cpu_count() or 2) - 1,
    ) -> Dict[str, List[Dict]]:
        """Run grid calculation for all receptors in parallel.

        Args:
            save_to: Root output directory.
            cpu: Total CPU cores to distribute across workers.

        Returns:
            Dict mapping absolute receptor file path → ``list[dict]``.
            Each dict has keys ``rank``, ``probability``, ``center``,
            ``grid_size``.  Failed receptors are omitted.
        """
        save_to = Path(save_to).resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        n_files        = len(self.files)
        max_workers    = min(cpu, n_files)
        cpu_per_worker = max(1, cpu // max_workers)

        print(f"Starting batch grid calculation ({self.mode} mode)")
        print(f"Using {max_workers} parallel workers, {cpu_per_worker} CPUs per worker")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_one,
                    file,
                    save_to,
                    cpu_per_worker,
                    self.mode,
                    self.padding,
                ): file
                for file in self.files
            }

            for future in as_completed(futures):
                res       = future.result()
                file_path = res["file"]

                if res["status"] == "success":
                    result = res["result"]

                    if self.mode == "p2rank" and self.max_pockets is not None:
                        result = result[: self.max_pockets]

                    self.results[file_path] = result

                    print(
                        f"✅ {Path(file_path).name}  →  "
                        f"{len(result)} pocket(s)"
                    )
                    pprint.pprint(result)
                else:
                    print(
                        f"❌ {Path(file_path).name} failed: "
                        f"{res.get('error', 'Unknown error')}"
                    )

        print("✅ Batch grid calculation completed.")
        return self.results
