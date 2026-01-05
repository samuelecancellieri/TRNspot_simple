"""
I/O utility functions for TRNspot
=================================

Functions for directory management, file tracking, and checkpointing.
"""

import os
import shutil
import hashlib
import json
from datetime import datetime
from typing import Optional, Dict, Any


def setup_directories(
    output_dir: str,
    figures_dir: Optional[str] = None,
    clean: bool = False,
) -> Dict[str, str]:
    """
    Create necessary output directories for the pipeline.

    Parameters
    ----------
    output_dir : str
        Main output directory
    figures_dir : str, optional
        Figures directory. If None, created as {output_dir}/figures
    clean : bool
        If True, remove existing directories before creating

    Returns
    -------
    dict
        Dictionary with paths to created directories
    """
    if figures_dir is None:
        figures_dir = os.path.join(output_dir, "figures")

    directories = {
        "output": output_dir,
        "figures": figures_dir,
        "logs": os.path.join(output_dir, "logs"),
        "celloracle": os.path.join(output_dir, "celloracle"),
        "hotspot": os.path.join(output_dir, "hotspot"),
        "figures_qc": os.path.join(figures_dir, "qc"),
        "figures_grn": os.path.join(figures_dir, "grn"),
        "figures_hotspot": os.path.join(figures_dir, "hotspot"),
    }

    if clean:
        for directory in directories.values():
            shutil.rmtree(directory, ignore_errors=True)

    for directory in directories.values():
        os.makedirs(directory, exist_ok=True)

    print("✓ Created output directories")
    return directories


def track_files(
    output_dir: str,
    tracked_file: str = "tracked_files.txt",
) -> None:
    """
    Create a list of all output files in the output directory.

    Parameters
    ----------
    output_dir : str
        Directory to scan for files
    tracked_file : str
        Name of the output tracking file
    """
    tracked_files_path = os.path.join(output_dir, tracked_file)

    with open(tracked_files_path, "w") as f:
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if file != tracked_file:
                    f.write(os.path.join(root, file) + "\n")

    print(f"✓ Tracked files saved to: {tracked_files_path}")


def compute_input_hash(
    input_path: Optional[str] = None,
    **params: Any,
) -> str:
    """
    Compute hash of input file and parameters for checkpoint verification.

    Parameters
    ----------
    input_path : str, optional
        Path to input file
    **params
        Additional parameters to include in hash

    Returns
    -------
    str
        MD5 hash string
    """
    hash_obj = hashlib.md5()

    # Hash input file path and modification time
    if input_path:
        hash_obj.update(str(input_path).encode())
        if os.path.exists(input_path):
            mtime = os.path.getmtime(input_path)
            hash_obj.update(str(mtime).encode())

    # Hash relevant parameters
    for key, value in sorted(params.items()):
        hash_obj.update(f"{key}:{value}".encode())

    return hash_obj.hexdigest()


def write_checkpoint(
    log_dir: str,
    step_name: str,
    input_hash: str,
    **metadata: Any,
) -> str:
    """
    Write checkpoint file for a completed step.

    Parameters
    ----------
    log_dir : str
        Directory for log files
    step_name : str
        Name of the pipeline step
    input_hash : str
        Hash of input data/parameters
    **metadata
        Additional metadata to store

    Returns
    -------
    str
        Path to checkpoint file
    """
    os.makedirs(log_dir, exist_ok=True)
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    checkpoint_data = {
        "step_name": step_name,
        "input_hash": input_hash,
        "timestamp": datetime.now().isoformat(),
        "metadata": metadata,
    }

    with open(checkpoint_file, "w") as f:
        json.dump(checkpoint_data, f, indent=2)

    return checkpoint_file


def check_checkpoint(
    log_dir: str,
    step_name: str,
    input_hash: str,
) -> bool:
    """
    Check if a step has already been completed with the same input.

    Parameters
    ----------
    log_dir : str
        Directory containing checkpoint files
    step_name : str
        Name of the pipeline step
    input_hash : str
        Hash of current input data/parameters

    Returns
    -------
    bool
        True if valid checkpoint exists, False otherwise
    """
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    if not os.path.exists(checkpoint_file):
        return False

    try:
        with open(checkpoint_file, "r") as f:
            checkpoint_data = json.load(f)

        if checkpoint_data.get("input_hash") == input_hash:
            print(f"  ⏭ Checkpoint found - skipping {step_name}")
            print(f"     (completed at {checkpoint_data.get('timestamp')})")
            return True
    except Exception as e:
        print(f"  ⚠ Error reading checkpoint: {e}")

    return False


def get_checkpoint_metadata(
    log_dir: str,
    step_name: str,
) -> Optional[Dict[str, Any]]:
    """
    Get metadata from a checkpoint file.

    Parameters
    ----------
    log_dir : str
        Directory containing checkpoint files
    step_name : str
        Name of the pipeline step

    Returns
    -------
    dict or None
        Checkpoint metadata if exists, None otherwise
    """
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    if not os.path.exists(checkpoint_file):
        return None

    try:
        with open(checkpoint_file, "r") as f:
            return json.load(f)
    except Exception:
        return None
