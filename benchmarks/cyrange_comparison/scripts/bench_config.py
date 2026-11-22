"""Shared settings for the RESUPER vs CYRANGE-reference benchmark.

Every script imports its paths from here, so a single file controls where
things are read from and written to. Nothing in this project writes into the
RESUPER repository (REPO_DIR); it is only read.
"""
import logging
import os
import sys
from datetime import datetime

# --- Locations -------------------------------------------------------------
BENCH_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Default: this folder lives at <repository>/benchmarks/cyrange_comparison. Override with RESUPER_REPO.
REPO_DIR = os.environ.get("RESUPER_REPO", os.path.abspath(os.path.join(BENCH_DIR, "..", "..")))
REPO_SRC = os.path.join(REPO_DIR, "src")
REPO_TESTDATA = os.path.join(REPO_DIR, "testdata")

DATA_DIR = os.path.join(BENCH_DIR, "data")
RAW_PDB_DIR = os.path.join(DATA_DIR, "raw_pdb")          # downloaded PDB files
META_DIR = os.path.join(RAW_PDB_DIR, "metadata")         # RCSB entry metadata (JSON)
MANIFEST_DRAFT = os.path.join(DATA_DIR, "dataset_manifest_DRAFT.csv")
MANIFEST = os.path.join(DATA_DIR, "dataset_manifest.csv")
INVENTORY = os.path.join(DATA_DIR, "structure_inventory.csv")
RESULTS_DIR = os.path.join(BENCH_DIR, "results")
LOG_DIR = os.path.join(BENCH_DIR, "logs")

# --- RCSB download endpoints (public, no credentials) -----------------------
RCSB_FILE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
USER_AGENT = "resuper-benchmark/1.0 (academic benchmark script)"

# Do not let Python write __pycache__ folders into the RESUPER repository when we
# import its code.
sys.dont_write_bytecode = True


def import_resuper():
    """Make the RESUPER package in REPO_SRC importable (read-only use)."""
    if REPO_SRC not in sys.path:
        sys.path.insert(0, REPO_SRC)
    import resuper  # noqa: F401
    return resuper


def setup_logging(name: str) -> logging.Logger:
    """Log to the console and to logs/<name>_<timestamp>.log."""
    os.makedirs(LOG_DIR, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(LOG_DIR, f"{name}_{stamp}.log")
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s")
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(path)):
        h.setFormatter(fmt)
        logger.addHandler(h)
    logger.info("log file: %s", path)
    return logger
