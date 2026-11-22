"""Unpack the bundled RESUPER/reference run outputs (kept as zip archives so the folder stays under
GitHub's 100-files-per-web-upload limit). Run once after checkout, before re-running the pipeline:
    python unpack_runs.py
Each archive is extracted next to itself; existing files are left unchanged."""
import os, zipfile
HERE = os.path.dirname(os.path.abspath(__file__))
ARCHIVES = ["results/phase4_runs.zip", "selection_study/results/phase4/runs.zip",
            "selection_study/results/confirm/runs.zip", "selection_study/results/expanded/runs.zip",
            "selection_study/results/phase9/verify16/resuper_outputs.zip"]
for a in ARCHIVES:
    path = os.path.join(HERE, a); dest = os.path.dirname(path)
    with zipfile.ZipFile(path) as zf:
        new = [m for m in zf.namelist() if not os.path.exists(os.path.join(dest, m))]
        for m in new:
            if os.path.isabs(m) or ".." in m.split("/"):
                raise SystemExit(f"unsafe path in {a}: {m}")
            zf.extract(m, dest)
    print(f"{a}: {len(new)} files extracted")
