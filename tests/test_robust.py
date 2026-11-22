"""Defective inputs never produce a traceback: an answer or a FATAL ERROR message."""
import io
import contextlib

import pytest

from resuper.cli import main

BASE = "testdata/final.pdb"


def _models():
    models, cur = [], None
    for line in open(BASE):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith("ENDMDL"):
            models.append(cur)
            cur = None
        elif cur is not None and line.startswith("ATOM"):
            cur.append(line)
    return models


def _write(path, models, bare=False):
    with open(path, "w") as f:
        for i, m in enumerate(models):
            if not bare:
                f.write(f"MODEL     {i + 1:4d}\n")
            f.writelines(m)
            if not bare:
                f.write("ENDMDL\n")
        f.write("END\n")


def _is_ca(line, num=40):
    return int(line[22:26]) == num and line[12:16] == " CA "


DEFECTS = {
    "miss_ca_all": lambda ms: [[l for l in m if not _is_ca(l)] for m in ms],
    "miss_ca_one": lambda ms: [ms[0]] + [[l for l in ms[1] if not _is_ca(l)]] + ms[2:],
    "nan": lambda ms: [ms[0]] + [[l[:30] + "     nan" + l[38:] if _is_ca(l) else l
                                   for l in ms[1]]] + ms[2:],
    "junk": lambda ms: [ms[0]] + [ms[1][:10] + ["NOT A PDB LINE\n", ms[1][10][:40] + "\n"]
                                  + ms[1][11:]] + ms[2:],
    "tiny": lambda ms: [[l for l in m if int(l[22:26]) <= 6] for m in ms],
}


@pytest.mark.parametrize("method", ["reference", "default"])
@pytest.mark.parametrize("kind", list(DEFECTS) + ["bare", "single", "empty"])
def test_defective_input(tmp_path, kind, method):
    ms = _models()
    path = str(tmp_path / f"{kind}.pdb")
    if kind == "empty":
        open(path, "w").close()
    elif kind == "single":
        _write(path, ms[:1])
    elif kind == "bare":
        _write(path, ms, bare=True)
    else:
        _write(path, DEFECTS[kind](ms))
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        rc = main(["--method", method, path])
    text = out.getvalue()
    if rc == 0:
        assert "Optimal range" in text or "WARNING" in text
    else:
        assert rc == 1 and "FATAL ERROR" in text
    if kind in ("miss_ca_all", "miss_ca_one", "nan", "junk", "bare"):
        assert rc == 0 and "Optimal range" in text
