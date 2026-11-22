"""Read RESUPER / CYRANGE text reports into Python data.

Both programs print the same line format, e.g.

    PDB coordinate file "1a67.pdb" read, 16 rigid conformers.
    Optimal range 10..70, 82..84, 91..116: RMSD 1.17 A, 2 gaps, 90 residues.
    Optimal range A4..A92, B3..B40: RMSD 0.52 A, 1 gap, 128 residues.

Each "Optimal range" line is one domain. A domain consists of one or more
segments. In a multi-chain structure each segment carries its chain letter.
Residue numbers may be negative (e.g. tags numbered -3..-1).

A residue is identified by (chain, number). In single-chain reports the chain
is "" (empty). Use `residue_set(domain, default_chain)` to compare reports
from single-chain files with chain IDs from the structure.
"""
import re
from dataclasses import dataclass, field
from typing import List, Optional, Set, Tuple

RANGE_LINE = re.compile(
    # RMSD may be negative: CYRANGE 2.0 prints "RMSD -1.00 A" for a degenerate
    # 1-residue range (seen for 1BFW). Such domains are flagged, not dropped.
    r"Optimal range (?P<ranges>.+?): RMSD (?P<rmsd>-?[0-9.]+) A, "
    r"(?P<gaps>\d+) gaps?, (?P<nres>\d+) residues?\.")
SEGMENT = re.compile(r"^(?P<c1>[A-Za-z]?)(?P<a>-?\d+)(?:\.\.(?P<c2>[A-Za-z]?)(?P<b>-?\d+))?$")
HEADER = re.compile(r'PDB coordinate file "(?P<file>[^"]*)" read, (?P<n>\d+) rigid conformers')


@dataclass
class Domain:
    segments: List[Tuple[str, int, int]]   # (chain, first, last), inclusive
    rmsd: float
    n_gaps: int
    n_residues: int
    text: str                              # the range text as printed

    def residue_set(self, default_chain: str = "") -> Set[Tuple[str, int]]:
        out = set()
        for ch, a, b in self.segments:
            ch = ch or default_chain
            out.update((ch, i) for i in range(a, b + 1))
        return out


@dataclass
class Report:
    domains: List[Domain] = field(default_factory=list)
    n_conformers: Optional[int] = None
    input_name: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    fatal: Optional[str] = None
    program: str = ""                      # "CYRANGE" / "RESUPER" header text

    @property
    def ok(self) -> bool:
        return self.fatal is None


def parse_segments(text: str) -> List[Tuple[str, int, int]]:
    segs = []
    for part in text.split(","):
        m = SEGMENT.match(part.strip())
        if not m:
            raise ValueError(f"cannot read range segment {part!r} in {text!r}")
        a = int(m.group("a"))
        b = int(m.group("b")) if m.group("b") is not None else a
        segs.append((m.group("c1") or "", a, b))
    return segs


def parse_report(text: str) -> Report:
    rep = Report()
    for line in text.splitlines():
        s = line.strip()
        if "CYRANGE (v." in s or "RESUPER (v." in s:
            rep.program = s
        m = HEADER.search(s)
        if m:
            rep.input_name, rep.n_conformers = m.group("file"), int(m.group("n"))
            continue
        m = RANGE_LINE.search(s)
        if m:
            segs = parse_segments(m.group("ranges"))
            rep.domains.append(Domain(segs, float(m.group("rmsd")), int(m.group("gaps")),
                                      int(m.group("nres")), m.group("ranges")))
            continue
        if "FATAL ERROR" in s:
            rep.fatal = s
        elif "WARNING" in s:
            rep.warnings.append(s)
    # Consistency check: the residue count printed must equal the segments' size.
    for d in rep.domains:
        if d.rmsd < 0:
            rep.warnings.append(f"degenerate domain {d.text}: negative RMSD {d.rmsd}")
        if len(d.residue_set()) != d.n_residues:
            rep.warnings.append(f"parser check: {d.text} has {len(d.residue_set())} "
                                f"residues but the report says {d.n_residues}")
    return rep


def pooled(report: Report, default_chain: str = "") -> Set[Tuple[str, int]]:
    """All residues selected in any domain."""
    out = set()
    for d in report.domains:
        out |= d.residue_set(default_chain)
    return out


def jaccard(a: set, b: set) -> Optional[float]:
    """|A ∩ B| / |A ∪ B|; None when both sets are empty (undefined)."""
    union = a | b
    return len(a & b) / len(union) if union else None
