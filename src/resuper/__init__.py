"""RESUPER — ordered regions and superposition ranges in protein ensembles.

RESUPER identifies rigid domains and the residues to superimpose in a
multi-conformer structure bundle. The default method is RESUPER's own
interpretable procedure; a reference-compatibility mode reproduces the
CYRANGE 2.0 report for side-by-side comparison.
"""
from .core import Params, Domain, evaluate
from .pdb import parse_pdb

__all__ = ["Params", "Domain", "evaluate", "parse_pdb"]
__version__ = "1.0.0"
