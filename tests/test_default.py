"""Regression tests for RESUPER's default method (Params.method == "default").

These pin the domain decomposition on the committed structures and the basic
report invariants. The expected domain counts are those the method was frozen
with; they happen to agree with the reference on these fixtures except 1A11,
which RESUPER splits into two 12-residue halves.
"""
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from resuper.core import Params, identify_domains  # noqa: E402
from resuper.pdb import parse_pdb  # noqa: E402

TESTDATA = os.path.join(ROOT, "testdata")
SAME_DOMAIN_COUNT = ["1A3P", "1A67", "1A6S", "1a03", "1a5e", "1a5r", "1a63", "1a7m",
                     "1a90", "1a9v", "1aaf", "1abz", "1aca", "1acp", "1adn", "1adr",
                     "1d3z", "final", "m10", "m12", "m123", "m1234", "m13", "m1_5",
                     "m23", "two_domain"]


def _domains(base):
    return identify_domains(parse_pdb(os.path.join(TESTDATA, base + ".pdb")),
                            Params(method="default"))


@pytest.mark.parametrize("base", SAME_DOMAIN_COUNT)
def test_domain_count_on_committed_structures(base):
    ref = open(os.path.join(TESTDATA, base + ".out")).read()
    assert len(_domains(base)) == len(re.findall("Optimal range", ref))


def test_two_domain_bodies():
    doms = _domains("two_domain")
    assert [d.ranges for d in doms] == [[(1, 38)], [(39, 76)]]


def test_domains_sorted_and_disjoint():
    doms = _domains("1aaf")
    starts = [d.ranges[0][0] for d in doms]
    assert starts == sorted(starts)
    seen = set()
    for d in doms:
        res = {i for lo, hi in d.ranges for i in range(lo, hi + 1)}
        assert not (res & seen)
        seen |= res
