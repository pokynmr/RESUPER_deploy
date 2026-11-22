"""Regression tests for the multiconformer state-detection feature.

Given a bundle of conformers, `detect_states` decides how many distinct
conformational states the ensemble adopts (see resuper/states.py) using a
rotation-invariant triadistance-PCA descriptor plus a Ward/silhouette rule.

These tests assert the *classification* (1 vs >=2 states) and, for the
canonical two-state case 7SA5, the exact interleaved 10/10 conformer
membership -- the states alternate through the bundle rather than occupying
consecutive model numbers. The last section covers the per-state analysis:
running the ordered-region method on one state's conformers alone.
"""
import os
import re

import numpy as np
import pytest

from resuper.cli import (format_per_state_report, format_states_report, main,
                         parse_args)
from resuper.core import Params, identify_domains, identify_domains_per_state
from resuper.pdb import parse_pdb, subset_conformers
from resuper.states import StateInfo, StatesReport, detect_states

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTDATA = os.path.join(ROOT, "testdata")
# The four Pin1 bundles the manuscript validates the detector on. 7SA5, 7SUQ and
# 7SUR carry "two-state solution NMR structure" in their deposition TITLE, so for
# those the expected state count is external evidence. 8SG2 (Pin1 with the PKC
# C-terminal tail) makes no such claim; it is kept as the marginal case, and it
# does score closest to the floor.
PIN1_PATHS = {
    "7SA5": os.path.join(TESTDATA, "7sa5.pdb"),   # TITLE: "TWO-STATE SOLUTION
                                                  #  NMR STRUCTURE OF APO PIN1"
    "7SUQ": os.path.join(TESTDATA, "7suq.pdb"),   # TITLE: two-state + FFpSPR
    "7SUR": os.path.join(TESTDATA, "7sur.pdb"),   # TITLE: two-state + pCDC25c
    "8SG2": os.path.join(TESTDATA, "8sg2.pdb"),   # PKC C-tail, no two-state title
}

# 7SA5's detected membership: state A is interleaved in model number
# {1,3,6,7,10,11,14,15,17,19}, not a contiguous block of model numbers. Frozen
# as a regression anchor -- the point is the interleaving: the two states
# alternate through the bundle, so conformer order carries no information about
# the states.
KNOWN_A = frozenset({1, 3, 6, 7, 10, 11, 14, 15, 17, 19})


def _parse_fixture(path):
    """Parse a test ensemble, skipping when it is one of the fixtures that are not
    committed (the AISAR bundles: scripts/aisar_bundles.py --fetch --testdata)."""
    if not os.path.exists(path):
        pytest.skip(f"{os.path.basename(path)} not present "
                    "(run scripts/aisar_bundles.py --fetch --testdata)")
    return parse_pdb(path)


def _detect(pdb_name, path=None):
    full = path or os.path.join(TESTDATA, pdb_name + ".pdb")
    return detect_states(_parse_fixture(full))


def _partition(report):
    """Return the two partition sets for a 2-state report (conformer indices)."""
    a, b = [], []
    for i, lab in enumerate(report.labels):
        (a if lab == 0 else b).append(i + 1)
    return a, b


def _detect_7sa5():
    path = PIN1_PATHS["7SA5"]
    if not os.path.exists(path):
        pytest.skip(f"7SA5 fixture not present at {path}")
    return detect_states(parse_pdb(path))


@pytest.mark.parametrize("pdb", [
    "gluc_s12", "cdk_s12", "cdk_c4d3",
])
def test_two_state_cases_report_two_states(pdb):
    report = _detect(pdb)
    assert report.n_states == 2


@pytest.mark.parametrize("pdb", [
    "gluc_nmr", "cdk_nmr", "1A3P",
])
def test_single_state_cases_report_single_state(pdb):
    report = _detect(pdb)
    assert report.n_states == 1
    # A single state holds every conformer.
    assert report.states[0].size == report.labels.size


def test_7sa5_is_two_states():
    report = _detect_7sa5()
    assert report.n_states == 2
    assert report.silhouette.get(2, 0.0) >= 0.60


def test_7sa5_membership_matches_known_interleaving():
    report = _detect_7sa5()
    if report.n_states != 2:
        pytest.fail(f"expected 2 states, got {report.n_states}")
    a, b = _partition(report)
    # The interleaved 10/10 ground-truth split must be one of the two parts.
    assert set(a) == KNOWN_A or set(b) == KNOWN_A
    # And together they cover every conformer 1..20.
    assert sorted(set(a) | set(b)) == list(range(1, 21))


def _detect_pin1(pdb_id):
    path = PIN1_PATHS[pdb_id]
    if not os.path.exists(path):
        pytest.skip(f"{pdb_id} fixture not present at {path}")
    return detect_states(parse_pdb(path))


@pytest.mark.parametrize("pdb_id", ["7SA5", "7SUQ", "7SUR", "8SG2"])
def test_pin1_two_state_analogs_report_two_states(pdb_id):
    report = _detect_pin1(pdb_id)
    # All four Pin1 bundles must come out as 2 states; three of them are titled
    # two-state at deposition, 8SG2 is the marginal case kept as a stress test.
    assert report.n_states == 2, f"{pdb_id}: expected 2 states"
    # Clean split: both states carry a real fraction of the conformers.
    assert min(s.population for s in report.states) >= 0.30
    # Balanced 10/10 split across the 20 conformers of every one of them.
    assert report.labels.size == 20
    a, b = _partition(report)
    assert sorted(set(a) | set(b)) == list(range(1, 21))


# --- how the detection is reported -------------------------------------------
# A two-state ensemble interleaves its conformers, so a state must never be
# printed as one first..last span: for 7SA5 that collapses both states to
# overlapping "1-19" / "2-20" ranges, which hides the split entirely.

def test_members_are_the_real_model_numbers():
    # Renumber the bundle so conformer positions and MODEL ids disagree; the
    # reported membership must follow the file, not the array position.
    st = _parse_fixture(os.path.join(TESTDATA, "gluc_s12.pdb"))
    st.model_ids = [i + 21 for i in range(st.n_models)]
    report = detect_states(st)
    assert report.n_states == 2
    for state in report.states:
        assert all(21 <= m <= 20 + st.n_models for m in state.members)
        # `indices` are the 0-based positions -- a writer selects with these.
        assert all(m == st.model_ids[i] for m, i in zip(state.members, state.indices))


def test_indices_cover_every_conformer_exactly_once():
    report = _detect("cdk_c4d3")
    seen = [i for s in report.states for i in s.indices]
    assert sorted(seen) == list(range(report.labels.size))


def test_report_compacts_only_consecutive_runs():
    members = [1, 3, 6, 7, 10, 11, 14, 15, 17, 19]
    other = [m for m in range(1, 21) if m not in members]
    rep = StatesReport(
        n_states=2, labels=np.zeros(20, dtype=int), silhouette={},
        states=[StateInfo(members, list(range(10)), 0.5),
                StateInfo(other, list(range(10, 20)), 0.5)])
    line = format_states_report(20, rep)
    assert "conformers 1, 3, 6-7, 10-11, 14-15, 17, 19 -> state 1" in line
    assert "conformers 2, 4-5, 8-9, 12-13, 16, 18, 20 -> state 2" in line
    # The old first..last span must not reappear.
    assert "1-19" not in line and "2-20" not in line


def test_7sa5_report_shows_the_interleaving():
    if not os.path.exists(PIN1_PATHS["7SA5"]):
        pytest.skip("7SA5 fixture not present")
    st = parse_pdb(PIN1_PATHS["7SA5"])
    report = detect_states(st)
    line = format_states_report(st.n_models, report)
    # An interleaved state has interior gaps, so it cannot print as one span.
    assert "6-7" in line, line
    assert not re.search(r"conformers \d+-\d+ -> state", line), line
    # The two states name every conformer exactly once.
    named = [m for s in report.states for m in s.members]
    assert sorted(named) == list(range(1, st.n_models + 1))


# --- the per-state analysis ----------------------------------------------------
# A state is an ensemble in its own right, so `identify_domains` runs again over
# that state's conformers alone (`identify_domains_per_state`). On cdk_c4d3 that
# is not cosmetic: the full bundle calls 71 residues ordered in 4 domains, while
# either state of 5 conformers calls 108 ordered in 1-2. The gaps between the
# bundle's domains are the hinge *between* the states, not disorder inside
# either one. A same-size random subsample stays near 71 residues in 4 domains,
# so the widening belongs to the states and not to the smaller conformer count.

CDK_C4D3 = os.path.join(TESTDATA, "cdk_c4d3.pdb")


def _signature(domains):
    return [(d.ranges, d.n_gaps, round(d.rmsd, 9), d.n_residues) for d in domains]


def test_subset_keeps_the_sequence_and_selects_the_conformers():
    st = _parse_fixture(CDK_C4D3)
    sub = subset_conformers(st, [3, 0, 1])
    assert sub.n_models == 3
    assert sub.residues == st.residues          # same sequence, same residue objects
    # Model numbers follow the selected conformers, in the order asked for.
    assert sub.model_ids == [st.model_ids[3], st.model_ids[0], st.model_ids[1]]
    # The full bundle's Stage A verdict does not carry over to the subset.
    assert sub.warnings == [] and sub.few_core is False
    # A residue complete in every conformer of the bundle is complete in a
    # subset of them, so the subset drops nothing extra.
    assert sub.dropped == st.dropped


def test_per_state_analyses_every_multi_conformer_state():
    st = _parse_fixture(CDK_C4D3)
    report = detect_states(st)
    per_state = identify_domains_per_state(st, Params(), report)
    assert len(per_state) == report.n_states
    for state, domains in per_state:
        assert len(state.indices) >= 2
        for d in domains:
            assert d.n_residues >= 1


def test_single_state_per_state_answer_is_the_ensemble_answer():
    # One state is every conformer, so analysing "per state" must reproduce the
    # ensemble answer exactly: the feature adds nothing when there is nothing
    # to split.
    st = _parse_fixture(os.path.join(TESTDATA, "1A3P.pdb"))
    report = detect_states(st)
    assert report.n_states == 1
    per_state = identify_domains_per_state(st, Params(), report)
    assert len(per_state) == 1
    assert _signature(per_state[0][1]) == _signature(identify_domains(st, Params()))


def test_two_state_bundle_reports_more_ordered_residues_within_each_state():
    st = _parse_fixture(CDK_C4D3)
    full = identify_domains(st, Params())
    # Frozen: what the untouched bundle analysis reports for cdk_c4d3.
    assert sum(d.n_residues for d in full) == 71
    assert len(full) == 4
    per_state = identify_domains_per_state(st, Params(), detect_states(st))
    for _state, domains in per_state:
        # Within one state more residues are ordered and fewer domains are
        # needed -- the inter-domain gaps were inter-state motion.
        assert sum(d.n_residues for d in domains) > 71
        assert len(domains) < len(full)


def test_per_state_report_lists_each_state_then_its_ranges():
    st = _parse_fixture(CDK_C4D3)
    per_state = identify_domains_per_state(st, Params(), detect_states(st))
    text = format_per_state_report(per_state, Params().min_cluster)
    lines = text.splitlines()
    assert lines[0] == "    State 1 (conformers 1-5):"
    assert "    State 2 (conformers 6-10):" in lines
    # One Optimal range line per domain, in the same wording the ensemble
    # report uses.
    assert sum(1 for ln in lines if ln.startswith("    Optimal range ")) == \
        sum(len(domains) for _s, domains in per_state)


def test_per_state_flag_implies_states():
    _params, ns = parse_args(["--per-state", CDK_C4D3])
    assert ns.per_state and ns.states
    _params, ns = parse_args([CDK_C4D3])
    assert not ns.per_state and not ns.states


def test_command_line_reports_the_states_and_their_ranges(capsys):
    st = _parse_fixture(CDK_C4D3)
    expected = 4 + sum(len(d) for _s, d in
                       identify_domains_per_state(st, Params(), detect_states(st)))
    assert main(["--per-state", CDK_C4D3]) == 0
    report = capsys.readouterr().out
    assert "Conformational states: 2" in report
    assert "State 1 (conformers 1-5):" in report and "State 2 (" in report
    # The ensemble-wide lines stay: the pair is what distinguishes an
    # inter-state hinge from disorder.
    assert report.count("Optimal range ") == expected
