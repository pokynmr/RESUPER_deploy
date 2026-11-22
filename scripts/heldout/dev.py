"""Development harness: score one sound-mode configuration on every benchmark at once.

Configuration overrides are module constants of resuper.core given as KEY=VALUE
(Python literals), applied in every worker process. Suites:
  truth   ubiquitin truth bench (13)             count / mixing / contamination / Jaccard
  set1    held-out set 1 (96 + 56 sweep)         same, plus sweep detections
  set2    held-out set 2 (96 + 56 sweep)
  set<k>  any further held-out set (setR: hinges on real ensembles, realhinge.py)
  real    134 real ensembles vs CYRANGE          same count, Jaccard, RMSD, catastrophes
  stab    perturbation stability                 ensembles whose domain count changes
The CYRANGE row of every suite is computed from the captured reference reports. The
criteria of PROTOCOL.md are checked at the end.

usage: /usr/bin/python3 scripts/heldout/dev.py [--suites truth,set1,set2,real]
                                               [KEY=VALUE ...]
"""
import ast, os, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.abspath("src"))
sys.path.insert(0, os.path.abspath("scripts"))
sys.path.insert(0, os.path.abspath("scripts/heldout"))
sys.path.insert(0, os.path.abspath("scripts/synthetic"))
import numpy as np
import heldout as H
import perturb as P
from eval_methods import parse_report
from resuper import core as C
from resuper.cli import format_report
from resuper.pdb import parse_pdb

CFG = {}


def init(cfg):
    for k, v in cfg.items():
        setattr(C, k, v)


def run(path):
    st = parse_pdb(path)
    return format_report("x", st.n_models, C.identify_domains(st, C.Params(method="default")))


def sets_of(text):
    return [{i for _, i in d[0]} for d in parse_report(text)]


def summary(text):
    doms = parse_report(text)
    sel = set().union(*[d[0] for d in doms]) if doms else set()
    n = sum(d[2] for d in doms)
    return len(doms), sel, (sum(d[1] * d[2] for d in doms) / n if n else float("nan"))


def truth_jobs():
    import truthbench as T
    jobs = []
    for name, spec in T.TRUTH.items():
        labels, ndom = T.truth_labels(spec)
        jobs.append((f"cywork/truth/{name}.pdb", f"testdata/oracle_truth/{name}.default",
                     (labels, ndom), False))
    return jobs


def heldout_jobs(s):
    if s == "R":
        import realhinge as RH
        return [(f"{RH.OUT}/{name}.pdb", f"{RH.REF}/{name}.out", RH.truth(case),
                 "sweep" in name) for name, case in RH.design().items()]
    H.use_set(s)
    jobs = []
    for name, case in H.design().items():
        labels, ndom = H.truth(case)
        jobs.append((f"{H.OUT}/{name}.pdb", f"{H.REF}/{name}.out", (labels, ndom),
                     "sweep" in name))
    return jobs


def score_truth(jobs, outs):
    res = {"ref": [], "default": []}
    for (pdb, ref, (labels, ndom), sweep), out in zip(jobs, outs):
        for who, text in (("ref", open(ref).read()), ("default", out)):
            res[who].append((H.score(labels, ndom, sets_of(text)), sweep))
    rows = {}
    for who, r in res.items():
        main = [x for x, sw in r if not sw]
        sw = [x for x, s in r if s]
        rows[who] = {"count": sum(x["count"] for x in main), "mixing": sum(x["mixing"] for x in main),
                     "contam": np.mean([x["contam"] for x in main]),
                     "jaccard": np.mean([x["jaccard"] for x in main]),
                     "sweep": sum(x["n"] == 2 for x in sw), "n": len(main), "nsweep": len(sw)}
    return rows


def score_real(bases, outs):
    same, jac, cat, rm_s, rm_r, diff = 0, [], [], [], [], []
    for b, out in zip(bases, outs):
        nr, sr, rr = summary(open(f"testdata/{b}.out").read())
        ns, ss, rs = summary(out)
        same += nr == ns
        if nr != ns:
            diff.append(f"{b}({ns}/{nr})")
        jac.append(len(sr & ss) / max(1, len(sr | ss)))
        rm_s.append(rs); rm_r.append(rr)
        if rs > rr + 1.0:
            cat.append(f"{b}({rs:.2f}/{rr:.2f})")
    return {"same": same, "n": len(bases), "jaccard": np.mean(jac), "rmsd": np.nanmean(rm_s),
            "ref_rmsd": np.nanmean(rm_r), "catastrophes": cat, "diff": diff}


def score_stab(bases, outs):
    rep = dict(outs)
    changed_s, changed_r = set(), set()
    for b in bases:
        for p in P.PERTS[1:]:
            if len(parse_report(rep[f"{b}__ident"])) != len(parse_report(rep[f"{b}__{p}"])):
                changed_s.add(b)
            a = parse_report(open(f"{P.REF}/{b}__ident.out").read())
            c = parse_report(open(f"{P.REF}/{b}__{p}.out").read())
            if len(a) != len(c):
                changed_r.add(b)
    return {"default": sorted(changed_s), "ref": sorted(changed_r)}


def main():
    suites = ["truth", "set1", "set2", "real"]
    cfg = {}
    for a in sys.argv[1:]:
        if a.startswith("--suites="):
            suites = a.split("=", 1)[1].split(",")
        elif "=" in a:
            k, v = a.split("=", 1)
            cfg[k] = ast.literal_eval(v)
    print("config:", cfg or "defaults")
    plan = {}
    for s in suites:
        if s == "truth":
            plan[s] = truth_jobs()
        elif s.startswith("set"):
            plan[s] = heldout_jobs(s[3:])
        elif s == "real":
            plan[s] = P.cases()
        elif s == "stab":
            plan[s] = [f"{b}__{p}" for b in P.cases() for p in P.PERTS]
    paths = []
    for s, jobs in plan.items():
        if s == "real":
            paths += [f"testdata/{b}.pdb" for b in jobs]
        elif s == "stab":
            paths += [f"{P.OUT}/{n}.pdb" for n in jobs]
        else:
            paths += [j[0] for j in jobs]
    with Pool(initializer=init, initargs=(cfg,)) as pool:
        outs = dict(zip(paths, pool.map(run, paths, chunksize=1)))
    verdict = []
    for s, jobs in plan.items():
        if s == "real":
            r = score_real(jobs, [outs[f"testdata/{b}.pdb"] for b in jobs])
            print(f"real   same n {r['same']}/{r['n']}  Jaccard {r['jaccard']:.3f}  "
                  f"RMSD {r['rmsd']:.3f} (ref {r['ref_rmsd']:.3f})  "
                  f"catastrophes {len(r['catastrophes'])} {' '.join(r['catastrophes'])}")
            print(f"       differ: {' '.join(r['diff'])}")
            verdict += [("6 real same n >= 130", r["same"] >= 130),
                        ("7 no catastrophe", not r["catastrophes"])]
        elif s == "stab":
            r = score_stab(P.cases(), [(n, outs[f"{P.OUT}/{n}.pdb"]) for n in jobs])
            print(f"stab   count changes: sound {len(r['sound'])} {' '.join(r['sound'])} | "
                  f"ref {len(r['ref'])}")
            verdict.append(("8 stability", len(r["default"]) <= len(r["ref"])))
        else:
            rows = score_truth(jobs, [outs[j[0]] for j in jobs])
            for who in ("ref", "default"):
                x = rows[who]
                sw = f"  sweep {x['sweep']}/{x['nsweep']}" if x["nsweep"] else ""
                print(f"{s:6s} {who:5s} count {x['count']:3d}/{x['n']}  mix {x['mixing']:3d}  "
                      f"contam {x['contam']:.3f}  Jaccard {x['jaccard']:.3f}{sw}")
            a, b = rows["default"], rows["ref"]
            verdict += [(f"1 {s} count", a["count"] >= b["count"]),
                        (f"2 {s} mixing", a["mixing"] <= b["mixing"]),
                        (f"3 {s} contamination", a["contam"] <= b["contam"] + 0.005),
                        (f"4 {s} Jaccard", a["jaccard"] >= b["jaccard"])]
            if a["nsweep"]:
                verdict.append((f"5 {s} sweep", a["sweep"] >= b["sweep"]))
    print("criteria: " + "  ".join(f"{k}:{'ok' if v else 'FAIL'}" for k, v in verdict))


if __name__ == "__main__":
    main()
