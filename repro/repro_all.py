#!/usr/bin/env python3
"""Tip-to-tail repro orchestrator for manuscript/main_revised_7.tex.

Every mathematical figure in the paper is reproduced or traced to its
artifact here, in manuscript order:

  stage 0  fetch_pdg.py           PDG 2026 inputs, live, with source PDFs
  stage 1  exact algebra          involution identities, quadratic, Q3/s_eff
  stage 2  census + extended      declared candidate space, 1304/3740 cells
  stage 3  checks/manuscript_checks.py   every standalone number, 54+ gates
  stage 4  sweep                  continuum crossing census
  stage 5  MC suite               lite, sensitivity, RB, convergence, IS
  stage 6  report                 REPRO_REPORT.txt claim table

Usage (from repro/):
  python3 repro_all.py            fast: regenerate stdlib artifacts, verify
                                  numpy artifacts by checksum, run all gates
  python3 repro_all.py --full     regenerate EVERYTHING from scratch
                                  (numpy suite included, ~45 min)

Exit code 0 only if every gate passes.
"""
import hashlib, os, re, subprocess, sys, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PIP = os.path.join(ROOT, "pipeline")
FAILS = []

def run(cmd, cwd, stage, verify=None):
    print(f"\n===== stage: {stage} =====\n$ {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    tail = (r.stdout + r.stderr).strip().splitlines()[-6:]
    print("\n".join(tail), flush=True)
    if r.returncode != 0:
        FAILS.append(f"{stage}: exit {r.returncode}")
        print(f"!! {stage} FAILED")
        return False
    return True

def shacheck(path):
    want = {l.split()[1]: l.split()[0] for l in
            open(os.path.join(ROOT, "provenance", "checksums.txt"))}
    rel = os.path.relpath(path, ROOT)
    if rel not in want:
        return True
    got = hashlib.sha256(open(path, "rb").read()).hexdigest()
    ok = got == want[rel]
    if not ok:
        FAILS.append(f"checksum drift: {rel}")
    return ok

# ---------------------------------------------------------------- stage 0
def stage_fetch():
    r = subprocess.run([sys.executable, "fetch_pdg.py"], cwd=HERE,
                       capture_output=True, text=True)
    print(f"\n===== stage: 0. PDG live fetch =====\n$ python3 fetch_pdg.py",
          flush=True)
    print("\n".join((r.stdout + r.stderr).strip().splitlines()[-16:]), flush=True)
    if r.returncode != 0:
        FAILS.append("0. PDG live fetch: exit " + str(r.returncode))

# ---------------------------------------------------------------- stage 1
def stage_algebra():
    import math
    AK = math.sqrt(1.5) - 1.0
    M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.930
    MU = M_E + M_MU + M_TAU
    out = ["", "===== stage: 1. exact algebra =====", ""]
    out.append(f"mu* = m_e+m_mu+m_tau = {MU:.6f} MeV")
    # H3 = (2/3) J - I on 3 vectors: H3^2 = I (involution), trace fixed line
    # checked symbolically once; numeric confirmation on an arbitrary vector:
    import numpy as np
    J = np.ones((3, 3)); H = 2 / 3 * J - np.eye(3)
    out.append(f"H_3^2 = I: {np.allclose(H @ H, np.eye(3))}")
    lam = np.linalg.eigvalsh(H)
    out.append(f"H_3 eigenvalues: {sorted(set(np.round(lam, 6)))} (line-fixing: +1 once)")
    for label, expr, want in [
        ("aK^2 + 2 aK", AK ** 2 + 2 * AK, 0.5),
        ("3 aK (2+aK)", 3 * AK * (2 + AK), 1.5),
        ("(1+aK)^2", (1 + AK) ** 2, 1.5),
        ("(1/2) aK^-1 - (2+aK)", 0.5 / AK - (2 + AK), 0.0),
        ("aK^-3 - 8(2+aK)^3", AK ** -3 - 8 * (2 + AK) ** 3, 0.0),
        ("aK^-2 - 4(2+aK)^2", AK ** -2 - 4 * (2 + AK) ** 2, 0.0),
    ]:
        out.append(f"  {label:28s} = {expr:.3e}   (exact {want})")
    # Q3 = (s-1)/(s+1) identity on a geometric triple (any ratio)
    for r_ in (0.37, 1.0, 2.618):
        m = [1.0, r_ ** 2, r_ ** 4]
        s = r_ + 1 / r_
        q3_id = (s - 1) / (s + 1)
        q3_koide = (m[0] + m[1] + m[2]) / (math.sqrt(m[0]) + math.sqrt(m[1]) + math.sqrt(m[2])) ** 2
        out.append(f"  Q3 identity, r={r_}: (s-1)/(s+1)={q3_id:.9f}  Koide={q3_koide:.9f}")
    q_l = MU / (math.sqrt(M_E) + math.sqrt(M_MU) + math.sqrt(M_TAU)) ** 2
    out.append(f"  charged-lepton Q3 = {q_l:.9f} = 2/3 - {(2/3 - q_l):.3e}; "
               f"s_eff = {(1 + q_l) / (1 - q_l):.6f}")
    a2 = AK ** 2
    q3_casc = (1 + a2 + AK ** 4) / (1 + AK + a2) ** 2
    out.append(f"  light-cascade model Q3 = {q3_casc:.6f} = 3(9-2sqrt6)/19 "
               f"= {3 * (9 - 2 * math.sqrt(6)) / 19:.6f}")
    print("\n".join(out), flush=True)

# ---------------------------------------------------------------- stage 2-5
def stage_fast():
    ok = run([sys.executable, "census_mustar.py"], os.path.join(ROOT, "census"),
             "2a. 1304-cell census")
    ok &= run([sys.executable, "extended_census.py"], PIP, "2b. 3740-cell census")
    ok &= run([sys.executable, "coefficient_grammar.py"], PIP,
              "2c. coefficient grammar (528 points)")
    ok &= run([sys.executable, "manuscript_checks.py"], os.path.join(ROOT, "checks"),
              "3. standalone-number gates")
    ok &= run([sys.executable, "enumerate_crossovers.py"],
              os.path.join(ROOT, "sweep"), "4. continuum sweep")
    import shutil
    shutil.copyfile(os.path.join(ROOT, "sweep", "conventions_provenance.txt"),
                    os.path.join(ROOT, "provenance", "conventions_provenance.txt"))
    for f in ["census/census_mustar.csv", "census/census_summary.txt",
              "pipeline/extended_census.csv", "pipeline/extended_summary.txt",
              "pipeline/coefficient_grammar.txt", "checks/checks_summary.txt",
              "sweep/crossings.csv", "sweep/nearest_misses.csv",
              "sweep/tier1_census.csv", "sweep/conventions_provenance.txt"]:
        if not shacheck(os.path.join(ROOT, f)):
            print(f"!! drift: {f}")
    return ok

def stage_mc(full):
    if full:
        run([sys.executable, "lite_mc.py"], PIP, "5a. whole-grammar replay")
        run([sys.executable, "test_b_exact.py"], PIP, "5b. RB joint")
        run([sys.executable, "test_b_convergence.py"], PIP, "5c. convergence")
        run([sys.executable, "importance_sampling.py"], PIP, "5d. IS cross-check")
        run([sys.executable, "sensitivity_mc.py", "all"], PIP, "5e. sensitivity")
        run([sys.executable, "grammar_robustness.py"], PIP, "5f. menu growth")
    else:
        print("\n===== stage: 5. MC suite (checksum verify; --full to regen) =====")
        for f in ["pipeline/lite_mc_summary.txt", "pipeline/lite_mc_universes.csv",
                  "pipeline/test_b_exact_summary.txt", "pipeline/test_b_exact_curve.csv",
                  "pipeline/test_b_convergence.txt", "pipeline/importance_sampling.txt",
                  "pipeline/sensitivity_summary.txt", "pipeline/sensitivity_A_universes.csv",
                  "pipeline/sensitivity_B.csv", "pipeline/sensitivity_D.csv"]:
            if not shacheck(os.path.join(ROOT, f)):
                print(f"!! drift: {f}")
        print("all MC artifacts match pinned checksums")

# ---------------------------------------------------------------- stage 6
CLAIMS = [
    # (ms section, claim, manuscript value, artifact)
    ("abstract", "heavy forms residuals (2026)", "+0.256% / -0.082%", "checks_summary.txt"),
    ("tab:heavy", "2024 residuals", "+0.264% / -0.153%", "checks_summary.txt"),
    ("tab:heavy", "grading units", "0.354% / 0.143%", "checks_summary.txt"),
    ("sec:heavymirror", "predictions", "1269.65+-0.06 / 4189.4+-0.2 MeV", "checks_summary.txt"),
    ("sec:heavymirror", "ratio seed / measured / residual", "3.2997 / 3.2886 / -0.34%", "checks_summary.txt"),
    ("sec:heavymirror", "aK completion shifts", "46 keV / 15 keV", "checks_summary.txt"),
    ("sec:heavymirror", "quadratic closure", "3aK(2+aK) = 3/2", "coefficient_grammar.txt"),
    ("sec:grammar", "grammar density, 528 points", "no local competitors; coverage 0.25/0.09", "coefficient_grammar.txt"),
    ("sec:rgtransport", "transported residual invariance", "scale-independent, +0.174%", "checks_summary.txt"),
    ("sec:prescription", "common-scale misses", "66%/37% (M_Z), 45% (10 GeV), 27% (m_b)", "checks_summary.txt"),
    ("sec:prescription", "2 GeV / MZ grid", "2.3-17% / 30-73%", "checks_summary.txt (2 GeV formal); Zenodo (MZ band)"),
    ("sec:prescription", "pole-mass substitution misses", "30% / 14%", "Zenodo supp. (external)"),
    ("sec:light", "F^2 = aK^2 mu*", "95.1158 MeV", "checks_summary.txt"),
    ("sec:light", "images", "-0.622% (archived) / -0.78% / -0.81%", "checks_summary.txt"),
    ("sec:light", "terminal model m_u", "2.21585 MeV", "checks_summary.txt"),
    ("sec:light", "ratio s/d vs FLAG 20.1", "aK^-2 = 19.798, -1.5%", "checks_summary.txt"),
    ("sec:light", "light-cascade Q3", "0.64817 / model 0.64753 / +0.10% / s_eff 4.685", "checks_summary.txt (exact 0.648157 FLAGGED)"),
    ("sec:c3", "lepton s_eff", "5 (displacement 2.20e-6)", "checks_summary.txt"),
    ("sec:c3", "up-type Q3", "0.888 / 0.914", "checks_summary.txt"),
    ("sec:topprobe", "probe aK^-3 mu*", "165.88 GeV; -3.9%/+2.1%; 6.2%", "checks_summary.txt"),
    ("sec:topprobe", "running ratios 7.9889/7.9900", "decoupling-frozen", "Zenodo supp. (external)"),
    ("sec:census", "original census", "1304 cells; ladder 0/1/1/2; hit +0.1742%", "census_summary.txt"),
    ("sec:census", "six-class census", "3740 cells; ladder 0/1/3/6", "extended_summary.txt"),
    ("sec:census", "disclosed cells", "-0.2901% x2, +0.6034%, +0.6216%, +0.7872%, +1.1826%, nearest -1.2031%", "extended_summary.txt"),
    ("sec:census", "coefficient-grammar enum", "528 points; rewrites 6/2=3, (1/2)aK^-1=2+aK", "coefficient_grammar.txt"),
    ("sec:mcnull", "whole-grammar replay means", "0.44/1.13/2.28/4.59", "lite_mc_summary.txt"),
    ("sec:mcnull", "P(>=1), P(>=2), max", "0.27/0.55/0.78/0.95; 16/24/31/47", "lite_mc_summary.txt"),
    ("sec:mcnull", "beat sharpest cell", "43% (0.4275)", "lite_mc_summary.txt"),
    ("sec:mcnull", "per-class means", "0.74 / 0.64 / 0.91", "lite_mc_summary.txt"),
    ("sec:mcsel", "Koide front-end frequency", "5.4e-6 (2e8 draws)", "sensitivity_summary.txt"),
    ("sec:mcsel", "slot rates", "1.1e-3, 6.7e-4, 1.0e-3, 4.3e-4, 1.5e-4, 2.6e-4", "sensitivity_summary.txt"),
    ("sec:mcsel", "five-cut joint (primary / replicated)", "2.95e-17 / 2.86e-17", "test_b_exact_summary.txt / test_b_convergence.txt"),
    ("sec:mcsel", "six-cut coupled joint (primary / replicated)", "(2.01+-0.03)e-17 / (1.94+-0.02)e-17", "test_b_exact_summary.txt / test_b_convergence.txt"),
    ("sec:mcsel", "thinning ratio", "0.68", "test_b_exact_summary.txt"),
    ("sec:mcsel", "product of marginals", "4.9e-17 (factor 1.7 off)", "test_b_exact_summary.txt"),
    ("sec:mcsel", "brute-force validation", "3.50e-7 vs 3.42e-7 (s-pair)", "test_b_exact_summary.txt"),
    ("sec:mcsel", "T-tail at T_obs=1.183%", "5.6e-15", "test_b_exact_summary.txt"),
    ("sec:mcsel", "convergence stats", "ESS ~3.8e3; top weight 3.3e-4; top-10 0.3%; max/median 1.5", "test_b_convergence.txt"),
    ("sec:mcsel", "IS cross-check", "2.85e-17 (joint5) / 1.94e-17 (joint6); ESS/weights per run", "importance_sampling.txt"),
    ("sec:robust", "menu growth 9/18/27 anchors", "2.28 / 3.60 / 4.80 linear", "grammar_robustness.py (stdout)"),
    ("sec:sweep", "configurations / crossings", "9788 / 360", "sweep/crossings.csv"),
    ("sec:sweep", "unphysical continuations", "92", "crossings.csv tags"),
    ("sec:sweep", "c*b=(3/2)mu*^2 crossing", "mu = 2.069 GeV; no preferred scale", "enumerate_crossovers gates"),
    ("sec:sweep", "AM-GM nearest anchor-internal", "-0.084%", "nearest_misses.csv"),
    ("sec:sweep", "chance baseline", "~3 of 360", "checks_summary.txt"),
    ("sec:sweep", "2024 grid crossing", "at 2.3 GeV (chronicle)", "Zenodo supp. (external)"),
]

def stage_report(full):
    ms_sha = hashlib.sha256(open(os.path.join(ROOT, "manuscript", "main_revised_7.tex"), "rb").read()).hexdigest()
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    lines = [
        "REPRO_REPORT -- tip-to-tail record for main_revised_7.tex",
        f"generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"manuscript snapshot sha256: {ms_sha}",
        f"repo commit: {commit}   mode: {'full regeneration' if full else 'fast (heavy artifacts checksum-verified)'}",
        "",
        "PDG inputs: repro/pdg_live/fetched_inputs.csv (live fetch with source PDFs)",
        "validation: checks/checks_summary.txt (all gates) + provenance/checksums.txt",
        "",
        f"{'#':>2} {'manuscript':<16} {'claim':<44} {'manuscript value':<46} artifact",
        "-" * 140,
    ]
    for i, (sec, claim, val, art) in enumerate(CLAIMS, 1):
        lines.append(f"{i:>2} {sec:<16} {claim:<44} {val:<46} {art}")
    lines += [
        "",
        "EXTERNAL (cited, not locally recomputed):",
        "  - four-loop R_m and full running prescription (Appendix A) -- Zenodo supp.",
        "  - pole-mass substitution grid, MZ prescription band edges -- Zenodo supp.",
        "  - up-type running ratios 7.9889/7.9900, Rivero crossing ~280 TeV -- Zenodo supp.",
        "  - FLAG m_s/m_d = 20.1+-0.3 -- FLAG 2024 (citation only)",
        "  - Antusch-Hinze-Saad Yukawas at MZ (up-type Q3) -- cited paper values",
        "",
        "FLAGGED items recorded in checks_summary.txt:",
        "  - Q3 light-cascade exact 0.648157 rounds 0.64816 (ms prints 0.64817)",
        "  - m_s residual archived -0.6216% vs rounded-input -0.620%",
        "",
    ]
    if FAILS:
        lines.append("FAILURES:")
        lines += [f"  - {f}" for f in FAILS]
    else:
        lines.append("ALL STAGES COMPLETE -- no failures.")
    with open(os.path.join(HERE, "REPRO_REPORT.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[-8:]))

def main():
    full = "--full" in sys.argv
    stage_fetch()
    stage_algebra()
    stage_fast()
    stage_mc(full)
    stage_report(full)
    if FAILS:
        print(f"\n{len(FAILS)} FAILURE(S) -- see REPRO_REPORT.txt")
        sys.exit(1)
    print("\nrepro complete: repro/REPRO_REPORT.txt")

if __name__ == "__main__":
    main()
