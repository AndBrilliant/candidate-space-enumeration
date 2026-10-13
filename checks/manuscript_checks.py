#!/usr/bin/env python3
"""Manuscript claim checks: every standalone number in main_revised_17.tex that
is not already reproduced by the census, pipeline, or sweep artifacts.

Run from the checks/ directory (stdlib only):
    python3 manuscript_checks.py
Prints PASS/FAIL gates and writes checks_summary.txt.

v17 (post-audit): gates synced to the post-audit manuscript. The historical
signed-amplitude ensemble gates are retained under "historical:" labels; the
corrected finite-band conditional law is gated from koide_production_v24/.
"""
import csv, json, math, os

AK = math.sqrt(1.5) - 1.0
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.930
MU = M_E + M_MU + M_TAU                       # 1883.099374 MeV
R_M = 1.01750
U_S, D_S, S_S = 2.16 * R_M, 4.70 * R_M, 92.9 * R_M   # at mu*
M_C, M_B = 1272.9, 4186.0                     # PDG 2026 self scale
M_C24, M_B24 = 1273.0, 4183.0                 # PDG 2024 self scale (hist.)
U_C_ERR, U_B_ERR = 4.5, 6.0
M_T, M_TT = 172600.0, 162500.0                # direct, m_t(m_t)
M_W, M_Z, M_H = 80362.5, 91187.9, 125130.0    # PDG 2026, MeV

gates, lines = [], []
def check(name, got, want, tol, unit=""):
    ok = abs(got - want) <= tol
    gates.append(ok)
    lines.append(f"[{'PASS' if ok else 'FAIL'}] {name}: {got:.6g}{unit} "
                 f"(manuscript: {want:.6g}{unit}, tol {tol:g})")
    return ok

# ---- heavy forms, predictions and residuals --------------------------------
pred_c, pred_b = 3 * AK * MU, (2 + AK) * MU
check("pred m_c = 3aK mu*", pred_c, 1269.65, 0.005, " MeV")
check("pred m_b = (2+aK) mu*", pred_b, 4189.4, 0.05, " MeV")
check("pred unc m_c (mu* only)", 3 * AK * 0.090, 0.06, 0.005, " MeV")
check("pred unc m_b (mu* only)", (2 + AK) * 0.090, 0.2, 0.02, " MeV")
check("2026 residual m_c", (M_C / pred_c - 1) * 100, 0.256, 0.0005, "%")
check("2026 residual m_b", (M_B / pred_b - 1) * 100, -0.082, 0.0005, "%")
check("2024 residual m_c", (M_C24 / pred_c - 1) * 100, 0.264, 0.0005, "%")
check("2024 residual m_b", (M_B24 / pred_b - 1) * 100, -0.153, 0.0005, "%")
check("grading unit charm", U_C_ERR / M_C * 100, 0.354, 0.0005, "%")
check("grading unit bottom", U_B_ERR / M_B * 100, 0.143, 0.0005, "%")

# ---- product anchor and ratio ----------------------------------------------
check("product anchor 3aK*(2+aK)", 3 * AK * (2 + AK), 1.5, 1e-12)
ratio_pred = (2 + AK) / (3 * AK)
check("ratio seed alphaK^-2/6", AK ** -2 / 6, ratio_pred, 1e-12)
check("ratio seed value", ratio_pred, 3.2997, 0.0001)
check("measured b/c ratio", M_B / M_C, 3.2886, 0.0001)
check("ratio residual", (M_B / M_C / ratio_pred - 1) * 100, -0.34, 0.005, "%")

# ---- alpha_K completion sensitivity (measured-lepton completion) -----------
def k4_minus(m1, m2, m3):
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    p = math.sqrt(m1 * m2) + math.sqrt(m1 * m3) + math.sqrt(m2 * m3)
    return s - 2 * math.sqrt(p)

K4SQ = k4_minus(M_E, M_MU, M_TAU) ** 2        # measured-lepton completion^2
F2 = AK ** 2 * MU
rel_daK2 = K4SQ / F2 - 1                      # -2.610e-5
check("aK^2 relative displacement (measured completion)", rel_daK2, -2.6e-5,
      5e-6)
check("keV shift m_c", abs(3 * MU * AK * rel_daK2 / 2) * 1e3, 17.0, 1.0,
      " keV")
check("keV shift m_b", abs(MU * AK * rel_daK2 / 2) * 1e3, 6.0, 1.0, " keV")

# ---- light sector -----------------------------------------------------------
check("F^2 = aK^2 mu*", F2, 95.1158, 0.0005, " MeV")
check("strange-anchor residual s/F^2", (S_S / F2 - 1) * 100, -0.6204, 0.0005,
      "%")
check("census cell k4(e,mu,tau,-)^2 vs m_s", (K4SQ / S_S - 1) * 100, 0.622,
      0.001, "%")
check("F^2 vs m_s slot residual", (F2 / S_S - 1) * 100, 0.624, 0.001, "%")
check("upper mirror image m_s^2/(mu* m_d)", (S_S ** 2 / (MU * D_S) - 1) * 100,
      -0.78, 0.005, "%")
md_model = AK ** 4 * MU
check("m_d model = aK^4 mu*", md_model, 4.8043, 0.0005, " MeV")
mu_model = math.sqrt(2 * M_E * md_model)
check("m_u model terminal", mu_model, 2.21585, 0.0005, " MeV")
check("terminal mass residual", (U_S / mu_model - 1) * 100, -0.81, 0.005, "%")
check("light ratio aK^-2", AK ** -2, 19.798, 0.001)
check("ratio vs FLAG 20.1", (AK ** -2 / 20.1 - 1) * 100, -1.5, 0.01, "%")

# ---- Koide coordinate s_eff --------------------------------------------------
q_l = (M_E + M_MU + M_TAU) \
      / (math.sqrt(M_E) + math.sqrt(M_MU) + math.sqrt(M_TAU)) ** 2
check("Q_l displacement from 2/3", (2 / 3 - q_l) * 1e6, 2.20, 0.05)
s_eff = (1 + q_l) / (1 - q_l)
check("s_eff charged leptons", s_eff, 4.99996, 1e-5)
check("s=5 deviation", (s_eff / 5 - 1) * 100, -0.0008, 0.0002, "%")

# ---- s-ladder rows (Table sladder) -------------------------------------------
def q3_of(m1, m2, m3):
    return (m1 + m2 + m3) / (math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)) ** 2

q_hd = q3_of(M_C, math.sqrt(M_C * M_B), M_B)
check("heavy doublet Q3 (measured)", q_hd, 0.3526, 5e-5)
check("heavy doublet s (measured)", (1 + q_hd) / (1 - q_hd), 2.0892, 5e-4)
q_ew = q3_of(M_W, M_Z, M_T)
check("EW triple Q3", q_ew, 0.3435, 5e-5)
check("EW triple s_eff", (1 + q_ew) / (1 - q_ew), 2.0466, 5e-4)

# ---- EW pointer (top section) -------------------------------------------------
h2z = (M_H / 1000.0) ** 2 / (M_Z / 1000.0)    # GeV
check("m_H^2/m_Z", h2z, 171.7, 0.05, " GeV")
check("m_H^2/m_Z vs m_t direct", (h2z / (M_T / 1000.0) - 1) * 100, -0.52,
      0.01, "%")

# ---- up-type triple (Antusch 2026 Yukawas) -----------------------------------
Y = (7.04e-6, 3.56e-3, 0.967)
v_dir = [math.sqrt(y) for y in Y]
v_inv = [1 / math.sqrt(y) for y in Y]
qf = lambda v: sum(x * x for x in v) / sum(v) ** 2
check("Q_uct direct", qf(v_dir), 0.888, 0.001)
check("Q_uct reciprocal", qf(v_inv), 0.914, 0.001)

# ---- top probe ----------------------------------------------------------------
probe = AK ** -3 * MU
check("top probe aK^-3 mu*", probe, 165880.0, 60.0, " MeV")
check("probe identity aK^-3 = 8(2+aK)^3", AK ** -3 - 8 * (2 + AK) ** 3, 0.0,
      1e-9)
check("probe vs direct m_t", (probe / M_T - 1) * 100, -3.9, 0.05, "%")
check("probe vs m_t(m_t)", (probe / M_TT - 1) * 100, 2.1, 0.05, "%")
check("scheme separation", (M_T / M_TT - 1) * 100, 6.2, 0.05, "%")
check("target m_t for 9Q=8", 166400.0 / M_TT * 100 - 100, 2.4, 0.05, "%")

# ---- light-cascade Q3 (measured vs closed form) ------------------------------
_md, _ms = 4.70 * R_M, 92.9 * R_M
q3_meas = (MU + _md + _ms) / (math.sqrt(MU) + math.sqrt(_md) + math.sqrt(_ms)) ** 2
_a2, _a4 = AK ** 2, AK ** 4
q3_model = (1 + _a2 + _a4) / (1 + AK + _a2) ** 2
check("Q3 light cascade measured", q3_meas, 0.64816, 2e-5)
lines.append(f"note: canonical-input Q3 = {q3_meas:.6f} rounds to 0.64816, "
             f"matching the v17 manuscript; the v7-era 0.64817 drift is "
             f"resolved.")
check("Q3 model closed form", q3_model, 0.64753, 1e-5)
check("Q3 closed radical 3(9-2sqrt6)/19", 3 * (9 - 2 * math.sqrt(6)) / 19,
      0.64753, 1e-5)
check("Q3 model deviation", (q3_meas / q3_model - 1) * 100, 0.10, 0.005, "%")
check("Q3 measured s_eff", (1 + q3_meas) / (1 - q3_meas), 4.685, 0.002)

# ---- historical ensemble: archived signed-amplitude replay gates -------------
# These gate the ARCHIVED artifacts of the historical (flawed-law) replay,
# retained for provenance. The manuscript quotes these values explicitly
# labeled as the historical ensemble; do not "fix" them. The corrected law
# is gated separately below.
import re as _re
here = os.path.dirname(os.path.abspath(__file__))
_pl = os.path.join(here, "..", "pipeline")
_exact = open(os.path.join(_pl, "test_b_exact_summary.txt")).read()
_m = _re.search(r"baseline: joint5 exact = ([\d.e+-]+).*joint6 \(coupled\) = "
                r"([\d.e+-]+)", _exact)
if _m:
    j5, j6 = float(_m.group(1)), float(_m.group(2))
    check("historical: six-cut joint6 < joint5", j6 - j5, -9.42e-18, 1e-20)
    _t1 = _re.search(r"slot thresholds \(matching-or-better\): ([\d.]+%)",
                     _exact)
    if _t1:
        check("historical: archived slot-1 threshold 0.625%",
              float(_t1.group(1)[:-1]), 0.625, 0.0005)
    check("historical: six-cut ratio joint6/joint5", j6 / j5, 0.681, 0.005)
else:
    lines.append("[SKIP] six-cut gate: coupled summary line not yet present")
_conv = open(os.path.join(_pl, "test_b_convergence.txt")).read()
_mc = _re.search(r"combined joint6 \(25e6 universes\): ([\d.e+-]+)", _conv)
if _mc:
    check("historical: six-cut replicated in coupled band",
          float(_mc.group(1)), 1.9e-17, 3.0e-18)
_m5 = _re.search(r"combined \(25e6 universes\): ([\d.e+-]+)", _conv)
if _m5:
    check("historical: five-cut regression vs archived",
          float(_m5.group(1)), 2.8565e-17, 1e-20)
_sens = open(os.path.join(_pl, "sensitivity_summary.txt")).read()
_mf = _re.search(r"D0_rejection_ref: s/F\^2=([\d.e+-]+)", _sens)
if _mf:
    check("historical: slot-1 procedural frequency rounds 1.1e-3",
          float(_mf.group(1)), 1.1e-3, 1e-5)

# ---- corrected finite-band conditional law (koide_production_v24) ------------
_corr_path = os.path.join(here, "..", "koide_production_v24",
                          "corrected_results.json")
if os.path.exists(_corr_path):
    _c = json.load(open(_corr_path))
    _law = _c["laws"]
    check("corrected: baseline five-cut", _law["baseline"]["five"]["mean"],
          3.39e-17, 5e-19)
    check("corrected: baseline six-cut", _law["baseline"]["six"]["mean"],
          2.31e-17, 5e-19)
    check("corrected: doubled-span five-cut", _law["double"]["five"]["mean"],
          3.33e-18, 1e-19)
    check("corrected: doubled-span six-cut", _law["double"]["six"]["mean"],
          2.27e-18, 1e-19)
    check("corrected: uniform five-cut", _law["uniform"]["five"]["mean"],
          1.38e-27, 1e-28)
    check("corrected: uniform six-cut", _law["uniform"]["six"]["mean"],
          9.4e-28, 1e-28)
    check("corrected: halved-span exactly zero (five)",
          _law["half"]["five"]["mean"], 0.0, 0.0)
    check("corrected: halved-span exactly zero (six)",
          _law["half"]["six"]["mean"], 0.0, 0.0)
    check("corrected: half-span support proof",
          0.0 if _law["half"]["support_proof"]["no_overlap_proven"] else 1.0,
          0.0, 0.0)
    _marg_want = [1.1043e-3, 6.9840e-4, 1.0452e-3, 4.5410e-4, 1.4467e-4,
                  2.6518e-4]
    for _i, (_g, _w) in enumerate(zip(_c["marginals"]["mean"], _marg_want)):
        check(f"corrected: slot-{_i+1} marginal", _g, _w,
              max(1e-7, 0.002 * abs(_w)))
    check("corrected: physical Koide max error",
          _c["checks"]["physical_koide_max_error"], 0.0, 2.2e-6)
else:
    lines.append("[SKIP] corrected-law gates: koide_production_v24/"
                 "corrected_results.json not present")

# ---- sweep-derived: common-scale failures and chance baseline ----------------
with open(os.path.join(here, "..", "sweep", "sweep_table_pdg2026.csv"),
          newline="", encoding="ascii") as f:
    grid = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]

def interp(mu, col):
    for a, b in zip(grid, grid[1:]):
        if a["mu_GeV"] <= mu <= b["mu_GeV"]:
            t = (math.log(mu) - math.log(a["mu_GeV"])) / \
                (math.log(b["mu_GeV"]) - math.log(a["mu_GeV"]))
            la, lb = math.log(a[col]), math.log(b[col])
            return math.exp(la + t * (lb - la))
    raise ValueError(mu)

ANCH = 1.5 * (MU / 1000.0) ** 2          # (3/2) mu*^2 in GeV^2
for mu, want_p, want_r, tag in ((91.1876, -66.0, 37.0, "M_Z"),
                                (10.0, -45.0, None, "10 GeV"),
                                (4.186, -27.0, None, "m_b")):
    c, b = interp(mu, "c_GeV"), interp(mu, "b_GeV")
    check(f"common-scale product miss at {tag}", (c * b / ANCH - 1) * 100,
          want_p, 0.75, "%")
    if want_r is not None:
        check(f"common-scale ratio miss at {tag}",
              (b / c / ratio_pred - 1) * 100, want_r, 0.75, "%")

# transported-product residual invariance (Sec. rgtransport): constant by
# construction under matched transport; verify numerically over the grid
res = [(lambda c, b: c * b * (1.2729 * 4.186) / (ANCH * c * b) - 1)(
    row["c_GeV"], row["b_GeV"]) for row in grid]
check("transported residual invariance", max(abs(r - res[0]) for r in res),
      0.0, 1e-12)
check("transported residual value", res[0] * 100, 0.174, 0.001, "%")

with open(os.path.join(here, "..", "sweep", "crossings.csv"),
          newline="", encoding="ascii") as f:
    _cross = list(csv.DictReader(f))
n_cross = len(_cross)
check("sweep crossing count", n_cross, 360, 0)
expected = n_cross * (5 * 2 * 5e-3) / math.log(400.0)
check("chance baseline (uniform-expected)", expected, 3.0, 0.1)

# distinguished-scale window census (v17 sweep claim): exactly one
# configuration family (m_s^2 = 3 m_e mu*, two roots x2 orientations) lands
# in the M_Z window; no crossing lands at m_c, m_b, m_t, or mu*.
def _in_window(mu0):
    sel = [float(r["crossing_mu_GeV"]) for r in _cross
           if abs(float(r["crossing_mu_GeV"]) - mu0) / mu0 < 5e-3]
    return sel
_mz = _in_window(91.1879)
check("M_Z window crossing rows (x2 orientations)", len(_mz), 4, 0)
check("M_Z window unique roots", len(set(round(x, 3) for x in _mz)), 2, 0)
check("M_Z window root 1", min(_mz), 91.627, 0.002, " GeV")
check("M_Z window root 2", max(_mz), 91.629, 0.002, " GeV")
for _nm, _mu0 in (("m_c", 1.2729), ("m_b", 4.186), ("m_t", 172.60),
                  ("mu*", MU / 1000.0)):
    check(f"no crossings in {_nm} window", len(_in_window(_mu0)), 0, 0)

lines.append("")
lines.append(f"note: strange-anchor residual from the frozen inputs of this "
             f"repository is {(S_S / F2 - 1) * 100:.4f}%, matching the v17 "
             f"manuscript (-0.6204%, reported as -0.620%).")
lines.append(f"note: m_t(m_t) quoted 162.5 +2.1/-1.5 GeV; probe sits 3.38 GeV above "
             f"central value (3.38 = {probe/1000 - M_TT/1000:.2f}).")
lines.append(f"note: up-type running ratios 7.9889 / 7.9900 and Rivero-2026 crossing "
             f"near 280 TeV are Zenodo-sourced (no local running pipeline).")
lines.append(f"note: historical gates reference the archived signed-amplitude "
             f"ensemble (flawed law, retained for provenance); the corrected "
             f"finite-band law is gated from koide_production_v24/.")
text = "\n".join(lines) + "\n"
with open(os.path.join(here, "checks_summary.txt"), "w") as f:
    f.write(text)
print(text)
if not all(gates):
    raise SystemExit(f"{gates.count(False)} gate(s) FAILED")
print(f"ALL {len(gates)} GATES PASS")
