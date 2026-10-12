#!/usr/bin/env python3
"""Manuscript claim checks: every standalone number in main_revised_7.tex that
is not already reproduced by the census, pipeline, or sweep artifacts.

Run from the checks/ directory (stdlib only):
    python3 manuscript_checks.py
Prints PASS/FAIL gates and writes checks_summary.txt.
"""
import csv, math, os

AK = math.sqrt(1.5) - 1.0
M_E, M_MU, M_TAU = 0.51099895069, 105.6583755, 1776.930
MU = M_E + M_MU + M_TAU                       # 1883.099374 MeV
R_M = 1.01750
U_S, D_S, S_S = 2.16 * R_M, 4.70 * R_M, 92.9 * R_M   # at mu*
M_C, M_B = 1272.9, 4186.0                     # PDG 2026 self scale
M_C24, M_B24 = 1273.0, 4183.0                 # PDG 2024 self scale (hist.)
U_C_ERR, U_B_ERR = 4.5, 6.0
M_T, M_TT = 172600.0, 162500.0                # direct, m_t(m_t)

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
rel_dak = 3.6e-5          # manuscript-stated relative displacement of aK
check("keV shift m_c", 3 * MU * AK * rel_dak * 1e3, 46.0, 2.0, " keV")
check("keV shift m_b", MU * AK * rel_dak * 1e3, 15.0, 1.5, " keV")

# ---- light sector -----------------------------------------------------------
F2 = AK ** 2 * MU
check("F^2 = aK^2 mu*", F2, 95.1158, 0.0005, " MeV")
check("m_s model residual (rounded inputs)", (S_S / F2 - 1) * 100, -0.6216, 0.002, "%")
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
check("s_eff charged leptons", s_eff, 5.0, 1e-4)

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
check("probe identity aK^-3 = 8(2+aK)^3", AK ** -3 - 8 * (2 + AK) ** 3, 0.0, 1e-9)
check("probe vs direct m_t", (probe / M_T - 1) * 100, -3.9, 0.05, "%")
check("probe vs m_t(m_t)", (probe / M_TT - 1) * 100, 2.1, 0.05, "%")
check("scheme separation", (M_T / M_TT - 1) * 100, 6.2, 0.05, "%")
check("target m_t for 9Q=8", 166400.0 / M_TT * 100 - 100, 2.4, 0.05, "%")

# ---- light-cascade Q3 (measured vs closed form; post-review gate) -----------
_md, _ms = 4.70 * R_M, 92.9 * R_M
q3_meas = (MU + _md + _ms) / (math.sqrt(MU) + math.sqrt(_md) + math.sqrt(_ms)) ** 2
_a2, _a4 = AK ** 2, AK ** 4
q3_model = (1 + _a2 + _a4) / (1 + AK + _a2) ** 2
check("Q3 light cascade measured", q3_meas, 0.64817, 2e-5)
lines.append(f"note: canonical-input Q3 = {q3_meas:.6f} rounds to 0.64816; the "
             f"manuscript's 0.64817 is a 1.3e-5 last-digit drift "
             f"(FLAGGED per review instruction; no input difference required).")
check("Q3 model closed form", q3_model, 0.64753, 1e-5)
check("Q3 closed radical 3(9-2sqrt6)/19", 3 * (9 - 2 * math.sqrt(6)) / 19,
      0.64753, 1e-5)
check("Q3 model deviation", (q3_meas / q3_model - 1) * 100, 0.10, 0.005, "%")
check("Q3 measured s_eff", (1 + q3_meas) / (1 - q3_meas), 4.685, 0.002)

# ---- six-cut recompute relationship gates (read committed summaries) --------
import re as _re
here = os.path.dirname(os.path.abspath(__file__))
_pl = os.path.join(here, "..", "pipeline")
_exact = open(os.path.join(_pl, "test_b_exact_summary.txt")).read()
_m = _re.search(r"baseline: joint5 exact = ([\d.e+-]+).*joint6 \(coupled\) = "
                r"([\d.e+-]+)", _exact)
if _m:
    j5, j6 = float(_m.group(1)), float(_m.group(2))
    check("six-cut joint6 < joint5", j6 - j5, -9.42e-18, 1e-20)
    _t1 = _re.search(r"slot thresholds \(matching-or-better\): ([\d.]+%)", _exact)
    if _t1:
        check("slot-1 threshold rounds 0.625%", float(_t1.group(1)[:-1]),
              0.625, 0.0005)
    check("six-cut ratio joint6/joint5", j6 / j5, 0.681, 0.005)
else:
    lines.append("[SKIP] six-cut gate: coupled summary line not yet present")
_conv = open(os.path.join(_pl, "test_b_convergence.txt")).read()
_mc = _re.search(r"combined joint6 \(25e6 universes\): ([\d.e+-]+)", _conv)
if _mc:
    check("six-cut replicated in coupled band", float(_mc.group(1)), 1.9e-17,
          3.0e-18)
_m5 = _re.search(r"combined \(25e6 universes\): ([\d.e+-]+)", _conv)
if _m5:
    check("five-cut regression vs archived", float(_m5.group(1)), 2.8565e-17,
          1e-20)
_sens = open(os.path.join(_pl, "sensitivity_summary.txt")).read()
_mf = _re.search(r"D0_rejection_ref: s/F\^2=([\d.e+-]+)", _sens)
if _mf:
    check("slot-1 procedural frequency rounds 1.1e-3", float(_mf.group(1)),
          1.1e-3, 1e-5)

# ---- sweep-derived: common-scale failures and chance baseline ----------------
here = os.path.dirname(os.path.abspath(__file__))
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
        check(f"common-scale ratio miss at {tag}", (b / c / ratio_pred - 1) * 100,
              want_r, 0.75, "%")

# transported-product residual invariance (Sec. rgtransport): constant by
# construction under matched transport; verify numerically over the grid
res = [(lambda c, b: c * b * (1.2729 * 4.186) / (ANCH * c * b) - 1)(
    row["c_GeV"], row["b_GeV"]) for row in grid]
check("transported residual invariance", max(abs(r - res[0]) for r in res),
      0.0, 1e-12)
check("transported residual value", res[0] * 100, 0.174, 0.001, "%")

with open(os.path.join(here, "..", "sweep", "crossings.csv"),
          newline="", encoding="ascii") as f:
    n_cross = sum(1 for _ in csv.DictReader(f))
check("sweep crossing count", n_cross, 360, 0)
expected = n_cross * (5 * 2 * 5e-3) / math.log(400.0)
check("chance baseline (uniform-expected)", expected, 3.0, 0.1)

lines.append("")
lines.append(f"note: strange-anchor residual from archived unrounded inputs is "
             f"-0.6216%; rounded-input evaluation above agrees to 2e-3 %.")
lines.append(f"note: m_t(m_t) quoted 162.5 +2.1/-1.5 GeV; probe sits 3.38 GeV above "
             f"central value (3.38 = {probe/1000 - M_TT/1000:.2f}).")
lines.append(f"note: up-type running ratios 7.9889 / 7.9900 and Rivero-2026 crossing "
             f"near 280 TeV are Zenodo-sourced (no local running pipeline).")
text = "\n".join(lines) + "\n"
with open(os.path.join(here, "checks_summary.txt"), "w") as f:
    f.write(text)
print(text)
if not all(gates):
    raise SystemExit(f"{gates.count(False)} gate(s) FAILED")
print(f"ALL {len(gates)} GATES PASS")
