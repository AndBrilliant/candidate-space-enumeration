#!/usr/bin/env python3
"""Coefficient-grammar density check (manuscript Sec. enum, paragraph 3).

The six-class anchor census cannot test the two one-mass heavy forms,
m_c = 3 aK mu* and m_b = (2+aK) mu*, because 3 aK mu* and (2+aK) mu* are not
anchors in the declared menu. Their local environment is instead enumerated
directly in the additive coefficient grammar used throughout the series:

  multiplicative points:  (p/q) * aK^k * mu*,  1 <= p,q <= 8, -2 <= k <= 4
                          (8*8*7 = 448 raw points)
  additive points:        (p + q*aK) * mu*,    0 <= p,q <= 8, not both zero
                          (81 - 1 = 80 raw points)
  total raw grid: 528 points.

Checks (manuscript claims):
  1. raw grid size == 528;
  2. the only points within 0.3% of m_c(m_c) or m_b(m_b) are the two declared
     forms themselves plus their algebraic rewritings under the quadratic
     aK^2 + 2 aK = 1/2 (6/2 = 3; (1/2) aK^-1 = 2 + aK): one distinct value
     near each mass;
  3. at the 0.1% rung only the bottom form remains;
  4. coverage: a mass drawn log-uniformly from 0.5--20 GeV falls within 0.3%
     of some grammar point with probability ~= 0.25, and within 0.1% with
     probability ~= 0.09 (union of log-window intervals / log span).

Output: coefficient_grammar.txt. Python 3 + stdlib only.
"""
import math

AK = math.sqrt(1.5) - 1.0
MU = 1883.099374  # MeV
M_C, M_B = 1272.9, 4186.0  # PDG 2026 self-scale (MeV)

points = []
for p in range(1, 9):
    for q in range(1, 9):
        for k in range(-2, 5):
            points.append(((p / q) * AK ** k * MU, f"{p}/{q}*aK^{k}"))
for p in range(0, 9):
    for q in range(0, 9):
        if p == 0 and q == 0:
            continue
        points.append(((p + q * AK) * MU, f"{p}+{q}*aK"))

assert len(points) == 528, len(points)
distinct = sorted({round(v, 6) for v, _ in points})

lines = []
lines.append("coefficient-grammar density check (528 raw points; %d distinct values)"
             % len(distinct))
lines.append("quadratic identities (machine residual):")
lines.append("  aK^2 + 2*aK - 1/2      = %+.3e" % (AK ** 2 + 2 * AK - 0.5))
lines.append("  3*aK*(2+aK) - 3/2      = %+.3e" % (3 * AK * (2 + AK) - 1.5))
lines.append("  (1+aK)^2 - 3/2         = %+.3e" % ((1 + AK) ** 2 - 1.5))
lines.append("  (1/2)*aK^-1 - (2+aK)   = %+.3e" % (0.5 / AK - (2 + AK)))
lines.append("  aK^-3 - 8*(2+aK)^3     = %+.3e" % (AK ** -3 - 8 * (2 + AK) ** 3))

for name, mobs in (("m_c(m_c)", M_C), ("m_b(m_b)", M_B)):
    near = sorted((abs(v / mobs - 1.0), v / mobs - 1.0, lab) for v, lab in points
                  if abs(v / mobs - 1.0) < 0.01)
    lines.append("")
    lines.append(f"points within 1% of {name} = {mobs} MeV:")
    for a, r, lab in near:
        rung = "0.1%" if a < 1e-3 else ("0.3%" if a < 3e-3 else "1%")
        lines.append(f"  {100 * r:+.4f}%  [{rung} rung]  {lab}")
    sub = {round(v, 6) for v, lab in points if abs(v / mobs - 1.0) < 3e-3}
    lines.append(f"  distinct values within 0.3%: {len(sub)}")
    sub1 = {round(v, 6) for v, lab in points if abs(v / mobs - 1.0) < 1e-3}
    lines.append(f"  distinct values within 0.1%: {len(sub1)}")

def coverage(w):
    ivs = sorted((math.log(v * (1 - w)), math.log(v * (1 + w)))
                 for v in distinct if 500.0 <= v <= 20000.0)
    total, cs, ce = 0.0, *ivs[0]
    for s, e in ivs[1:]:
        if s > ce:
            total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    total += ce - cs
    return total / (math.log(20000.0) - math.log(500.0))

c3, c1 = coverage(0.003), coverage(0.001)
lines.append("")
lines.append(f"log-uniform coverage of [0.5, 20] GeV by grammar windows:")
lines.append(f"  0.3% windows: P = {c3:.3f}   (manuscript: ~0.25)")
lines.append(f"  0.1% windows: P = {c1:.3f}   (manuscript: ~0.09)")

text = "\n".join(lines) + "\n"
with open("coefficient_grammar.txt", "w") as f:
    f.write(text)
print(text)
