#!/usr/bin/env python3
"""Replicate the fixed-mass census, then enumerate running-mass crossovers."""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from itertools import combinations
from pathlib import Path

QUARKS = ("u", "d", "s", "c", "b", "t_MS")
FORMAL_COLUMNS = {"c": "c_formal", "b": "b_formal", "t_MS": "t_formal"}
LEPTONS = {"e": 0.00051099895, "mu_lepton": 0.1056583755, "tau": 1.776930}
ANCHORS = {
    "mustar": 1.883099374,
    "P_2me": 0.0010219979,
    "sqrt32mu": 2.306316,
    "tau": 1.776930,
    "mustar_plus_tau": 3.660029,
    "two_tau": 3.553860,
    "mustar_over_2": 0.941550,
    "two_mustar": 3.766199,
    "one_point_five_mustar": 2.824649,
}
EXTRA_CONSTANTS = {"t_direct": 172.60, "sq32mu2": 5.319095}
SELF_QUARKS = {
    "u": 0.00216,
    "d": 0.00470,
    "s": 0.0929,
    "c": 1.2729,
    "b": 4.186,
    "t_MS": 162.5,
}
CENSUS_MASSES = {
    "u": 0.00216,
    "d": 0.00470,
    "s": 0.0929,
    "c": 1.2729,
    "b": 4.186,
    "t_direct": 172.60,
    **LEPTONS,
}
COLUMNS = (
    "mu_GeV", "nf_active", "alpha_s", "u_GeV", "d_GeV", "s_GeV",
    "c_GeV", "b_GeV", "t_GeV", "c_formal", "b_formal", "t_formal",
)
TOLERANCES = (0.001, 0.0025, 0.005, 0.01)
TAUTOLOGY_EPSILON = 1e-6
ANCHOR_INTERNAL_MEMBERS = frozenset((*LEPTONS, *ANCHORS))
SCALE_FAMILIES = {
    "mustar": ("mustar", 1.0),
    "mustar_over_2": ("mustar", 0.5),
    "two_mustar": ("mustar", 2.0),
    "one_point_five_mustar": ("mustar", 1.5),
    "tau": ("tau", 1.0),
    "two_tau": ("tau", 2.0),
    "e": ("e", 1.0),
    "P_2me": ("e", 2.0),
}


def load_table(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise SystemExit(f"Unexpected columns: {reader.fieldnames}")
        rows = [{key: float(value) for key, value in row.items()} for row in reader]
    if len(rows) != 300:
        raise SystemExit(f"Expected 300 rows, found {len(rows)}")
    mus = [row["mu_GeV"] for row in rows]
    if mus[0] != 1.0 or mus[-1] != 400.0 or any(b <= a for a, b in zip(mus, mus[1:])):
        raise SystemExit("mu grid must increase strictly from 1.0 to 400 GeV")
    return rows


def tier2_values(row: dict[str, float] | None = None) -> dict[str, float]:
    values = {**LEPTONS, **ANCHORS, **EXTRA_CONSTANTS}
    if row is None:
        values.update(SELF_QUARKS)
    else:
        values.update({name: row[f"{name.removesuffix('_MS')}_GeV"] for name in QUARKS})
    return values


def tier2_residual(config: tuple[str, str, str, str], values: dict[str, float]) -> float:
    role, p, q, target = config
    if role == "CENTER":
        return values[p] * values[q] / values[target] ** 2 - 1.0
    return values[p] ** 2 / (values[q] * values[target]) - 1.0


def tier2_prediction(config: tuple[str, str, str, str], values: dict[str, float]) -> float:
    role, p, q, _ = config
    return math.sqrt(values[p] * values[q]) if role == "CENTER" else values[p] ** 2 / values[q]


def is_symbolic_identity(config: tuple[str, str, str, str]) -> bool:
    role, p, q, target = config
    if not all(name in SCALE_FAMILIES for name in (p, q, target)):
        return False
    p_base, p_factor = SCALE_FAMILIES[p]
    q_base, q_factor = SCALE_FAMILIES[q]
    target_base, target_factor = SCALE_FAMILIES[target]
    if p_base != q_base or p_base != target_base:
        return False
    if role == "CENTER":
        return math.isclose(p_factor * q_factor, target_factor ** 2)
    return math.isclose(p_factor ** 2, q_factor * target_factor)


def tier2_configs() -> list[tuple[str, str, str, str]]:
    names = list(QUARKS) + list(LEPTONS) + list(ANCHORS) + ["t_direct"]
    names = list(dict.fromkeys(names))
    configs = []
    for role in ("CENTER", "EDGE"):
        for p in names:
            for q in names:
                if p == q:
                    continue
                for target in names:
                    if target in (p, q):
                        continue
                    config = (role, p, q, target)
                    if not is_symbolic_identity(config):
                        configs.append(config)
    return configs


def census_rows() -> list[dict[str, object]]:
    rows = []

    def add(class_name: str, expression: str, numerator: float, denominator: float,
            center: str, left: str, right: str) -> None:
        residual = numerator / denominator - 1.0
        rows.append({
            "class": class_name,
            "expression": expression,
            "center": center,
            "left_edge": left,
            "right_edge": right,
            "numerator_GeV2": f"{numerator:.12g}",
            "denominator_GeV2": f"{denominator:.12g}",
            "residual": f"{residual:+.12g}",
            "abs_residual": f"{abs(residual):.12g}",
            "within_0p1pct": abs(residual) < 0.001,
            "within_0p25pct": abs(residual) < 0.0025,
            "within_0p5pct": abs(residual) < 0.005,
            "within_1pct": abs(residual) < 0.01,
        })

    for a, b in combinations(CENSUS_MASSES, 2):
        for anchor in ANCHORS:
            add("P", f"{a}*{b}/{anchor}^2", CENSUS_MASSES[a] * CENSUS_MASSES[b],
                ANCHORS[anchor] ** 2, anchor, a, b)

    for center in CENSUS_MASSES:
        edges = [name for name in CENSUS_MASSES if name != center]
        for a, b in combinations(edges, 2):
            add("T", f"{center}^2/({a}*{b})", CENSUS_MASSES[center] ** 2,
                CENSUS_MASSES[a] * CENSUS_MASSES[b], center, a, b)

    for center in CENSUS_MASSES:
        for anchor in ANCHORS:
            for mass in CENSUS_MASSES:
                numerator = ANCHORS[anchor] * CENSUS_MASSES[mass]
                denominator = CENSUS_MASSES[center] ** 2
                if abs(numerator / denominator - 1.0) <= TAUTOLOGY_EPSILON:
                    continue
                add("A", f"({anchor}*{mass})/{center}^2", numerator, denominator,
                    center, anchor, mass)
    return rows


def validate_census(rows: list[dict[str, object]]) -> list[str]:
    if len(rows) != 1304:
        raise SystemExit(f"TIER 1 COUNT GATE FAILED: expected 1304, found {len(rows)}")
    counts = {name: sum(row["class"] == name for row in rows) for name in ("P", "T", "A")}
    if counts != {"P": 324, "T": 252, "A": 728}:
        raise SystemExit(f"TIER 1 CLASS COUNT GATE FAILED: {counts}")

    ranked = sorted(rows, key=lambda row: (float(row["abs_residual"]), str(row["expression"])))
    hits = [row for row in rows if float(row["abs_residual"]) < 0.005]
    expected_hit = next(row for row in rows if row["class"] == "P" and row["expression"] == "c*b/sqrt32mu^2")
    if len(hits) != 1 or hits[0] is not expected_hit:
        raise SystemExit(f"TIER 1 HIT FAILED: {[(r['expression'], r['residual']) for r in hits]}")
    hit_residual = float(expected_hit["residual"])
    if not math.isclose(hit_residual, 0.001742, rel_tol=0, abs_tol=5e-7):
        raise SystemExit(f"TIER 1 G1 VALUE FAILED: {hit_residual}")

    non_hit = ranked[1]
    if non_hit["expression"] != "(mustar*d)/s^2":
        raise SystemExit(f"TIER 1 SEP FAILED: {non_hit['expression']} = {non_hit['residual']}")
    nearest_residual = float(non_hit["residual"])
    if not (0.024 < nearest_residual < 0.026):
        raise SystemExit(f"TIER 1 SEP VALUE FAILED: {nearest_residual}")
    ratio = abs(nearest_residual / hit_residual)
    if ratio < 14.5:
        raise SystemExit(f"TIER 1 SEP SEPARATION FAILED: ratio={ratio}")

    sweep = [sum(float(row["abs_residual"]) < tolerance for row in rows) for tolerance in TOLERANCES]
    if sweep != [0, 1, 1, 1]:
        raise SystemExit(f"TIER 1 TOLERANCE GATE FAILED: {sweep}")
    return [
        "PASS Tier 1 count: 1304 comparisons (P=324, T=252, A=728).",
        f"PASS hit gate: one sub-0.5% hit, c*b/sqrt32mu^2 - 1 = {hit_residual:+.9%}.",
        f"PASS sep gate: nearest non-hit, mustar*d/s^2 - 1 = {nearest_residual:+.9%}; "
        f"separation={ratio:.3f}x.",
        "PASS tolerance sweep 0.1/0.25/0.5/1%: 0/1/1/1 hits.",
    ]


def interpolate(rows: list[dict[str, float]], mu: float, column: str) -> float:
    for left, right in zip(rows, rows[1:]):
        if left["mu_GeV"] <= mu <= right["mu_GeV"]:
            fraction = ((math.log(mu) - math.log(left["mu_GeV"])) /
                        (math.log(right["mu_GeV"]) - math.log(left["mu_GeV"])))
            return left[column] + fraction * (right[column] - left[column])
    raise ValueError(f"mu={mu} lies outside the table")


def crossing_intervals(residuals: list[float]) -> list[int]:
    return [i for i, (a, b) in enumerate(zip(residuals, residuals[1:]))
            if a == 0.0 or b == 0.0 or a * b < 0.0]


def crossing_mu(left: dict[str, float], right: dict[str, float], r_left: float, r_right: float) -> float:
    fraction = 0.0 if r_left == 0.0 else -r_left / (r_right - r_left)
    return math.exp(math.log(left["mu_GeV"]) +
                    fraction * (math.log(right["mu_GeV"]) - math.log(left["mu_GeV"])))


def formal_members(config: tuple[str, str, str, str], row: dict[str, float]) -> list[str]:
    return sorted({name for name in config[1:]
                   if name in FORMAL_COLUMNS and int(row[FORMAL_COLUMNS[name]])})


def base_tags(config: tuple[str, str, str, str]) -> list[str]:
    members = set(config[1:])
    tags = []
    if members <= ANCHOR_INTERNAL_MEMBERS:
        tags.append("ANCHOR-INTERNAL")
    if {"t_direct", "t_MS"} <= members:
        tags.append("MIXED-T")
    return tags


def tags(config: tuple[str, str, str, str], row: dict[str, float] | None = None) -> tuple[list[str], list[str]]:
    result = base_tags(config)
    formal = formal_members(config, row) if row is not None else []
    if formal:
        result.append("FORMAL")
    return result, formal


def validate_grid(rows: list[dict[str, float]]) -> list[str]:
    alpha = interpolate(rows, 91.1876, "alpha_s")
    charm = interpolate(rows, 91.1876, "c_GeV")
    bottom = interpolate(rows, 91.1876, "b_GeV")
    expected = (0.118000, 0.6327, 2.8681)
    limits = (5e-6, 5e-4, 5e-4)
    if any(abs(actual - target) > limit
           for actual, target, limit in zip((alpha, charm, bottom), expected, limits)):
        raise SystemExit(f"TABLE INTEGRITY GATE FAILED: alpha_s={alpha}, c={charm}, b={bottom}")

    config = ("CENTER", "c", "b", "sqrt32mu")
    residuals = [tier2_residual(config, tier2_values(row)) for row in rows]
    candidates = []
    for i in crossing_intervals(residuals):
        if formal_members(config, rows[i]) or formal_members(config, rows[i + 1]):
            candidates.append(crossing_mu(rows[i], rows[i + 1], residuals[i], residuals[i + 1]))
    if not candidates or not math.isclose(candidates[0], 2.07, abs_tol=0.03):
        raise SystemExit(f"FORMAL CROSSING GATE FAILED: {candidates}")
    return [
        f"PASS table integrity at 91.1876 GeV: alpha_s={alpha:.6f}, c={charm:.4f}, b={bottom:.4f} GeV.",
        f"PASS formal c*b=(3/2)mu*^2 crossing: mu={candidates[0]:.6f} GeV.",
    ]


def role_names(config: tuple[str, str, str, str]) -> tuple[str, str, str]:
    if config[0] == "CENTER":
        return "LEFT_EDGE", "RIGHT_EDGE", "CENTER"
    return "CENTER", "KNOWN_EDGE", "PREDICTED_EDGE"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_crossings(path: Path, rows: list[dict[str, float]],
                    configs: list[tuple[str, str, str, str]]) -> tuple[int, dict[str, int]]:
    output = []
    sections = {"MASS-INVOLVING": 0, "ANCHOR-INTERNAL": 0}
    for config in configs:
        role, p, q, target = config
        p_role, q_role, target_role = role_names(config)
        residuals = [tier2_residual(config, tier2_values(row)) for row in rows]
        section = "ANCHOR-INTERNAL" if "ANCHOR-INTERNAL" in base_tags(config) else "MASS-INVOLVING"
        for i in crossing_intervals(residuals):
            left, right = rows[i], rows[i + 1]
            left_tags, left_formal = tags(config, left)
            right_tags, right_formal = tags(config, right)
            all_tags = list(dict.fromkeys(left_tags + right_tags))
            all_formal = sorted(set(left_formal + right_formal))
            output.append({
                "section": section,
                "config_id": ":".join(config),
                "relation_role": role,
                "p": p,
                "p_role": p_role,
                "q": q,
                "q_role": q_role,
                "target": target,
                "target_role": target_role,
                "crossing_mu_GeV": f"{crossing_mu(left, right, residuals[i], residuals[i + 1]):.10g}",
                "local_slope_dr_dlnmu": f"{(residuals[i + 1] - residuals[i]) / (math.log(right['mu_GeV']) - math.log(left['mu_GeV'])):.10g}",
                "left_mu_GeV": f"{left['mu_GeV']:.10g}",
                "right_mu_GeV": f"{right['mu_GeV']:.10g}",
                "tags": ";".join(all_tags) if all_tags else "NONE",
                "formal_members": ";".join(all_formal),
            })
            sections[section] += 1
    write_csv(path, output)
    return len(output), sections


def nearest_record(section: str, scope: str, rank: int, config: tuple[str, str, str, str],
                   values: dict[str, float], signed: float, mu: str,
                   row: dict[str, float] | None) -> dict[str, object]:
    role, p, q, target = config
    p_role, q_role, target_role = role_names(config)
    row_tags, formal = tags(config, row)
    magnitude = abs(signed)
    return {
        "section": section,
        "scope": scope,
        "rank": rank,
        "config_id": ":".join(config),
        "relation_role": role,
        "p": p,
        "p_role": p_role,
        "q": q,
        "q_role": q_role,
        "target": target,
        "target_role": target_role,
        "mu_GeV": mu,
        "prediction_GeV": f"{tier2_prediction(config, values):.12g}",
        "target_GeV": f"{values[target]:.12g}",
        "residual": f"{signed:+.12g}",
        "abs_residual": f"{magnitude:.12g}",
        "within_0p1pct": magnitude <= 0.001,
        "within_0p25pct": magnitude <= 0.0025,
        "within_0p5pct": magnitude <= 0.005,
        "within_1pct": magnitude <= 0.01,
        "tags": ";".join(row_tags) if row_tags else "NONE",
        "formal_members": ";".join(formal),
    }


def write_nearest(path: Path, rows: list[dict[str, float]],
                  configs: list[tuple[str, str, str, str]]) -> dict[str, int]:
    output = []
    counts = {}
    self_values = tier2_values()
    for section, want_internal in (("MASS-INVOLVING", False), ("ANCHOR-INTERNAL", True)):
        selected = [config for config in configs
                    if ("ANCHOR-INTERNAL" in base_tags(config)) == want_internal]
        self_ranked = sorted((abs(tier2_residual(config, self_values)), config,
                              tier2_residual(config, self_values)) for config in selected)[:25]
        for rank, (_, config, signed) in enumerate(self_ranked, 1):
            output.append(nearest_record(section, "SELF_SCALE", rank, config, self_values,
                                         signed, "SELF_SCALE", None))

        grid_ranked = []
        for config in selected:
            best = min((abs(tier2_residual(config, tier2_values(row))), index,
                        tier2_residual(config, tier2_values(row)))
                       for index, row in enumerate(rows))
            grid_ranked.append((best[0], config, best[1], best[2]))
        grid_ranked.sort(key=lambda item: (item[0], item[1]))
        for rank, (_, config, index, signed) in enumerate(grid_ranked[:25], 1):
            row = rows[index]
            output.append(nearest_record(section, "GRID", rank, config, tier2_values(row),
                                         signed, f"{row['mu_GeV']:.10g}", row))
        counts[f"{section}:SELF_SCALE"] = len(self_ranked)
        counts[f"{section}:GRID"] = min(25, len(grid_ranked))
    write_csv(path, output)
    return counts


def write_log(path: Path, source: Path, census_gates: list[str], grid_gates: list[str],
              config_count: int, crossing_count: int, crossing_sections: dict[str, int],
              nearest_sections: dict[str, int]) -> None:
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    anchors = "\n".join(f"- {name} = {value}" for name, value in ANCHORS.items())
    text = f"""GEOMETRIC-MEAN CENSUS AND CROSSOVER ENUMERATION

INPUT AND PROVENANCE
- Input: {source.name}
- SHA-256: {digest}
- 300 supplied log-spaced rows from 1.0 through 400 GeV.
- Supplied provenance only: PDG 2026 inputs; alpha_s(M_Z)=0.1180;
  m_c(m_c)=1.2729, m_b(m_b)=4.186, m_t(m_t)=162.5 GeV; u/d/s at
  2 GeV=2.16/4.70/92.9 MeV; four-loop running; thresholds decoupled
  at m_c, m_b, and m_t.
- No web data, RunDec, or other mass-running code was fetched, built, or run.

FIXED INPUTS (GeV)
- e=0.00051099895; mu_lepton=0.1056583755; tau=1.776930.
- t_direct=172.60; self-scale t_MS=162.5.
- Self-scale u/d/s/c/b=0.00216/0.00470/0.0929/1.2729/4.186.
{anchors}
- sq32mu2=5.319095 is retained as the separately supplied squared-anchor check value, (sqrt(3/2)mu*)^2; Tier 2 uses sqrt32mu itself as the universe member.

TIER 1: CENSUS REPLICATION
- Nine charged-fermion masses use u/d/s at 2 GeV, c(m_c), b(m_b), direct top,
  and pole leptons.
- Class P: C(9,2)*9 = 324 mass products divided by squared anchors.
- Class T: 9*C(8,2) = 252 center-squared mass triples; edge pairs unordered.
- Class A: 9*9*9 minus tau^2=tau*tau = 728 comparisons. Residual is
  anchor*mass/center^2-1, matching the positive squared-anchor convention.
- Algebraic identities satisfy |residual| <= {TAUTOLOGY_EPSILON:g} using supplied rounded values.
- Hit tests use strict |residual| < tolerance.

TIER 1 VALIDATION GATES
{chr(10).join('- ' + line for line in census_gates)}

TIER 2: OPEN EXPLORATION
- Universe: six running/self-scale quarks, t_direct, three leptons, and all anchors.
- Every ordered distinct (p,q), every distinct target, and CENTER/EDGE role is evaluated.
- CENTER residual: p*q/target^2-1; roles LEFT_EDGE, RIGHT_EDGE, CENTER.
- EDGE residual: p^2/(q*target)-1; roles CENTER, KNOWN_EDGE, PREDICTED_EDGE.
- Every adjacent-grid sign change is interpolated linearly in ln(mu). Local slope is
  the bracket secant dr/dln(mu). Grid nearest misses use sampled rows, not crossing zeros.
- FORMAL tags any participating c, b, or t_MS below its supplied threshold.
- MIXED-T tags configurations containing both t_direct and t_MS.
- ANCHOR-INTERNAL means all three members are from anchors/leptons.
  Anchor-internal results are separated from mass-involving rankings. Their near-misses
  are expected near-identities, not empirical relations: mustar/tau=1.0598 creates
  AM-GM-type proximities at approximately (0.06)^2/8, about 0.04%.
- Ordered/symmetric configurations remain distinct as requested. Symbolic scale-family
  identities (for example e, 2e, or mustar, mustar/2, 2*mustar) are excluded;
  arbitrary numeric near-matches are retained.
- Tier 2 configurations after distinct-member/self-match and identity exclusions: {config_count}.
- Crossing rows: {crossing_count}; sections: {crossing_sections}.
- Nearest-miss sections (25 each requested): {nearest_sections}.

GRID VALIDATION GATES
{chr(10).join('- ' + line for line in grid_gates)}
"""
    path.write_text(text, encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", nargs="?", default="sweep_table_pdg2026.csv")
    parser.add_argument("--output-dir", default=".")
    args = parser.parse_args()
    source = Path(args.input_csv).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)

    table = load_table(source)
    census = census_rows()
    census_gates = validate_census(census)
    grid_gates = validate_grid(table)

    configs = tier2_configs()
    write_csv(output / "tier1_census.csv", census)
    crossing_count, crossing_sections = write_crossings(output / "crossings.csv", table, configs)
    nearest_sections = write_nearest(output / "nearest_misses.csv", table, configs)
    write_log(output / "conventions_provenance.txt", source, census_gates, grid_gates,
              len(configs), crossing_count, crossing_sections, nearest_sections)

    print("\n".join(census_gates + grid_gates))
    print(f"Tier 2: {len(configs)} configurations, {crossing_count} crossings.")
    print(f"Crossing sections: {crossing_sections}")
    print(f"Nearest sections: {nearest_sections}")


if __name__ == "__main__":
    main()
