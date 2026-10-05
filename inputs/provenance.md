# Provenance

- Lepton masses, top mass: PDG 2026 (Int. J. Mod. Phys. A 41, 2630011).
- u, d, s at 2 GeV MSbar: PDG 2026 quark mass listings; c, b: MSbar self-consistent
  values 1272.9 / 4186.0 MeV.
- Light quarks are evaluated at mu* = sum of charged-lepton masses by multiplying the
  2 GeV values by R_m = 1.01750 (frozen running factor, CRunDec-validated in the
  companion running module).
- alpha_K = sqrt(3/2) - 1 is an exact surd; no measurement input.
- All Monte Carlo thresholds are *derived* from this file at runtime by
  analysis/deterministic_relations.py. No threshold is hand-typed into any MC.
- Historical pipeline/ and candidate_space_enumeration.py are preserved unchanged
  (commits e67211b, 4fc061c, 5414e31, fec9592). The manuscript-facing lineage is
  analysis/; marcus_repro/ is the independent checker.
