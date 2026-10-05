#!/usr/bin/env bash
# One-command clean-checkout reproduction (Marcus §XI task 11).
# Regenerates every manuscript-facing numerical table into outputs/.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

echo "== 1. deterministic relations =="
(cd analysis && $PY deterministic_relations.py)

echo "== 2. census =="
(cd analysis && $PY census.py)

echo "== 3. whole-grammar Monte Carlo (Table IV) =="
(cd analysis && $PY whole_grammar_mc.py)

echo "== 4. corrected conditional joint (finite-window Koide law) =="
(cd analysis && $PY conditional_joint.py)

echo "== 5. continuum sweep =="
(cd analysis && $PY continuum_sweep.py)

echo "== 6. running factor (documentation) =="
(cd analysis && $PY running.py)

echo "== 7. regression tests =="
$PY -m pytest tests/ -q

echo "== done. outputs/ regenerated. =="
