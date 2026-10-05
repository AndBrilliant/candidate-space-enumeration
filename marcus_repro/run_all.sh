#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python deterministic_checks.py
python grammar.py
python validate.py
python whole_grammar_mc.py --n 10000 --seed 20261001
python conditional_replay.py --n 1000000 --seed 731994
python multi_seed.py --n 1000000
python importance_replay.py --n 500000
