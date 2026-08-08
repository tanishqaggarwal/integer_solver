#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frontier as F

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'core_baseline.json')
# core frontier: 4 base solvers, s=4..8, all four effort tiers.
F.sweep(sizes=[4, 5, 6, 7, 8],
        solvers=['sa', 'pt', 'tabu', 'sb'],
        efforts=['low', 'mid', 'high', 'extreme'],
        encoder='baseline', mode='wallace',
        n_caps=dict(sa=30, pt=16, tabu=350, sb=500),
        time_cap=30.0, out_json=OUT, label='core')
print("CORE DONE")
