# Agent A LOG

## t0 — setup
- s9/*.pkl caches were absent (gitignored). Rebuilt: atomize/poly/gates/fwd (~60 s).
- Verified baseline: best/new_instance_partial_39026.json -> 39026/39033. CONFIRMED.

## Structure re-derivation (probe1..probe5, region.py, ahandles.py)
- 7 nonzero atoms; 12 equations E; 33 atoms in E; 39 equations in region; 11 zero-cost knobs.
- Private-handle census (ALL variables, not just free inputs): 1,562 atoms have a private
  variable; 326 of them with granularity 1. Prior lab census (s10/handles.py) restricted to
  free inputs -> 1,249, all p-quantised. My census is a strict superset. None in E.
- Found four knobs the prior generator list missed: x1613, x1844, x21574, x29305
  (all with 0 atoms outside the 33-atom region).
