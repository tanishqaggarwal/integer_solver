# RESUME_O — agent O (representative sweep / channel tension).  Live notes.

Best verified so far: **39,026 / 39,033** (the existing deliverable, re-verified with
`solve_lab/checker.py`: satisfied 39026/39033, failing [12231,12270,12350,14584,18673,22044,29125]).
Nothing of mine beats it yet.

Key finding so far (see final report): the 39,026 witness's entire residual lives in 8 atoms
touched by exactly 12 equations, and SEVEN variables are private to that region
(x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864 — they occur in no atom outside it).
Exhaustive over all 4,095 subsets of those 12 equations: exactly ONE is integrally satisfiable,
the witness's own 5.  So 39,026 is exactly optimal for that knob set at that configuration.
Next: grow the region (regiongrow.py).
