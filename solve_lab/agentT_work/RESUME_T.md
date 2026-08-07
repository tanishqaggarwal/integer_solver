# RESUME_T — agent T (auditor).  Status: in progress.

## Confirmed by my own re-runs
1. **Deliverable 39,026/39,033** — `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
   -> satisfied 39026/39033, failing `[12231,12270,12350,14584,18673,22044,29125]`. CONFIRMED.
2. **F's peeling certificate** — `python3 agentF_work/peel_cert.py` -> `certificate verified: True,
   39033 of 39033, rank(M)=39033, dim ker(M)=0`. Reproduced from cold. CONFIRMED.
3. **Pivot magnitudes** — F claims "all pivots +-1 or +-2".  `peel_cert.py` DOES NOT CHECK THIS
   (it only tests pivot != 0).  I measured them: 37,889 pivots = 1, 1,144 = 2, none other.
   Claim CONFIRMED, but it had never been verified by the script cited as verifying it.
   (Note: pivot magnitude is irrelevant to ker(M)=0 over Q/Z — any nonzero pivot suffices.
   It matters only for the "char != 2" adornment.)
4. **M is FAITHFUL** — nobody had checked this, and it is what ker(M)=0 actually rests on.
   `agentT_work/t_faithful.py`, `t_faithful2.py`: evaluate all 39,033 atoms at a given
   assignment (no forward re-derivation) and compare rows with `M*a != 0` against
   `checker.py`'s own `load_equations`/`evaluate_all` failing list, at 10 points:
   all-zeros (11,684), 4 saved partials (7 / 9 / 12 / 20), 3 random small, 2 random big
   (39,033).  **Exact list equality at every point.**  So `eq_e = 0 <=> (M a)_e = 0` really
   does hold for the instance, and ker(M)=0 => all-atoms-zero is an equivalence. CONFIRMED.
   Bonus datum: at the deliverable exactly **7 atoms** are nonzero and they produce exactly
   **7** nonzero rows.
