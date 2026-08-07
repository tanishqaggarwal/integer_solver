# The fleet's own closures, re-scored by me with `checker.py`

Provenance read out of the first line of each `agentT_work/t_close2wj_*.log` /
`t_close2w_*.log`; scores are mine (`ah_score.py`, = `checker.py`'s loader + evaluator).

| file | \|S\| | ON-set draw | score | footprint |
|---|---|---|---|---|
| `close_T2ctl.json` | 2 | historical `[24601,2081]` | 39,018 | THE-15 |
| `close_T3.json` | 3 | `Random(7)` | 39,018 | THE-15 |
| `close_T5.json` | 5 | `Random(7)` | 39,018 | THE-15 |
| `close_T6.json` | 6 | `Random(7)` | 39,018 | THE-15 |
| `close_T7.json` | 7 | `Random(7)` | 39,018 | THE-15 |
| `close_T8.json` | 8 | `Random(7)` | 39,002 | 31 eqs (did not close) |
| `close_T8w.json`, `close_T8pair.json` | 8 | `Random(7)` | 39,018 | THE-15 |
| `close_T17w.json` | 17 | `Random(7)` | 39,003 | 30 eqs (did not close) |
| `close_T17g.json`, `close_T17j.json` | 17 | `Random(7)` | 39,018 | THE-15 |
| `close_T32g.json` | 32 | `Random(7)` **(prefix chain)** | 38,996 | 37 eqs (4 atoms) |
| `close_T32h.json`, `close_M32.json`, `close_T32mix.json` | 32 | `Random(7)` **(prefix chain)** | 39,005 | 28 eqs (3 atoms) |
| `close_T32f.json`, `close_T32z.json` | 32 | `Random(7)` **(prefix chain)** | 39,018 | THE-15 |
| `close_T64.json` | 64 | `Random(7)` **(prefix chain)** | 39,018 | THE-15 |
| `close_T128s7fix.json` | 128 | `Random(7)` **(prefix chain)** | 39,018 | THE-15 |
| `close_T128s59.json` | 128 | **`Random(59)` — independent** | 39,018 | THE-15 |
| `close_S2/S3/S5/S8/S17.json` (agent L) | 2,3,5,8,17 | agent L | 39,001 / 38,962 / 38,961 / 38,938 / 38,852 | did not close |

Two things this table says that the summary documents do not:

1. **The `|S| = 32/64/128` seed-7 rows really are one sample.**  Their log heads are the
   identical prefix `[19745, 12134, 35062, 23262, 17710, …]`.  The `|S| = 8` and `|S| = 17`
   seed-7 rows share only the first element, so they are *not* prefixes of that chain.
2. **`|S| = 128` closes to 39,018 with the identical footprint on an INDEPENDENT seed (59).**
   `UPPER_BOUND_MAP.md` §S5 records `|S| = 128` as "stalled and gave up"; that was
   `t_close2wj_T128.log`, and agent T's own T36 transposition fix resolved it.  The stalled
   run's log ends with `NO JOINT ROOT mod 116507 (sampled)` — a *sampled* miss, which is a
   statement about 400 draws out of 116,507, not about the instance.
3. **Several rows below 39,018 are at LOW `|S|` (8, 17, 32).**  A score below the ceiling is
   therefore not a high-`|S|` phenomenon at all; it is which pass happened to fire.
