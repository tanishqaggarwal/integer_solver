# FLEET — multi-agent solver campaign: roster, stop policy, restart procedure

This file is written and owned by the **coordinator session**, not by any solver agent.
Solver agents must not modify it. It exists so the campaign can be stopped at a session
limit and restarted afterwards without losing work or re-deciding the plan.

## Roster (10 independent solvers)

Every agent: reads `PROMPT.txt` as its mission, re-verifies the 39,026 partial itself,
treats prior lab conclusions as hypotheses, writes only into its own directory, runs no
git commands, and verifies every claimed improvement with `solve_lab/checker.py`.

| Agent | Work dir | Angle |
|-------|----------|-------|
| A | `agentA_work/` | Exact integer linear algebra: HNF / Smith normal form, LLL / BKZ on the lift lattice; retest "39,026 is optimal" under a larger movable-variable set |
| B | `agentB_work/` | Independent re-parse of the raw file; own model, ignoring the lab's atom/gate decomposition |
| C | `agentC_work/` | SAT / SMT / CP / MIP and algebraic solver encodings of the reduced core |
| D | `agentD_work/` | Large-scale stochastic search; fast exact incremental evaluator, parallel tempering / LNS across multiple restart partials |
| E | `agentE_work/` | Structural attack on the 7 failing equations: exact dependency neighborhood, computed free deformation space, simultaneous repair |
| F | `agentF_work/` | p-adic and multi-modular lifting: solve mod many primes, Hensel-lift, CRT-reconstruct |
| G | `agentG_work/` | Relaxation and rounding: homotopy, SDP / moment, MIP relaxations as a compass, exact integer snapping |
| H | `agentH_work/` | Decomposition: dependency graph, separators, treewidth, block DP over a small separator |
| I | `agentI_work/` | Build from scratch rather than repair: complete CDCL-style search, hard variables decided at the top of the tree |
| J | `agentJ_work/` | Attack the reduced parameterization (the claimed thirteen 296-bit numbers); verify it first, then solve directly in that space |

**Forbidden for all agents:** reverse-engineering how the instance file was generated — no
PRNG or coefficient-template forensics, no emission-order or index-ordering forensics aimed
at recovering a generator. Agent E originally held that angle and was redirected; artifacts
it produced under it are withdrawn, not built upon. The equations are to be analyzed as
mathematical objects only.

## Binding limits (measured, not assumed)

| Resource | Measured | Status |
|----------|----------|--------|
| CPU | 4 cores, 10 agents | **BINDING** — ~2.5x oversubscribed |
| Memory | 15 GB total, ~12 GB available | ample |
| Disk | ~30 GB available; `solve_lab/` ~350 MB | ample |
| Coordinator context | grows with each agent report | the real session limit |

CPU is the near constraint, so agents are compute-starved rather than blocked. That is
acceptable for search-style work but means wall-clock results arrive slower than a 10x
fan-out suggests.

## Stop policy

A stop is *appropriate* only when the agent's work is already durable. Order matters:

1. Send each agent a flush request (`SendMessage`): write its findings to
   `agent<X>_work/LOG.md` and `agent<X>_work/RESUME_<X>.md`, save any verified assignment
   JSON, and report back — within one tool round.
2. Wait for the reports, or for a short grace period.
3. Only then `TaskStop` anything still running.
4. Coordinator commits and pushes all agent directories.
5. Update the "Campaign state" section below.

Never `TaskStop` first: agents hold their best results in context, and an un-flushed stop
discards exactly the work the campaign is for.

Under CPU pressure specifically, prefer **thinning the fleet** over stopping it — stop the
agents whose angle has plateaued and let the rest have the cores.

## Restart procedure (after a limit resets)

1. Read this file, then `RESUME.md` and `STATE.json`.
2. Re-verify the deliverable: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`.
3. For each agent to resume: read its `agent<X>_work/RESUME_<X>.md`, then relaunch a fresh
   agent with the same angle from the roster above, pointing it at that file so it continues
   rather than restarts. Agent work dirs are committed, so a fresh container recovers them
   from the branch.
4. Do not relaunch an angle that reported a definitive dead end — record it under "Retired
   angles" below and spend the core elsewhere.

## Campaign state

- Best verified partial: **39,026 / 39,033**, `best/new_instance_partial_39026.json`,
  failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.
- Agents A–J launched; none has reported yet. No agent-produced improvement to date.
- Retired angles: generator inversion (dropped by user instruction, agent E redirected).
