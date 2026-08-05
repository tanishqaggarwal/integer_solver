#!/usr/bin/env python3
"""Generic partial reconstruction for a fresh instance of the trapdoor.

The intended partial = "all free inputs 0, forward-evaluate the gate orientation"
(verified on the old instance: every one of its free inputs is 0, all nonzero vars are
gate-defined and forced by the pin atoms). So:
  1. gates.jsonl gives  x_t = rhs  definitions (target t, rhs string).
  2. Greedy topological orientation: a var is ready when it is a free input, a pin, or its
     chosen gate's rhs vars are all ready. Free inputs default to 0.
  3. Seed pins (single-variable atoms force a value).
  4. Evaluate gate targets in topological order (exact integer).
Writes rebuilt_partial.json and reports checker score."""
import json, re, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

def main():
    t0 = time.time()
    gates = []
    with open('atoms/gates.jsonl') as f:
        for line in f:
            d = json.loads(line)
            gates.append((d['t'], d['rhs'], tuple(d['vids'])))
    A = load_atoms()

    # ---- pins: single-variable atoms force a value ----
    val = [0] * NVARS
    pinned = [False] * NVARS
    for p in A:
        vs = atom_vars(p)
        if len(vs) == 1:
            v = next(iter(vs))
            c0 = p.get((), 0); c1 = p.get((v,), 0); c2 = p.get((v, v), 0)
            if c2 == 0 and c1 != 0 and (-c0) % c1 == 0:
                x = (-c0) // c1
                if not pinned[v]:
                    val[v] = x; pinned[v] = True

    # ---- candidate definers ----
    cand = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates):
        cand[t].append(gi)
    targets = set(cand)
    all_vars = set(range(NVARS))
    # ready = free inputs (never a target) + pins
    ready = [False] * NVARS
    order = []
    for v in range(NVARS):
        if v not in targets or pinned[v]:
            ready[v] = True

    # greedy topological definition
    gate_unready = [0] * len(gates)
    gates_using = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates):
        u = 0
        for v in vids:
            if not ready[v]:
                u += 1
            gates_using[v].append(gi)
        gate_unready[gi] = u
    definer = {}
    q = deque(gi for gi in range(len(gates)) if gate_unready[gi] == 0)
    while q:
        gi = q.popleft()
        t, rhs, vids = gates[gi]
        if ready[t]:
            continue
        definer[t] = gi; order.append(t); ready[t] = True
        for gj in gates_using[t]:
            gate_unready[gj] -= 1
            if gate_unready[gj] == 0:
                q.append(gj)
    ndef = len(definer)
    nfree = NVARS - ndef - sum(1 for v in range(NVARS) if pinned[v] and v not in definer)
    print(f"oriented: {ndef} gate-defined, {sum(pinned)} pinned, "
          f"{NVARS-ndef-sum(1 for v in range(NVARS) if pinned[v])} free(=0)  ({time.time()-t0:.1f}s)", flush=True)

    # ---- evaluate in topo order ----
    VAR = re.compile(r'x_(\d+)')
    code = {}
    for t in order:
        gi = definer[t]
        code[t] = compile(VAR.sub(r'v[\1]', gates[gi][1]), '<r>', 'eval')
    ns = {'__builtins__': {}}
    for t in order:
        ns['v'] = val
        val[t] = eval(code[t], ns)

    out = {f"x_{i}": val[i] for i in range(NVARS)}
    json.dump(out, open('rebuilt_partial.json', 'w'))
    print(f"nonzero vars: {sum(1 for x in val if x)}  ({time.time()-t0:.1f}s)", flush=True)

if __name__ == '__main__':
    main()
