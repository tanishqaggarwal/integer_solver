#!/usr/bin/env python3
"""Solve the twist by circuit inversion (not local repair).

The two verifier checks are linear combos of a 17-gate sub-circuit; they vanish iff those
gates hold. The gaps route into factorable product gates: x_38045=x_15298*x_22162,
x_10156=x_15298*x_30213. Setting x_15298=1 (via boolean OR gates: free inputs x_3475=1,
x_9429=1) makes the products pass values through, so x_37892=x_30213 and x_13682=x_22162.
Set x_30213=BIGCONST, x_22162=HUGE2 and the gaps 602/1465 close with all chain multipliers 0.
Forward-evaluate and verify."""
import json, re
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

def main():
    A = load_atoms()
    base = {int(k[2:]): x for k, x in json.load(open('rebuilt_partial.json')).items()}
    V0 = base[23917]; H1 = abs(A[602].get((), 0)); H2 = abs(A[1465].get((), 0))
    BIGCONST = H1 // 8863713
    gates = []
    with open('atoms/gates.jsonl') as f:
        for line in f:
            d = json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))

    # pins (single-var atoms)
    val = [0] * NVARS; pinned = [False] * NVARS
    for p in A:
        vs = atom_vars(p)
        if len(vs) == 1:
            v = next(iter(vs)); c0 = p.get((), 0); c1 = p.get((v,), 0); c2 = p.get((v, v), 0)
            if c2 == 0 and c1 != 0 and (-c0) % c1 == 0 and not pinned[v]:
                val[v] = (-c0) // c1; pinned[v] = True

    # OVERRIDES: the input settings that invert the sub-circuit
    override = {3475: 1, 9429: 1, 22162: H2, 30213: BIGCONST, 18956: BIGCONST, 24468: H2}
    for v, x in override.items():
        val[v] = x; pinned[v] = True

    # greedy topological orientation (same as rebuild_partial)
    cand = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates): cand[t].append(gi)
    targets = set(cand)
    ready = [False] * NVARS
    for v in range(NVARS):
        if v not in targets or pinned[v]: ready[v] = True
    gate_unready = [0] * len(gates); gates_using = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates):
        u = 0
        for v in vids:
            if not ready[v]: u += 1
            gates_using[v].append(gi)
        gate_unready[gi] = u
    definer = {}; order = []
    q = deque(gi for gi in range(len(gates)) if gate_unready[gi] == 0)
    while q:
        gi = q.popleft(); t, rhs, vids = gates[gi]
        if ready[t]: continue
        definer[t] = gi; order.append(t); ready[t] = True
        for gj in gates_using[t]:
            gate_unready[gj] -= 1
            if gate_unready[gj] == 0: q.append(gj)

    VAR = re.compile(r'x_(\d+)')
    code = {t: compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in order}
    ns = {'__builtins__': {}}
    for t in order:
        ns['v'] = val; val[t] = eval(code[t], ns)

    out = {f"x_{i}": val[i] for i in range(NVARS)}
    json.dump(out, open('solved_v2.json', 'w'))
    # quick self-report on the twist atoms
    def ev(poly):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        return s
    nz = [i for i in range(len(A)) if ev(A[i]) != 0]
    print(f"nonzero atoms: {len(nz)} -> {nz[:20]}", flush=True)
    print(f"nonzero vars: {sum(1 for x in val if x)}", flush=True)

if __name__ == '__main__':
    main()
