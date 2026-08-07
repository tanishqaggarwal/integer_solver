#!/usr/bin/env python3
"""Validate jmodel2.pkl: reconstruct each eq from (mult, kind, terms, atoms)
and compare with the raw LHS, evaluated at a random point modulo a big prime."""
import pickle, re, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
VAR = re.compile(r'x_(\d+)')
Q = (1 << 127) - 1
NV = 38748

M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
eqs, atoms = M['eqs'], M['atoms']

random.seed(12345)
v = [random.randrange(Q) for _ in range(NV)]
ns = {'v': v, '__builtins__': {}}

acodes = [compile(VAR.sub(r'v[\1]', a), '<a>', 'eval') for a in atoms]
print(f"compiled {len(acodes)} atoms", file=sys.stderr)
aval = [eval(c, ns) % Q for c in acodes]
print("atoms evaluated", file=sys.stderr)

bad = 0
t0 = time.time()
with open(EQ) as f:
    for idx, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        lhs = line.rsplit('=', 1)[0]
        raw = eval(compile(VAR.sub(r'v[\1]', lhs), '<e>', 'eval'), ns) % Q
        e = eqs[idx]
        s = 0
        for c, aid in e['terms']:
            s += c * aval[aid]
        s %= Q
        rec = pow(s, e['kind'], Q)
        rec = (rec * e['mult']) % Q
        if rec != raw:
            bad += 1
            if bad <= 5:
                print("MISMATCH eq", idx, e['kind'], e['mult'], len(e['terms']))
print(f"mismatches: {bad} / {len(eqs)}   ({time.time()-t0:.1f}s)")
