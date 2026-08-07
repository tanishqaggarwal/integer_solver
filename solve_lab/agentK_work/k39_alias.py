#!/usr/bin/env python3
"""K39: re-run TEST 1 against the ALIAS form.

TEST 1 searched all 38,748 variables for the predicted composition as a literal (raw and
shifted) and found nothing.  Q measures an additive/aliasing layer between a slot's mux output
and its parent's input, shape  x24468 = x13682 + 12354891*x34243.  If compositions only ever
appear aliased, a literal search comes back empty even though the value is present.

Two things are measured here.

(A) ARE THE ALIASING TERMS EVEN NONZERO IN MY CLOSURE?  I seed handles to 0.  If the additive
    terms of the hand-off layer are handles, they vanish and the literal search was already
    valid.  Check the actual example and every atom of that shape.

(B) THE ALIAS SEARCH ITSELF.  For every pair of variables that CO-OCCUR IN AN ATOM (so the
    circuit really does relate them) and every coefficient appearing in that atom, test
        v[w] +- c*v[t]  ==  V
    for V = the predicted composition's X (raw and shifted) and Y.  Co-occurring pairs, not
    all pairs, is what makes this both cheap and meaningful."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from k26_drive import FORBID

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
handles, leafsel, otherbools, wires = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
ORDER = handles + leafsel + otherbools + defvars + wires
S = FD.SHIFT
ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}


def close(on):
    seed = {u: 0 for u in handles}
    for u in leafsel: seed[u] = 1 if u in on else 0
    for u in otherbools: seed[u] = 0
    v, _ = C.close(seed, ORDER, forbid=FORBID)
    return v


def comp(es):
    R = FD.INF
    for e in es: R = FD.add(R, bypow[e])
    return R


ON_EXP = [0, 1]
v = close(set(exp2sel[e] for e in ON_EXP))
pr = comp(ON_EXP)
TARGETS = {'X_shifted': (pr[0] - S) % P, 'X_raw': pr[0] % P, 'Y': pr[1] % P}
print('ON exponents', ON_EXP, ' predicted composition X=%d... Y=%d...' % (pr[0] % 10 ** 9, pr[1] % 10 ** 9))

# ---------- (A) are the aliasing terms nonzero in this closure? -----------------------
print('\n(A) the hand-off / aliasing layer in THIS closure')
print('    x24468 =', v[24468] % 10 ** 12, '  x13682 =', v[13682] % 10 ** 12,
      '  x34243 =', v[34243], '  12354891*x34243 =', 12354891 * v[34243] % P)
print('    x24468 - x13682 =', (v[24468] - v[13682]) % P)
alias = re.compile(r'^\(x(\d+)-\(x(\d+)\+\((\d+)\*x(\d+)\)\)\)$')
alias2 = re.compile(r'^\(\(x(\d+)-x(\d+)\)-\((\d+)\*x(\d+)\)\)$')
nz_alias = 0; tot_alias = 0
for nm in C.names:
    m = alias.match(nm) or alias2.match(nm)
    if m:
        tot_alias += 1
        t = int(m.group(4))
        if v[t] % P: nz_alias += 1
print('    atoms of alias shape: %d ; of those with a NONZERO additive term here: %d'
      % (tot_alias, nz_alias))

# ---------- (B) alias search over co-occurring pairs ----------------------------------
print('\n(B) alias search: v[w] +- c*v[t] == V over variables that co-occur in an atom')
coeff_of_atom = []
for i, nm in enumerate(C.names):
    cs = set(int(x) for x in re.findall(r'(?<![\dx])(\d{1,9})(?=\*)', nm))
    cs.add(1)
    coeff_of_atom.append(cs)

pairs = set()
for i, vs in enumerate(C.avars):
    if len(vs) > 8: continue
    for a in vs:
        for b in vs:
            if a != b:
                for c in coeff_of_atom[i]:
                    pairs.add((a, b, c))
print('    (wire, wire, coefficient) triples from real atoms:', len(pairs))

found = collections.defaultdict(list)
for (w, t, c) in pairs:
    base = v[w]
    d1 = (base - c * v[t]) % P
    d2 = (base + c * v[t]) % P
    for lab, V in TARGETS.items():
        if d1 == V: found[lab].append(('x%d - %d*x%d' % (w, c, t)))
        elif d2 == V: found[lab].append(('x%d + %d*x%d' % (w, c, t)))
for lab in TARGETS:
    hits = found[lab]
    print('    %-10s : %s' % (lab, (hits[:6] if hits else 'NOT FOUND in any alias form')))

# control: the same search for a value that IS present, to prove the search works
ctrl = v[12186]
hits = [1 for (w, t, c) in list(pairs)[:200000] if (v[w] - c * v[t]) % P == ctrl % P]
print('\n    CONTROL (searching for a value known to be on a wire): %d hits -> search works'
      % len(hits))
