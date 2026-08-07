"""W: extract ALL 383 law blocks and the EXACT INTEGER form of their congruences.
Everything recomputed from EQUATIONS.txt through my own model; nothing is taken on faith
from another agent's symbolic expansion.  Output: w_blocks.json
"""
import sys, os, re, json, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model

d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)

P_VAR, Q_VAR = 26064, 24453
# ---- alias union-find over atoms  x_a - x_b  -------------------------------
par = {}
def find(x):
    par.setdefault(x, x)
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
def uni(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb
ALIAS = re.compile(r'^x_(\d+) - x_(\d+)$')
for s in A:
    m = ALIAS.match(s)
    if m: uni(int(m.group(1)), int(m.group(2)))
PCLS = find(P_VAR)
PALIAS = set(v for v in par if find(v) == PCLS)
print('P-alias class size:', len(PALIAS))

def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]
def defs_using(v, pat):
    out = []
    for a in short(v):
        m = re.fullmatch(pat, A[a])
        if m: out.append((a, m))
    return out

blocks = []
Ecands = [(i, s) for i, s in enumerate(A) if re.fullmatch(r'x_\d+ - \(x_\d+ \+ x_%d\)' % Q_VAR, s)]
fails = []
for i0, s0 in Ecands:
    g = re.findall(r'x_(\d+)', s0); Ev, Dv = int(g[0]), int(g[1])
    b = {'Eatom': i0, 'E': Ev, 'D': Dv}
    # D = i1i2i5 sum tree:  D = (Dp + i_x) ; Dp = (i_a + i_b)
    m = re.fullmatch(r'x_%d - \(x_(\d+) \+ x_(\d+)\)' % Dv, A[model_def(Dv)] if False else '')
    # E * A^2
    c = defs_using(Ev, r'x_(\d+) - x_%d \* x_(\d+)' % Ev)
    if len(c) != 1: fails.append((i0, 'EA2', len(c))); continue
    Mv, Sv = int(c[0][1].group(1)), int(c[0][1].group(2))
    c = defs_using(Sv, r'x_%d - x_(\d+) \* x_(\d+)' % Sv)
    if len(c) != 1 or c[0][1].group(1) != c[0][1].group(2): fails.append((i0, 'A2', 0)); continue
    Av = int(c[0][1].group(1))
    # N1 = M - B^2
    c = defs_using(Mv, r'x_(\d+) - \(x_%d - x_(\d+)\)' % Mv)
    if len(c) != 1: fails.append((i0, 'N1', len(c))); continue
    N1v, Tv = int(c[0][1].group(1)), int(c[0][1].group(2))
    c = defs_using(Tv, r'x_%d - x_(\d+) \* x_(\d+)' % Tv)
    if len(c) != 1 or c[0][1].group(1) != c[0][1].group(2): fails.append((i0, 'B2', 0)); continue
    Bv = int(c[0][1].group(1))
    # A = i1 - i2 ,  B = i4 - i3
    c = defs_using(Av, r'x_%d - \(x_(\d+) - x_(\d+)\)' % Av)
    if len(c) != 1: fails.append((i0, 'Adef', len(c))); continue
    i1, i2 = int(c[0][1].group(1)), int(c[0][1].group(2))
    c = defs_using(Bv, r'x_%d - \(x_(\d+) - x_(\d+)\)' % Bv)
    if len(c) != 1: fails.append((i0, 'Bdef', len(c))); continue
    i4, i3 = int(c[0][1].group(1)), int(c[0][1].group(2))
    # N2:  A*(i3+i6) - B*(i2-i5)
    c = [t for t in defs_using(Av, r'x_(\d+) - x_(\d+) \* x_%d' % Av)
         if int(t[1].group(2)) != Av]
    if len(c) != 1: fails.append((i0, 'AH', len(c))); continue
    Wv, Gv = int(c[0][1].group(1)), int(c[0][1].group(2))
    c = defs_using(Gv, r'x_%d - \(x_(\d+) \+ x_(\d+)\)' % Gv)
    if len(c) != 1: fails.append((i0, 'Gdef', len(c))); continue
    ga, gb = int(c[0][1].group(1)), int(c[0][1].group(2))
    i6 = ga if gb == i3 else gb
    c = defs_using(Wv, r'x_(\d+) - \(x_%d - x_(\d+)\)' % Wv)
    if len(c) != 1: fails.append((i0, 'N2', len(c))); continue
    N2v, Vv = int(c[0][1].group(1)), int(c[0][1].group(2))
    c = defs_using(Vv, r'x_%d - x_%d \* x_(\d+)' % (Vv, Bv))
    if len(c) != 1: fails.append((i0, 'BJ', len(c))); continue
    Jv = int(c[0][1].group(1))
    c = defs_using(Jv, r'x_%d - \(x_(\d+) - x_(\d+)\)' % Jv)
    if len(c) != 1: fails.append((i0, 'Jdef', len(c))); continue
    ja, i5 = int(c[0][1].group(1)), int(c[0][1].group(2))
    b.update(dict(i1=i1, i2=i2, i3=i3, i4=i4, i5=i5, i6=i6, A=Av, B=Bv, N1=N1v, N2=N2v,
                  Jx=ja, G=Gv))
    # ---- the congruence layer: consumers  x_Y - c * N1  and  x_Y - c * N2
    def mults(nv):
        out = []
        for a in short(nv):
            m = re.fullmatch(r'x_(\d+) - (-?\d+) \* x_%d' % nv, A[a])
            if m: out.append((int(m.group(1)), int(m.group(2)), a))
        return out
    m1, m2 = mults(N1v), mults(N2v)
    b['m1'] = m1; b['m2'] = m2
    blocks.append(b)
print('blocks parsed:', len(blocks), 'failures:', len(fails))
from collections import Counter
print('failure kinds:', Counter(f[1] for f in fails))
json.dump({'blocks': blocks, 'fails': fails, 'palias': sorted(PALIAS)},
          open('w_blocks.json', 'w'))
