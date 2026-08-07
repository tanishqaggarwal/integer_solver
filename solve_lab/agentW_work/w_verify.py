"""W stage 9: DIRECT RECOMPUTATION.  For every one of the 383 blocks, expand each of its 3
congruence atoms and 2 off-pin atoms all the way down through the definition DAG and check
the exact polynomial identity
    cong_k  ==  a_k * L * (c_k1*N1 + c_k2*N2)  -  c_k * P * u_k
    offpin  ==  a'  * (1-L) * i               -  c' * P * u'
with N1 = E*A^2 - B^2, N2 = A*(i3+i6) - B*(i2-i5), A=i1-i2, B=i4-i3, E=i1+i2+i5+Q.
Nothing is assumed from any other agent's expansion."""
import sys, os, re, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model, sympy as sp
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
PALIAS = set(json.load(open('w_blocks.json'))['palias'])
Q_VAR = 24453
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]
DEF = {}
for v in range(38748):
    for a in short(v):
        m = re.fullmatch(r'x_%d - (.*)' % v, A[a])
        if m: DEF.setdefault(v, m.group(1))
SY = {}
def sym(v):
    if v not in SY: SY[v] = sp.Symbol('x%d' % v)
    return SY[v]
def parse(s):
    ns = {('x_%s' % m): sym(int(m)) for m in set(re.findall(r'x_(\d+)', s))}
    return sp.sympify(s, locals=ns)
def expand_to(expr, stop, rounds=60):
    for _ in range(rounds):
        sub = {}
        for s_ in expr.free_symbols:
            v = int(str(s_)[1:])
            if v in stop or v not in DEF: continue
            sub[s_] = parse('(' + DEF[v] + ')')
        if not sub: break
        expr = expr.xreplace(sub)
    return sp.expand(expr)

PS = sp.Symbol('P'); QS = sp.Symbol('Q')
ok = Counter(); bad = []
t0 = time.time()
for bi, b in enumerate(blocks):
    L = b['L']
    stop = {b['i1'], b['i2'], b['i3'], b['i4'], b['i5'], b['i6'], L, Q_VAR} | PALIAS
    i1, i2, i3, i4, i5, i6 = (sym(b[k]) for k in ('i1','i2','i3','i4','i5','i6'))
    Ls = sym(L)
    Aa = i1 - i2; Bb = i4 - i3; Ee = i1 + i2 + i5 + sym(Q_VAR)
    N1 = sp.expand(Ee*Aa**2 - Bb**2); N2 = sp.expand(Aa*(i3+i6) - Bb*(i2-i5))
    for cg in b['congs']:
        r = cg['ring']
        atom = A[r['atom']]
        e = expand_to(parse(atom), stop | {r['u']})
        # substitute every P-alias symbol by PS
        e = e.xreplace({sym(v): PS for v in PALIAS if sym(v) in e.free_symbols})
        want = sp.expand(r['aout']*Ls*(r['c1']*N1 + r['c2']*N2) - r['chand']*PS*sym(r['u']))
        if sp.expand(e - want) == 0: ok['cong'] += 1
        elif sp.expand(e + want) == 0: ok['cong(sign)'] += 1
        else: bad.append(('cong', b['E'], r['atom'], sp.srepr(sp.expand(e-want))[:120]))
    if bi % 50 == 0: print('  block %d/%d  %.0fs  ok=%s bad=%d' % (bi, len(blocks), time.time()-t0, dict(ok), len(bad)), flush=True)
print('DONE %.0fs' % (time.time()-t0))
print('congruence identities verified:', dict(ok), 'mismatches:', len(bad))
for t in bad[:6]: print('  ', t)
json.dump({'ok': dict(ok), 'bad': len(bad)}, open('w_verify.json', 'w'))
