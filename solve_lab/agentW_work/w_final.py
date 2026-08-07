"""W stage 12: (i) symbolic verification of the 766 OFF-PINS, (ii) handle privacy,
(iii) identify the deliverable's single degenerate block and the 7 failing atoms."""
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
def expand_to(expr, stop, rounds=40):
    for _ in range(rounds):
        sub = {}
        for s_ in expr.free_symbols:
            v = int(str(s_)[1:])
            if v in stop or v not in DEF: continue
            sub[s_] = parse('(' + DEF[v] + ')')
        if not sub: break
        expr = expr.xreplace(sub)
    return sp.expand(expr)
PS = sp.Symbol('P')
PATS = [
 (re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'),                lambda m: (1, m.group(1), m.group(2), 1, m.group(3))),
 (re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),               lambda m: (1, m.group(1), m.group(2), -1, m.group(3))),
 (re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), lambda m: (m.group(1), m.group(2), m.group(3), 1, m.group(4))),
 (re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$'),     lambda m: (1, m.group(1), m.group(2), m.group(3), m.group(4))),
]
ok = Counter(); bad = []; t0 = time.time()
handles = []
for b in blocks:
    L = b['L']
    NL = [int(re.match(r'x_(\d+)', A[a]).group(1)) for a in short(L)
          if re.fullmatch(r'x_\d+ - \(1 - x_%d\)' % L, A[a])][0]
    for iv in (b['i5'], b['i6']):
        for a in short(iv):
            s = A[a]
            for pat, f in PATS:
                m = pat.fullmatch(s)
                if not m: continue
                aout, v1, v2, ch, H = f(m); v1, v2, H = int(v1), int(v2), int(H)
                if {v1, v2} != {NL, iv}: continue
                u = None
                m2 = re.fullmatch(r'x_%d - x_(\d+) \* x_(\d+)' % H, DEF.get(H, '') and ('x_%d - %s' % (H, DEF[H])))
                if m2:
                    p, q = int(m2.group(1)), int(m2.group(2))
                    u = q if p in PALIAS else (p if q in PALIAS else None)
                if u is None: bad.append(('nohandle', b['E'], a)); break
                handles.append(u)
                stop = {iv, L, u} | PALIAS
                e = expand_to(parse(s), stop)
                e = e.xreplace({sym(vv): PS for vv in PALIAS if sym(vv) in e.free_symbols})
                want = sp.expand(int(aout)*(1-sym(L))*sym(iv) - int(ch)*PS*sym(u))
                if sp.expand(e - want) == 0: ok['offpin'] += 1
                elif sp.expand(e + want) == 0: ok['offpin(sign)'] += 1
                else: bad.append(('mismatch', b['E'], a, str(sp.expand(e-want))[:100]))
                break
print('off-pin identities verified: %s  mismatches: %d   (%.0fs)' % (dict(ok), len(bad), time.time()-t0))
for t in bad[:5]: print('  ', t)

# --- handle privacy: every handle u must occur in exactly ONE atom (its P*u product) ------
allh = list(handles) + [cg['ring']['u'] for b in blocks for cg in b['congs']]
print('total handles:', len(allh), 'distinct:', len(set(allh)))
occ = Counter()
for u in set(allh):
    sh = [a for a in byvar[u] if len(A[a]) < 200]
    occ[len(sh)] += 1
print('handles by #short atoms they occur in:', occ.most_common())
shared=[u for u in set(allh) if len([a for a in byvar[u] if len(A[a])<200])>1]
print('handles in >1 short atom:', len(shared))
for u in shared[:6]:
    print('   x_%d :'%u, [A[a][:90] for a in byvar[u] if len(A[a])<200])
# are the H = P*u vars used anywhere but the congruence/off-pin atom?
print()
dl = json.load(open('w_deliv.json'))
deg = [r for r in dl if r['tag'].startswith('DEGENERACY')]
print('deliverable degenerate block:', deg)
off = [r['E'] for r in dl if r['tag'] == 'gate off']
print('gate-off-with-live-inputs blocks (E vars):', off)
