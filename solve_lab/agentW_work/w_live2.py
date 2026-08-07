"""W stage 6b: exact value-set abstraction over the liveness cone.
Each cone var gets the SET of integer values it can take, computed bottom-up.  If every gate
lands in {0,1} the 'gate off' route is exactly L=0, which the mux alignment already kills."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter, deque
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]
def defs(v): return [A[a] for a in short(v) if re.fullmatch(r'x_%d - .*' % v, A[a])]
BOOL = set()
for s in A:
    m = re.fullmatch(r'x_(\d+) \* x_(\d+) - x_(\d+)', s)
    if m and m.group(1) == m.group(2) == m.group(3): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'x_(\d+) \* \(x_(\d+) - 1\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'2 \* x_(\d+) \* \(1 - x_(\d+)\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))

RULES = [('prod', r'x_%d - x_(\d+) \* x_(\d+)'), ('not', r'x_%d - \(1 - x_(\d+)\)'),
         ('alias', r'x_%d - x_(\d+)'), ('sum', r'x_%d - \(x_(\d+) \+ x_(\d+)\)'),
         ('diff', r'x_%d - \(x_(\d+) - x_(\d+)\)'), ('pin', r'x_%d - (-?\d+)')]
node = {}
gates = set(b['L'] for b in blocks)
Q = deque(gates); seen = set()
while Q:
    v = Q.popleft()
    if v in seen: continue
    seen.add(v)
    got = None
    for k, pat in RULES:
        for s in defs(v):
            m = re.fullmatch(pat % v, s)
            if m:
                if k == 'pin': got = (k, [int(m.group(1))])
                else: got = (k, [int(g) for g in m.groups()])
                break
        if got: break
    if got is None: node[v] = ('leaf', [])
    else:
        node[v] = got
        if got[0] != 'pin': Q += got[1]
print('cone size:', len(seen))
print('node kinds:', Counter(node[v][0] for v in seen).most_common())
leaves = [v for v in seen if node[v][0] == 'leaf']
print('leaves:', len(leaves), ' of which boolean-constrained:', sum(1 for v in leaves if v in BOOL))
print('leaves NOT boolean-constrained:', [v for v in leaves if v not in BOOL][:20])

CAP = 8
val = {}
def evalv(v, stack=()):
    if v in val: return val[v]
    if v in stack: return None          # cycle
    k, ar = node[v]
    if k == 'pin': r = frozenset(ar)
    elif k == 'leaf': r = frozenset({0, 1}) if v in BOOL else None
    else:
        ps = [evalv(p, stack + (v,)) for p in ar]
        if any(p is None for p in ps) or any(len(p) > CAP for p in ps): r = None
        else:
            if k == 'prod': r = frozenset(a * b for a in ps[0] for b in ps[1])
            elif k == 'not': r = frozenset(1 - a for a in ps[0])
            elif k == 'alias': r = ps[0]
            elif k == 'sum': r = frozenset(a + b for a in ps[0] for b in ps[1])
            elif k == 'diff': r = frozenset(a - b for a in ps[0] for b in ps[1])
    if r is not None and len(r) > CAP: r = None
    val[v] = r
    return r
sys.setrecursionlimit(20000)
res = Counter(); bad = []
for b in blocks:
    r = evalv(b['L'])
    if r is None: res['UNBOUNDED']+= 1; bad.append(b['E'])
    else:
        res[tuple(sorted(r))] += 1
        if not set(r) <= {0, 1}: bad.append(b['E'])
print('GATE value sets over all 383 blocks:', res.most_common())
print('gates not provably in {0,1}:', len(bad))
allv = Counter()
for v in seen:
    r = evalv(v)
    allv['UNBOUNDED' if r is None else ('sub{0,1}' if set(r) <= {0,1} else tuple(sorted(r)))] += 1
print('whole cone value-set census:', allv.most_common(8))
json.dump({'bad': bad}, open('w_live2.json', 'w'))
