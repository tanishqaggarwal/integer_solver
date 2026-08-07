"""W stage 6: the liveness layer.  Is the gate L forced into {0,1}?  If L can be a nonzero
multiple of P the law goes vacuous WITHOUT the output being killed -> a third family."""
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

# --- global census of booleanity atoms -------------------------------------
BOOL = set()
for s in A:
    m = re.fullmatch(r'x_(\d+) \* x_(\d+) - x_(\d+)', s)
    if m and m.group(1) == m.group(2) == m.group(3): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'x_(\d+) \* \(x_(\d+) - 1\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'2 \* x_(\d+) \* \(1 - x_(\d+)\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))
print('vars with an explicit booleanity atom:', len(BOOL))

# --- close the liveness cone downward from every gate ------------------------
gates = set(b['L'] for b in blocks)
print('distinct gate vars:', len(gates))
seen = set(); Q = deque(gates); kind = Counter(); leafkind = Counter(); leaves = set()
parents = {}
while Q:
    v = Q.popleft()
    if v in seen: continue
    seen.add(v)
    ds = defs(v)
    if not ds: leafkind['nodef'] += 1; leaves.add(v); continue
    ks = set()
    for s in ds:
        m = re.fullmatch(r'x_%d - x_(\d+) \* x_(\d+)' % v, s)
        if m: ks.add('prod'); parents[v] = [int(m.group(1)), int(m.group(2))]; Q += parents[v]; continue
        m = re.fullmatch(r'x_%d - \(1 - x_(\d+)\)' % v, s)
        if m: ks.add('not'); parents[v] = [int(m.group(1))]; Q += parents[v]; continue
        m = re.fullmatch(r'x_%d - x_(\d+)' % v, s)
        if m: ks.add('alias'); parents[v] = [int(m.group(1))]; Q += parents[v]; continue
        m = re.fullmatch(r'x_%d - (-?\d+)' % v, s)
        if m: ks.add('pin=' + m.group(1)); leaves.add(v); continue
        m = re.fullmatch(r'x_%d - \(x_(\d+) \+ x_(\d+)\)' % v, s)
        if m: ks.add('sum'); parents[v] = [int(m.group(1)), int(m.group(2))]; Q += parents[v]; continue
        ks.add('OTHER:' + s[:70]); leaves.add(v)
    kind[tuple(sorted(ks))] += 1
print('liveness cone size:', len(seen))
print('node kinds:', kind.most_common(10))
print('leaves:', len(leaves))
nb = [v for v in seen if v not in BOOL]
print('liveness-cone vars WITHOUT an explicit booleanity atom:', len(nb))
# a var is provably in {0,1} if it is boolean-pinned, or a product/1-x/alias of such
prov = set(v for v in seen if v in BOOL)
for _ in range(60):
    add = set()
    for v in seen:
        if v in prov: continue
        ds = defs(v)
        for s in ds:
            m = re.fullmatch(r'x_%d - (0|1)' % v, s)
            if m: add.add(v)
        p = parents.get(v)
        if p and all(q in prov for q in p):
            if any(re.fullmatch(r'x_%d - x_\d+ \* x_\d+' % v, s) or
                   re.fullmatch(r'x_%d - \(1 - x_\d+\)' % v, s) or
                   re.fullmatch(r'x_%d - x_\d+' % v, s) for s in ds): add.add(v)
    if not add: break
    prov |= add
print('liveness-cone vars PROVABLY in {0,1}:', len(prov), 'of', len(seen))
notprov = sorted(seen - prov)
print('not provably boolean:', len(notprov), notprov[:12])
for v in notprov[:8]:
    print('   x_%d  defs=%s  boolatom=%s' % (v, [s[:70] for s in defs(v)], v in BOOL))
gates_np = [b['E'] for b in blocks if b['L'] not in prov]
print('GATES not provably boolean:', len(gates_np))
json.dump({'prov': sorted(prov), 'cone': sorted(seen), 'notprov': notprov,
           'gates': sorted(gates)}, open('w_live.json', 'w'))
