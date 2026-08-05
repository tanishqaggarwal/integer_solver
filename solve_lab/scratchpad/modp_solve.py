import sys, os, json, pickle, time
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import atomlib as A
import heal_harness as H
from collections import defaultdict
p = A.p
N = A.NVARS
atoms = A.ATOMS

# --- wire set ---
par0 = list(range(N))
def find0(x):
    while par0[x] != x:
        par0[x] = par0[par0[x]]; x = par0[x]
    return x
for poly in atoms:
    if all(len(vs) <= 1 for vs, c in poly):
        lin = [(vs[0], c) for vs, c in poly if len(vs) == 1]
        const = sum(c for vs, c in poly if len(vs) == 0)
        if len(lin) == 2 and const == 0 and abs(lin[0][1]) == 1 and abs(lin[1][1]) == 1:
            a, b = find0(lin[0][0]), find0(lin[1][0])
            if a != b: par0[b] = a
rw = find0(26064)
wire = set(v for v in range(N) if find0(v) == rw)

# --- collect linear-mod-p relations from atoms (drop wire terms) ---
lin_rels = []
extra_pins = []   # (var, 0) pins for wire members and wire-product outputs
for v in wire:
    extra_pins.append((v, 0))

# product-gate outputs with a wire factor -> pin to 0
# gate defs from harness
for t, gi in H.definer.items():
    _, rhs, vids = H.gates[gi]
    # a product gate: rhs contains '*' between two x_ vars; detect via poly of the atom? simpler: use vids
    # H gate vids are the referenced vars. A pure product gate has exactly form x_a*x_b.
    pass

# Instead: scan atoms of form 'x_t - (coef)*x_a*x_b' (gate-def products) => but easier:
# For every atom that is a single product term equfrom gate, check. We'll pin any gate output
# whose definition is a product with a wire factor. Use gates.jsonl rhs strings.
import re
prod_pin = set()
gate_rhs = {}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line)
        gate_rhs.setdefault(d['t'], d['rhs'])
for t, rhs in gate_rhs.items():
    # find var factors in a product
    # product gate looks like "x_a * x_b" possibly with coeff; detect two x_ multiplied
    if '*' in rhs:
        vs = [int(m) for m in re.findall(r'x_(\d+)', rhs)]
        # if it's a simple product of exactly the referenced vars and any is wire
        if any(v in wire for v in vs) and t not in wire:
            prod_pin.add(t)
for t in prod_pin:
    extra_pins.append((t, 0))
print(f"wire members: {len(wire)}; product-with-wire outputs pinned to 0: {len(prod_pin)}")

# linear relations from atoms
for ai, poly in enumerate(atoms):
    rel = defaultdict(int); const = 0; ok = True
    for vs, c in poly:
        if len(vs) == 0:
            const = (const + c) % p
        elif len(vs) == 1:
            v = vs[0]
            if v in wire:
                pass  # ≡0
            else:
                rel[v] = (rel[v] + c) % p
        else:
            # product term; if any factor wire -> 0; else genuine quadratic -> skip this atom (nonlinear)
            if any(v in wire for v in vs):
                pass
            else:
                ok = False; break
    if ok:
        rel = {k: v for k, v in rel.items() if v % p}
        lin_rels.append((rel, const % p))
print(f"linear-mod-p relations: {len(lin_rels)}")

# --- union-find + gaussian solve mod p (like lin_solve/gauss_solve) ---
CONST = N
parU = list(range(N+1)); mult = [1]*(N+1); off = [0]*(N+1)
def find2(x):
    chain = []
    r = x
    while parU[r] != r:
        chain.append(r); r = parU[r]
    m, o = 1, 0
    for v in reversed(chain):
        m = (mult[v]*m) % p; o = (mult[v]*o + off[v]) % p
    for v in chain:
        parU[v] = r; mult[v] = (mult[v]* (1)) % p  # will recompute below properly
    # recompute properly per node
    # simpler: recompute mult/off for each node in chain relative to root
    # do second pass
    accm, acco = {}, {}
    m2, o2 = 1, 0
    for v in reversed(chain):
        m2 = (mult[v]*m2) % p; o2 = (mult[v]*o2 + off[v]) % p
    # We only need (r, m, o) for x itself:
    # recompute from scratch to be safe
    return r
# The above find2 is messy; implement clean recursive with memo via iterative compose.
sys.setrecursionlimit(3000000)
parU = list(range(N+1)); mult = [1]*(N+1); off = [0]*(N+1)
def find3(x):
    if parU[x] == x:
        return x, 1, 0
    r, m, o = find3(parU[x])
    M = (mult[x]*m) % p; O = (mult[x]*o + off[x]) % p
    parU[x] = r; mult[x] = M; off[x] = O
    return r, M, O
contra = []
def pin(v, k):
    rv, mv, ov = find3(v)
    if rv == CONST:
        if (mv + ov - k) % p != 0: contra.append(('pin', v, k))
        return
    parU[rv] = CONST; mult[rv] = 0; off[rv] = ((k - ov) * pow(mv, -1, p)) % p
def union(a, b, m, o):
    ra, ma, oa = find3(a); rb, mb, ob = find3(b)
    if ra == rb:
        if (ma - m*mb) % p != 0 or (oa - (m*ob + o)) % p != 0: contra.append(('u', a, b))
        return
    if ra == CONST:
        K = (ma + oa) % p; pin(b, ((K - o)*pow(m, -1, p)) % p); return
    if rb == CONST:
        L = (mb + ob) % p; pin(a, (m*L + o) % p); return
    inv = pow((m*mb) % p, -1, p)
    parU[rb] = ra; mult[rb] = (ma*inv) % p; off[rb] = ((oa - m*ob - o)*inv) % p
# process pins first
for v, k in extra_pins:
    pin(v, k)
# process 1- and 2-var linear rels via union-find, collect bigger ones
red = []
for rel, const in lin_rels:
    vs = list(rel)
    if len(vs) == 0:
        if const % p != 0: contra.append(('empty', const))
    elif len(vs) == 1:
        v = vs[0]; pin(v, (-const*pow(rel[v], -1, p)) % p)
    elif len(vs) == 2:
        va, vb = vs; ca, cb = rel[va], rel[vb]
        union(va, vb, (-cb*pow(ca, -1, p)) % p, (-const*pow(ca, -1, p)) % p)
    else:
        red.append((rel, const))
print(f"after UF: contradictions={len(contra)}, bigger rels deferred={len(red)}")

# reduce deferred rels through UF, then gaussian
rows = []; consts = []
for rel, const in red:
    r = defaultdict(int); cc = const % p
    for v, coef in rel.items():
        rv, m, o = find3(v)
        if rv == CONST:
            cc = (cc + coef*(m + o)) % p
        else:
            r[rv] = (r[rv] + coef*m) % p; cc = (cc + coef*o) % p
    r = {k: v % p for k, v in r.items() if v % p}
    if not r:
        if cc % p != 0: contra.append(('red-empty', cc))
    else:
        rows.append(r); consts.append(cc)
# ADD CORE CONDITIONS x_14853 ≡ x_12186, x_24908 ≡ x_16742
for a, b in [(14853, 12186), (24908, 16742)]:
    r = defaultdict(int); cc = 0
    for v, coef in [(a, 1), (b, p-1)]:
        rv, m, o = find3(v)
        if rv == CONST: cc = (cc + coef*(m+o)) % p
        else: r[rv] = (r[rv] + coef*m) % p; cc = (cc + coef*o) % p
    r = {k: v % p for k, v in r.items() if v % p}
    if not r:
        print(f"  CORE {a},{b} reduces to const {cc} -> {'OK' if cc%p==0 else 'CONTRADICTION'}")
        if cc % p != 0: contra.append(('core', a, b))
    else:
        rows.append(r); consts.append(cc)
print(f"gaussian rows: {len(rows)}, pre-gauss contradictions: {len(contra)}")

# sparse gaussian
col_rows = defaultdict(set)
for i, row in enumerate(rows):
    for v in row: col_rows[v].add(i)
alive = set(range(len(rows))); pivots = {}
gcontra = 0
order = sorted(alive, key=lambda i: len(rows[i]))
t0 = time.time()
for i in order:
    if i not in alive: continue
    row = rows[i]
    if not row:
        if consts[i] % p != 0: gcontra += 1
        alive.discard(i); continue
    pc = min(row, key=lambda v: len(col_rows[v] & alive))
    pivots[pc] = i; alive.discard(i)
    inv = pow(row[pc], -1, p)
    rows[i] = {v: (c*inv) % p for v, c in row.items()}; consts[i] = (consts[i]*inv) % p
    row = rows[i]
    for j in list(col_rows[pc]):
        if j == i or j not in alive: continue
        f = rows[j].get(pc, 0)
        if not f: continue
        for v, c in row.items():
            nv = (rows[j].get(v, 0) - f*c) % p
            if nv: rows[j][v] = nv; col_rows[v].add(j)
            elif v in rows[j]: del rows[j][v]; col_rows[v].discard(j)
        consts[j] = (consts[j] - f*consts[i]) % p
        if not rows[j] and consts[j] % p != 0: gcontra += 1; alive.discard(j)
print(f"pivots={len(pivots)}, gaussian contradictions={gcontra}, total contra={len(contra)}, time={time.time()-t0:.0f}s")
print("RESULT:", "CONSISTENT (mod-p witness with wire-lock exists in linear span)" if (gcontra==0 and len(contra)==0) else f"INCONSISTENT")

# ============ EXTRACT residue assignment and check ALL atoms mod p ============
print("\n=== extracting residue assignment ===")
val39022 = H.loadd('best_agentA_39022.json')
def rootval_free(r):
    return val39022.get(r, 0) % p
# resolve pivot roots: rows[i] with pivot pc: pc = consts[i] - sum(coef*v for v!=pc)
# iterate to fixpoint (rows should be reduced)
root_val = {}
# free roots first: all roots that are not pivots and not CONST
allroots = set()
for v in range(N):
    r,_,_ = find3(v); allroots.add(r)
allroots.discard(CONST)
pivot_cols = set(pivots.keys())
free_roots = [r for r in allroots if r not in pivot_cols]
for r in free_roots:
    root_val[r] = rootval_free(r)
# resolve pivots iteratively
remaining = dict(pivots)
for _ in range(50):
    progressed = False
    for pc, i in list(remaining.items()):
        row = rows[i]
        deps = [v for v in row if v != pc]
        if all((v in root_val or v == CONST) for v in deps):
            s = consts[i]
            for v in deps:
                vv = 1 if v == CONST else root_val[v]
                s = (s - row[v]*vv) % p
            root_val[pc] = s % p
            del remaining[pc]; progressed = True
    if not remaining or not progressed: break
print(f"resolved pivots: {len(pivots)-len(remaining)}/{len(pivots)}, unresolved={len(remaining)}")
for pc in remaining:  # any leftover -> set free
    root_val[pc] = 0
root_val[CONST] = 1

# full var residues
resid = [0]*N
for v in range(N):
    r, m, o = find3(v)
    rv = 1 if r == CONST else root_val.get(r, 0)
    resid[v] = (m*rv + o) % p

# check all atoms mod p
bad_lin = bad_quad = bad_higher = 0
bad_examples = []
import json as _json
mc = _json.load(open('/home/user/integer_solver/solve_lab/scratchpad/modp_class.json'))
quadset = set(mc['quad_atoms'])
for ai, poly in enumerate(atoms):
    val = 0
    for vs, c in poly:
        term = c % p
        for v in vs: term = (term*resid[v]) % p
        val = (val + term) % p
    if val % p != 0:
        if ai in quadset: bad_quad += 1
        else: bad_lin += 1
        if len(bad_examples) < 20: bad_examples.append((ai, A.ATOM_REPR[ai][:60]))
print(f"\nATOM CHECK MOD P: bad_linear={bad_lin}, bad_quadratic={bad_quad}")
print("sample bad atoms:")
for ai, r in bad_examples: print(f"  atom {ai} ({'quad' if ai in quadset else 'lin'}): {r}")
