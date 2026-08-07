"""Slot-structure oracle on the CORRECTED engine.

Mechanism: turning ON a leaf inside subtree S saturates S and splits S's class into
S's two slot supports.  So we probe many bases, each saturating a different part of
the tree, and take the maximal common refinement of the observed leaf partitions.

Signature at each base = delta of that base's OWN bad-atom set (the natural residual
cluster there), which avoids the private-pin blow-up of a full-support signature.

Two leaves stay together iff their signatures agree in EVERY base where NEITHER is ON
(the correction from refine2.py: a leaf that is ON in a base must not be compared there).
"""
import sys, os, json, time, pickle, collections, random
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2, mcore2 as M

P = M.P
B = M.bools()
NB = len(B)
vd = M.load_vec()
seed_d = E2.seed_of(vd)
DELIV_ON = [f for f in B if vd[f] != 0]
neutral = dict(seed_d)
for f in DELIV_ON:
    neutral.pop(f, None)

blocks = json.load(open('blocks8.json'))['blocks']
tree = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/tree96.json'))
gsup = {k: set(v['gsup']) for k, v in tree.items()}

# ---- build the base list: one base per tree stage (saturating that stage), plus extras
rnd = random.Random(7)
bases = [('neutral', {}), ('deliv', {f: 1 for f in DELIV_ON})]
seen = set()
for k, g in sorted(gsup.items(), key=lambda kv: -len(kv[1])):
    g2 = sorted(g & set(B))
    if len(g2) < 2:
        continue
    for f in g2[:2]:                      # two different saturating leaves per stage
        key = (k, f)
        if f in seen:
            continue
        seen.add(f)
        bases.append((f'stage{k}_leaf{f}', {f: 1}))
# a few random pairs to break ties
for i in range(20):
    a, b = rnd.sample(B, 2)
    bases.append((f'rnd{i}', {a: 1, b: 1}))

print(f'{len(bases)} bases, {NB} leaves', flush=True)

sigs = {}          # basename -> {leaf: sig or 'ON'}
t0 = time.time()
for n, (name, extra) in enumerate(bases):
    s = dict(neutral); s.update(extra)
    v0 = E2.forward(s)
    bad0 = E2.badatoms(v0)
    coords = sorted(bad0)
    if not coords:
        continue
    d = {}
    for f in B:
        if v0[f] != 0:
            d[f] = 'ON'; continue
        b1, _ = fast2.resid_delta(v0, bad0, {f: 1})
        d[f] = tuple((b1.get(a, 0) - bad0.get(a, 0)) % P for a in coords)
    sigs[name] = d
    if n % 10 == 0:
        print(f'  [{n}/{len(bases)}] {time.time()-t0:.0f}s  ({len(coords)} coords)', flush=True)

pickle.dump(sigs, open('orefine_sigs.pkl', 'wb'))

# ---- maximal common refinement, ON-corrected ----
def together(f, g):
    for name, d in sigs.items():
        sf, sg = d[f], d[g]
        if sf == 'ON' or sg == 'ON':
            continue                     # not comparable at this base
        if sf != sg:
            return False
    return True

parent = {f: f for f in B}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

for i in range(NB):
    for j in range(i + 1, NB):
        f, g = B[i], B[j]
        if find(f) != find(g) and together(f, g):
            union(f, g)

grp = collections.defaultdict(list)
for f in B:
    grp[find(f)].append(f)
newblocks = sorted((sorted(v) for v in grp.values()), key=lambda s: (-len(s), s[0]))
print(f'\nCOMMON REFINEMENT: {len(newblocks)} blocks, sizes {[len(x) for x in newblocks]}', flush=True)
print('previous (E engine, 297 cfgs): 8 blocks, sizes [178, 41, 21, 6, 3, 3, 3, 1]')

json.dump({'blocks': newblocks}, open('blocks_corrected.json', 'w'))

# ---- how many tree96 stages does this refinement now resolve? ----
bset = [set(x) for x in newblocks]
resolved = {}
for k, g in gsup.items():
    gl = g & set(B)
    if len(gl) < 2:
        continue
    parts = [sorted(gl & bb) for bb in bset if gl & bb]
    if len(parts) >= 2:
        resolved[k] = [len(p) for p in parts]
print(f'\ntree96 stages SPLIT by the refinement: {len(resolved)}')
for k, v in sorted(resolved.items(), key=lambda kv: -sum(kv[1]))[:30]:
    print(f'  stage {k}: |gsup∩leaves|={sum(v)} -> parts {v}')
json.dump({'split_stages': {k: v for k, v in resolved.items()}},
          open('oracle_stage_splits.json', 'w'))
print('done', flush=True)
