import sys, os, json, re
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import atomlib as A
p = A.p

# 1) find atom(s) that are exactly "x_26064 - p" or similar single-var pins to p
print("=== atoms mentioning x_26064 ===")
for ai in A.VAR_ATOMS[26064]:
    print(f"  atom {ai}: {A.ATOM_REPR[ai][:120]!r}  vars={sorted(A.ATOM_VARS[ai])} in {len(A.ATOM_EQS[ai])} eqs {A.ATOM_EQS[ai][:6]}")

# 2) look for the equation that forces x_26064.  Find equations that contain x_26064
lines = [L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR = re.compile(r'x_(\d+)')
eq_with = [i for i,L in enumerate(lines) if 26064 in [int(m) for m in VAR.findall(L)]]
print(f"\n=== equations containing x_26064: {eq_with} ===")
for i in eq_with[:3]:
    print(f"  eq {i} (len {len(lines[i])}): {lines[i][:300]}")

# 3) The wire: which vars are 'copies' forming the union-find class of x_26064.
# Build union-find from identity gates: target = single var (copy) or target = -var etc.
gate_def = {}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line)
        gate_def.setdefault(d['t'], []).append((d['rhs'], tuple(d['vids'])))

# find gates that are pure copies: rhs is 'x_k' (single var, coeff1) => t == k
copy_edges = []
for t, defs in gate_def.items():
    for rhs, vids in defs:
        s = rhs.strip()
        m = re.fullmatch(r'x_(\d+)', s)
        if m:
            copy_edges.append((t, int(m.group(1))))
# union-find
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb
for a,b in copy_edges:
    union(a,b)
root = find(26064)
wire = [x for x in list(parent) if find(x)==root]
print(f"\n=== wire (copy-class of x_26064): {len(wire)} members ===")
print(f" contains x_28599? {28599 in wire}; x_17499? {17499 in wire}; x_26874? {26874 in wire}; x_13859? {13859 in wire}; x_15616? {15616 in wire}")
print(f" sample: {sorted(wire)[:30]}")
