import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

# gate defs: target -> vids  (for computed vars)
gate_out = H.gate_out
freeinp = H.freeinp

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()

def decompose(ai):
    """Return (primary_free, coeff, other_vars) for a small gap atom.
    primary = free var appearing in a single-var term. Prefer the one whose
    other-term partner is a gate output (the matched computed var)."""
    poly = A.ATOMS[ai]
    singles = [(vl[0], c) for vl, c in poly if len(vl) == 1]
    free_singles = [(v, c) for v, c in singles if v in freeinp]
    return free_singles, singles, poly

def zero_atom(ai):
    """Set a free primary var so atom -> 0. Returns primary var or None."""
    poly = A.ATOMS[ai]
    val = A.eval_atom(ai, H.val)
    # find free var with single-var term
    for vl, c in poly:
        if len(vl) == 1 and vl[0] in freeinp:
            f = vl[0]
            # atom = c*f + rest ; set f so total 0
            if val % c == 0:
                H.val[f] -= val // c
                return f
    return None

# Run cascade with proper fixer, collect all gaps + the fail set at each step
seen = {}       # ai -> repr
gap_primary = {}  # ai -> primary free var
for step in range(200):
    # small nonzero atoms
    sn = [ai for ai in range(A.NATOM) if len(A.ATOM_VARS[ai]) <= 6 and A.eval_atom(ai, H.val) != 0]
    if not sn:
        print(f"step {step}: ALL small gaps zero. fails={len(H.fails())}")
        break
    for ai in sn:
        seen.setdefault(ai, A.ATOM_REPR[ai])
    fixed_any = False
    for ai in sn:
        f = zero_atom(ai)
        if f is not None:
            gap_primary[ai] = f
            fixed_any = True
    H.forward()
    if not fixed_any:
        print(f"step {step}: stuck, small gaps not linearly fixable: {sn}")
        break
else:
    print("hit 200 steps")

print(f"\nTotal distinct gaps: {len(seen)}, fails now={len(H.fails())}")
# Build dependency: for each gap, primary free var f; the 'computed' vars are the
# gate-output vars in the atom. Edge j->i if primary_j in anc[computed_var_i].
primaries = set(gap_primary.values())
print(f"primaries ({len(primaries)}): {sorted(primaries)}")

edges = []
gapinfo = {}
for ai, r in seen.items():
    f = gap_primary.get(ai)
    computed_vars = [v for v in A.ATOM_VARS[ai] if v in gate_out]
    # which primaries does this gap depend on (via computed vars' ancestors)?
    deps = set()
    for cv in computed_vars:
        deps |= (H.anc[cv] & primaries)
    deps.discard(f)
    gapinfo[ai] = (f, computed_vars, deps)
    for d in deps:
        edges.append((d, f, ai))
    print(f"  gap {ai} [{r[:55]}] primary=x_{f} computed_gates={computed_vars} depends_on_primaries={sorted(deps)}")

# cycle check among primaries
import collections
adj = collections.defaultdict(set)
for d, f, ai in edges:
    adj[d].add(f)
# detect cycle via DFS
WHITE, GREY, BLACK = 0, 1, 2
color = {}
cyc = [False]
def dfs(u):
    color[u] = GREY
    for w in adj[u]:
        if color.get(w, WHITE) == GREY:
            cyc[0] = True
        elif color.get(w, WHITE) == WHITE:
            dfs(w)
    color[u] = BLACK
for u in list(primaries):
    if color.get(u, WHITE) == WHITE:
        dfs(u)
print(f"\nCYCLE among gap primaries: {cyc[0]}")
print(f"edges (dep -> primary): {[(d,f) for d,f,ai in edges]}")
