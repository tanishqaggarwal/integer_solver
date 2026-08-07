"""Translate P's merge law (N1 = E*A^2 - B^2, N2 = A(i3+i6) - B(i2-i5), A = i1-i2, B = i4-i3)
   into MY parse, then test degeneracy (A = B = 0) directly against my own measurements.

   Atom indices are NOT comparable across agent directories, so nothing here is imported from
   another agent's numbering -- the law is re-found structurally in my own atom table."""
import sys, json, collections, pickle, re
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H, engine as E, fast
P = C.P


def deftext(var, depth=0, maxd=3):
    """Recursively unfold the atom that defines x_var."""
    a = (E.definer[var] if var<len(E.definer) else None)
    if a is None:
        return "x_%d(FREE)" % var
    t = H.atoms[a]
    if depth >= maxd:
        return "a%d[%s]" % (a, t[:60])
    def sub(m):
        v = int(m.group(1))
        if v == var:
            return "x_%d" % v
        return "(" + deftext(v, depth + 1, maxd) + ")"
    return re.sub(r'x_(\d+)', sub, t)


print("=== structure of the two cluster rows in MY parse ===", flush=True)
for a in (20215, 28647):
    print("a%d = %s" % (a, H.atoms[a]))
    for m in re.finditer(r'x_(\d+)', H.atoms[a]):
        v = int(m.group(1))
        d = (E.definer[v] if v<len(E.definer) else None)
        print("    x_%d <- %s" % (v, ("a%d: %s" % (d, H.atoms[d][:110])) if d is not None else "FREE"))
    print()

print("=== unfolded 2 levels ===", flush=True)
for a in (20215, 28647):
    print("a%d = %s\n" % (a, deftext_top := re.sub(r'x_(\d+)',
          lambda m: "(" + deftext(int(m.group(1)), 1, 2) + ")", H.atoms[a])[:1200]), flush=True)

# ---- direct degeneracy test -------------------------------------------------
# Degeneracy in P's sense makes the block's residual atoms vanish IDENTICALLY, i.e. the atom
# is exactly 0 AND stops responding to every knob.  Test both over my whole BFS image.
seen = pickle.load(open('bfs_image.pkl', 'rb'))
print("=== exact (not mod p) values of the cluster rows over the converged BFS image ===", flush=True)
zero_exact = collections.Counter()
nresp = collections.Counter()
for k, a in sorted(seen.items(), key=lambda kv: len(kv[1]))[:48]:
    seed = dict(C.BASE); seed.update(a)
    v0 = E.forward(seed)
    ns = {'v': v0, '__builtins__': {}}
    for r in C.ROWS:
        val = eval(H.acodes[r], ns)
        if val == 0:
            zero_exact[r] += 1
print("configs (of %d) where each row is EXACTLY 0:" % len(seen),
      {("a%d" % r): zero_exact.get(r, 0) for r in C.ROWS}, flush=True)
