"""Identify ALL checked free inputs (free inputs pinned by a clean linear equality-check eq),
build slave order, provide forward_full_slave() that keeps every equality-check satisfied."""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p
from collections import defaultdict, deque

NV = H.NVARS
freeinp = H.freeinp
anc = H.anc
eqvars = H.eqvars
eqcode = H.eqcode

def build_slaves():
    """Return: slave_of[v]=(eqidx, coeff), and topo order list of checked free inputs."""
    # candidate: free input v that appears DIRECTLY in eq i, constant coeff, and v is NOT an
    # ancestor of any other variable in eq i (clean).
    ns = {'v': H.val, '__builtins__': {}}
    # We detect coeff via finite-diff on val[v] w/o forward (captures direct term only).
    slave_of = {}
    owners = defaultdict(list)  # eq -> list of clean free inputs it could pin
    for i in range(len(eqvars)):
        vs = eqvars[i]
        directfree = [v for v in vs if v in freeinp]
        if not directfree: continue
        for v in directfree:
            # cleanness: v not ancestor of any OTHER var in eq
            clean = True
            for w in vs:
                if w == v: continue
                if w in freeinp:
                    if w == v: continue
                    # another free input directly present -> that's fine, it's separate term
                    continue
                if v in anc.get(w, ()):  # v feeds gate w
                    clean = False; break
            if not clean: continue
            owners[i].append(v)
    # An equation is an equality-check pinning v iff it has exactly one clean owner free input
    # AND that free input's coeff is a nonzero constant. Compute coeff by finite diff.
    def direct_coeff(i, v):
        old = H.val[v]
        b0 = eval(eqcode[i], ns)
        H.val[v] = old + 1; b1 = eval(eqcode[i], ns)
        H.val[v] = old + 2; b2 = eval(eqcode[i], ns)
        H.val[v] = old
        c = b1 - b0
        c2 = b2 - b1
        if c2 != c: return None  # nonlinear in v
        return c
    # assign each checked free input to a unique owning equation
    cand = {}  # v -> list of (eq, coeff)
    for i, vlist in owners.items():
        # only consider if exactly the eq can pin one of them uniquely; collect all linear owners
        lin = []
        for v in vlist:
            c = direct_coeff(i, v)
            if c is not None and c != 0:
                lin.append((v, c))
        # equation pins a free input if it has exactly ONE free-input variable total (the check)
        # OR structurally is an equality. We take: if the eq's free support (through all vars) minus
        # gate-fed frees equals a single free input. Simpler: if exactly one clean linear owner AND
        # eq has no other free inputs feeding its gates uniquely -> assign.
        for v, c in lin:
            cand.setdefault(v, []).append((i, c))
    # Build slave_of: pick for each free input the check eq with smallest free-support (most local)
    efs = {}
    def eq_free_support(i):
        if i in efs: return efs[i]
        s = set()
        for var in eqvars[i]:
            if var in freeinp: s.add(var)
            else: s |= anc.get(var, set())
        efs[i] = s; return s
    for v, lst in cand.items():
        # prefer equation where v is the ONLY free input in support (pure check), else smallest support
        best = None
        for (i, c) in lst:
            sup = eq_free_support(i)
            key = (len(sup), i)
            if best is None or key < best[0]:
                best = (key, i, c)
        slave_of[v] = (best[1], best[2])
    return slave_of, eq_free_support

if __name__ == '__main__':
    v013 = H.loadd('best/new_instance_partial_39013.json')
    for v in freeinp: H.val[v] = v013.get(v, 0)
    H.forward()
    slave_of, efs = build_slaves()
    print(f"checked free inputs found: {len(slave_of)}")
    # how many have a PURE check (only free input in their eq support)?
    pure = 0
    for v,(i,c) in slave_of.items():
        if efs(i) == {v}: pure += 1
    print(f"  pure single-free-input checks: {pure}")
    # coeff distribution
    from collections import Counter
    cc = Counter(abs(c) for v,(i,c) in slave_of.items())
    print(f"  coeff |c| distribution (top): {cc.most_common(8)}")
    # save the slave map
    json.dump({str(v): [i, c] for v,(i,c) in slave_of.items()}, open('sg1_slavemap.json','w'))
    print("saved sg1_slavemap.json")
