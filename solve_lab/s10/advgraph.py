"""S11 step 87: the advice constraint graph, and a Gauss-Seidel solve of it.

§85-86 reduced the instance to thirteen 296-bit advice values, and §86 showed the
k*p part is irrelevant -- every k scores identically, so only the residues matter.
It also showed the failures MOVE rather than clear: setting x24548 correctly makes
a7930 and a41512 vanish and a21617 and a37662 appear.  That is not conservation, it
is a CYCLE: x24548's target depends on x14623 and x14623's target depends back.

So build the graph.  Each advice value is constrained by one of two shapes:

    two-sided   c*(x_i - y_i) - p*h            ->  x_i ≡ y_i (mod p),  y_i computed
    constant    w*(x_i - C)   - ...            ->  x_i ≡ C  (mod p),  C written in
                                                    the instance as a literal

The constant pins are ground truth -- nothing to solve, just read them off.  The
two-sided ones form a dependency graph among the thirteen; if it is a DAG, one
Gauss-Seidel sweep in topological order sets every advice value correctly and the
instance is solved.  If it has cycles, the cycles are the whole remaining problem
and are now small enough to name.

Usage: advgraph.py [state.json] [SWEEPS]
"""
import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'B7_39004.json'
SWEEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 8
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
FREE = [t for t in range(L.NVARS) if t not in L.definer]
ADV = [t for t in FREE if 280 <= v[t].bit_length() <= 330]
ADVS = set(ADV)
print('%s: score %d; %d advice values'
      % (src, L.NEQ - len(L.failing_eqs(av)), len(ADV)), flush=True)

CON = []          # (advice, kind, payload, check)
for t in ADV:
    for a in sorted(L.var_atoms[t]):
        if a in L.atom_out:
            continue
        pl = L.polys[a]
        if len(pl) > 4:
            continue
        lin = {m[0]: c for m, c in pl.items() if len(m) == 1}
        quad = {m: c for m, c in pl.items() if len(m) == 2}
        # two-sided gadget c*(x_t - y) - handle: the OPERAND is the variable whose
        # coefficient is exactly minus x_t's; the leftover +-1 term is p*handle.
        # Filtering the operand by size was wrong -- y is often zero at a given state.
        if t in lin and len(lin) >= 2:
            tgt = [w for w, c in lin.items() if w != t and c == -lin[t]]
            if tgt:
                CON.append((t, 'two', (tgt[0], lin[t]), a))
                continue
        # w*(x_t - C): the same variable multiplies x_t and the literal constant
        for m, c in quad.items():
            if t in m:
                w = m[0] if m[1] == t else m[1]
                cm = [(mm, cc) for mm, cc in pl.items()
                      if len(mm) == 1 and mm[0] == w]
                if cm and w != t:
                    CON.append((t, 'const', (w, (-cm[0][1]) * pow(c, -1, P) % P), a))
                    break
print('\n--- constraints found ---')
for t, kind, pay, a in CON:
    if kind == 'two':
        y, co = pay
        d = (v[t] - v[y]) % P
        print('  x%-6d ≡ x%-6d (mod p)   via a%-6d  %s'
              % (t, y, a, 'HOLDS' if d == 0 else 'FAILS'))
    else:
        w, C = pay
        d = (v[t] - C) % P
        print('  x%-6d ≡ CONSTANT  (mod p)  via a%-6d (multiplier x%d)  %s'
              % (t, a, w, 'HOLDS' if d == 0 else 'FAILS'))

_, freelist, SVS = suppfree.build(v, modp=None)


def fsupp(t):
    m = SVS[t] if t < len(SVS) else 0
    return {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}


print('\n--- dependency graph over the advice set ---')
dep = {}
for t, kind, pay, a in CON:
    if kind != 'two':
        dep.setdefault(t, set())
        continue
    y = pay[0]
    d = fsupp(y) & ADVS - {t}
    dep.setdefault(t, set()).update(d)
    print('  x%-6d target x%-6d depends on advice %s' % (t, y, sorted(d)))
order, temp, perm, cyc = [], set(), set(), []


def visit(n, stack):
    if n in perm:
        return
    if n in temp:
        cyc.append(stack[stack.index(n):] + [n])
        return
    temp.add(n)
    for m in sorted(dep.get(n, ())):
        visit(m, stack + [m])
    temp.discard(n)
    perm.add(n)
    order.append(n)


for t in sorted(dep):
    visit(t, [t])
print('  topological order: %s' % order)
if cyc:
    print('  CYCLES: %s' % cyc[:6])
else:
    print('  the advice graph is a DAG -- one sweep should solve it')

print('\n--- Gauss-Seidel sweeps ---', flush=True)
best = L.NEQ - len(L.failing_eqs(av))


def handles(c, v):
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    out = []
    for i in range(len(freelist)):
        if (m >> i) & 1:
            u = freelist[i]
            g = jacZ(u, v, [c]).get(c, 0)
            if g and g % P == 0:
                out.append((u, g))
    return out


for sw in range(SWEEPS):
    changed = 0
    for t in order:
        cs = [c for c in CON if c[0] == t]
        for _, kind, pay, a in cs:
            tgt = v[pay[0]] % P if kind == 'two' else pay[1]
            k = v[t] // P
            new = k * P + tgt % P
            if new != v[t]:
                v[t] = new
                changed += 1
        ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    for a in [x for x in range(L.NA) if x not in L.atom_out and av[x]]:
        if av[a] % P:
            continue
        for h, g in handles(a, v):
            if av[a] % g == 0:
                v[h] = v[h] - av[a] // g
                ad.fwd(v, rounds=6)
                av = L.all_atom_values(v)
                break
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('  sweep %-2d changed %-3d  score %-6d  nonzero checks %s'
          % (sw, changed, s, nz[:10]), flush=True)
    if s > best:
        best = s
        T.save(v, os.path.join(HERE, 'AG_%d.json' % s))
        print('    *** saved AG_%d.json' % s, flush=True)
    if not changed and not nz:
        print('*** FIXED POINT WITH NO RESIDUAL')
        break
print('\nbest %d' % best)
