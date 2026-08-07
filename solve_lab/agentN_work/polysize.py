"""Size the witness region as an integer POLYNOMIAL system in the 7 zero-collateral knobs.

Everything I have measured so far linearises: the model rows are finite differences of the knobs, so
only affine solutions are ever findable.  The atoms carry genuine products (`x_26874 * x_6947` and
the like), so the true response surface is polynomial.  Here the products are carried rather than
probed away.

Set Y_j = base_j + t_j for the 7 knobs and hold every other free input at its witness value.  Each
DAG variable is then a polynomial in t_1..t_7 over Z, and each region row is the linear core of its
equation, evaluated on those polynomials.  This measures the system BEFORE attempting to solve it:
per-variable total degree by exact propagation through the expression trees, the number of variables
that actually depend on the knobs, and the resulting monomial-count bound.
"""
import ast, json, os, sys
from collections import defaultdict
import ev, model
import optN
from optN import make, build, WIT, POOL, fr, FREE, FR0
from sqaudit import square_base

HERE = os.path.dirname(os.path.abspath(__file__))
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']
eq_terms = d['eq_terms']
definer = ev.F['definer']
order = ev.F['order']

st = make(WIT)
b0 = build(st)
Rl, KNOBS = b0['R'], b0['knobs']
print('witness region |R| = %d, zero-collateral knobs = %d: %s' % (len(Rl), len(KNOBS), KNOBS))

KN = set(KNOBS)


def bits(x):
    o = []
    while x:
        q = x & -x
        o.append(q.bit_length() - 1)
        x ^= q
    return o


# which variables depend on at least one knob?  fr.csup[atom] is a free-input bitmask; build the
# analogous thing for variables by propagating through the definition DAG.
dep = bytearray(38748)
for Y in KNOBS:
    dep[Y] = 1
deg = [0] * 38748
for Y in KNOBS:
    deg[Y] = 1


def expr_degree(node, degof):
    if isinstance(node, ast.Constant):
        return 0
    if isinstance(node, ast.Name):
        return degof(int(node.id[2:]))
    if isinstance(node, ast.UnaryOp):
        return expr_degree(node.operand, degof)
    if isinstance(node, ast.BinOp):
        L = expr_degree(node.left, degof)
        R = expr_degree(node.right, degof)
        if isinstance(node.op, (ast.Add, ast.Sub)):
            return max(L, R)
        if isinstance(node.op, ast.Mult):
            return L + R
        if isinstance(node.op, ast.Pow):
            try:
                e = ast.literal_eval(node.right)
            except Exception:
                e = 2
            return L * e
        raise ValueError(ast.dump(node.op))
    raise ValueError(type(node).__name__)


TREE = {}


def tree(a):
    t = TREE.get(a)
    if t is None:
        t = ast.parse(atom_src[a], mode='eval').body
        TREE[a] = t
    return t


print('\npropagating dependence and exact total degree through the definition DAG ...', flush=True)
ndep = len(KNOBS)
for v in order:
    a = definer[v]
    if a < 0:
        continue
    t = tree(a)
    # definition is x_v = rest, where atom = x_v - rest ; the degree of x_v is the degree of rest
    if not (isinstance(t, ast.BinOp) and isinstance(t.op, ast.Sub) and isinstance(t.left, ast.Name)):
        continue
    rhs = t.right
    vs = [u for u in atom_vars[a] if u != v]
    if not any(dep[u] for u in vs):
        continue
    dep[v] = 1
    ndep += 1
    deg[v] = expr_degree(rhs, lambda u: deg[u])
print('variables depending on at least one of the 7 knobs: %d' % ndep, flush=True)
if ndep:
    dd = [deg[v] for v in range(38748) if dep[v]]
    print('   their total degrees: min %d, max %d, mean %.1f' % (min(dd), max(dd), sum(dd) / len(dd)),
          flush=True)

print('\nper-row degree of the witness region (linear core of each equation):', flush=True)
rows = []
maxdeg = 0
for e in Rl:
    m, sq, tl = eq_terms[e]
    dmax = 0
    natoms = 0
    for c, a in tl:
        if not any(dep[u] for u in atom_vars[a]):
            continue
        natoms += 1
        sb = square_base(a)
        node = tree(a)
        if sb is not None:
            node = ast.parse(sb, mode='eval').body    # strip to the linear core
        dmax = max(dmax, expr_degree(node, lambda u: deg[u]))
    rows.append((e, dmax, natoms))
    maxdeg = max(maxdeg, dmax)
    print('   eq %-6d  degree %-6d  knob-dependent atoms %d' % (e, dmax, natoms), flush=True)

import math
print('\n=== system size ===', flush=True)
print('unknowns: 7 (the zero-collateral knobs), over Z')
print('equations: %d' % len(Rl))
print('max total degree: %d' % maxdeg)
nm = math.comb(7 + maxdeg, 7) if maxdeg < 2000 else None
print('dense monomial count in 7 vars at degree %d: C(%d,7) = %s'
      % (maxdeg, 7 + maxdeg, ('%.3e' % nm) if nm else 'astronomically large'))
print('Bezout bound (product of degrees) for the %d rows: %d^%d = 10^%.1f'
      % (len(Rl), maxdeg, len(Rl), len(Rl) * math.log10(maxdeg) if maxdeg > 1 else 0))
json.dump(dict(knobs=KNOBS, R=Rl, rows=[[e, dg, na] for e, dg, na in rows],
               maxdeg=maxdeg, ndep=ndep),
          open(os.path.join(HERE, 'runs', 'polysize.json'), 'w'), indent=1)
