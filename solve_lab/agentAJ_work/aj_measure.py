#!/usr/bin/env python3
"""AJ: measure the structure of EQUATIONS.txt for a QUBO cost model.

Single pass, no materialisation of any QUBO.  Per equation we record:
  nvars_distinct, total multiplicative degree, max |integer literal| bitlength,
  number of monomials in the fully expanded form (exact, via a dict-based
  expansion that we abort if it exceeds a cap), and a shape signature.
"""
import ast, re, sys, json, collections

EQ = '/home/user/integer_solver/EQUATIONS.txt'
CAP = 200000  # abort expansion beyond this many monomials


def degree_and_lits(node):
    """Return (degree, maxlit_bits, varset). Degree = total poly degree."""
    if isinstance(node, ast.Expression):
        return degree_and_lits(node.body)
    if isinstance(node, ast.Constant):
        v = abs(int(node.value))
        return 0, v.bit_length(), frozenset()
    if isinstance(node, ast.Name):
        return 1, 0, frozenset([node.id])
    if isinstance(node, ast.UnaryOp):
        return degree_and_lits(node.operand)
    if isinstance(node, ast.BinOp):
        dl, ll, vl = degree_and_lits(node.left)
        dr, lr, vr = degree_and_lits(node.right)
        if isinstance(node.op, ast.Mult):
            return dl + dr, max(ll, lr), vl | vr
        if isinstance(node.op, ast.Pow):
            return dl * int(node.right.value), max(ll, lr), vl | vr
        return max(dl, dr), max(ll, lr), vl | vr
    raise TypeError(ast.dump(node)[:80])


# --- exact expansion into monomials (dict: sorted var tuple -> coeff) ---
def expand(node):
    if isinstance(node, ast.Expression):
        return expand(node.body)
    if isinstance(node, ast.Constant):
        c = int(node.value)
        return {(): c} if c else {}
    if isinstance(node, ast.Name):
        return {(node.id,): 1}
    if isinstance(node, ast.UnaryOp):
        e = expand(node.operand)
        if isinstance(node.op, ast.USub):
            return {k: -v for k, v in e.items()}
        return e
    if isinstance(node, ast.BinOp):
        L = expand(node.left)
        R = expand(node.right)
        out = {}
        if isinstance(node.op, ast.Mult):
            if len(L) * len(R) > CAP:
                raise OverflowError
            for a, ca in L.items():
                for b, cb in R.items():
                    k = tuple(sorted(a + b))
                    out[k] = out.get(k, 0) + ca * cb
        elif isinstance(node.op, (ast.Add, ast.Sub)):
            s = -1 if isinstance(node.op, ast.Sub) else 1
            out = dict(L)
            for b, cb in R.items():
                out[b] = out.get(b, 0) + s * cb
        else:
            raise TypeError('op')
        return {k: v for k, v in out.items() if v}
    raise TypeError('node')


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
    rows = []
    degh = collections.Counter()
    lith = collections.Counter()
    nvh = collections.Counter()
    monh = collections.Counter()
    # per-equation max monomial coefficient bitlength (after expansion)
    coefh = collections.Counter()
    aborted = 0
    with open(EQ) as f:
        for i, line in enumerate(f):
            if i >= limit:
                break
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0]
            t = ast.parse(lhs, mode='eval')
            d, lb, vs = degree_and_lits(t)
            degh[d] += 1
            lith[lb] += 1
            nvh[len(vs)] += 1
            try:
                mons = expand(t)
                nm = len(mons)
                mc = max((abs(c).bit_length() for c in mons.values()), default=0)
                # degree of expanded form
                monh[nm] += 1
                coefh[mc] += 1
                rows.append((i, d, len(vs), lb, nm, mc))
            except OverflowError:
                aborted += 1
                rows.append((i, d, len(vs), lb, -1, -1))
    out = {
        'n_equations': len(rows),
        'aborted_expansions': aborted,
        'degree_hist': dict(degh),
        'litbits_hist': dict(lith),
        'nvars_hist': dict(nvh),
        'monomials_hist': dict(sorted(monh.items())),
        'expanded_coefbits_hist': dict(sorted(coefh.items())),
    }
    json.dump(out, open('/home/user/integer_solver/solve_lab/agentAJ_work/measure_summary.json', 'w'), indent=1)
    with open('/home/user/integer_solver/solve_lab/agentAJ_work/measure_rows.tsv', 'w') as g:
        g.write('idx\tdeg\tnvars\tlitbits\tnmon\tcoefbits\n')
        for r in rows:
            g.write('\t'.join(map(str, r)) + '\n')
    print(json.dumps({k: v for k, v in out.items() if 'hist' not in k}, indent=1))
    print('degree_hist', dict(sorted(degh.items())))
    print('nvars_hist(top)', dict(sorted(nvh.items())[:15]))
    print('litbits_hist(top)', dict(sorted(lith.items())[:8]), '... max', max(lith))
    print('coefbits max', max(coefh) if coefh else None)
    print('monomials max', max(monh) if monh else None)


main()
