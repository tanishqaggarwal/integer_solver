#!/usr/bin/env python3
"""Decompose each equation into gate 'atoms'.

Each equation's LHS, after stripping outer scalar multipliers and squares,
is a top-level Add-chain of terms.  Each term is (coef)*(atom) where a clean
atom has the form  x_target - rhs  (a gate definition).  Terms that are not
of that shape mark the equation as a 'constraint' equation.

Outputs solve_lab/atoms/gates.jsonl and solve_lab/atoms/summary.json
"""
import ast, re, json, time, sys
from collections import defaultdict, Counter

EQ_PATH = __file__.rsplit('/', 1)[0] + '/../EQUATIONS.txt'
OUT = __file__.rsplit('/', 1)[0] + '/atoms'

def unparse(node):
    return ast.unparse(node)

def const_val(node):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
        return -node.operand.value
    return None

def strip_outer(node):
    """Reduce s*E -> E, E*E -> E (square). Return core node whose ==0 is the constraint."""
    while True:
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            a, b = node.left, node.right
            ca, cb = const_val(a), const_val(b)
            if ca is not None and cb is not None:
                # constant*constant: whole thing constant; leave as-is
                return node
            if ca is not None:
                node = b; continue
            if cb is not None:
                node = a; continue
            # both non-constant
            sa, sb = unparse(a), unparse(b)
            if sa == sb:
                node = a; continue  # square E*E -> E
            return node  # genuine product P*Q
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            node = node.operand; continue
        return node

def flatten_add(node):
    """Flatten a top-level Add chain into terms. Do NOT descend into Sub (that is
    the atom's own `x_t - rhs`) or into Mult. Returns list of (sign, term_node)."""
    terms = []
    def rec(n, sign):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            rec(n.left, sign); rec(n.right, sign)
        else:
            terms.append((sign, n))
    rec(node, 1)
    return terms

def strip_coef(node):
    """term = coef*atom or atom. Return (coef_or_None, atom_node)."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        ca, cb = const_val(node.left), const_val(node.right)
        if ca is not None and cb is None:
            return ca, node.right
        if cb is not None and ca is None:
            return cb, node.left
    return None, node

def classify_atom(node):
    """If node is  x_t - rhs  return (target_id, rhs_node). Else None."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub):
        left = node.left
        if isinstance(left, ast.Name) and left.id.startswith('x_'):
            return int(left.id[2:]), node.right
    return None

VAR_RE = re.compile(r'x_(\d+)')

def rhs_info(rhs_node):
    s = unparse(rhs_node)
    vids = tuple(sorted(set(int(m) for m in VAR_RE.findall(s))))
    return s, vids

def main():
    t0 = time.time()
    gates = []               # (target, rhs_str, rhs_vids, eq_idx)
    constraint_eqs = []      # eq_idx with non-atom terms
    pure_wiring = 0
    term_shape = Counter()
    n = 0
    with open(EQ_PATH) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            n += 1
            lhs = line.rsplit('=', 1)[0]
            tree = ast.parse(lhs, mode='eval').body
            core = strip_outer(tree)
            terms = flatten_add(core)
            all_clean = True
            eq_gates = []
            for sign, term in terms:
                coef, atom = strip_coef(term)
                res = classify_atom(atom)
                if res is None:
                    all_clean = False
                    term_shape[type(atom).__name__ + ('/' + type(atom.op).__name__ if isinstance(atom, ast.BinOp) else '')] += 1
                    continue
                tgt, rhs = res
                s, vids = rhs_info(rhs)
                eq_gates.append((tgt, s, vids, i))
            if all_clean:
                pure_wiring += 1
            else:
                constraint_eqs.append(i)
            gates.extend(eq_gates)
    print(f"parsed {n} equations in {time.time()-t0:.1f}s")
    print(f"pure-wiring equations: {pure_wiring}")
    print(f"constraint equations (have non-atom terms): {len(constraint_eqs)}")
    print(f"total gate atoms collected: {len(gates)}")
    print(f"non-atom term shapes: {term_shape.most_common(15)}")

    # dedup gates by (target, rhs_str)
    seen = set()
    uniq = []
    for g in gates:
        key = (g[0], g[1])
        if key not in seen:
            seen.add(key)
            uniq.append(g)
    print(f"distinct (target,rhs) gates: {len(uniq)}")

    # targets: how many distinct, multiplicity
    by_target = defaultdict(list)
    for tgt, s, vids, i in uniq:
        by_target[tgt].append(s)
    multi = {t: rs for t, rs in by_target.items() if len(rs) > 1}
    print(f"distinct targets: {len(by_target)}   targets with >1 distinct rhs: {len(multi)}")

    # inputs = vars that appear but are never a target
    all_vars = set()
    for tgt, s, vids, i in uniq:
        all_vars.add(tgt); all_vars.update(vids)
    inputs = sorted(all_vars - set(by_target.keys()))
    print(f"total vars referenced in gates: {len(all_vars)}   input (never-target) vars: {len(inputs)}")
    print(f"sample inputs: {inputs[:30]}")

    import os
    os.makedirs(OUT, exist_ok=True)
    with open(OUT + '/gates.jsonl', 'w') as g:
        for tgt, s, vids, i in uniq:
            g.write(json.dumps({"t": tgt, "rhs": s, "vids": vids, "eq": i}) + "\n")
    with open(OUT + '/summary.json', 'w') as g:
        json.dump({
            "n_eqs": n, "pure_wiring": pure_wiring,
            "n_constraint_eqs": len(constraint_eqs),
            "n_gates_uniq": len(uniq), "n_targets": len(by_target),
            "n_multi_target": len(multi), "n_inputs": len(inputs),
            "inputs": inputs, "constraint_eqs": constraint_eqs,
        }, g)
    with open(OUT + '/multi_targets.json', 'w') as g:
        json.dump({str(t): rs for t, rs in multi.items()}, g, indent=1)
    print("wrote atoms/gates.jsonl, summary.json, multi_targets.json")

if __name__ == '__main__':
    main()
