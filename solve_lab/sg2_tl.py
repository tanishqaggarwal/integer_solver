#!/usr/bin/env python3
"""Mod-p tangent-linear (derivative) forward evaluator.
Base point = agentA. Propagate a perturbation dval over free inputs through the gate DAG mod p,
and compute the induced derivative of every atom. Used to find TRUE mod-p sensitivity."""
import heal_harness as H
import sg2_lib as L
import ast, re, pickle
p = H.p
atoms = L.load_atoms_full()
A = {a['idx']: a for a in atoms}
freeinp = H.freeinp
order = H.order
gates = H.gates
definer = H.definer

# --- parse each gate rhs into an AST; build value+deriv evaluator ---
VAR = re.compile(r'x_(\d+)')
gate_ast = {}
for t in order:
    _, rhs, _ = gates[definer[t]]
    gate_ast[t] = ast.parse(rhs, mode='eval').body

# base values mod p
val = [0]*L.NVARS
def set_base(assignment):
    for v in freeinp:
        val[v] = assignment.get(v, 0) % p
    for t in order:
        val[t] = evast(gate_ast[t]) % p

def evast(node):
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name): return val[int(node.id[2:])]
    if isinstance(node, ast.UnaryOp): return (-evast(node.operand)) % p
    if isinstance(node, ast.BinOp):
        a = evast(node.left); b = evast(node.right)
        if isinstance(node.op, ast.Add): return (a+b) % p
        if isinstance(node.op, ast.Sub): return (a-b) % p
        if isinstance(node.op, ast.Mult): return (a*b) % p
    raise ValueError(ast.dump(node))

# derivative propagation: given dval[] for free inputs, compute dval[] for gates in order
dval = [0]*L.NVARS
def devast(node):
    """derivative of node w.r.t. current dval, values from val[]."""
    if isinstance(node, ast.Constant): return 0
    if isinstance(node, ast.Name): return dval[int(node.id[2:])]
    if isinstance(node, ast.UnaryOp): return (-devast(node.operand)) % p
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Add): return (devast(node.left)+devast(node.right)) % p
        if isinstance(node.op, ast.Sub): return (devast(node.left)-devast(node.right)) % p
        if isinstance(node.op, ast.Mult):
            a = evast(node.left); b = evast(node.right)
            da = devast(node.left); db = devast(node.right)
            return (a*db + b*da) % p
    raise ValueError

def tl_forward():
    for t in order:
        dval[t] = devast(gate_ast[t]) % p

# derivative of an atom (poly) at current val/dval
def datom(poly):
    s = 0
    for m, c in poly.items():
        # d(prod) = sum_j (prod of others) * dval[x_j]
        for j in range(len(m)):
            term = c
            for k, u in enumerate(m):
                if k == j: term = (term * dval[u]) % p
                else: term = (term * val[u]) % p
            s = (s + term) % p
    return s % p

# atom value at current val (mod p)
def vatom(poly):
    s = 0
    for m, c in poly.items():
        t = c
        for u in m: t = (t*val[u]) % p
        s = (s+t) % p
    return s % p

if __name__ == '__main__':
    vA = H.loadd('best_agentA_39022.json')
    set_base(vA)
    # confirm G1,G2 residues nonzero, everything else zero mod p
    nzmod = [a['idx'] for a in atoms if vatom(a['poly']) != 0]
    print(f"atoms nonzero mod p at agentA: {nzmod}")
