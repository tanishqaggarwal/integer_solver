"""Classify atoms as gates (define a var) vs checks; build the gate DAG over variables."""
import model, ast, re, sys, pickle, os
from collections import Counter, defaultdict

VAR_RE = re.compile(r'x_(\d+)')
d = model.get()
atom_src = d['atom_src']; atom_vars = d['atom_vars']; eq_terms = d['eq_terms']
NA = len(atom_src)

# parse each atom into a polynomial dict: monomial(tuple sorted vars) -> coeff
def poly_of(node):
    """return dict monomial->coeff ; monomial = tuple of var ids (sorted, with multiplicity)"""
    if isinstance(node, ast.Constant):
        return {(): node.value} if node.value else {}
    if isinstance(node, ast.Name):
        return {(int(node.id[2:]),): 1}
    if isinstance(node, ast.UnaryOp):
        p = poly_of(node.operand)
        if isinstance(node.op, ast.USub): return {k:-v for k,v in p.items()}
        return p
    if isinstance(node, ast.BinOp):
        L = poly_of(node.left); R = poly_of(node.right)
        if isinstance(node.op, ast.Add):
            o = dict(L)
            for k,v in R.items():
                o[k] = o.get(k,0)+v
                if o[k]==0: del o[k]
            return o
        if isinstance(node.op, ast.Sub):
            o = dict(L)
            for k,v in R.items():
                o[k] = o.get(k,0)-v
                if o[k]==0: del o[k]
            return o
        if isinstance(node.op, ast.Mult):
            o = {}
            for k1,v1 in L.items():
                for k2,v2 in R.items():
                    k = tuple(sorted(k1+k2))
                    o[k] = o.get(k,0)+v1*v2
                    if o[k]==0: del o[k]
            return o
    raise ValueError(ast.dump(node))

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'polys.pkl')
if os.path.exists(CACHE):
    polys = pickle.load(open(CACHE,'rb'))
else:
    polys = [poly_of(ast.parse(s, mode='eval').body) for s in atom_src]
    pickle.dump(polys, open(CACHE,'wb'))

deg = [max((len(k) for k in p), default=0) for p in polys]
print('atom degree hist:', Counter(deg).most_common())
print('atom monomial count hist:', Counter(len(p) for p in polys).most_common(10))

# a candidate definition: variable v with monomial (v,) coeff +-1 and v not in any other monomial
def defcands(p):
    out = []
    for k,c in p.items():
        if len(k)==1 and abs(c)==1:
            v = k[0]
            if not any(v in k2 for k2 in p if k2!=k):
                out.append(v)
    return out

cands = [defcands(p) for p in polys]
print('atoms with >=1 def candidate:', sum(1 for c in cands if c))
print('def-candidate count hist:', Counter(len(c) for c in cands).most_common(10))
pickle.dump(cands, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'defcands.pkl'),'wb'))
