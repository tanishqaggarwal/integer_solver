#!/usr/bin/env python3
"""Independent parser: EQUATIONS.txt -> structural model (atoms, gates)."""
import ast, re, sys, pickle, time, collections

EQ = '/home/user/integer_solver/EQUATIONS.txt'

def to_py(line):
    return line.rsplit('=',1)[0].strip()

# ---- structural extraction ----
def flatten_add(node):
    """Return list of (sign, node) summands of a +/- chain."""
    out=[]
    def rec(n, sgn):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            rec(n.left, sgn); rec(n.right, sgn)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            rec(n.left, sgn); rec(n.right, -sgn)
        else:
            out.append((sgn,n))
    rec(node,1)
    return out

def const_of(n):
    if isinstance(n, ast.Constant) and isinstance(n.value,int): return n.value
    if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
        c=const_of(n.left if False else n.operand)
        return None if c is None else -c
    return None

def src(n):
    return ast.unparse(n)

def main():
    t0=time.time()
    eqs=[]
    with open(EQ) as f:
        lines=[l.strip() for l in f if l.strip()]
    print(f"lines {len(lines)} read {time.time()-t0:.1f}s", file=sys.stderr)
    atom_id={}   # canonical source -> id
    atoms=[]     # id -> source
    eq_terms=[]  # per eq: (kind, outer_scalar, [(coef, atom_id)])
    for li,line in enumerate(lines):
        e=to_py(line)
        tree=ast.parse(e, mode='eval').body
        # peel structure
        kind='lin'
        outer=1
        node=tree
        # square?  A*A with identical unparse
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
            L,R=node.left,node.right
            cl,cr=const_of(L),const_of(R)
            if cl is not None:
                outer*=cl; node=R
            elif cr is not None:
                outer*=cr; node=L
            else:
                if src(L)==src(R):
                    kind='sq'; node=L
        # possibly nested again: C * ((-1)*(X))
        changed=True
        while changed:
            changed=False
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
                L,R=node.left,node.right
                cl,cr=const_of(L),const_of(R)
                if cl is not None: outer*=cl; node=R; changed=True
                elif cr is not None: outer*=cr; node=L; changed=True
                elif kind=='lin' and src(L)==src(R):
                    kind='sq'; node=L; changed=True
            elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
                outer=-outer; node=node.operand; changed=True
        # now node is either the inner sum, or sum of  C1*(A)+C2*(A)  form
        terms=[]
        for sgn,s in flatten_add(node):
            c=sgn; sub=s
            # peel leading constant multipliers
            while isinstance(sub, ast.BinOp) and isinstance(sub.op, ast.Mult):
                cl=const_of(sub.left); cr=const_of(sub.right)
                if cl is not None: c*=cl; sub=sub.right
                elif cr is not None: c*=cr; sub=sub.left
                else: break
            while isinstance(sub, ast.UnaryOp) and isinstance(sub.op, ast.USub):
                c=-c; sub=sub.operand
            key=src(sub)
            aid=atom_id.get(key)
            if aid is None:
                aid=len(atoms); atom_id[key]=aid; atoms.append(key)
            terms.append((c,aid))
        eq_terms.append((kind,outer,terms))
        if li%5000==0: print(f"  {li} {time.time()-t0:.0f}s", file=sys.stderr)
    print(f"parsed {time.time()-t0:.1f}s  distinct atoms {len(atoms)}", file=sys.stderr)
    with open('model.pkl','wb') as f:
        pickle.dump({'atoms':atoms,'eq_terms':eq_terms}, f)
    # stats
    kinds=collections.Counter(k for k,_,_ in eq_terms)
    print("kinds", kinds)
    print("terms per eq histogram", collections.Counter(len(t) for _,_,t in eq_terms).most_common(20))

main()
