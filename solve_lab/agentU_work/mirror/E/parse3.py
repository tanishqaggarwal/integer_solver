#!/usr/bin/env python3
"""Gate-level parse mirroring the emitter's left-nested sum structure."""
import ast, sys, pickle, time, collections
EQ='/home/user/integer_solver/EQUATIONS.txt'
def const_of(n):
    if isinstance(n, ast.Constant) and isinstance(n.value,int): return n.value
    if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
        c=const_of(n.operand); return None if c is None else -c
    return None
def peel(n, c):
    """peel constant multipliers / unary minus"""
    while True:
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            c=-c; n=n.operand; continue
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
            cl=const_of(n.left); cr=const_of(n.right)
            if cl is not None: c*=cl; n=n.right; continue
            if cr is not None: c*=cr; n=n.left; continue
        return n,c
def main():
    t0=time.time()
    lines=[l.strip() for l in open(EQ) if l.strip()]
    atom_id={}; atoms=[]; eq_terms=[]; sq_flag=[]
    out=[]
    def spine(node, coef):
        parts=[]
        cur=node
        while isinstance(cur, ast.BinOp) and isinstance(cur.op, ast.Add):
            parts.append(cur.right); cur=cur.left
        parts.append(cur); parts.reverse()
        for p in parts:
            sub,c=peel(p, coef)
            if isinstance(sub, ast.BinOp) and isinstance(sub.op, ast.Add):
                spine(sub, c)
            else:
                out.append((c, sub))
    nsq=0
    for li,line in enumerate(lines):
        e=line.rsplit('=',1)[0].strip()
        node=ast.parse(e, mode='eval').body
        outer=1; issq=False
        while True:
            node,outer=peel(node,outer)
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult) \
               and ast.unparse(node.left)==ast.unparse(node.right):
                issq=True; node=node.left; continue
            break
        if issq: nsq+=1
        out.clear(); spine(node,1)
        terms=[]
        for c,sub in out:
            k=const_of(sub)
            if k is not None: terms.append((c*k,-1)); continue
            key=ast.unparse(sub)
            aid=atom_id.get(key)
            if aid is None:
                aid=len(atoms); atom_id[key]=aid; atoms.append(key)
            terms.append((c,aid))
        eq_terms.append((issq,outer,terms))
        if li%10000==0: print(f"  {li} {time.time()-t0:.0f}s natoms={len(atoms)}", file=sys.stderr)
    print(f"parsed {time.time()-t0:.1f}s atoms={len(atoms)} sq={nsq}", file=sys.stderr)
    pickle.dump({'atoms':atoms,'eq_terms':eq_terms}, open('model3.pkl','wb'))
    import re
    pat=collections.Counter(re.sub(r'x_\d+','V',a) for a in atoms)
    print("distinct shapes", len(pat))
    for k,v in pat.most_common(30): print(v, repr(k[:150]))
main()
