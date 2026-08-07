#!/usr/bin/env python3
import ast, sys, pickle, time, collections
EQ='/home/user/integer_solver/EQUATIONS.txt'
def const_of(n):
    if isinstance(n, ast.Constant) and isinstance(n.value,int): return n.value
    if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
        c=const_of(n.operand); return None if c is None else -c
    return None
def main():
    t0=time.time()
    lines=[l.strip() for l in open(EQ) if l.strip()]
    atom_id={}; atoms=[]; eq_terms=[]
    out=[]
    def expand(n, c):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            expand(n.left,c); expand(n.right,c); return
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            expand(n.left,c); expand(n.right,-c); return
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            expand(n.operand,-c); return
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
            cl=const_of(n.left); cr=const_of(n.right)
            if cl is not None: expand(n.right, c*cl); return
            if cr is not None: expand(n.left, c*cr); return
        k=const_of(n)
        if k is not None:
            out.append((c*k, None)); return
        out.append((c, ast.unparse(n)))
    for li,line in enumerate(lines):
        e=line.rsplit('=',1)[0].strip()
        node=ast.parse(e, mode='eval').body
        kind='lin'; outer=1
        changed=True
        while changed:
            changed=False
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
                cl=const_of(node.left); cr=const_of(node.right)
                if cl is not None: outer*=cl; node=node.right; changed=True
                elif cr is not None: outer*=cr; node=node.left; changed=True
                elif kind=='lin' and ast.unparse(node.left)==ast.unparse(node.right):
                    kind='sq'; node=node.left; changed=True
            elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
                outer=-outer; node=node.operand; changed=True
        out.clear(); expand(node,1)
        terms=[]
        for c,key in out:
            if key is None:
                terms.append((c,-1)); continue
            aid=atom_id.get(key)
            if aid is None:
                aid=len(atoms); atom_id[key]=aid; atoms.append(key)
            terms.append((c,aid))
        eq_terms.append((kind,outer,terms))
        if li%10000==0: print(f"  {li} {time.time()-t0:.0f}s natoms={len(atoms)}", file=sys.stderr)
    print(f"parsed {time.time()-t0:.1f}s distinct atoms {len(atoms)}", file=sys.stderr)
    pickle.dump({'atoms':atoms,'eq_terms':eq_terms}, open('model2.pkl','wb'))
    print("kinds", collections.Counter(k for k,_,_ in eq_terms))
    print("nterms hist", sorted(collections.Counter(len(t) for _,_,t in eq_terms).items())[:25])
    # atom usage
    use=collections.Counter()
    for _,_,t in eq_terms:
        for c,a in t:
            if a>=0: use[a]+=1
    print("atom usage hist", sorted(collections.Counter(use.values()).items())[:20])
main()
