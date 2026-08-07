import sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
v=forward([0]*L.NVARS)
def expr(x):
    if x not in outs: return None
    c,rest=LN[x]
    return {m:(-cc)//c for cc,m in rest if cc}
def as_add(x):
    d=expr(x)
    if d and len(d)==2 and all(len(m)==1 and d[m]==1 for m in d):
        return tuple(sorted(m[0] for m in d))
    return None
def as_mul(x):
    d=expr(x)
    if d and len(d)==1:
        m=next(iter(d))
        if len(m)==2 and d[m]==1: return tuple(sorted(m))
    return None
def as_copy(x):
    d=expr(x)
    if d and len(d)==1:
        m=next(iter(d))
        if len(m)==1 and d[m]==1: return m[0]
    return None
def as_or(x):
    d=expr(x)
    if not (d and len(d)==2): return None
    it=list(d.items())
    pos=[m for m,c in it if c==1 and len(m)==1]; neg=[m for m,c in it if c==-1 and len(m)==1]
    if len(pos)!=1 or len(neg)!=1: return None
    A=pos[0][0]; B=neg[0][0]
    a=as_add(A); b=as_mul(B)
    if a and b and a==b: return a
    return None
def leaves(root):
    out=[]; st=[root]; seen=set()
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        o=as_or(x)
        if o: st.extend(o); continue
        cp=as_copy(x)
        if cp is not None: st.append(cp); continue
        out.append(x)
    return out
allv=[]
for r in [7715,34554]:
    lv=leaves(r); allv+=lv
    print('x_%d OR-tree leaves: %d'%(r,len(lv)))
    for x in sorted(lv):
        print('   x_%-6d val=%-3s def=%s'%(x,v[x],(L.atom_src[outs[x]][:150] if x in outs else 'FREE')))
