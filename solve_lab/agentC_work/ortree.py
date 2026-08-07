import sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from lib2 import *
v=forward([0]*L.NVARS)
# classify gate: returns ('copy',u) ('add',u,w) ('mul',u,w) ('sub',u,w) ('or',u,w) etc
def gate(x):
    a=outs.get(x)
    if a is None: return ('free',)
    Pp=L.polys[a]; c,rest=LN[x]
    # rest is list of (coef, monomial) with x removed => x = -(1/c)*sum
    terms=[( -cc//c if cc% c==0 else None, m) for cc,m in rest]
    return ('gate',a,rest,c)
def expr(x):
    """return normalized dict monomial->coef for x's definition in terms of its immediate inputs"""
    a=outs.get(x)
    if a is None: return None
    c,rest=LN[x]
    return {m:(-cc)//c for cc,m in rest if cc}, a
def is_or(x):
    e=expr(x)
    if not e: return None
    d,a=e
    # x = u + w - u*w
    if len(d)!=3: return None
    lin=[m for m in d if len(m)==1]; quad=[m for m in d if len(m)==2]
    if len(lin)!=2 or len(quad)!=1: return None
    u,w=lin[0][0],lin[1][0]
    if d[lin[0]]!=1 or d[lin[1]]!=1: return None
    if set(quad[0])!={u,w} or d[quad[0]]!=-1: return None
    return (u,w)
def is_copy(x):
    e=expr(x)
    if not e: return None
    d,a=e
    if len(d)==1:
        m=next(iter(d))
        if len(m)==1 and d[m]==1: return m[0]
    return None
def leaves(root):
    out=[]; st=[root]; seen=set()
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        o=is_or(x)
        if o: st.extend(o); continue
        cp=is_copy(x)
        if cp is not None: st.append(cp); continue
        out.append(x)
    return out
for r in [7715,34554]:
    lv=leaves(r)
    print('x_%d OR-tree leaves: %d'%(r,len(lv)))
    for x in lv:
        e=expr(x)
        print('   x_%-6d val=%-3s def=%s'%(x,v[x],(L.atom_src[outs[x]][:120] if x in outs else 'FREE')))
