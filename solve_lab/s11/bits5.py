"""All 256 single-bit messages, scored exactly in GF(p).

Every ON bit pins two dedicated variables to fixed residues; every OFF bit leaves its two free.
So LOW-weight messages have the most continuous freedom -- the opposite of where the search has
been looking (it has only ever perturbed the 2-bit checkpoint pattern).  A 1-bit message is the
minimum allowed by the OR gate, and there are only 256 of them.
"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
D=json.load(open(os.path.join(HERE,'data','bits1.json')))
PIN={int(k):[(a,x,int(C)) for a,x,C in v] for k,v in D['pins'].items()}
real=sorted(PIN)
tr=json.load(open(os.path.join(HERE,'data','bits_trees.json')))
tree={}
for k,v in tr.items():
    for b in v: tree[b]=k
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
def msg(S):
    """canonical state for message set S: bits on, their loads set, everything else off/zero"""
    v=list(base)
    for b in real:
        v[b]=1 if b in S else 0
        for a,x,C in PIN[b]:
            v[x]=(C%P) if b in S else 0
            Pp=L.polys[a]
            hs=[m[0] for m,c in Pp.items() if len(m)==1 and m[0]!=b]
            if hs and hs[0] in FREE: v[hs[0]]=0
    forwardp(v); return v
def fails(v): return [a for a in CHK if evalp(L.polys[a],v)]
t0=time.time(); res=[]
for b in real:
    v=msg({b})
    F=fails(v)
    res.append((len(F), b, tree.get(b,'?'), v[7715], v[34554], F[:6]))
res.sort()
print(f"scanned 256 single-bit messages ({time.time()-t0:.0f}s)")
print("  distribution:", dict(sorted(collections.Counter(r[0] for r in res).items())))
print("  best:")
for n,b,tk,U,V,F in res[:20]:
    print(f"    x{b:6d} tree {tk} U={U} V={V}  failing checks mod p = {n:3d}  {F}")
json.dump([[r[0],r[1],r[2]] for r in res], open(os.path.join(HERE,'data','bits5.json'),'w'))
