"""Buy the missing knob by breaking one cheap gate.

The continuous deficit is exactly 1.  Breaking a gate atom frees its output variable -- that is
precisely the trick the 39,026 checkpoint uses -- at a cost of the equations containing that
atom.  Many gate atoms sit in very few equations, so if one of them frees a variable that moves
the deficit rows (a25676 / a42245) independently of x14623 and x31339, the whole system closes
for a handful of failing equations.
"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
import fw
from gmp1 import evalp, solvep
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]

def forwardp_frozen(v, frozen):
    for comp in fw.ORDER:
        if len(comp)==1:
            u=comp[0]
            if u in frozen: continue
            x=solvep(L.definer[u],u,v)
            if x is not None: v[u]=x
        else:
            for _ in range(60):
                ch=False
                for u in comp:
                    if u in frozen: continue
                    x=solvep(L.definer[u],u,v)
                    if x is not None and x!=v[u]: v[u]=x; ch=True
                if not ch: break
    return v

base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
FAIL=[a for a in CHK if bd[a]]
KEY=[7930,21617,25676,42245,29539,40826,41512,33792,40562,40623]
GATES=[a for a in range(L.NA) if L.atom_out.get(a) is not None]
cost={a:len(L.atom2eq.get(a,{})) for a in GATES}
dist=collections.Counter(cost.values())
print("gate atoms by #equations:", dict(sorted(dist.items())[:10]))
cheap=[a for a in GATES if cost[a]<=6]
print(f"gate atoms in <=6 equations: {len(cheap)}")
t0=time.time(); hits=[]
for i,g in enumerate(cheap):
    t=L.atom_out[g][1]
    v=list(base); v[t]=(v[t]+1)%P
    forwardp_frozen(v,{t})
    d={a:(evalp(L.polys[a],v)-bd[a])%P for a in CHK}
    d={a:x for a,x in d.items() if x}
    if not d: continue
    touch=[a for a in KEY if a in d]
    if touch: hits.append((cost[g],g,t,touch,len(d)))
    if i%400==0: print(f"   {i}/{len(cheap)} hits={len(hits)} ({time.time()-t0:.0f}s)", flush=True)
hits.sort()
print(f"gates whose freed output touches a KEY row: {len(hits)}  ({time.time()-t0:.0f}s)")
for c,g,t,touch,nd in hits[:30]:
    print(f"   a{g} (cost {c} eqs) frees x{t}: touches {touch}  total checks moved {nd}")
json.dump([[c,g,t,touch,nd] for c,g,t,touch,nd in hits], open(os.path.join(HERE,'data','gmp26.json'),'w'))
