"""Direct local analysis of the 39,026 deliverable's own assignment (no re-propagation)."""
import sys, json, collections, time, math
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD=[0]*E.NV
for k,val in d.items(): vD[int(k.split('_')[1])]=int(val)
ns={'v':vD,'__builtins__':{}}
bad={i:eval(H.acodes[i],ns) for i in range(len(H.atoms))}
bad={i:r for i,r in bad.items() if r}
print("bad atoms:",sorted(bad))
ff=E.eqfails(bad); print("fails",len(ff),"score",39033-len(ff))
for a in sorted(bad):
    vs=H.avars[a]
    print(f"\na{a} nf={NF[a]} : {H.atoms[a][:160]}")
    print(f"   vars={len(vs)}  occ-counts={[(u,len(H.occ[u])) for u in vs]}")
    priv=[u for u in vs if len(H.occ[u])==1]
    print(f"   PRIVATE vars (occ==1): {priv}")
