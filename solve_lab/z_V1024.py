import os,sys,json,time,itertools
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
# local cone: gates that are ancestors-side of the 27 eqs (non-free vars appearing in their var closure)
LOCAL=FAILS11+RIPPLE16
need=set()
for i in LOCAL:
    need|=H.eqvars[i]
# expand to all gate-ancestors: a var's gate def uses vids; recompute in H.order restricted
localgates=[k for k,t in enumerate(H.order) if t in need or (H.anc[t] & set())]
# Simpler: recompute ALL gates whose target is in the transitive cone of the 27 eqs.
# transitive: vars in eqs -> their gate vids -> ...
cone=set(need)
changed=True
gdef_vids={t:H.gates[H.definer[t]][2] for t in H.order}
while changed:
    changed=False
    add=set()
    for t in list(cone):
        if t in gdef_vids:
            for u in gdef_vids[t]:
                if u not in cone: add.add(u)
    if add: cone|=add; changed=True
localorder=[k for k,t in enumerate(H.order) if t in cone]
ns={'v':H.val,'__builtins__':{}}
def fwd_local():
    for k in localorder: H.val[H.order[k]]=eval(H.gcode[k],ns)
def rip_resid():
    return tuple(eval(H.eqcode[i],ns)%p for i in RIPPLE16)
print(f"local cone gates: {len(localorder)}")
t0=time.time()
results=[]
allzero=[]
for pat in itertools.product([0,1],repeat=10):
    for b,val in zip(bits10,pat): H.val[b]=val
    fwd_local()
    # close gaps mod p
    H.val[7068]=H.val[2099]%p; H.val[4432]=H.val[19964]%p
    fwd_local()  # propagate x_7068,x_4432 into ripple gates
    r=rip_resid()
    nz=sum(1 for x in r if x!=0)
    results.append((nz,pat))
    if nz==0: allzero.append(pat)
print(f"scanned 1024 in {time.time()-t0:.0f}s")
results.sort()
print("min nonzero-count patterns:")
for nz,pat in results[:12]:
    print(f"  nz={nz} pat={pat}  bits={[bits10[i] for i in range(10) if pat[i]]}")
print(f"\npatterns with ALL 16 ripple ≡0 mod p: {len(allzero)}")
for pat in allzero[:20]: print("   ",pat, [bits10[i] for i in range(10) if pat[i]])
