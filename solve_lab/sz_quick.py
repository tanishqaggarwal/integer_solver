import heal_harness as H, sz_engine as E
import time
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
t=time.time(); F=H.fails(); print(f"fails()={ (time.time()-t)*1000:.0f}ms  n={len(F)}")
# 39022 coset sweeps: x_17325 (steps of 7376877p in x_642->G1), x_9413 (steps of p in x_28730->G2)
base=len(F); best=base
for k in range(-5,6):
    H.val[17325]=k; H.forward(); n=len(H.fails())
    if n!=base: print("x_17325=",k,"->",n)
    best=min(best,n)
H.val[17325]=0
for k in range(-5,6):
    H.val[9413]=k; H.forward(); n=len(H.fails())
    if n!=base: print("x_9413=",k,"->",n)
    best=min(best,n)
H.val[9413]=0
# joint small grid
for a in range(-2,3):
  for b in range(-2,3):
    H.val[17325]=a; H.val[9413]=b; H.forward(); n=len(H.fails())
    if n<base: print("joint",a,b,"->",n)
    best=min(best,n)
H.val[17325]=0;H.val[9413]=0;H.forward()
print("39022 coset-sweep best:",best,"(baseline",base,")")
