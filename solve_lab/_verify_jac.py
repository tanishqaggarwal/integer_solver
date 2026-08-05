import heal_harness as H, jac_lib as J
from collections import defaultdict
import random
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F0=H.fails()
print("baseline fails",len(F0))
vd=J.build_duals()
# sanity: dual values match forward values
mism=sum(1 for t in H.order if vd[t].v != H.val[t]%p)
print("value mismatches (should be 0):",mism)
# pick a few equations, compare dual grad to central finite-diff over a few free vars
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
inv2=pow(2,p-2,p)
random.seed(1)
test_eqs=[2554,6816,8124,12231,29125]+random.sample(range(len(H.eqcode)),5)
ok=True
for i in test_eqs:
    rv,grad=J.eq_jac_row(i,vd)
    # test against FD on the free vars in grad (a few)
    cols=list(grad.keys())[:6]
    for c in cols:
        v=J.freelist[c]
        o=base[v]
        setfree(v,o+1); rp=eval(H.eqcode[i],ns)%p
        setfree(v,o-1); rm=eval(H.eqcode[i],ns)%p
        setfree(v,o)
        fd=((rp-rm)*inv2)%p
        if fd!=grad[c]%p:
            print(f"MISMATCH eq{i} col{c}(x_{v}): dual={grad[c]%p} fd={fd}"); ok=False
    # also residual matches
    rr=eval(H.eqcode[i],ns)%p
    if rr!=rv: print(f"resid mismatch eq{i}: dual={rv} real={rr}"); ok=False
print("all grads match FD:",ok)
# report the G1/G2 residual mod p from the 11 fails
for i in F0:
    rv,grad=J.eq_jac_row(i,vd)
    print(f"eq {i}: resid mod p nonzero? {rv!=0}  #grad nonzero cols {len(grad)}")
