import json
import heal_harness as H
from collections import defaultdict
p=H.p
pinrec=json.load(open('pinrec.json'))
bysel=defaultdict(list)
for i,sel,tgt,const,coef,handle in pinrec: bysel[sel].append((tgt,const%p,handle))
bits=sorted(bysel)
# baseline: fc_partial
d=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
S,T=35389,6671
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
baseS,baseT=V[S]%p,V[T]%p
print(f"baseline S={baseS}, T={baseT}")
# trace S's definition
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
print(f"S=x_35389 def: {gdef.get(35389)}")
print(f"T=x_6671 def: {gdef.get(6671)}")
# activate a bit properly: set bit=1, targets=const
def activate(b,on=True):
    saved={b:V[b]}; V[b]=1 if on else 0
    for tgt,const,handle in bysel[b]:
        saved[tgt]=V[tgt]
        if on: V[tgt]=const
    knobs=[b]+[t for t,_,_ in bysel[b]]
    fwd_from(knobs)
    return saved
def restore(saved):
    for v,x in saved.items(): V[v]=x
    fwd_from(list(saved))
# measure marginal S,T for each core bit
core_bits=[]
for b in bits:
    sv=activate(b); dS=(V[S]-baseS*0)%p; dT=(V[T])%p
    dSm=(V[S]%p); dTm=(V[T]%p)
    restore(sv)
    if dSm!=baseS or dTm!=baseT: core_bits.append((b,dSm,dTm))
print(f"\ncore-affecting bits: {len(core_bits)}")
# test additivity: activate first 2 core bits together
if len(core_bits)>=2:
    b1,b2=core_bits[0][0],core_bits[1][0]
    s1=activate(b1); S1,T1=V[S]%p,V[T]%p; restore(s1)
    s2=activate(b2); S2,T2=V[S]%p,V[T]%p; restore(s2)
    sa=activate(b1); sb=activate(b2); S12,T12=V[S]%p,V[T]%p; restore(sb); restore(sa)
    # additive prediction: S12 == S1+S2-baseS
    pred=(S1+S2-baseS)%p
    print(f"bit {b1}: S={S1}")
    print(f"bit {b2}: S={S2}")
    print(f"both: S={S12}")
    print(f"additive pred (S1+S2-base): {pred}  MATCH={S12==pred}")
    print(f"OR-like (same as single): {S12==S1}")
