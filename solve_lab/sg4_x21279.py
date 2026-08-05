import heal_harness as H
import json
p=H.p
gates={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.setdefault(d['t'],[]).append((d['rhs'],tuple(d['vids'])))
def gdef(t): 
    return gates.get(t,[('FREE',())])
# structural gate chain for the key vars
for t in [21279,4287,9062,20434,36760,19892,25297,19964,20492,2099,37158,8731,9118,28730,9413]:
    for rhs,vids in gdef(t):
        print(f"x_{t} = {rhs[:80]}   vids={vids}")
print("\n=== values at agentA baseline (core solved) ===")
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
for t in [21279,4287,9062,20434,36760,19892,25297,19964,20492,2099,37158,8731,9118,4432,7068,28730,9413,19892]:
    vv=V[t]
    print(f"x_{t} = {vv if abs(vv)<10**12 else str(vv)[:8]+'...('+str(len(str(abs(vv))))+'d)'}  (mod p={vv%p if abs(vv)<10**90 else (vv%p)})")
# gaps
g2 = V[4432]-V[20492]
g1 = V[7068]-V[37158]
print("\ng2 = x_4432-x_20492 =", g2)
print("g1 = x_7068-x_37158 =", g1)
print("G2 gap (x_4432-x_19964) =", V[4432]-V[19964])
print("G1 gap (x_7068-x_2099) =", V[7068]-V[2099])
