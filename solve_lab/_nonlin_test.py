import heal_harness as H, json, random
from collections import defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
bysel=defaultdict(list)
for rec in pins:
    _,sel,tgt,const,coef,handle=rec; bysel[sel].append((tgt,const%p,handle))
sels=set(bysel)
tree2=sorted(H.anc[11150]&sels)   # 78 bits
LOADS=[11150,25739,37758]
def setbits(active):
    for v in H.freeinp: H.val[v]=0
    for b in active:
        H.val[b]=1
        for tgt,const,handle in bysel[b]: H.val[tgt]=const
    H.forward()
    return tuple(H.val[n]%p for n in LOADS)
def coef_at(baseline_bits, testbits):
    base=setbits(baseline_bits)
    out={}
    for b in testbits:
        if b in baseline_bits:
            newset=[x for x in baseline_bits if x!=b]   # turn OFF
        else:
            newset=list(baseline_bits)+[b]              # turn ON
        v=setbits(newset)
        out[b]=tuple((v[k]-base[k])%p for k in range(3))
    return base,out
# Baseline A: all-zero
baseA=[]
bA,cA=coef_at(baseA, tree2)
respA=[b for b in tree2 if any(cA[b][k] for k in range(3))]
print("all-zero baseline: bits changing load:",respA)
# Baseline B: bit 32669 + 3 random others
random.seed(1)
extra=random.sample([b for b in tree2 if b!=32669],4)
baseB=[32669]+extra
bB,cB=coef_at(baseB, tree2)
print("baseline B active bits:",baseB,"loads:",[str(x)[:12] for x in bB])
# compare coefficient of each bit A vs B
diff=[b for b in tree2 if cA[b]!=cB[b]]
print("bits whose single-flip effect DIFFERS between baseline A and B:",len(diff),"/",len(tree2))
print("  -> if >0, loads are NONLINEAR (interaction terms) in the bits")
# show a few
for b in diff[:8]:
    print("  bit",b,"cA=",[str(x)[:10] for x in cA[b]],"cB=",[str(x)[:10] for x in cB[b]])
# Also: pair-additivity test from all-zero for pairs including 32669
print("\nPair-additivity test (all-zero baseline):")
base0=setbits([])
c_single={}
for b in tree2:
    v=setbits([b]); c_single[b]=tuple((v[k]-base0[k])%p for k in range(3))
random.seed(2)
testpairs=[(32669,b) for b in random.sample([x for x in tree2 if x!=32669],5)]
testpairs+=[tuple(random.sample(tree2,2)) for _ in range(5)]
for i,j in testpairs:
    v=setbits([i,j])
    actual=tuple((v[k]-base0[k])%p for k in range(3))
    pred=tuple((c_single[i][k]+c_single[j][k])%p for k in range(3))
    print(f"  pair({i},{j}): additive={actual==pred}  actualL1={str(actual[0])[:14]} predL1={str(pred[0])[:14]}")
