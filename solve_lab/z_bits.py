import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import Counter,defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
print(f"total pins: {len(pins)}")
sels=[r[1] for r in pins]
uniq_sels=sorted(set(sels))
print(f"unique selectors: {len(uniq_sels)}")
# how many selectors are free inputs?
free_sels=[s for s in uniq_sels if s in H.freeinp]
print(f"selectors that are FREE inputs: {len(free_sels)}")
gate_sels=[s for s in uniq_sels if s not in H.freeinp]
print(f"selectors that are GATE outputs: {len(gate_sels)}")
# the 10 "local verifier" bits from the mission
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
print(f"\n10 local-verifier bits status:")
for b in bits10:
    print(f"  x_{b}: {'FREE' if b in H.freeinp else 'gate'} , is_selector={b in uniq_sels}")
# how many pins does each of the 10 bits control?
selcount=Counter(sels)
print("\npins controlled by each of 10 bits:", {b:selcount.get(b,0) for b in bits10})
# value of each selector at agentA (are they 0/1?)
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
print("\nselector values at agentA (nonzero ones):")
nz=[(s,H.val[s]%p) for s in uniq_sels if H.val[s]%p!=0]
print(f"  {len(nz)} of {len(uniq_sels)} selectors are nonzero mod p")
print("  sample values:", [(s,(v if v<1000 else 'big')) for s,v in nz[:15]])
# how many are exactly 0 or 1?
v01=Counter('0' if H.val[s]%p==0 else '1' if H.val[s]%p==1 else 'other' for s in uniq_sels)
print("  selector value distribution:",dict(v01))
