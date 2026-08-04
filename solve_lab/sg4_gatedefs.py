import heal_harness as H
import json
p=H.p
# gate definitions
gates={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.setdefault(d['t'],[]).append((d['rhs'],tuple(d['vids'])))
for t in [2099,19964,642,28730,24908,29322,3558,4432,7068,14853,12186,16742,6947,26874,23754]:
    if t in gates:
        for rhs,vids in gates[t]:
            print(f"x_{t} = {rhs[:150]}   vids={vids}")
    else:
        print(f"x_{t}: FREE INPUT (no gate)")
