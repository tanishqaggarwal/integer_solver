import heal_harness as H
import json
p=H.p
# ALL gate defs (there may be multiple per target)
allg=defaultdict=__import__('collections').defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); allg[d['t']].append((d['rhs'],tuple(d['vids'])))
for t in [36065,31033,20434,31822,12553,16738,6788,9062,5814,14865]:
    print(f"x_{t}: {len(allg[t])} defs")
    for rhs,vids in allg[t]:
        print(f"    = {rhs[:90]}  vids={vids}")
# which def does the harness use? check H.order / definer
# H doesn't expose definer directly, but gcode is compiled from definer. Let me find via H.anc
print("\nharness anc[36065] (free ancestors):", sorted(a for a in H.anc.get(36065,set()))[:20])
print("harness anc[31033]:", sorted(H.anc.get(31033,set())))
