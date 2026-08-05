import json
import heal_harness as H
p=H.p
# check x_21279 definition and x_20434
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
print(f"x_21279 def: {gdef.get(21279)}")
print(f"x_20434: {'FREE' if 20434 in H.freeinp else 'gate'}, def={gdef.get(20434)}")
print(f"x_9062: {'FREE' if 9062 in H.freeinp else 'gate'}")
# load sy config, check current
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
print(f"\ncurrent: x_9062={V[9062]}, x_20434={V[20434]}, x_21279={V[21279]}")
F0=set(H.fails())
print(f"baseline fails: {len(F0)}")
# try x_20434=0
if 20434 in H.freeinp:
    V[20434]=0
    H.forward()
    F=set(H.fails())
    print(f"\nwith x_20434=0: x_21279={V[21279]}, fails={len(F)}")
    print(f"  fixed: {sorted(F0-F)}")
    print(f"  broken: {sorted(F-F0)}")
