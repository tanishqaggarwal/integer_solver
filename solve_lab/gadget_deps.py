import json
import heal_harness as H
p=H.p
pinrec=json.load(open('pinrec.json'))
bittargets={}  # target -> selector bit
for i,sel,tgt,const,coef,handle in pinrec: bittargets[tgt]=sel
sels=set(s for _,s,_,_,_,_ in pinrec)
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
for n in [7181,2239,15963,3177,27177,4306,9106,31731,17925,27019,3349,34310,8731,9118,6418,31861,14865,12553]:
    fr='FREE' if n in H.freeinp else 'gate'
    df=gdef.get(n,('',()))[0][:60]
    pinbit = bittargets.get(n)
    print(f"x_{n}: {fr} pin_by=x_{pinbit} def={df}")
# free ancestors of gadget-critical gates, and which are pin targets (loaded by bits)
print("\n=== controlling bits of gadget quantities ===")
for n in [27177,4306,2239,7181]:
    anc=H.anc.get(n,set())&H.freeinp
    ancbits={a:bittargets.get(a) for a in anc}
    loaded=[(a,bittargets[a]) for a in anc if a in bittargets]
    plainfree=[a for a in anc if a not in bittargets]
    print(f"x_{n}: {len(anc)} free-anc")
    print(f"   loaded-by-bits: {loaded}")
    print(f"   plain-free (not pinned): {plainfree}")
