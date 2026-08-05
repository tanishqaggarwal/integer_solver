import json
import heal_harness as H
p=H.p
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
# trace x_15298, x_9062, x_20434, x_21279 control
def show(n,depth=0,seen=None):
    if seen is None: seen=set()
    if n in seen or depth>4: return
    seen.add(n)
    fr='FREE' if n in H.freeinp else 'gate'
    df=gdef.get(n,('',()))
    print(f"{'  '*depth}x_{n} [{fr}]: {df[0][:60]}")
    if n not in H.freeinp:
        for c in df[1]:
            show(c,depth+1,seen)
for root in [15298,9062,20434,21279]:
    print(f"\n===== {root} =====")
    show(root)
# is x_2081 required for x_15298? check pins that x_2081 controls and whether x_15298 depends on x_2081
pinrec=json.load(open('pinrec.json'))
x2081_loads=[(tgt,const) for i,sel,tgt,const,coef,handle in pinrec if sel==2081]
x4287_loads=[(tgt,const) for i,sel,tgt,const,coef,handle in pinrec if sel==4287]
print(f"\nx_2081 loads: {[(t) for t,c in x2081_loads]}")
print(f"x_4287 loads: {[(t) for t,c in x4287_loads]}")
