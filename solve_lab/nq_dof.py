import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
F=set(H.fails())
print(f"fails: {len(F)}")

# equation footprint of each free input (structural, quadrant-independent)
eq_frees=[]
for i in range(len(H.lines)):
    s=set()
    for v in H.eqvars[i]:
        s|=H.anc.get(v, {v} if v in H.freeinp else set())
    eq_frees.append(s & H.freeinp)
foot=defaultdict(set)
for i,fs in enumerate(eq_frees):
    for f in fs: foot[f].add(i)

# free inputs feeding the 15 fails
brk_frees=set()
for i in F: brk_frees|=eq_frees[i]
private=[f for f in brk_frees if not (foot[f]-F)]
leak={f:(len(foot[f]-F)) for f in brk_frees if (foot[f]-F)}
print(f"free inputs feeding 15 fails: {len(brk_frees)}")
print(f"PRIVATE (feed only fails): {len(private)}: {sorted(private)}")
# low-leak knobs (feed few outside)
lowleak=sorted(leak.items(), key=lambda kv:kv[1])[:25]
print(f"low-leak knobs (free: #outside-eqs): {lowleak}")
# key target vars and their free ancestors
for name in [27177,4306,23754,33168,9106,31731,35619,9629,2239,26874,6947]:
    a=H.anc.get(name,set())&H.freeinp
    isfree = name in H.freeinp
    print(f"x_{name}: {'FREE' if isfree else 'gate'}, #free-anc={len(a)}, anc(private)={sorted(a&set(private))[:10]}")
