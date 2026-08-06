import heal_harness as H
import json, pickle
from collections import defaultdict
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
d=H.loadd('best/new_instance_partial_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def incr(w,nv):
    H.val[w]=nv
    for k in desc_of[w]: H.val[H.order[k]]=eval(H.gcode[k],{'v':H.val,'__builtins__':{}})
def av(i):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
def free_anc_atom(i):
    s=set()
    for m,c in ATOMS[i]:
        for v in m:
            if v in H.freeinp: s.add(v)
            else: s|=H.anc.get(v,set())
    return s
# rigorously re-test linearity of the 4 nonlinear cert atoms (ALL free ancestors)
print("=== rigorous linearity of the 4 crux atoms ===")
for a in [30614,18081,37034,37229]:
    fa=sorted(free_anc_atom(a)); nonl=[]
    for f in fa:
        b0=av(a); incr(f,base[f]+1); b1=av(a); incr(f,base[f]+2); b2=av(a); incr(f,base[f])
        if (b1-b0)%p!=(b2-b1)%p: nonl.append(f)
    print(f"  atom {a} ({reprs[a][:45]}): #free_anc={len(fa)}, nonlinear-in free: {nonl[:8]}")
    # which var makes it nonlinear
    vs=sorted(set(x for m,c in ATOMS[a] for x in m))
    for v in vs:
        isgate = v not in H.freeinp
        print(f"      x_{v} free={v in H.freeinp} #anc={len(H.anc.get(v,{v}))}")
