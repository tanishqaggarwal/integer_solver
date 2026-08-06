import heal_harness as H, sz_engine as E
import json,re,time
from collections import defaultdict
p=H.p
RIP=E.RIP
E.classify()
r7,r4=E.setup()
F=set(H.fails()); assert F==set(RIP), (len(F),sorted(F)[:5])
print("target state OK: 16 ripple fail")

# ripple cone: all vars ripple eqs transitively depend on (via definers)
t0=time.time()
cone=set()
stack=[]
for e in RIP:
    for w in H.eqvars[e]:
        if w not in cone: cone.add(w); stack.append(w)
while stack:
    w=stack.pop()
    gi=H.definer.get(w)
    if gi is None: continue
    for u in H.gates[gi][2]:
        if u not in cone: cone.add(u); stack.append(u)
print(f"ripple cone size {len(cone)} ({time.time()-t0:.1f}s)")

# parse bilinear product gates x_c = x_a * x_b  (a!=b), both free-input, both currently 0, t in cone
VAR=re.compile(r'x_(\d+)')
prod=[]
for line in open('atoms/gates.jsonl'):
    d=json.loads(line); rhs=d['rhs']; vids=d['vids']; t=d['t']
    if t not in cone: continue
    m=re.fullmatch(r'\s*x_(\d+)\s*\*\s*x_(\d+)\s*',rhs)
    if not m: continue
    a,b=int(m.group(1)),int(m.group(2))
    if a==b: continue
    prod.append((t,a,b))
print(f"bilinear product gates in cone: {len(prod)}")
bothzero=[(t,a,b) for (t,a,b) in prod if H.val[a]==0 and H.val[b]==0]
print(f"  both-factor-zero: {len(bothzero)}")
freeboth=[(t,a,b) for (t,a,b) in bothzero if a in H.freeinp and b in H.freeinp]
print(f"  both factors FREE inputs: {len(freeboth)}")
# fanout of each free var: how many gates reference it (as vids) -> 'clean' if a,b feed few
fan=defaultdict(int)
for gi,(tt,rhs,vids) in enumerate(H.gates):
    for u in vids: fan[u]+=1
clean=[(t,a,b) for (t,a,b) in freeboth if fan[a]<=1 and fan[b]<=1]
print(f"  clean (each factor feeds <=1 gate): {len(clean)}")
semiclean=[(t,a,b) for (t,a,b) in freeboth if fan[a]<=2 and fan[b]<=2]
print(f"  semiclean (<=2): {len(semiclean)}")
json.dump([[t,a,b] for t,a,b in freeboth], open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/sz_freeboth.json','w'))
# also non-both-zero products in cone with at least one free zero factor & other nonzero (single-activate)
print("\nsample clean slacks (t,a,b):",clean[:8])
