import sys, os, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
p=H.p
def cone(v, seen=None):
    if seen is None: seen=set()
    if v in seen: return seen
    seen.add(v); k=H.definer.get(v)
    if k is None: return seen
    for iv in H.gates[k][2]: cone(iv,seen)
    return seen
treebits=sorted(v for v in cone(15298) if v in H.freeinp)
vA=H.loadd('best_agentA_39022.json')
def reset():
    for v in H.freeinp: H.val[v]=vA.get(v,0)
reset(); H.forward()
g1_0=(H.val[7068]-H.val[2099])%p; g2_0=(H.val[4432]-H.val[19964])%p
base=set(H.fails())
on=[b for b in treebits if vA.get(b,0)==1]; off=[b for b in treebits if vA.get(b,0)==0]
print(f'tree bits: {len(treebits)} total, {len(on)} ON, {len(off)} OFF in agentA')
live=0; sample=on[:8]+off[:8]
for b in sample:
    reset(); H.val[b]=1-vA.get(b,0); H.forward()
    g1=(H.val[7068]-H.val[2099])%p; g2=(H.val[4432]-H.val[19964])%p
    F=set(H.fails())
    dq=(g1!=g1_0) or (g2!=g2_0)
    live+=dq
    print(f'  flip bit {b} ({vA.get(b,0)}->{1-vA.get(b,0)}): dQ={"LIVE" if dq else "inert"} fails={len(F)}')
print(f'\n{live}/{len(sample)} sampled selector flips move the obstruction mod p')
