import sys, os, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
p=H.p
core=json.load(open('core_relevant_bits.json'))['core_relevant']
vA=H.loadd('best_agentA_39022.json')
def reset():
    for v in H.freeinp: H.val[v]=vA.get(v,0)
reset(); H.forward()
g1_0=(H.val[7068]-H.val[2099])%p; g2_0=(H.val[4432]-H.val[19964])%p
base=set(H.fails())
print(f'baseline: G1={str(g1_0)[:18]} G2={str(g2_0)[:18]} fails={len(base)}')
print()
print('flip each core bit ON, measure dG1,dG2 mod p and collateral fails:')
live=0; inert=0
for b in core:
    reset(); H.val[b]=1; H.forward()
    g1=(H.val[7068]-H.val[2099])%p; g2=(H.val[4432]-H.val[19964])%p
    F=set(H.fails())
    dq = (g1!=g1_0) or (g2!=g2_0)
    if dq: live+=1
    else: inert+=1
    tag='LIVE' if dq else 'inert(modp)'
    # only print live ones or a few
    if dq or b==core[0]:
        print(f'  bit {b}: dG1={"Y" if g1!=g1_0 else "n"} dG2={"Y" if g2!=g2_0 else "n"} fails={len(F)} (new {len(F-base)}) [{tag}]')
print()
print(f'summary: LIVE bits (move obstruction mod p) = {live}, inert = {inert}, of {len(core)}')
