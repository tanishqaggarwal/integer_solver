import os,sys,json,time
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict,deque
p=H.p
pins=json.load(open('pinrec.json'))
selectors=set(r[1] for r in pins)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
var_eqs=defaultdict(list)
for i,vs in enumerate(H.eqvars):
    for v in vs: var_eqs[v].append(i)
# seed: local verifier eqs
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
# closure over CONTINUOUS free inputs only (bits fixed)
E=set(FAILS11+RIPPLE16)
K=set()
for _ in range(30):
    # free inputs (non-bit) feeding E
    newK=set()
    for i in E:
        for v in H.eqvars[i]:
            if v in H.freeinp and v not in selectors: newK.add(v)
    # eqs downstream of K
    newE=set()
    for v in newK:
        newE.update(var_eqs[v])
        for k in desc_of[v]: newE.update(var_eqs[H.order[k]])
    if newK<=K and newE<=E:
        break
    K|=newK; E|=newE
E=sorted(E); K=sorted(K)
print(f"closure: {len(E)} eqs, {len(K)} continuous free inputs")
# how many bits appear in E's ancestry
Eanc=set()
for i in E:
    for v in H.eqvars[i]: Eanc|=H.anc.get(v,{v})
Ebits=sorted(Eanc&selectors)
print(f"bits in closure ancestry: {len(Ebits)}: {Ebits}")
json.dump({'E':E,'K':K,'bits':Ebits}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/lc.json','w'))
