#!/usr/bin/env python3
"""AUDIT T3b: dormancy at CLOSED states and at the deliverable's own free-input values."""
import os,sys,json,random
H=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentH_work'))
sys.path.insert(0,H); os.chdir(H)
import ev, fast
from fast import St, chk
from chain import close_trace
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
FREE=ev.F['free0']; solo=[X for X in FREE if len(chk[X])==1]
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
def census(st,nm,base_dorm=None):
    g={'dormant':0,'p':0,'other':0}; dorm=set()
    for X in solo:
        a=chk[X][0]
        q=st.clone().set_free({X:st.fv.get(X,0)+1})
        d=q.av[a]-st.av[a]
        if d==0: g['dormant']+=1; dorm.add(X); continue
        if d%p==0: g['p']+=1
        else: g['other']+=1
    extra='' if base_dorm is None else '  newly-live vs all-zero: %d'%len(base_dorm-dorm)
    print('%-46s %s%s'%(nm,g,extra),flush=True)
    return dorm
b=census(St({}),'all-zero (uncosed)')
# closed 1-selector state, exactly hsweep's base
st0=St({}); bsel=BITS['A'][0]
base,ok,tr,frz=close_trace(st0.clone().set_free({bsel:1}),frozen=set(ALL))
print('hsweep base closed state score',base.score(),flush=True)
census(base,'CLOSED 1-sel x_%d (hsweep actual base)'%bsel,b)
# deliverable's own free values
d=json.load(open('../best/new_instance_partial_39026.json'))
full={int(k[2:]):int(v) for k,v in d.items()}
fv={u:full.get(u,0) for u in FREE}
w=St(fv); print('witness-derived state score',w.score(),flush=True)
census(w,'DELIVERABLE free-input values',b)
# random multi-selector configs
random.seed(3)
for t in range(3):
    on=random.sample(sorted(ALL),8)
    s=St({u:1 for u in on})
    census(s,'random 8-selector cfg #%d'%t,b)
