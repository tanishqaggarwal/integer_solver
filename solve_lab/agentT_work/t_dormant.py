#!/usr/bin/env python3
"""AUDIT T3: agent H's handle-carrier census excluded 722 of 1,865 one-check-atom free
inputs as 'dormant'.  handles.py measures dormancy by a unit bump at the ALL-ZERO state
St({}); hsweep.py then prices the survivors at a DIFFERENT state (the closed 1-selector
state, selector = BITS['A'][0] = x_47).  Re-measure dormancy at several configurations."""
import os,sys,json,time
H=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentH_work'))
sys.path.insert(0,H); os.chdir(H)
import ev, fast
from fast import St, chk
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
FREE=ev.F['free0']
solo=[X for X in FREE if len(chk[X])==1]
print('one-check-atom free inputs (H\'s "solo" class):',len(solo))
BITS=json.load(open('bits.json'))
cfgs=[('all-zero (handles.py base)',{}),
      ("1-sel x_%d (hsweep base)"%BITS['A'][0],{BITS['A'][0]:1}),
      ('1-sel x_%d (other A bit)'%BITS['A'][7],{BITS['A'][7]:1}),
      ('1-sel x_%d (a B bit)'%BITS['B'][0],{BITS['B'][0]:1}),
      ('witness 2-sel {x_24601,x_2081}',{24601:1,2081:1})]
base_dormant=None
for nm,sel in cfgs:
    st=St({})
    if sel: st=st.clone().set_free({k:v for k,v in sel.items()})
    gran={'dormant':0,'p':0,'other':0}; dorm=set()
    for X in solo:
        a=chk[X][0]
        g=st.clone().set_free({X:st.fv.get(X,0)+1})
        d=g.av[a]-st.av[a]
        if d==0: gran['dormant']+=1; dorm.add(X); continue
        if d%p==0: gran['p']+=1
        else: gran['other']+=1
    if base_dormant is None: base_dormant=dorm
    newly_live=len(base_dormant-dorm)
    print('%-38s  %s   dormant-at-all-zero that are LIVE here: %d'%(nm,gran,newly_live),flush=True)
