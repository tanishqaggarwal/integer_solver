#!/usr/bin/env python3
"""AUDIT T27b -- M's engine is exact AT THE WITNESS.  A scorer can be exact at its calibration
point and wrong elsewhere, and the enumeration's value is its verdict on the other 4,095.
Spot-check non-witness subsets of H12 against checker.py."""
import os,sys,json,random,time
LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0,os.path.join(LAB,'agentM_work')); sys.path.insert(0,os.path.join(LAB,'agentE_work'))
sys.set_int_max_str_digits(20_000_000)
os.chdir(os.path.join(LAB,'agentM_work'))
import ieng
PF=json.load(open('pfamily.json'))
H12=sorted({v['h'] for v in PF['incident_7'].values()})
sys.path.insert(0,LAB)
import checker as CK
codes,_=CK.load_equations()
D4=[642,28730,29854,31864]
rnd=random.Random(3)
tests=[tuple(D4),(28730,),(642,29854),tuple(sorted(D4+[H12[0]])),tuple(sorted(D4+[h for h in H12 if h not in D4][:2]))]
seen=set(tests)
while len(tests)<9:
    k=rnd.randint(1,6); S=tuple(sorted(rnd.sample(H12,k)))
    if S not in seen: seen.add(S); tests.append(S)
print('H12 = %s\n'%H12)
print('%-42s %-9s %-9s %s'%('subset','M engine','checker','agree'))
ok=bad=0
for S in tests:
    try: r=ieng.tune(list(S),want=True)
    except Exception as e:
        print('%-42s ENGINE ERROR %s'%(str(S)[:42],e)); continue
    if not r.get('ok'):
        print('%-42s engine: not ok (%s)'%(str(S)[:42],r.get('why'))); continue
    if r.get('changes') is None:
        v=ieng.V_UNC
    else:
        v,_=ieng.apply_delta(ieng.V_UNC,r['changes'],r['pin'])
    sc=len(codes)-len(CK.evaluate_all(codes,v))
    agree = (sc==r['score'])
    ok+=agree; bad+=(not agree)
    print('%-42s %-9d %-9d %s'%(str(S)[:42],r['score'],sc,'YES' if agree else '*** NO ***'))
print('\nagree %d / disagree %d'%(ok,bad))
