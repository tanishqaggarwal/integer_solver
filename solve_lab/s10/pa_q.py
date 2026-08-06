"""Can a37887 = Q^2 be held at 0 while a22231 is used as a free compensator?

Q = a22231 + (linear form in the compensator-family variables).  If ANY variable in
that linear form is movable at zero collateral cost, Q = 0 can be met without spending
a22231's freedom, and the placement becomes |E|=12, |S|=8, c=2 -> 6 failing.
"""
import os, sys, collections, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
print('a37887 src:'); print(L.atom_src[37887])
print('\nvars of a37887:',sorted(L.avars[37887]))
QV=[18253,37720,30108,34600,23754,7945,23642,23822,37254,15324,35619,9629,4432,19964,28730]
print('\nvalues + who pins them:')
for u in QV:
    d=L.definer.get(u)
    cons=[a for a in L.var_atoms[u]]
    print(f' x_{u:<6} val={v[u]} definer=a{d}')
    for a in cons:
        print(f'      a{a} fp={len(L.atom2eq[a])} gate={L.atom_out.get(a)} : {L.atom_src[a][:95]}')
# census lookup: harmless variables (fail stays 7 with unchanged support)
cen={}
for f in ('pa_census_a.json','pa_census_b.json'):
    cen.update(json.load(open(os.path.join(HERE,f))))
BASE={22229,22230,35758,35759,35760,35761,35762}
harmless=[int(t) for t,r in cen.items() if r['fail']==7 and set(r['nz'])==BASE]
print(f'\nharmless perturbations (support unchanged, failing stays 7): {len(harmless)}')
# which harmless perturbations move any Q variable?
import random
random.seed(3)
D=random.randrange(1,2**60)
hit=[]
t0=time.time()
QS=set(QV)
for i,t in enumerate(harmless):
    w=list(v)
    ch,_=L.ripple(w,{t:v[t]+D},maxsteps=60000,block=set())
    mv=[u for u in QS if u in ch]
    if mv: hit.append((t,mv))
    if i%400==0: print(' ',i,len(hit),f'{time.time()-t0:.0f}s',flush=True)
print('harmless perturbations that MOVE a Q-variable:',len(hit))
for t,mv in hit[:40]: print(f'   x_{t} moves {mv}')
json.dump(hit,open(os.path.join(HERE,'pa_q_hits.json'),'w'))
