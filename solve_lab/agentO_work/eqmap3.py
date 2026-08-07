import sys, json, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
bad=E.badatoms(vd); ff=E.eqfails(bad); BAD=sorted(bad)
print('failing eqs',sorted(ff))
# equations containing any bad atom
touch=collections.defaultdict(dict)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a in bad: touch[e][a]=c
print('equations touching a bad atom:',len(touch))
for e in sorted(touch):
    s=sum(c*bad[a] for a,c in touch[e].items())
    print(f'eq{e}: {touch[e]}  sum{"!=0 FAIL" if s else "=0 ok"}')
# per-atom: how many failing eqs contain it
for a in BAD:
    ineq=[e for e in ff if a in touch[e]]
    print(f'atom {a}: in {len(ineq)} failing eqs {ineq}')
