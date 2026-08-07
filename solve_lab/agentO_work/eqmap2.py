import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
bad=E.badatoms(vd); ff=E.eqfails(bad)
print('failing eq indices',sorted(ff))
BAD=sorted(bad)
eqof=collections.defaultdict(list)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for t in terms:
        a=t[0] if isinstance(t,(tuple,list)) else t
        if a in bad: eqof[a].append(e)
for a in BAD:
    print(f'atom {a}: eqs {eqof[a]}  "{H.atoms[a][:60]}"')
# which equations would survive if atom X were zero
for a in BAD:
    rem={x for x in BAD if x!=a}
    survive=set()
    for e,(issq,outer,terms) in enumerate(H.eqt):
        ats={(t[0] if isinstance(t,(tuple,list)) else t) for t in terms}
        if ats & rem: survive.add(e)
    print(f'  if atom {a} were 0 (others unchanged): failing eqs would be <= {len(survive & set(ff))}')
