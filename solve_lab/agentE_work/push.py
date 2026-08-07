import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(5000000)
import engine as E, iterfix, harness as H
s={int(k):int(v) for k,v in json.load(open(sys.argv[1])).items()}
frozen={18956,1530,1603}
v=E.forward(s); av=E.badatoms(v)
print("start bad",sorted(av),"fails",len(E.eqfails(av)),flush=True)
# selector MUX for the (1,1) branch: x_22162 = x_13682, x_30213 = x_18956 - x_32237
for _ in range(6):
    s[22162]=v[13682]; s[30213]=v[18956]-v[32237]; v=E.forward(s)
av=E.badatoms(v); print("after MUX fixpoint bad",sorted(av),"fails",len(E.eqfails(av)),flush=True)
ns,hist,ok=iterfix.iterate(s,frozen|{22162,30213},iters=8,exclude=set(),log=sys.stdout)
v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
print("FINAL fails",len(ff),"score",39033-len(ff),"bad",sorted(av),flush=True)
json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('push_%d.json'%(39033-len(ff)),'w'))
json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('push_seed.json','w'))
for a in sorted(av):
    print("  ATOM",a,H.atoms[a][:150],flush=True)
