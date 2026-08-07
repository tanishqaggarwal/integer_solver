"""Flip a boolean selector, repair ITS pins with the general iterated closure+solve, then run
   the cluster solve.  Composition of the two things that each work in isolation."""
import sys, json, math, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, iterfix, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
res=pickle.load(open('boolknob.pkl','rb'))
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
# bits that move a10187 (scarce: 23) first, then a20212
b10=[f for f,(d,c,cur) in res.items() if d[10187]%P]
b20=[f for f,(d,c,cur) in res.items() if d[20212]%P]
order=sorted(b10)+sorted(set(b20)-set(b10))
print(f"a10187 movers {len(b10)}, a20212 movers {len(b20)}",flush=True)
best=None
for f in order[:60]:
    t0=time.time()
    s=dict(base); s[f]=1
    frozen={18956,1530,1603,f}
    ns,hist,ok=iterfix.iterate(s,frozen,iters=5,exclude=set(),log=open('/dev/null','w'))
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"x_{f}: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:10]} ({time.time()-t0:.0f}s)",flush=True)
    if best is None or len(ff)<best[0]:
        best=(len(ff),dict(ns),sorted(av),f)
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('boolfix_%d.json'%(39033-len(ff)),'w'))
        json.dump({str(a):str(int(b)) for a,b in ns.items()}, open('boolfix_seed.json','w'))
print("BEST",best[0],"bit",best[3],"bad",best[2],flush=True)
