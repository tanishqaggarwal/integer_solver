"""usage: sweep_run.py <outjsonl> <maxr> <maxv> <combo1> [<combo2> ...]
   combo = comma-separated bit list, or 'empty'."""
import sys, json, time, os
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import sweep, engine as E
out=sys.argv[1]; maxr=int(sys.argv[2]); maxv=int(sys.argv[3])
combos=[[] if c=='empty' else [int(x) for x in c.split(',')] for c in sys.argv[4:]]
BEST=os.environ.get('SAVEBEST')
f=open(out,'a',buffering=1)
best=10**9
for combo in combos:
    s=dict(sweep.C.base)
    for b in combo: s[b]=1
    r,msg,t=sweep.run(s,maxr=maxr,maxv=maxv)
    rec={'combo':combo,'msg':msg,'t':round(t,1),'maxr':maxr,'maxv':maxv}
    if r:
        n,ns,av,v=r
        rec.update(fails=n,score=39033-n,bad=av)
        if BEST and n<best:
            best=n
            json.dump({f"x_{i}":str(int(v[i])) for i in range(E.NV) if v[i]!=0},
                      open(f'{BEST}_{39033-n}.json','w'))
            json.dump({str(a):str(int(b)) for a,b in ns.items()}, open(f'{BEST}_{39033-n}_seed.json','w'))
    f.write(json.dumps(rec)+'\n')
    print(rec, flush=True)
