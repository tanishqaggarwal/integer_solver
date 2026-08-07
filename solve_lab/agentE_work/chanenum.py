"""Exact channel-set x representative enumeration with the LOG 16 simultaneous solve."""
import sys, json, itertools, time, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import channels as C, engine as E
CFG=[({}, 'cfg0'), ({1530:0},'cfg5'), ({1530:0,1603:0},'cfg7')]
NREP=int(sys.argv[1]) if len(sys.argv)>1 else 3
best=(10**9,None,None,None)
for extra,tag in CFG:
    seed0=dict(C.base); seed0.update(extra)
    v0,bad0,cls=C.channels(seed0)
    keys=sorted(cls,key=lambda k:-len(cls[k]))
    print(f"=== {tag}: {len(keys)} channels, sizes {[len(cls[k]) for k in keys]}, base bad={sorted(bad0)}",flush=True)
    rnd=random.Random(11)
    for mask in itertools.product([0,1],repeat=len(keys)):
        reps=[]
        for i,k in enumerate(keys):
            if mask[i]: reps.append(cls[k][:NREP])
            else: reps.append([None])
        for combo in itertools.product(*reps):
            s=dict(seed0)
            for f in combo:
                if f is not None: s[f]=1
            t0=time.time()
            try: r=C.simsolve(s)
            except Exception as e:
                print(f"  mask{mask} {combo}: ERR {type(e).__name__}",flush=True); continue
            if r is None:
                print(f"  mask{mask} {combo}: no solution ({time.time()-t0:.0f}s)",flush=True); continue
            n,ns,av,v=r
            print(f"  mask{mask} {combo}: fails={n} score={39033-n} bad={av[:8]} ({time.time()-t0:.0f}s)",flush=True)
            if n<best[0]:
                best=(n,dict(ns),av,tag)
                json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('chan_%d.json'%(39033-n),'w'))
                json.dump({str(a):str(int(b)) for a,b in ns.items()}, open('chan_seed.json','w'))
                print(f"    *** NEW BEST {39033-n} ***",flush=True)
            if n==0: sys.exit(0)
print("BEST",best[0],39033-best[0],best[3],best[2],flush=True)
