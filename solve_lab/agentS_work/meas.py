"""Exact delta table: every cluster-cone knob x every cluster row, mod p and over Z.
   NO filtering: booleans (0<->1) and non-booleans (+1,+2,+7 affinity test) both included."""
import sys, json, collections, time, pickle
sys.path.insert(0,'.')
import common as C
import engine as E, fast
P=C.P; ROWS=C.ROWS

def measure(seed, tag, rows=None, extra_knobs=()):
    rows = rows if rows is not None else ROWS
    v0=E.forward(seed); bad0=E.badatoms(v0)
    cand=sorted(set(C.cluster_cone())|set(C.CLUSTERKN)|set(extra_knobs))
    out={}
    for f in cand:
        o=v0[f]
        b1,_=fast.resid_delta(v0,bad0,{f:o+1})
        d1={a:(b1.get(a,0)-bad0.get(a,0)) for a in rows}
        kind='bool' if C.isbool(f) else 'int'
        if kind=='bool':
            # boolean: only meaningful flip is 0<->1
            if o not in (0,1):
                kind='bool?'
            b=b1 if o==0 else None
            if o==1:
                bb,_=fast.resid_delta(v0,bad0,{f:0}); d1={a:(bb.get(a,0)-bad0.get(a,0)) for a in rows}
            out[f]=(kind,d1,None)
        else:
            b2,_=fast.resid_delta(v0,bad0,{f:o+2})
            d2={a:(b2.get(a,0)-bad0.get(a,0)) for a in rows}
            aff = all(d2[a]==2*d1[a] for a in rows)
            out[f]=('aff' if aff else 'nonaff', d1, d2)
    return v0,bad0,out

if __name__=='__main__':
    CFGS=[({}, 'cfg0'), ({1530:0},'cfg5'), ({1603:0},'cfg6'), ({1530:0,1603:0},'cfg7')]
    res={}
    for extra,tag in CFGS:
        s=dict(C.BASE); s.update(extra)
        t0=time.time()
        v0,bad0,out=measure(s,tag)
        R={a:bad0.get(a,0) for a in ROWS}
        print(f"=== {tag} extra={extra}  bad={sorted(bad0)}  ({time.time()-t0:.0f}s)",flush=True)
        for a in ROWS: print(f"    R[{a}] = {'0' if R[a]==0 else str(R[a])[:40]+'...('+str(len(str(abs(R[a]))))+'d)'}  mod p = {R[a]%P}")
        # classify by mod-p delta vector on ROWS
        cls=collections.defaultdict(list); kinds=collections.Counter()
        for f,(k,d1,d2) in out.items():
            kinds[k]+=1
            key=tuple(d1[a]%P for a in ROWS)
            if any(key): cls[key].append((f,k))
        print(f"    knob kinds: {dict(kinds)}")
        print(f"    {sum(len(v) for v in cls.values())} knobs move some row mod p; {len(cls)} distinct mod-p delta classes")
        for key,mem in sorted(cls.items(), key=lambda kv:-len(kv[1])):
            nb=sum(1 for _,k in mem if k.startswith('bool')); ni=len(mem)-nb
            sup=[ROWS[i] for i,x in enumerate(key) if x]
            print(f"      x{len(mem):4d} (b{nb}/i{ni}) support={sup}  reps={[f for f,_ in mem[:4]]}")
        res[tag]={'bad':{str(a):str(v) for a,v in bad0.items()},'out':{str(f):(k,{str(a):str(x) for a,x in d1.items()}) for f,(k,d1,d2) in out.items()}}
    pickle.dump(res, open('meas.pkl','wb'))
