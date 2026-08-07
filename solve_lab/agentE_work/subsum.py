"""Two-congruence 0/1 subset-sum over the selector bits, by enumerating class MULTIPLICITIES.
   Coefficients and targets are re-measured at every configuration (both are config-dependent)."""
import sys, json, math, re, collections, itertools, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
ROWS=[7389,10187,20212,20215,28647]
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
def isb(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False
CAND=sorted(set().union(*[set(E.cone(a)[1]) for a in ROWS]))
BOOLS=[f for f in CAND if isb(f)]
print(f"cluster cone: {len(CAND)} free vars, {len(BOOLS)} boolean",flush=True)

def measure(seed, tag):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    R={a:bad0.get(a,0) for a in ROWS}
    # combined coefficients per bit
    co={}
    for f in BOOLS:
        o=v0[f]
        if o==1: continue
        b1,_=fast.resid_delta(v0,bad0,{f:1})
        d={a:(b1.get(a,0)-bad0.get(a,0)) for a in ROWS}
        c1=(d[20212]+d[28647])%P
        c2=(d[20215]+d[10187])%P
        if c1 or c2: co[f]=(c1,c2)
    cls=collections.defaultdict(list)
    for f,(c1,c2) in co.items(): cls[(c1,c2)].append(f)
    T1=(R[28647]-R[20212])%P      # target for congruence 1 (see LOG 22 derivation)
    T2=(R[20215]-R[10187])%P
    print(f"--- {tag}: bad={sorted(bad0)}; {len(co)} moving bits, {len(cls)} classes",flush=True)
    for k,v in sorted(cls.items(), key=lambda kv:-len(kv[1]))[:9]:
        print(f"    x{len(v)}: c1={str(k[0])[:24]}.. c2={str(k[1])[:24]}..",flush=True)
    print(f"    targets T1={str(T1)[:26]}.. T2={str(T2)[:26]}..",flush=True)
    return v0,bad0,cls,T1,T2,R

def enumerate_mult(cls,T1,T2,cap=400000):
    keys=sorted(cls, key=lambda k:-len(cls[k]))
    ns=[len(cls[k]) for k in keys]
    tot=1
    for n in ns: tot*= (n+1)
    print(f"    multiplicity space size {tot}",flush=True)
    if tot>cap:
        print("    too large for direct enumeration",flush=True); return None
    hits=[]
    for m in itertools.product(*[range(n+1) for n in ns]):
        s1=0; s2=0
        for mi,k in zip(m,keys):
            if mi: s1=(s1+mi*k[0])%P; s2=(s2+mi*k[1])%P
        if s1==T1%P and s2==T2%P: hits.append(m)
    print(f"    multiplicity solutions: {len(hits)}  {hits[:5]}",flush=True)
    return hits,keys

CFG=[({},'cfg0 baseline (x_1530,x_1603)'),
     ({490:1},'cfg1 +x_490'),
     ({1530:0},'cfg5 only x_1603'),
     ({1530:0,1603:0},'cfg7 no selectors')]
for extra,tag in CFG:
    s=dict(base); s.update(extra)
    try:
        v0,bad0,cls,T1,T2,R=measure(s,tag)
        r=enumerate_mult(cls,T1,T2)
        if r and r[0]:
            hits,keys=r
            for m in hits[:3]:
                ns=dict(s)
                for mi,k in zip(m,keys):
                    for f in cls[k][:mi]: ns[f]=1
                v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
                print(f"    APPLIED {m}: fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:10]}",flush=True)
                if len(ff)<28:
                    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('subsum_%d.json'%(39033-len(ff)),'w'))
                    json.dump({str(a):str(int(b)) for a,b in ns.items()}, open('subsum_seed.json','w'))
    except Exception as e:
        print(f"--- {tag}: ERR {type(e).__name__} {e}",flush=True)
