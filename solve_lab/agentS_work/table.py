"""Full exact move table at a configuration: every cone knob -> (support over ALL atoms, exact deltas)."""
import sys, json, collections, pickle, time
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P

FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}

def build(seed, knobs=None):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    if knobs is None:
        knobs=sorted(set(C.cluster_cone())|set(C.CLUSTERKN))
    T={}
    for f in knobs:
        o=v0[f]
        b1,_=fast.resid_delta(v0,bad0,{f:o+1})
        keys=set(b1)|set(bad0)
        d1={a:b1.get(a,0)-bad0.get(a,0) for a in keys}
        d1={a:x for a,x in d1.items() if x}
        b2,_=fast.resid_delta(v0,bad0,{f:o+2})
        d2={a:b2.get(a,0)-bad0.get(a,0) for a in keys}
        aff=all(d2.get(a,0)==2*d1.get(a,0) for a in keys)
        isb=C.isbool(f)
        T[f]=dict(cur=o,aff=aff,bool=isb,d=d1)
    return v0,bad0,T

if __name__=='__main__':
    seed=dict(C.BASE)
    t0=time.time(); v0,bad0,T=build(seed); print("built %.0fs"%(time.time()-t0))
    print("bad0:",{a:NF[a] for a in sorted(bad0)},"total eqfails",len(E.eqfails(bad0)))
    # group knobs by exact support
    grp=collections.defaultdict(list)
    for f,r in T.items(): grp[tuple(sorted(r['d']))].append(f)
    print(f"\n{len(T)} knobs, {len(grp)} distinct exact supports")
    rows=[]
    for sup,fs in sorted(grp.items(), key=lambda kv:(len(kv[0]),-len(kv[1]))):
        if not sup: continue
        cost=len(set().union(*[FOOT[a] for a in sup]))
        nb=sum(1 for f in fs if T[f]['bool'])
        na=sum(1 for f in fs if T[f]['aff'])
        rows.append((len(sup),cost,len(fs),nb,na,sup,fs[:5]))
    for L,cost,n,nb,na,sup,fs in rows:
        if L<=4:
            print(f"  |sup|={L} eqcost={cost:3d} x{n:4d} (bool {nb}, aff {na}) sup={sup} nf={[NF[a] for a in sup]} reps={fs}")
    pickle.dump({'bad0':bad0,'T':T,'NF':NF},open('table_cfg0.pkl','wb'))
    # which knobs move a20215 or a28647
    print("\n--- knobs touching a20215 or a28647 ---")
    for f,r in sorted(T.items()):
        if 20215 in r['d'] or 28647 in r['d']:
            print(f"  x_{f} bool={r['bool']} aff={r['aff']} sup={sorted(r['d'])} nf={[NF[a] for a in sorted(r['d'])]}")
