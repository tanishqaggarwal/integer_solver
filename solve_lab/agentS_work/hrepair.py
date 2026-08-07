"""Handle-repaired scoring: every atom with a pure single-atom handle whose step divides its
   residual can be zeroed independently.  Score = eqfails AFTER all such repairs."""
import sys, collections, pickle, time, math, json, os
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}

_HC='handle_cache.pkl'
HAND=pickle.load(open(_HC,'rb')) if os.path.exists(_HC) else {}
_base=None
def _basestate():
    global _base
    if _base is None:
        v0=E.forward(dict(C.BASE)); _base=(v0,E.badatoms(v0))
    return _base

def find_handle(a):
    """Pure handle for atom a: free var whose exact delta support is {a}. Cached."""
    if a in HAND: return HAND[a]
    v0,bad0=_basestate()
    try: order,fr,seen=E.cone(a)
    except Exception: HAND[a]=None; return None
    best=None
    for f in fr:
        o=v0[f]
        try:
            b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
        except Exception: continue
        keys=set(b1)|set(bad0)
        d1={x:b1.get(x,0)-bad0.get(x,0) for x in keys}; d1={x:y for x,y in d1.items() if y}
        if len(d1)!=1 or a not in d1: continue
        if b2.get(a,0)-bad0.get(a,0)!=2*d1[a]: continue     # need affine
        s=d1[a]
        if best is None or abs(s)<abs(best[1]): best=(f,s)
    HAND[a]=best
    return best

def save(): pickle.dump(HAND,open(_HC,'wb'))

def repair(seed, verbose=False):
    """Return (nfails, remaining_bad, full_assignment_vector)."""
    v=E.forward(seed); bad=E.badatoms(v)
    ns=dict(seed); rem={}
    for a,R in bad.items():
        h=find_handle(a)
        if h is None: rem[a]=R; continue
        f,s=h
        if R % s: rem[a]=R; continue
        ns[f]=v[f] - R//s
    if ns!=seed:
        v=E.forward(ns); bad2=E.badatoms(v)
    else: bad2=bad
    ff=E.eqfails(bad2)
    return len(ff), bad2, v, ns

if __name__=='__main__':
    for name,seedfile in [('triple8(E best)','triple8_seed.json')]:
        seed={int(k):int(v) for k,v in json.load(open(seedfile)).items()}
        t0=time.time(); n,bad,v,ns=repair(seed)
        print(f"{name}: raw eqfails vs repaired -> {n} fails, score {39033-n}, bad={sorted(bad)} ({time.time()-t0:.0f}s)")
        for a in sorted(bad): print(f"    a{a} nf={NF[a]} handle={find_handle(a)[0] if find_handle(a) else None}")
    save()
    print("handle cache size",len(HAND))
