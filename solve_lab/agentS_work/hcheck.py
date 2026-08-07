"""For a list of atoms, find their pure handles (free vars whose exact delta support = {a}) and steps."""
import sys, collections, pickle, time, math, json
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
seed=dict(C.BASE); v0=E.forward(seed); bad0=E.badatoms(v0)
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}
TARG=[30787,20649,20652,32148,28033,28035,28037,26958,40306,747,20215,28647,20212,7389,10187,
      2849,2850,33297,33298,36152,36153,24177,24178,17061,17062]
cache={}
def probe(f):
    if f in cache: return cache[f]
    o=v0[f]
    b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
    keys=set(b1)|set(bad0)
    d1={a:b1.get(a,0)-bad0.get(a,0) for a in keys}; d1={a:x for a,x in d1.items() if x}
    d2={a:b2.get(a,0)-bad0.get(a,0) for a in keys}
    aff=all(d2.get(a,0)==2*d1.get(a,0) for a in keys)
    cache[f]=(d1,aff); return cache[f]
out={}
for a in TARG:
    try: order,fr,seen=E.cone(a)
    except Exception as e: print(f"a{a}: cone ERR {e}"); continue
    hs=[]
    t0=time.time()
    for f in fr:
        d1,aff=probe(f)
        if len(d1)==1 and a in d1: hs.append((f,d1[a],aff))
    def cls(x):
        if abs(x)==P: return 'p'
        if x%P==0: return '%d*p'%abs(x//P)
        if abs(x)==1: return 'UNIT(+-1)'
        return 'other(%dbits,gcd_p=%d)'%(x.bit_length(),math.gcd(abs(x),P))
    print(f"a{a}: nf={NF.get(a)} conefree={len(fr)} pure handles={len(hs)}: "
          f"{[(f,cls(x),aff) for f,x,aff in hs[:6]]}  ({time.time()-t0:.0f}s)",flush=True)
    out[a]=[(f,int(x),aff) for f,x,aff in hs]
json.dump({str(k):[(f,str(x),aff) for f,x,aff in v] for k,v in out.items()}, open('hcheck.json','w'))
