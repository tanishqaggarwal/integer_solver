"""Exact behaviour of x_11559: the only knob linking the cluster to a footprint-1 atom."""
import sys, json, collections, time, math
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}
base=dict(C.BASE); v0=E.forward(base); bad0=E.badatoms(v0)
print("a40306 footprint:",NF[40306],"a26958:",NF[26958],"a10187:",NF[10187])
print("a26958:",H.atoms[26958])
print("a40306 (first 300):",H.atoms[40306][:300])
print()
prev=None
for k in [0,1,2,3,4,5,10,100,10**20]:
    b,_=fast.resid_delta(v0,bad0,{11559:k})
    d={a:b.get(a,0)-bad0.get(a,0) for a in set(b)|set(bad0)}
    d={a:x for a,x in d.items() if x}
    s=[]
    for a in sorted(d):
        x=d[a]
        s.append(f"a{a}:{'p' if abs(x)==P else ('%d*p'%(x//P) if x%P==0 else '%db'%x.bit_length())}")
    print(f"  x_11559={k}: {', '.join(s)}")
# is the map k -> deltas linear?
b1,_=fast.resid_delta(v0,bad0,{11559:1}); b2,_=fast.resid_delta(v0,bad0,{11559:2})
for a in (10187,26958,40306):
    d1=b1.get(a,0)-bad0.get(a,0); d2=b2.get(a,0)-bad0.get(a,0)
    print(f"  a{a}: d(1)={d1.bit_length()}b  d(2)==2*d(1)? {d2==2*d1}  d(2)/d(1)={d2/d1 if d1 else None:.6f}" if d1 else f"  a{a}: d1=0")
