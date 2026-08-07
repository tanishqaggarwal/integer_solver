"""Exact per-unit deltas of the trade knobs, then the exact combined congruence(s)."""
import sys, json, collections, pickle, math
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
base=dict(C.BASE); v0=E.forward(base); bad0=E.badatoms(v0)
def dl(f,step=1):
    o=v0[f]; b1,_=fast.resid_delta(v0,bad0,{f:o+step})
    keys=set(b1)|set(bad0)
    return {a:b1.get(a,0)-bad0.get(a,0) for a in keys if b1.get(a,0)-bad0.get(a,0)}
TR={}
for f in (6083,14853,31339,18956,22820,11436,26489,37012,30468,33169,7497):
    TR[f]=dl(f)
    print(f"x_{f}: "+", ".join(f"a{a}:{('p' if abs(x)==P else ('%d*p'%(x//P) if x%P==0 else str(x) if abs(x)<10**12 else '%dbits'%x.bit_length()))}" for a,x in sorted(TR[f].items())))
print()
# handle moduli
MOD={7389:P,10187:7942211*P,20212:P,20215:P,747:P}
NOH={28647,30787,26958,40306}
print("handle moduli:",{a:('p' if m==P else '%d*p'%(m//P)) for a,m in MOD.items()},"no-handle:",NOH)
# derive combined congruence for a28647 (no handle, traded via x_14853 into a20212 and x_6083 into a7389)
s14=TR[14853]; s60=TR[6083]
print("\nx_14853: a20212 %+d , a28647 %+d per unit"%(s14[20212],s14[28647]))
print("x_6083 : a7389  %+d , a28647 %+d per unit"%(s60[7389],s60[28647]))
img=pickle.load(open('image.pkl','rb'))['img']
print("\nover the %d sampled image points, evaluate the combined invariants:"%len(img))
ROWS=C.ROWS  # [7389,10187,20212,20215,28647]
hits=collections.Counter()
for key,c in img.items():
    R7,R10,R20212,R20215,R28647=key
    # zero a28647 exactly using n14 (Δ20212=s14[20212], Δ28647=s14[28647]) and n60
    # need: R28647 + s14[28647]*n14 + s60[28647]*n60 == 0
    #       R20212 + s14[20212]*n14 == 0 mod p    -> n14 == -R20212 * inv(s14[20212]) mod p
    #       R7389  + s60[7389]*n60  == 0 mod p    -> n60 == -R7389  * inv(s60[7389])  mod p
    n14=(-R20212*pow(s14[20212]%P,-1,P))%P
    n60=(-R7*pow(s60[7389]%P,-1,P))%P
    inv=(R28647 + s14[28647]*n14 + s60[28647]*n60)%P
    # a20215: handle step p -> need R20215 == 0 mod p, or trade via x_31339 into a10187 (mod 7942211p)
    hits[(inv==0, R20215==0)]+=c
    if inv==0: print("  *** combined invariant for a28647 == 0 for image point",str(key[4])[:20])
print("counts (a28647-combined==0, a20215==0):",dict(hits))
# what IS the invariant value set?
vals=collections.Counter()
for key,c in img.items():
    R7,R10,R20212,R20215,R28647=key
    n14=(-R20212*pow(s14[20212]%P,-1,P))%P
    n60=(-R7*pow(s60[7389]%P,-1,P))%P
    vals[(R28647 + s14[28647]*n14 + s60[28647]*n60)%P]+=c
print("\ndistinct combined-invariant values (%d):"%len(vals))
for v,c in vals.most_common(): print(f"   x{c:4d} {v}")
print("\ndistinct a20215 mod p values:")
for v,c in collections.Counter({k[3]:0 for k in img}).items(): pass
s=collections.Counter()
for key,c in img.items(): s[key[3]]+=c
for v,c in s.most_common(): print(f"   x{c:4d} {v}")
