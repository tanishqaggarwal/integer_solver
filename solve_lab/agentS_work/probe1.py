import sys, json, collections
sys.path.insert(0,'.')
import common as C
import engine as E, fast
P=C.P
seed=dict(C.BASE)
v0=E.forward(seed); bad0=E.badatoms(v0)
print("bad0:", {a:len(str(abs(x))) for a,x in bad0.items()}, "(digits)")
KN=[30163,18956,14853,6083,31339,11559,30468,33169,11436,22820,26489,37012,4393,14393]
for f in KN:
    if f>=E.NV: continue
    o=v0[f]
    try:
        b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
    except Exception as ex:
        print(f"x_{f}: ERR {ex}"); continue
    keys=set(b1)|set(b2)|set(bad0)
    d1={a:b1.get(a,0)-bad0.get(a,0) for a in keys}
    d2={a:b2.get(a,0)-bad0.get(a,0) for a in keys}
    d1={a:x for a,x in d1.items() if x}
    aff=all(d2.get(a,0)==2*d1.get(a,0) for a in keys)
    print(f"\nx_{f}: cur={o} affine={aff} moves {len(d1)} atoms: {sorted(d1)}")
    for a in sorted(d1):
        R=bad0.get(a,0); D=d1[a]
        div = (R% D==0) if D else None
        print(f"    a{a}: R={len(str(abs(R)))}d  D={len(str(abs(D)))}d  R%D==0? {div}  need n={-R//D if (D and R%D==0) else 'n/a'}")
