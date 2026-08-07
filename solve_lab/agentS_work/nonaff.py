import sys, json, collections, time, math
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
base=dict(C.BASE); v0=E.forward(base); bad0=E.badatoms(v0)
for a in (30787,26958,40306,726,28647,20215):
    try: order,fr,seen=E.cone(a)
    except Exception as e: fr=[]
    print(f"a{a}: {H.atoms[a][:180]}")
    print(f"    cone free vars ({len(fr)}): {fr[:12]}")
print()
for f in (30163,11559):
    print(f"--- sweep x_{f} (currently {v0[f]}) ---")
    for k in [1,2,3,5,7,11,-1,-2,100,12345]:
        b,_=fast.resid_delta(v0,bad0,{f:v0[f]+k})
        keys=set(b)|set(bad0)
        d={a:b.get(a,0)-bad0.get(a,0) for a in keys}
        d={a:x for a,x in d.items() if x}
        newvals={a:b.get(a,0) for a in sorted(set(d)|set(bad0))}
        print(f"  step {k:6d}: support={sorted(d)}  a28647%p={newvals.get(28647,0)%P if 28647 in newvals or 28647 in bad0 else '-'}")
