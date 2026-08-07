import sys, json, math
sys.path.insert(0,'.')
import common as C
import engine as E, fast
P=C.P
seed=dict(C.BASE)
v0=E.forward(seed); bad0=E.badatoms(v0)
def delta(f,step=1):
    o=v0[f]; b1,_=fast.resid_delta(v0,bad0,{f:o+step})
    keys=set(b1)|set(bad0)
    return {a:b1.get(a,0)-bad0.get(a,0) for a in keys if b1.get(a,0)-bad0.get(a,0)}
PURE={22820:20215, 11436:20212, 14393:20212, 26489:7389, 37012:10187}
print("p =",P,"(%d bits)"%P.bit_length())
D={}
for f,a in PURE.items():
    d=delta(f); D[f]=d[a]
    x=d[a]
    print(f"\nx_{f} -> a{a}: D = {x}")
    print(f"   bits={x.bit_length()}  D%p={x%P}  D/p={x//P if x%P==0 else 'n/a'}")
    if x%P==0: print(f"   D = p * {x//P}   ({(x//P).bit_length()} bits)")
# other knobs
print("\n--- handles for a28647 ---")
for f in (30163,14853,6083):
    d=delta(f)
    for a,x in sorted(d.items()):
        print(f"x_{f} -> a{a}: bits={x.bit_length()} D%p={x%P if x else 0} div_by_p={x%P==0}"+(f" D/p bits={(x//P).bit_length()}" if x%P==0 else ""))
print("\n--- residuals mod p ---")
for a,r in sorted(bad0.items()):
    print(f"R[{a}]: bits={r.bit_length()}  mod p = {r%P}")
