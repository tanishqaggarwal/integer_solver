"""Bit flips as DISCRETE knobs.

900 of the 1,156 free bits leave the failing SET at the same six checks.  That is not the same
as leaving the RESIDUES alone -- if the residues move, bit subsets are a discrete move set on top
of the (inconsistent) continuous one, and the question becomes whether some subset drives the six
residues into the continuous system's reachable set.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BITS=[]
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.append(ks[0][0])
BITS=sorted(set(u for u in BITS if u in FREE))
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
forwardp(base)
TGT=[7930,29539,35759,35760,40826,41512]
b0=[evalp(L.polys[a],base) for a in TGT]
def fails(v): return tuple(a for a in CHK if evalp(L.polys[a],v))
f0=fails(base)
print("base failing:",f0)
same=[]; t0=time.time()
for u in BITS:
    v=list(base); v[u]=(1-base[u])%P; forwardp(v)
    if fails(v)!=f0: continue
    r=[evalp(L.polys[a],v) for a in TGT]
    same.append((u,[ (r[j]-b0[j])%P for j in range(6)]))
print(f"bits preserving the failing set: {len(same)}  ({time.time()-t0:.0f}s)")
moving=[(u,d) for u,d in same if any(d)]
print(f"  of those, bits that MOVE the six residues: {len(moving)}")
for u,d in moving[:15]:
    print(f"    x{u}: nonzero components {[j for j in range(6) if d[j]]}")
json.dump([[u]+[str(x) for x in d] for u,d in same], open(os.path.join(HERE,'data','gmp15.json'),'w'))
