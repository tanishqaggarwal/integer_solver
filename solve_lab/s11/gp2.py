import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v=load_raw(src)
WIRE=set(u for u in range(L.NVARS) if v[u]==P)
AV=[atomval(a,v) for a in range(L.NA)]
nz=[a for a in range(L.NA) if AV[a]!=0]
nzp=[a for a in range(L.NA) if AV[a]%P!=0]
print(os.path.basename(src))
print(" atoms nonzero exactly:",len(nz), nz)
print(" atoms nonzero mod p  :",len(nzp), nzp)
for a in nzp:
    print(f"   a{a} out={L.atom_out.get(a)} val%p={AV[a]%P}")
# ABSORBER census: atoms containing a monomial  h*w  with w a wire and h free & used nowhere else
free=set(u for u in range(L.NVARS) if u not in L.definer)
absorb=0; noabs=[]
for a in range(L.NA):
    ok=False
    for m,c in L.polys[a].items():
        w=[u for u in m if u in WIRE]; o=[u for u in m if u not in WIRE]
        if len(w)==1 and len(o)==1 and o[0] in free and len(L.var_atoms[o[0]])==1 and abs(c)==1:
            ok=True; break
    if ok: absorb+=1
    else: noabs.append(a)
print(f"\natoms with a PRIVATE free p-handle (can absorb any multiple of p): {absorb} of {L.NA}")
print(f"atoms WITHOUT one: {len(noabs)}")
gn=[a for a in noabs if L.atom_out.get(a) is not None]
cn=[a for a in noabs if L.atom_out.get(a) is None]
print(f"   gates without handle: {len(gn)}  checks without handle: {len(cn)}")
