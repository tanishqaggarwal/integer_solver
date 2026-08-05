import sys, os, json, re
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

# build signed wire class (members = ±root)
par=list(range(H.NVARS)); sgn=[1]*H.NVARS
def find(x):
    if par[x]==x: return x,1
    r,s=find(par[x]); par[x]=r; sgn[x]*=s; return r,sgn[x]
def uni(a,b,s):
    ra,sa=find(a); rb,sb=find(b)
    if ra==rb: return
    par[rb]=ra; sgn[rb]=sa*s*sb
for poly in A.ATOMS:
    if all(len(vs)<=1 for vs,c in poly):
        lin=[(vs[0],c) for vs,c in poly if len(vs)==1]
        const=sum(c for vs,c in poly if len(vs)==0)
        if len(lin)==2 and const==0 and abs(lin[0][1])==1 and abs(lin[1][1])==1:
            (va,ca),(vb,cb)=lin
            uni(va,vb,-(ca*cb))
r0,_=find(26064)
wire={v:find(v)[1] for v in range(H.NVARS) if find(v)[0]==r0}
print(f"wire size {len(wire)}; x_26064 free? {26064 in H.freeinp}; is root? {r0==26064}")

# is the wire root a free input or gate output?
print(f"root {r0}: free={r0 in H.freeinp}")
wire_free=[v for v in wire if v in H.freeinp]
print(f"free wire members: {len(wire_free)}: {sorted(wire_free)[:20]}")

# Test wire=1: load 39022, set all free wire members to sgn*1, then also override gate-output wire
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F0=H.fails()
print(f"\nwire=p baseline fails: {len(F0)}")
# now set wire=1: override every wire member to sgn (after forward, force them)
# do forward then overwrite wire members and re-eval equations (won't be self-consistent for gates,
# but tells us how many eqs break from wire value alone)
for v,s in wire.items():
    H.val[v]= s*1
# recompute non-wire gates? forward recomputes all gates incl wire copies from root.
# Instead: set root's free control. Find the free input controlling the wire.
# The wire members are gate outputs (copies). The root r0 - is it computed from a free input?
# Check anc of r0:
print(f"anc(root {r0}) free inputs: {sorted(H.anc.get(r0,set()))[:10]}")
print(f"anc(x_26064) free inputs: {sorted(H.anc.get(26064,set()))[:10]}")

# just brute: overwrite all wire vals to sgn*1 and count how many equations fail (ignoring gate self-consistency)
F1=H.fails()
print(f"after overwriting wire={1}: fails={len(F1)} (delta {len(F1)-len(F0)})")
# which atoms nonzero now
nz=A.nonzero_atoms(H.val)
print(f"nonzero atoms: {len(nz)}")
# categorize: how many are the '37110-type' (x_26064-p) unpacking checks
small=[ai for ai,val in nz if len(A.ATOM_VARS[ai])<=3]
print(f"small nonzero atoms: {len(small)}: {small[:20]}")
