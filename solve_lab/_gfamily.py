import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('gadget_handled.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=H.fails()
print(f"gadget_handled fails: {len(F)}: {F}")
# The 21 gadget atoms are satisfied here. G1/G2 (11 eqs) open.
# Find all free vars that, when perturbed, affect the CURRENTLY-SATISFIED gadget-related eqs.
# Gadget-related satisfied eqs = those that were failing in my branch-B attempt but satisfied here.
gadget_eqs=[287, 1531, 2043, 3081, 6494, 7425, 8273, 8470, 8687, 9708, 13790, 18030, 22563, 23083, 27706, 29926, 30383, 32138, 35561, 37297, 38740]
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
# candidate knobs: free vars in gadget eqs + G1/G2 knobs
knobs=sorted(set(v for i in gadget_eqs for v in H.eqvars[i] if v in H.freeinp) | {9118,8731,17325,9413,7068,4432})
touch={v:sorted(set(desc_of[v])) for v in knobs}
def setf(v,x):
    H.val[v]=x
    for k in touch[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def gadresid():
    return [eval(H.eqcode[i],ns)%p for i in gadget_eqs]
def g1g2():
    return (7376877*H.val[642]+H.val[2099]-H.val[7068], H.val[4432]-H.val[19964]-H.val[28730])
inv2=pow(2,p-2,p)
# Jacobian of gadget residuals over knobs (mod p)
J=[[0]*len(knobs) for _ in gadget_eqs]
g10,g20=g1g2()
gvec=[]  # (dG1,dG2) per knob
for j,v in enumerate(knobs):
    o=base[v]
    setf(v,o+1); rp=gadresid(); gp1=g1g2()
    setf(v,o-1); rm=gadresid(); gm1=g1g2()
    setf(v,o)
    for k in range(len(gadget_eqs)): J[k][j]=((rp[k]-rm[k])*inv2)%p
    gvec.append(((gp1[0]-gm1[0]), (gp1[1]-gm1[1])))
# null space of J mod p = gadget-preserving directions. Among those, is there one with (dG1,dG2)!=0?
# Row reduce J to find pivot columns; free columns give null space basis.
def rref_track(M,ncol):
    M=[row[:] for row in M]; rows=len(M); r=0; piv={}
    for c in range(ncol):
        pr=None
        for i in range(r,rows):
            if M[i][c]%p: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; iv=pow(M[r][c]%p,p-2,p)
        M[r]=[(x*iv)%p for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][t]-f*M[r][t])%p for t in range(ncol)]
        piv[c]=r; r+=1
    return M,piv
Mr,piv=rref_track(J,len(knobs))
pivcols=set(piv); freecols=[c for c in range(len(knobs)) if c not in pivcols]
print(f"gadget Jacobian: {len(gadget_eqs)} eqs x {len(knobs)} knobs, rank={len(piv)}, nullity={len(freecols)}")
# For each null-space basis vector, compute (dG1,dG2)
found=False
for fc in freecols:
    # null vector: x_fc=1, x_pivcol = -Mr[piv[pc]][fc]
    nv=[0]*len(knobs); nv[fc]=1
    for pc,rr in piv.items(): nv[pc]=(-Mr[rr][fc])%p
    dG1=sum(nv[j]*gvec[j][0] for j in range(len(knobs)))%p
    dG2=sum(nv[j]*gvec[j][1] for j in range(len(knobs)))%p
    if dG1%p or dG2%p:
        found=True
print(f"Exists gadget-preserving direction changing G1 or G2 (mod p)? {found}")
print(f"free knobs (nullspace dim): {len(freecols)}")
