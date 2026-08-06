import heal_harness as H
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
# message residue constants from pin atoms 3269,3271,3277
C_a=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110 # x_31861
C_b=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706 # x_14865
C_c=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319  # x_6418
print("current residues vs message residues:")
for name,C in [(31861,C_a),(14865,C_b),(6418,C_c)]:
    print(f"  x_{name}%p={V[name]%p}")
    print(f"    C%p ={C%p}  match={V[name]%p==C%p}")
# check x_12553 - is there a pin? current residue
print(f"  x_12553%p={V[12553]%p}")
# Now: what are C_c-C_a (x_17925 residue), and does gadget vanish at message residues?
print(f"\nx_17925 residue = (C_c-C_a)%p = {(C_c-C_a)%p}")
print(f"current x_17925%p = {V[17925]%p}")
