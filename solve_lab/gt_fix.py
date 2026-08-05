import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F0=set(H.fails())
print(f"baseline fails ({len(F0)}): {sorted(F0)}")

# G1: want x_7068 = x_2099 + 7376877*p*k. gap = (x_7068-x_2099) mod (7376877*p); subtract it.
M1=7376877*p
gap1=(H.val[7068]-H.val[2099])%M1
H.val[7068]-=gap1
k1=(H.val[7068]-H.val[2099])//M1
H.val[17325]=k1   # so x_642 = p*k1, 7376877*x_642 = 7376877*p*k1 = x_7068-x_2099
# G2: want x_28730 = x_4432 - x_19964, x_28730 = p*x_9413. So x_4432-x_19964 ≡0 mod p.
gap2=(H.val[4432]-H.val[19964])%p
H.val[4432]-=gap2
k2=(H.val[4432]-H.val[19964])//p
H.val[9413]=k2
H.forward()
F1=set(H.fails())
print(f"after leaf-fix fails ({len(F1)}): {sorted(F1)}")
print(f"  fixed: {sorted(F0-F1)}")
print(f"  newly broken: {sorted(F1-F0)}")
# check G1,G2 now
G1 = 7376877*H.val[642]+H.val[2099]-H.val[7068]
G2 = H.val[4432]-H.val[19964]-H.val[28730]
print(f"  G1={G1}  G2={G2}")
