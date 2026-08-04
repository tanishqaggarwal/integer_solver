import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F=H.fails()
print(f"39022 fails ({len(F)}): {F}")
# Which of the key vars are free inputs vs gates?
for name in [642,2099,7068,4432,19964,28730,17325,9413,28599,17499]:
    print(f"x_{name}: {'FREE' if name in H.freeinp else 'gate'}  val%p={H.val[name]%p}")
# The gaps
G1 = 7376877*H.val[642]+H.val[2099]-H.val[7068]
G2 = H.val[4432]-H.val[19964]-H.val[28730]
print(f"\nG1 = 7376877*x_642 + x_2099 - x_7068 = {G1}")
print(f"   G1%p = {G1%p},  G1//p = {G1//p}")
print(f"G2 = x_4432 - x_19964 - x_28730 = {G2}")
print(f"   G2%p = {G2%p},  G2//p = {G2//p}")
# residue gaps that need to be zero
print(f"\n(x_7068 - x_2099) % p = {(H.val[7068]-H.val[2099])%p}")
print(f"(x_7068 - x_2099) % (7376877*p) = {(H.val[7068]-H.val[2099])%(7376877*p)}")
print(f"(x_4432 - x_19964) % p = {(H.val[4432]-H.val[19964])%p}")
