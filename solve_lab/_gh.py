import heal_harness as H
p=H.p
d=H.loadd('gadget_handled.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); V=H.val
F=H.fails()
print(f"gadget_handled.json: {len(F)} fails: {F}")
print(f"x_4287={V[4287]}, x_2081={V[2081]}, x_9062={V[9062]}, x_21279={V[21279]}")
print(f"x_2099={V[2099]%p==V[9118]%p and 'equals x_9118' or V[2099]}, x_9118={V[9118]}")
print(f"x_19964 equals x_8731? {V[19964]==V[8731]}, x_8731={V[8731]}")
# gaps
g1=7376877*V[642]+V[2099]-V[7068]
g2=V[4432]-V[19964]-V[28730]
print(f"G1={g1} (mod p={g1%p}), bits={abs(g1).bit_length()}")
print(f"G2={g2} (mod p={g2%p}), bits={abs(g2).bit_length()}")
print(f"x_7068={V[7068]}, x_4432={V[4432]}, x_642={V[642]}, x_28730={V[28730]}")
