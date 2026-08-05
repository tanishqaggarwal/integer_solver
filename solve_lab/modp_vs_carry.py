import heal_harness as H
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
F0=set(H.fails())
print(f"39013 baseline fails: {len(F0)}: {sorted(F0)[:25]}")
# status of core control vars
for n in [14853,12186,24908,16742,29322,3558]:
    print(f"  x_{n}: {'FREE' if n in H.freeinp else 'gate'}")
# core residues
print(f"x_29322 = x_14853-x_12186, %p = {(H.val[14853]-H.val[12186])%p}")
print(f"x_3558  = x_24908-x_16742, %p = {(H.val[24908]-H.val[16742])%p}")
# apply linear residue fix: make x_14853 ≡ x_12186, x_24908 ≡ x_16742 mod p
if 14853 in H.freeinp:
    H.val[14853]-=(H.val[14853]-H.val[12186])%p
if 16742 in H.freeinp:
    H.val[16742]+=(H.val[24908]-H.val[16742])%p
H.forward()
F=set(H.fails())
newbroken=sorted(F-F0)
print(f"\nafter residue fix: {len(F)} fails; newly broken: {len(newbroken)}: {newbroken[:25]}")
print(f"x_29322%p now = {(H.val[14853]-H.val[12186])%p}, x_3558%p now = {(H.val[24908]-H.val[16742])%p}")
# for each newly broken eq, does it fail mod p or only in Z (carry)?
ns={'v':H.val,'__builtins__':{}}
modp_fail=0; carry_only=0
for i in newbroken:
    r=eval(H.eqcode[i],ns)
    if r%p==0: carry_only+=1
    else: modp_fail+=1
print(f"\nOF THE NEWLY BROKEN: {modp_fail} fail MOD P, {carry_only} fail ONLY in Z (carry/quotient)")
