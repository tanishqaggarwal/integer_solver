import heal_harness as H
p=H.p
# Load the 39013 baseline into free inputs
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=H.fails()
print("baseline fails:", len(F))
print(sorted(F))
# core vars
for name,idx in [('x_14853',14853),('x_12186',12186),('x_24908',24908),('x_16742',16742),
                 ('x_29322',29322),('x_3558',3558),('x_35389',35389),('x_6671',6671)]:
    print(name, H.val[idx]%p, "  isfree", idx in H.freeinp)

print("\n=== Apply core fix ===")
# does x_24908 depend on x_16742 or x_14853?
print("anc x_24908 contains 16742?", 16742 in H.anc[24908], " 14853?", 14853 in H.anc[24908], " 12186?", 12186 in H.anc[24908])
print("anc x_29322:", 29322 in H.anc, "->", sorted(H.anc[29322])[:10] if 29322 in H.anc else None)
print("anc x_3558 ->", sorted(H.anc[3558])[:10] if 3558 in H.anc else None)
# apply fix
H.val[14853] -= (H.val[14853]-H.val[12186])%p
H.val[16742] += (H.val[24908]-H.val[16742])%p
H.forward()
F2=H.fails()
print("after fix fails:", len(F2), sorted(F2))
# residues
print("x_29322 modp", H.val[29322]%p, " x_3558 modp", H.val[3558]%p)
print("x_35389 modp", H.val[35389]%p, " x_6671 modp", H.val[6671]%p)
