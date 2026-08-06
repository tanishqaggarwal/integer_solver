import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=H.fails()
print("nfree:", len(H.freeinp))
print("nvars:", H.NVARS, "ngates(order):", len(H.order), "neq:", len(H.eqcode))
print("baseline fails (20):", sorted(F))

# The 23 that break on core fix
break23=[3408,3841,4134,4526,5069,7276,15440,15724,15927,21600,22139,22825,27289,27999,28718,29305,31134,31269,32463,33195,36387,36390,38888]
allrel = sorted(set(F)|set(break23))
print("relevant eqs:", len(allrel))

# free ancestors of relevant equations
def eq_free_anc(i):
    s=set()
    for v in H.eqvars[i]:
        s |= (H.anc[v] if v in H.anc else ({v} if v in H.freeinp else set()))
    return s

relfree=set()
for i in allrel:
    relfree |= eq_free_anc(i)
print("free inputs touching relevant eqs:", len(relfree))

# Which free inputs touch the 20 verifier eqs vs 23 break eqs
vfree=set()
for i in F: vfree |= eq_free_anc(i)
bfree=set()
for i in break23: bfree |= eq_free_anc(i)
print("free touching 20 verifier:", len(vfree))
print("free touching 23 break:", len(bfree))
print("overlap:", len(vfree&bfree))
print("core vars free anc of x_29322:", 29322 in H.anc and sorted(H.anc[29322]))
