import heal_harness as H, re
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
# branch B activation + close G1,G2 via free x_9118,x_8731 (x_17325=x_9413=0 already)
base[4287]=1; base[2081]=1
base[9118]=base[7068]   # x_2099=x_9118 -> G1: x_2099-x_7068=0
base[8731]=base[4432]   # x_19964=x_8731 -> G2: x_4432-x_19964=0
for v in H.freeinp: H.val[v]=base[v]
H.forward(); V=H.val
F=H.fails()
print(f"branch B with x_9118=x_7068, x_8731=x_4432: {len(F)} fails")
# gaps now
g1=7376877*V[642]+V[2099]-V[7068]; g2=V[4432]-V[19964]-V[28730]
print(f"G1={g1} G2={g2}  (x_2099={V[2099]==V[7068] and '=x_7068 OK' or V[2099]}, x_19964={V[19964]==V[4432] and '=x_4432 OK' or V[19964]})")
print(f"first 30 fails: {F[:30]}")
# which of these fails are NEW (not in original 11)?
orig11={2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125}
newf=[i for i in F if i not in orig11]; fixed=[i for i in orig11 if i not in F]
print(f"NEW fails (gadget etc): {len(newf)}: {newf}")
print(f"of orig 11, now fixed: {sorted(fixed)}")
print(f"of orig 11, still failing: {sorted(set(F)&orig11)}")
