import heal_harness as H, re
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1; base[9118]=base[7068]; base[8731]=base[4432]
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=H.fails()
orig11={2554, 6816, 8124, 8680, 9421, 12231, 12350, 14584, 29125,12270,22044}
newf=[i for i in F if i not in orig11]
# which new fails involve x_9118 or x_8731?
g9118=[i for i in newf if 9118 in H.eqvars[i]]
g8731=[i for i in newf if 8731 in H.eqvars[i]]
print(f"new gadget fails: {len(newf)}")
print(f"  involving x_9118: {g9118}")
print(f"  involving x_8731: {g8731}")
# handles per ENDGAME: x_9629,x_6947,x_33168. which fails involve them?
for h in [9629,6947,33168]:
    print(f"  x_{h} in fails: {[i for i in F if h in H.eqvars[i]]}, free? {h in H.freeinp}")
# Look at smallest gadget fail involving x_9118
cands=sorted(newf, key=lambda i: len(H.eqvars[i]))
print("\nsmallest new fails (by #vars):", [(i,len(H.eqvars[i])) for i in cands[:6]])
