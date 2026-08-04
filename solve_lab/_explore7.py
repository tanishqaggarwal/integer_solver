import heal_harness as H, json
p=H.p
# For each control free, which gates does it define, and which equations are "checks"
# definer maps gate-output-var -> gate index. free inputs define nothing.
# Let's find, for each var v, the equations it appears in, and classify.
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
import re
VAR=re.compile(r'x_(\d+)')
eq_vars=[set(int(m) for m in VAR.findall(L)) for L in lines]

# which vars are gate-defined
gatedef=set(H.definer.keys())
print("controls free status:", {v:(v in H.freeinp) for v in [14853,12186,16742,24908]})
print("x_24908 definer gate rhs:", H.gates[H.definer[24908]][1] if 24908 in H.definer else "NONE")

# The gates that use x_14853 as input (candidates to be defined by it):
for cv in [14853,12186,16742]:
    print(f"\n=== x_{cv} ===")
    # gates that have cv as an input
    usedin=[]
    for gi,(t,rhs,vids) in enumerate(H.gates):
        if cv in vids:
            usedin.append((t,rhs,gi==H.definer.get(t,-1)))
    for t,rhs,isdef in usedin:
        print(f"  gate t=x_{t} rhs={rhs[:60]}  {'<-DEFINER(forward computes this)' if isdef else '(check/not-definer-for-this-gate)'}")
    print(f"  # equations containing x_{cv}: {sum(1 for ev in eq_vars if cv in ev)}")
