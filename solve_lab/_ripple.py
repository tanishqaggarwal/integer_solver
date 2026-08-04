import heal_harness as H, re
p=H.p
VAR=re.compile(r'x_(\d+)')
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
# find eqs containing x_7068 and x_4432
for tgt in [7068,4432]:
    eqs=[i for i in range(len(lines)) if tgt in H.eqvars[i]]
    F=set([2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125])
    inF=[e for e in eqs if e in F]; notF=[e for e in eqs if e not in F]
    print(f"x_{tgt}: appears in {len(eqs)} eqs; {len(inF)} failing, {len(notF)} currently-satisfied: {notF}")
# For x_7068's currently-satisfied eqs, what free vars do they share, and are they also G1/G2-like?
