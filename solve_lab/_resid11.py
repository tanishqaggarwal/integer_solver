import heal_harness as H, json, re
p=H.p
VAR=re.compile(r'x_(\d+)')
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
ns={'v':H.val,'__builtins__':{}}
# for each failing eq: integer LHS, variables, which free, mod-p value
for i in F:
    L=lines[i]; lhs=L.rsplit('=',1)[0]
    val=eval(H.eqcode[i],ns)
    vs=sorted(set(int(m) for m in VAR.findall(lhs)))
    freev=[v for v in vs if v in H.freeinp]
    print(f"eq {i}: LHS={val}")
    print(f"   |LHS| bits={val.bit_length()}, LHS mod p={val%p}, LHS/p={'int:'+str(val//p) if val%p==0 else 'NOT div by p'}")
    print(f"   #vars={len(vs)} #free={len(freev)} free={freev[:20]}")
