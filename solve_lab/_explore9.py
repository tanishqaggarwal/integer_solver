import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
def atomvars(ai):
    S=set()
    for vl,c in atoms[ai]['poly']: S|=set(vl)
    return S
gatedef=set(H.definer.keys())
for ai in [25170,27902,42851,43834,44270]:
    a=atoms[ai]
    print(f"\n=== atom#{ai} n_eq={a['n_eq']} eqs={a['eqs']} ===")
    print("repr:", a['repr'])
    vs=atomvars(ai)
    print("vars:", sorted(vs))
    print("free vars in atom:", sorted(v for v in vs if v in H.freeinp))
    print("gate vars in atom:", sorted(v for v in vs if v in gatedef))
