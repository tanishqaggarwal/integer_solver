"""Scan all 2^10 combos of the 10 residual-feeding selectors (raw, continuous held at agentA)."""
import _bitlab as L, heal_harness as H, itertools, json, time
SEL10=[2081, 4287, 5910, 11368, 13195, 17406, 18022, 22562, 23751, 28005]
base=dict(L.AGENTA_BITS)
res=[]
t=time.time()
for combo in itertools.product([0,1],repeat=10):
    b=dict(base)
    for s,v in zip(SEL10,combo): b[s]=v
    F=L.apply_pattern(b,twopass=False)
    res.append((len(F),combo,F))
res.sort(key=lambda x:x[0])
print('scanned 1024 in %.1fs'%(time.time()-t))
print('agentA combo (1,0,0,0,0,0,0,0,0,0)?  best 15:')
for n,combo,F in res[:15]:
    on=[SEL10[i] for i in range(10) if combo[i]]
    print(f'  {n} fails  ON={on}')
from collections import Counter
print('distribution:',sorted(Counter(n for n,_,_ in res).items())[:12])
json.dump([(n,list(c),f) for n,c,f in res], open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/scan10.json','w'))
