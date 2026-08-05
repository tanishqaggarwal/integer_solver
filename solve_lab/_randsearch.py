"""Broad randomized bit-pattern search (many seeds, various densities) + optional close.
Directly probes for any basin beating 11 fails."""
import _bitlab as L, heal_harness as H, random, time, json
best=(11,None)
t=time.time(); N=0
random.seed(12345)
log=[]
# Focus random flips on the 10 residual selectors + small random subsets of others
SEL10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
allsel=list(L.sels)
while time.time()-t<95:
    N+=1
    b=dict(L.AGENTA_BITS)
    mode=random.random()
    if mode<0.5:
        # random subset of the 10 residual selectors
        for s in SEL10: b[s]=random.randint(0,1)
    elif mode<0.8:
        # agentA + few random extra bits
        k=random.randint(1,4)
        for s in random.sample(allsel,k): b[s]=1-b[s]
    else:
        # random density pattern
        dens=random.choice([0.02,0.05,0.1,0.5])
        for s in allsel: b[s]=1 if random.random()<dens else 0
    F=L.apply_pattern(b,twopass=False)
    if len(F)<=best[0]:
        if len(F)<best[0] or best[1] is None:
            best=(len(F),{s:b[s] for s in allsel if b[s]})
            log.append((len(F),N))
print(f'tried {N} patterns in {time.time()-t:.0f}s. best fails={best[0]}')
print('best on-bits:',sorted(best[1].keys()) if best[1] else None)
print('improvement log (fails,iter):',log[-5:])
