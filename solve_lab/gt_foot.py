import heal_harness as H
from collections import defaultdict
p=H.p
# equation footprint of each free input: which equations depend on it (via anc)
# eqvars[i] -> for each var, its free ancestors
eq_frees=[]
for i in range(len(H.lines)):
    s=set()
    for v in H.eqvars[i]:
        s|=H.anc.get(v, {v} if v in H.freeinp else set())
    eq_frees.append(s & H.freeinp)
foot=defaultdict(set)  # free -> set of eqs
for i,fs in enumerate(eq_frees):
    for f in fs: foot[f].add(i)

broken={697, 1985, 5225, 10815, 16048, 17784, 17801, 22402, 23667, 24721, 27124, 28737, 29638, 29959, 35935, 37431}
core={2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125}
# free inputs feeding the broken eqs
brk_frees=set()
for i in broken: brk_frees|=eq_frees[i]
print(f"free inputs feeding 16 broken eqs: {len(brk_frees)}")
# classify: which of these are PRIVATE to broken (footprint ⊆ broken), and which leak
private=[]; leak=[]
for f in brk_frees:
    outside=foot[f]-broken
    if not outside: private.append(f)
    else: leak.append((f,len(outside),len(foot[f]&broken)))
print(f"PRIVATE frees (only feed broken, safe to move): {len(private)}: {sorted(private)[:40]}")
print(f"leaking frees (feed broken + others): {len(leak)}")
# how much does each broken eq depend on private frees?
for i in sorted(broken):
    pf=eq_frees[i]&set(private)
    print(f"  eq{i}: {len(eq_frees[i])} frees, {len(pf)} private: {sorted(pf)[:8]}")
