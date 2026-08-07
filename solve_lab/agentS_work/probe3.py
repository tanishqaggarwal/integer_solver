import sys, json, collections
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
# atom -> equations footprint
foot=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: foot[a].add(e)
cnt=collections.Counter({a:len(s) for a,s in foot.items()})
print("atoms with equation-footprint k:", collections.Counter(cnt.values()).most_common(12))
for a in (20215,28647,747,30787,20212,7389,10187,26958,40306,722,724,726):
    print(f"  a{a}: appears in {len(foot[a])} equations -> {sorted(foot[a])[:8]}")
# deliverable's 8 atoms
import checker
seed=dict(C.BASE); v0=E.forward(seed); bad0=E.badatoms(v0)
print("\ncfg0 bad:",sorted(bad0), "eqfails", len(E.eqfails(bad0)))
print("union footprint of {20215,28647}:", len(foot[20215]|foot[28647]))
