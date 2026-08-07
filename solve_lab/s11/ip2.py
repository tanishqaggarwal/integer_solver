"""IP #2 -- a GLOBAL LOWER BOUND on the cost of any defect placement.

A check atom that owns a PRIVATE handle (a free variable occurring in exactly one atom, this
one, linearly) can always be closed by itself, so it can never be the thing you are forced to
violate.  Absorbers must come from the checks WITHOUT such a handle.  Hence

    min cost of any single absorber  >=  min |equations|  over checks with no private handle
    min cost of a 2-absorber set     >=  min over pairs of the union of their equations

computed over the whole instance, independent of channel.
"""
import sys, os, json, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
HERE = os.path.dirname(os.path.abspath(__file__))

CHECKS = [a for a in range(L.NA) if L.atom_out.get(a) is None]
# count of atoms each variable occurs in
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}


def private_handle(a):
    """free var occurring ONLY in atom a, and linearly -> a is self-closable"""
    for u in L.avars[a]:
        if L.definer.get(u) is not None:
            continue
        if NAT[u] != 1:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        return u
    return None


selfclose = []
absorber = []
for a in CHECKS:
    (selfclose if private_handle(a) is not None else absorber).append(a)
print(f"checks total           : {len(CHECKS)}")
print(f"  self-closable (private handle): {len(selfclose)}")
print(f"  possible ABSORBERS            : {len(absorber)}")

cost = {a: len(L.atom2eq.get(a, {})) for a in CHECKS}
hist = collections.Counter(cost[a] for a in absorber)
print("\nabsorber cost histogram (equations -> #checks):")
for k in sorted(hist)[:14]:
    print(f"   {k:3d} eqs : {hist[k]}")
cheap = sorted(absorber, key=lambda a: cost[a])[:40]
print("\ncheapest absorbers:", [(a, cost[a]) for a in cheap[:14]])

# also: how cheap are the SELF-CLOSABLE ones (these can never be forced)?
h2 = collections.Counter(cost[a] for a in selfclose)
print("\nself-closable cost histogram (first few):")
for k in sorted(h2)[:8]:
    print(f"   {k:3d} eqs : {h2[k]}")

lo = cost[cheap[0]]
print(f"\nLOWER BOUND, single absorber : {lo} equations")
best = None
for a, b in itertools.combinations(cheap[:40], 2):
    u = len(set(L.atom2eq.get(a, {})) | set(L.atom2eq.get(b, {})))
    if best is None or u < best[0]:
        best = (u, a, b)
print(f"LOWER BOUND, 2-absorber set  : {best[0]} equations  (atoms a{best[1]}, a{best[2]})")
print(f"\n  => no defect placement anywhere in the instance can cost fewer than "
      f"{best[0]} failing equations, i.e. score <= {L.NEQ - best[0]}")
# where does the checkpoint's absorber sit?
CK = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
e = set()
for a in CK:
    e |= set(L.atom2eq.get(a, {}))
print(f"\ncheckpoint absorber {CK}: {len(e)} equations; "
      f"{[a for a in CK if a in set(absorber)]} are non-self-closable")
