"""The 15 failing equations are recoverable in ATOM space.  Are the required atom values
   REACHABLE?  Parametrise by the actual free handles and solve the integer system on the
   equation values directly."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
KNOBS = [26719, 26720, 26721, 26722, 26723, 28437]

v = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'finish3.json'))).items():
    v[int(k)] = int(x)
fw.forward(v)
av = L.all_atom_values(v)
S = sorted(L.failing_eqs(av))
print(f"failing equations: {len(S)}  nonzero atoms: {[a for a,x in enumerate(av) if x]}")

# candidate handles: free vars with an exact linear effect on any knob atom
cands = set()
for a in KNOBS:
    try:
        hs, base = deep.handles(v, a, locked=set())
        for t, d in hs:
            cands.add(t)
    except Exception:
        pass
cands = sorted(cands)
print(f"candidate handles from the knob atoms: {len(cands)} -> {cands[:25]}")


def eqvals(vv):
    a2 = L.all_atom_values(vv)
    return [L.eq_value(e, a2) for e in S]


base = eqvals(v)
print("base equation values (bit lengths):", [x.bit_length() for x in base][:15])

cols = []
used = []
for t in cands:
    old = v[t]
    v[t] = old + 1
    fw.forward(v)
    e1 = eqvals(v)
    v[t] = old + 2
    fw.forward(v)
    e2 = eqvals(v)
    v[t] = old
    fw.forward(v)
    d1 = [e1[i] - base[i] for i in range(len(S))]
    d2 = [e2[i] - base[i] for i in range(len(S))]
    if all(d2[i] == 2 * d1[i] for i in range(len(S))) and any(d1):
        cols.append(d1)
        used.append(t)
print(f"handles with EXACT LINEAR effect on the 15 equations: {len(used)} -> {used}")
if used:
    M = [[cols[j][i] for j in range(len(used))] for i in range(len(S))]
    rhs = [-base[i] for i in range(len(S))]
    x = solve_int(M, rhs)
    print("integer solution:", "FOUND" if x else "NONE")
    if x:
        for j, t in enumerate(used):
            v[t] += x[j]
        fw.forward(v)
        f = L.failing_eqs(L.all_atom_values(v))
        b = fw.bad_checks(v)
        print(f"AFTER: failing={len(f)} score={L.NEQ-len(f)} bad_checks={len(b)}")
        json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
                  open(os.path.join(HERE, 'data', 'realise_named.json'), 'w'))
    else:
        # how many equations can be recovered?
        import itertools
        for drop in range(1, 10):
            hit = None
            for combo in itertools.combinations(range(len(S)), drop):
                keep = [i for i in range(len(S)) if i not in combo]
                if solve_int([M[i] for i in keep], [rhs[i] for i in keep]) is not None:
                    hit = combo
                    break
            if hit:
                print(f"  best: recover {len(S)-drop} of {len(S)} (drop {drop}: {[S[i] for i in hit]}) "
                      f"-> failing {drop}, score {L.NEQ-drop}")
                break
