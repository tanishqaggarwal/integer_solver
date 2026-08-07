"""Frame-B collateral measurement for the two open carriers K1 and L.

agentH_work is imported READ-ONLY; nothing is written there.  H's model has its own atom
numbering (42,267 atoms) but the SAME equation numbering as checker.py, so all prices below
are quoted in equation space, which is what the score is anyway.

K1 = x_7068 - x_2099   is the external part of the atom  (x7068 - x2099) - 7376877*x_642
L  = x_4432 - x_19964  is the external part of the atom  (x4432 - x19964) - x_28730
Both x_642 and x_28730 are DETACHED (free) in frame B, so the private knob absorbs part of the
needed shift: x_642 enters with coefficient -7376877 and x_28730 with coefficient -1.
"""
import sys, os, json, time, re
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_probe.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


t0 = time.time()
import frameB
say('frameB imported %.0fs' % (time.time() - t0))
DET = [642, 28730, 29854, 31864]
t0 = time.time()
fr = frameB.Frame(DET)
say('frame built %.0fs  free=%d checks=%d' % (time.time() - t0, len(fr.free), len(fr.checks)))

W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v = [0] * frameB.NV
for k, val in W.items():
    v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v[u] for u in fr.free if v[u] != 0}
t0 = time.time()
st = frameB.State(fr, fv)
say('witness state built %.0fs' % (time.time() - t0))
say('score = %d   failing = %s' % (st.score(), sorted(st.fails)))
dif = [i for i in range(frameB.NV) if st.v[i] != v[i]]
say('variables differing from the witness: %d  %s' % (len(dif), dif[:8]))
assert st.score() == 39026 and not dif, 'frame B does not reproduce the witness'
say('frame B reproduces the witness bit-for-bit\n')

# locate my region's atoms in H's numbering by source text
src = frameB.atom_src
NORM = lambda s: re.sub(r'\s+', '', s)
WANT = {
    'a23616': 'x_7068 - x_2099 - 7376877 * x_642',
    'a23617': 'x_28730 - x_17499 * x_9413',
    'a23618': 'x_4432 - x_19964 - x_28730',
    'a36659': 'x_29854 - x_22665 * x_1329',
    'a36660': '5113045 * (x_7075 * x_9118) - x_29854',
    'a36661': 'x_31864 - x_28961 * x_10903',
    'a36662': 'x_7075 * x_8731',
    'a36663': 'x_31864',
    'a36664': 'x_642 - x_28599 * x_17325',
}
idx = {}
bysrc = {}
for i, s in enumerate(src):
    bysrc.setdefault(NORM(s), i)
for tag, s in WANT.items():
    j = bysrc.get(NORM(s))
    idx[tag] = j
    say('  %s -> H atom %s   check: %s' % (tag, j, (j in set(fr.checks)) if j is not None else '-'))
REG = set(j for j in idx.values() if j is not None)

say('\n--- zero-collateral census in frame B (reproduce H\'s nine)')
zc = []
for u in fr.free:
    ck = fr.chk.get(u)
    if ck and set(ck) <= REG:
        zc.append(u)
say('  free inputs whose checks lie entirely inside the region: %s' % sorted(zc))

say('\n--- who can move the two open carriers?')
for tag in ('a23616', 'a23618'):
    a = idx[tag]
    if a is None:
        say('  %s not found' % tag)
        continue
    sup = fr.SUPV.get(a, [])
    say('  %s (H atom %d): %d free inputs support it' % (tag, a, len(sup)))
    rows = []
    for u in sup:
        ck = set(fr.chk.get(u, []))
        out = ck - REG
        eqs = set()
        for b in out:
            eqs.update(fr.eq_of[b])
        rows.append((len(eqs), len(out), u))
    rows.sort()
    for n_eq, n_ck, u in rows[:15]:
        say('      x_%-6d outside-checks %4d  outside-equations %5d %s'
            % (u, n_ck, n_eq, '  <<< ZERO COLLATERAL' if n_eq == 0 else ''))
    say('      minimum outside-equation count: %d' % (rows[0][0] if rows else -1))
json.dump({'atom_index': idx, 'zero_collateral': sorted(zc)},
          open(OD + '/fb_probe.json', 'w'))
say('DONE')
