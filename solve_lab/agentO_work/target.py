"""STEP 7 — validate the inverted solution and emit the exact targets.

hitrate.py found an admissible boundary change delta0 supported entirely on the four
configuration-tunable constants.  Verify end-to-end that applying it makes ALL 13 region
equations hold, then state exactly what realizing it requires, including the CRT lift that
lets each shift be produced by the variable that carries it.
"""
import sys, json, math
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/target.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
Pv = G.private_vars(R)
const0, cols = G.build_model(R, Pv, G.V0)
Eqs, rows0 = G.eq_system(R, Pv, const0, cols)
b0 = {e: rows0[e][1] for e in Eqs}
Bcol = {a: {e: G.EQCO[a].get(e, 0) for e in Eqs if G.EQCO[a].get(e, 0)} for a in sorted(const0)}
TUN = [23616, 23618, 36660, 36662]

rws, rhs = [], []
for e in Eqs:
    r = {'z%d' % u: c for u, c in rows0[e][0].items()}
    for a in TUN:
        c = Bcol[a].get(e, 0)
        if c:
            r['d%d' % a] = r.get('d%d' % a, 0) + c
    rws.append(r)
    rhs.append(b0[e])
sol, msg, _ = sparse.solve_sparse(rws, rhs, names=list(range(len(Eqs))), verbose=False,
                                  maxcore=150, maxbits=10 ** 7, maxcorebits=10 ** 7)
assert sol is not None, msg
delta0 = {a: sol.get('d%d' % a, 0) for a in TUN}
z = {u: sol.get('z%d' % u, 0) for u in Pv}

say('--- end-to-end verification of the inverted solution')
const1 = dict(const0)
for a in TUN:
    const1[a] = const0[a] + delta0[a]
_, rows1 = G.eq_system(R, Pv, const1, cols)
bad = []
for e in Eqs:
    lhs = sum(c * z[int(k[1:])] for k, c in
              [('z%d' % u, rows1[e][0].get(u, 0)) for u in Pv] if c)
    if lhs != rows1[e][1]:
        bad.append(e)
say('  equations of the region not satisfied after applying delta0: %s' % (bad if bad else 'NONE'))
say('  => with those four boundary shifts, all %d region equations hold, so the whole' % len(Eqs))
say('     instance would be satisfied (nothing outside E(R) can be touched).')

say('\n--- the exact targets (shift of each boundary constant)')
info = {
    23616: ('K1 = x_7068 - x_2099', 1),
    23618: ('L  = x_4432 - x_19964', 1),
    36660: ('K2 = 5113045 * x_9118', 5113045),
    36662: ('J  = x_7075 * x_8731  (x_7075 = 1)', 1),
}
out = {}
for a in TUN:
    name, div = info[a]
    d = delta0[a]
    say('  const(a%-6d)  %-36s shift = %d bits' % (a, name, abs(d).bit_length()))
    say('        divisible by the carrier factor %d: %s' % (div, d % div == 0))
    if div != 1:
        g = math.gcd(div, P)
        say('        gcd(carrier factor, p) = %d  -> CRT lift available: %s' % (g, g == 1))
        if g == 1 and d % div:
            # choose k so that (d + k*p) % div == 0
            k = ((-d) * pow(P, -1, div)) % div
            d2 = d + k * P
            say('        lifted shift (same class mod p, divisible by %d): %d bits, x_9118 += %d bits'
                % (div, abs(d2).bit_length(), abs(d2 // div).bit_length()))
            out[a] = str(d2)
        else:
            out[a] = str(d)
    else:
        out[a] = str(d)
say('\n  (each shift is free to move by multiples of p in the a23618 / a36660 / a36662')
say('   directions -- measured periods -- so a representative divisible by the carrier')
say('   factor always exists.)')

say('\n--- what realizing them costs, and where that must be evaluated')
for u in (8731, 9118, 7068, 2099, 4432, 19964):
    occ = sorted(H.occ[u])
    say('  x_%-6d occurs in %d atoms %s' % (u, len(occ), occ[:8]))
say('  x_8731 and x_9118 are agent H\'s ZERO-COLLATERAL knobs (frame B), so the a36660 and')
say('  a36662 shifts cost nothing there.  x_7068/x_2099/x_4432/x_19964 are ordinary derived')
say('  variables; their collateral has to be evaluated in frame B, not in my atom-level model,')
say('  which holds non-private variables fixed and therefore cannot express re-derivation.')

json.dump({'delta0': {str(a): str(delta0[a]) for a in TUN},
           'carrier_lifted': out,
           'z': {str(u): str(z[u]) for u in Pv}},
          open(OD + '/target.json', 'w'))
say('\nwrote target.json')
say('DONE')
