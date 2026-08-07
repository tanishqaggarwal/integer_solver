"""Emit delta0 in a frame-independent, consumable form.

Atom NUMBERING differs between models (E 40,727 / H 42,267 / M its own), so everything here is
keyed by atom SOURCE TEXT and by variable index, both of which are shared.  Equation indices
are shared too (they match checker.py's failing-line indices).

Writes agentO_work/DELTA0_FOR_M.json  (machine) and  agentO_work/DELTA0_FOR_M.md  (human).
"""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H, sparse

OD = '/home/user/integer_solver/solve_lab/agentO_work'
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663

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

SRC = {a: H.atoms[a] for a in R}
CARRIER = {
    23616: {'expr': 'x_7068 - x_2099', 'vars': [7068, 2099], 'factor': 1,
            'note': 'external part of the atom that would define handle x_642; the private '
                    'knob x_642 enters with coefficient -7376877, so this shift only matters '
                    'modulo 7376877 = 3*2458959'},
    23618: {'expr': 'x_4432 - x_19964', 'vars': [4432, 19964], 'factor': 1,
            'note': 'external part of the atom that would define handle x_28730; the private '
                    'knob x_28730 enters with coefficient -1 here and +1 in the x_28730 - p*x_9413 '
                    'atom, so only the DIFFERENCE of the two directions is a new degree of freedom'},
    36660: {'expr': '5113045 * x_9118', 'vars': [9118], 'factor': 5113045,
            'note': 'carried by the free variable x_9118 (zero-collateral in frame B); the shift '
                    'may be moved by multiples of p, and gcd(5113045, p) = 1, so a representative '
                    'divisible by 5113045 always exists'},
    36662: {'expr': 'x_7075 * x_8731  (x_7075 = 1)', 'vars': [8731], 'factor': 1,
            'note': 'carried by the free variable x_8731 (zero-collateral in frame B)'},
}


def lift(a):
    d = delta0[a]
    f = CARRIER[a]['factor']
    if f == 1 or d % f == 0:
        return d
    k = ((-d) * pow(P, -1, f)) % f
    return d + k * P


out = {
    'what': 'delta0: an integral boundary shift that makes ALL 13 equations of the 39,026 '
            'witness residual region hold simultaneously.  Verified end-to-end in agent O\'s '
            'region model.  NOT verified as realisable - the collateral of the two derived-'
            'variable carriers is the open question.',
    'p': str(P),
    'region_atoms_by_source': {str(a): SRC[a] for a in R},
    'region_equations': Eqs,
    'private_variables': Pv,
    'private_solution_z': {str(u): str(z[u]) for u in Pv},
    'boundary_shifts': {},
    'blocking_moduli_without_the_shift': {
        'x_1329': 'p', 'x_9413': 'p', 'x_10903': 'p',
        'x_17325': '2458959 * p', 'x_642': '2458959'},
    'scan_verdict': 'a configuration scan is useless here: the four boundary quantities are '
                    'identically 0 across 35 configurations in E\'s frame (1 distinct value '
                    'out of 35 each), and the admissible lattice has index >= 2^768.',
    'how_to_price': 'apply the four shifts via the listed carrier variables, re-propagate, and '
                    'count failing equations.  Equation indices here match checker.py line '
                    'indices.  The witness fails [12231,12270,12350,14584,18673,22044,29125].',
}
for a in TUN:
    c = CARRIER[a]
    d = delta0[a]
    dl = lift(a)
    out['boundary_shifts'][str(a)] = {
        'atom_source': SRC[a],
        'external_expression': c['expr'],
        'carrier_variables': c['vars'],
        'carrier_factor': c['factor'],
        'shift': str(d),
        'shift_bits': abs(d).bit_length(),
        'shift_lifted_divisible_by_carrier_factor': str(dl),
        'carrier_increment': str(dl // c['factor']),
        'free_to_move_by': 'multiples of p' if a != 23616 else 'period exceeds 2458959*p (unmeasured)',
        'note': c['note'],
    }
json.dump(out, open(OD + '/DELTA0_FOR_M.json', 'w'), indent=1)

md = ['# delta0 — an exact lattice target for the 39,026 residual region',
      '',
      'Agent O.  Everything is keyed by atom **source text** and **variable index** (shared',
      'across models); atom *numbering* is not (E 40,727 / H 42,267 / M its own).  Equation',
      'indices ARE shared and match `checker.py` line indices.',
      '',
      '## What this is',
      'The 39,026 witness residual lives in 9 atoms touched by exactly 13 equations, with 8',
      'variables private to that region.  Over Q the system has a unique solution satisfying all',
      '13; over Z five coordinates are blocked, by moduli `p, p, p, 2458959, 2458959*p`.',
      'Solving `A z + B d = b0` over Z gives an integral boundary shift `d = delta0` supported',
      'on exactly the four constants that are NOT p-multipliers.  Applying it makes all 13',
      'region equations hold — **verified end-to-end in the region model**.',
      '',
      '**It is NOT verified as realisable.**  Two of the four carriers are free variables and',
      'cost nothing; the other two are derived, and their collateral is the open question.',
      'That is what needs pricing.',
      '',
      '## The four shifts', '']
for a in TUN:
    s = out['boundary_shifts'][str(a)]
    md += ['### atom `%s`' % s['atom_source'],
           '- external expression: `%s`' % s['external_expression'],
           '- carrier variables: %s   (carrier factor %d)' % (s['carrier_variables'], s['carrier_factor']),
           '- **shift = %d bits**, free to move by %s' % (s['shift_bits'], s['free_to_move_by']),
           '- carrier increment (already divided by the carrier factor): %d bits'
           % abs(int(s['carrier_increment'])).bit_length(),
           '- %s' % s['note'], '']
md += ['## Which are free',
       '- `x_8731` (atom `x_7075 * x_8731`) and `x_9118` (atom `5113045 * (x_7075 * x_9118) - x_29854`)',
       '  are agent H\'s zero-collateral knobs in frame B — these two shifts cost nothing there.',
       '- `x_7068 - x_2099` and `x_4432 - x_19964` are derived; **these two need pricing**.',
       '  They are the external parts of the atoms that would define handles `x_642` and `x_28730`,',
       '  i.e. two of the four handles the deliverable itself corrupts.',
       '',
       '## Two simplifications worth using',
       '- `x_642` enters `x_7068 - x_2099 - 7376877*x_642` with coefficient `-7376877`, and it is',
       '  private, so the `x_7068 - x_2099` shift **only matters modulo 7376877 = 3 x 2458959**.',
       '  That is a 23-bit condition, not a 2440-bit one.',
       '- `x_28730` enters `x_4432 - x_19964 - x_28730` with coefficient `-1` and the',
       '  `x_28730 - p*x_9413` atom with `+1`, both private, so only the DIFFERENCE of those two',
       '  directions is a genuinely new degree of freedom.',
       '',
       '## Do not scan configurations',
       'Measured: the four boundary quantities are identically 0 across 35 configurations in E\'s',
       'frame (1 distinct value out of 35 each) — a scan measures one point repeatedly — and the',
       'admissible lattice has index >= 2^768 (hit rate ~2^-767).',
       '',
       'Machine-readable copy: `DELTA0_FOR_M.json` (exact integers).', '']
open(OD + '/DELTA0_FOR_M.md', 'w').write('\n'.join(md))
print('wrote DELTA0_FOR_M.json and DELTA0_FOR_M.md')
for a in TUN:
    print('  a%d shift %d bits, carrier increment %d bits' %
          (a, abs(delta0[a]).bit_length(),
           abs(int(out['boundary_shifts'][str(a)]['carrier_increment'])).bit_length()))
