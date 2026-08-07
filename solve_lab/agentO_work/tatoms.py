"""T's 20 atoms: where does each one live, and what would making it nonzero cost?

T = 0 forces  delta(a23618) = -sum c_i delta(a_i)  over T's other 19 atoms.  If any of those 19
appears ONLY in eq8680, it is a FREE compensator: moving it keeps T = 0 achievable while
breaking nothing, which would restore the L direction delta0 needs at zero cost.  If they all
appear in many equations, every compensation is budgeted and the search below is the whole
question.
"""
import sys, collections, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/tatoms.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


issq, outer, tl = H.eqt[8680]
say('eq8680: issq=%s outer=%s terms=%d' % (issq, outer, len(tl)))
TAT = [(c, a) for c, a in tl]

eq_of = collections.defaultdict(list)
for e, (sq, o, terms) in enumerate(H.eqt):
    for c, a in terms:
        eq_of[a].append(e)

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * E.NV
for k, x in d.items():
    vw[int(k.split('_')[1])] = int(x)
badw = E.badatoms(vw)

say('\n%-6s %-5s %-7s %-9s %s' % ('atom', 'coef', '#eqs', 'value', 'source'))
free_comp = []
for c, a in TAT:
    eqs = eq_of[a]
    val = badw.get(a, 0)
    say('%-6d %-5s %-7d %-9s %s' % (a, c, len(eqs), abs(val).bit_length() if val else 0,
                                    H.atoms[a][:70]))
    if len(eqs) == 1:
        free_comp.append(a)
        say('        *** appears ONLY in eq8680 -> FREE COMPENSATOR')
say('\natoms of T appearing only in eq8680: %s' % free_comp)

say('\n--- which equations would each compensator disturb (beyond eq8680)?')
for c, a in TAT:
    eqs = [e for e in eq_of[a] if e != 8680]
    say('  a%-6d coef %-5s -> %d other equations %s' % (a, c, len(eqs), eqs[:12]))

# how many are E-frame definitions (identically zero) vs genuinely movable
say('\n--- is each atom a DEFINER in E\'s frame (identically satisfied) ?')
definer_atoms = {E.definer[u][0] for u in range(E.NV) if E.definer[u] is not None}
mov = []
for c, a in TAT:
    isdef = a in definer_atoms
    say('  a%-6d definer=%s%s' % (a, isdef, '' if isdef else '   <- a CHECK, can be nonzero'))
    if not isdef:
        mov.append(a)
say('atoms of T that are CHECKS in E\'s frame: %s' % mov)
json.dump({'T_terms': [[str(c), a] for c, a in TAT],
           'free_compensators': free_comp,
           'checks_in_E': mov},
          open(OD + '/tatoms.json', 'w'))
say('DONE')
