"""Agent AD -- the OFF-PIN family, measured separately.

Real instance: 766 off-pin conditions  a'_j*(1-L)*i_j = c'_j*P*u'_j, two per
block, i.e. one per private gadget output slot (i5, i6).

Claim to be tested, not assumed:
  (1) the number of ACTIVE off-pin conditions (those whose gate 1-L is 1) is
      2*(#blocks with L=0) = 2*(n-1-(|S|-1)) = 2*(n-|S|)  --  so it DECREASES
      with |S|, the opposite direction from the congruence family;
  (2) each active off-pin's variable pair (i5,i6 of a dead block) enters no other
      condition, because the only consumer of (i5,i6) is the mux multiplier
      cC = a*b = L, which is 0 exactly when the off-pin fires.  Hence w=0
      satisfies it and the off-pin family closes for EVERY subset at EVERY |S|.

(2) is checked numerically and adversarially: the private outputs of every dead
block are driven to arbitrary nonzero multiples of P and every other condition in
the system is re-evaluated to see whether anything moves.
"""
import json, os, random, sys
from collections import defaultdict
from ad_model import Curve, Instance, Coeffs, CANCEL, DEGEN

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    cs = {}
    for fn in ('ad_curves_partial.json', 'ad_curves.json', 'ad_curves2.json'):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            for k, v in json.load(open(p)).items():
                cs.setdefault(k, []).extend(v)
    return cs


def full_system_values(inst, S, priv, lifts, gate_mode='and'):
    """Evaluate every condition of the small system at an explicit assignment.

    priv[nid] = (i5, i6) integer values of the block's PRIVATE gadget outputs
    lifts[nid] = (t, s) lift of the node's OUTPUT slot
    Returns a list of (tag, value) that must all be 0 / satisfy the divisibility.
    The mux is arithmetic exactly as in the real encoding:
        X_v = cA*X_l + cB*X_r + cC*i5 ,  cA=a(1-b), cB=b(1-a), cC=ab
    """
    p = inst.cv.p
    info, _ = inst.eval_subset(S)
    slot = {}
    out = []
    for v in inst.nodes:
        if v.idx is not None:
            live, _k, pt = info[v.nid]
            slot[v.nid] = (pt[0], pt[1]) if live else (0, 0)
            continue
        la, Pl = info[v.left.nid][0], info[v.left.nid][2]
        lb, Pr = info[v.right.nid][0], info[v.right.nid][2]
        a, b = int(la), int(lb)
        cA, cB, cC = a * (1 - b), b * (1 - a), a * b
        i5, i6 = priv[v.nid]
        Xl, Yl = slot[v.left.nid]
        Xr, Yr = slot[v.right.nid]
        X = cA * Xl + cB * Xr + cC * i5
        Y = cA * Yl + cB * Yr + cC * i6
        slot[v.nid] = (X, Y)
        L = cC if gate_mode == 'and' else (1 if (a or b) else 0)
        if L:
            i1, i2, i3, i4 = Xr, Xl, Yl, Yr
            A = i1 - i2
            B = i4 - i3
            E = i1 + i2 + i5
            N1 = E * A * A - B * B
            N2 = A * (i3 + i6) - B * (i2 - i5)
            out.append(('cong-%d' % v.nid, (N1, N2)))
        else:
            out.append(('offpin-%d' % v.nid, (i5, i6)))
    out.append(('root', slot[inst.root.nid]))
    return out


def main():
    cs = load()
    rng = random.Random(31)
    tot_perturb, moved = 0, 0
    for nb in ('8', '12'):
        for cd in cs[nb][:2]:
            cv = Curve(cd['p'], cd['a'], cd['b'], cd['N'], cd['G'])
            n = int(nb)
            for mode in ('balanced', 'skew'):
                inst = Instance(cv, n, tree_mode=mode)
                cnt = defaultdict(lambda: [0, 0])
                bad_count = 0
                for S in range(1 << n):
                    info, blocks = inst.eval_subset(S)
                    if any(c == CANCEL for _v, c, _d in blocks):
                        continue
                    w = bin(S).count('1')
                    nact = sum(1 for _v, c, _d in blocks if c not in (0, DEGEN)) * 2
                    expect = 2 * (n - 1) if w == 0 else 2 * (n - w)
                    if nact != expect:
                        bad_count += 1
                    cnt[w][0] += nact
                    cnt[w][1] += 1
                # adversarial independence check on a sample of subsets
                for _ in range(300):
                    S = rng.randrange(1 << n)
                    info, blocks = inst.eval_subset(S)
                    if any(c == CANCEL for _v, c, _d in blocks):
                        continue
                    base_priv, base_lift = {}, {}
                    for v in inst.internal:
                        la, lb = info[v.left.nid][0], info[v.right.nid][0]
                        if la and lb:
                            pt = info[v.nid][2]
                            base_priv[v.nid] = (pt[0], pt[1])
                        else:
                            base_priv[v.nid] = (0, 0)
                        base_lift[v.nid] = (0, 0)
                    ref = full_system_values(inst, S, base_priv, base_lift)
                    pert = dict(base_priv)
                    ndead = 0
                    for v in inst.internal:
                        la, lb = info[v.left.nid][0], info[v.right.nid][0]
                        if not (la and lb):
                            ndead += 1
                            pert[v.nid] = (cv.p * rng.randrange(1, 50),
                                           cv.p * rng.randrange(1, 50))
                    if ndead == 0:
                        continue
                    got = full_system_values(inst, S, pert, base_lift)
                    for (t1, v1), (t2, v2) in zip(ref, got):
                        assert t1 == t2
                        if t1.startswith('offpin'):
                            continue        # this IS the perturbed variable
                        tot_perturb += 1
                        if v1 != v2:
                            moved += 1
                print('p=%-6d n=%-3d %-8s  active off-pins vs closed form '
                      '2*(n-|S|): mismatches=%d' % (cd['p'], n, mode, bad_count))
                print('    |S| : mean active off-pin conditions')
                print('    ' + '  '.join('%d:%d' % (w, cnt[w][0] // cnt[w][1])
                                         for w in sorted(cnt)))
    print('\nadversarial independence check: %d non-off-pin condition values '
          'compared after driving every dead block private output to a random '
          'nonzero multiple of P; %d moved.' % (tot_perturb, moved))
    print('OFF-PIN VERDICT: %s' % ('closes for every subset at every |S| '
                                   '(private variables, uncoupled)' if moved == 0
                                   else 'COUPLED -- claim refuted'))


if __name__ == '__main__':
    main()
