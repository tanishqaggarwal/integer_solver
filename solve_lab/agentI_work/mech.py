#!/usr/bin/env python3
"""MECHANISM TEST (corrected).

The accumulator is ADVICE, not a derived value: switching on extra selectors
RELEASES (x1,y1) instead of computing it.  So the correct test is:

 A. with exactly two selectors i,j on, the circuit's own (x1,y1),(x2,y2) are the
    two ladder points P_i,P_j, its X35389/X6671 equal the addition-law residuals
    of `P_i (+) P_j  vs  T`, and they vanish IFF P_i (+) P_j == T.
 B. with three on, (x1,y1) is released; pre-assigning it to the EC partial sum
    closes the newly created gadget exactly.

If the circuit computed anything other than the ladder, A fails.
"""
import os, re, sys, collections, random, json
from boolscore import Fast
from fp import P
HERE = os.path.dirname(os.path.abspath(__file__))
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
c3 = K * pow(3, -1, P) % P
BC = 64019533680030876408443198762210829058751700634554282185987325820393598524794


def eadd(A, B):
    if A is None:
        return B
    if B is None:
        return A
    xa, ya = A; xb, yb = B
    if xa == xb and (ya + yb) % P == 0:
        return None
    l = (3 * xa * xa % P * pow(2 * ya % P, -1, P)) % P if A == B \
        else ((yb - ya) * pow(xb - xa, -1, P)) % P
    xc = (l * l - xa - xb) % P
    return (xc, (l * (xa - xc) - ya) % P)


F = Fast(); M = F.M
sel = collections.defaultdict(list)
for i, s in enumerate(M.src):
    m = re.match(r'^X(\d+) \* \(X(\d+) - (\d+)\)', s)
    if m and int(m.group(3)) > 2**200:
        sel[int(m.group(1))].append((int(m.group(2)), int(m.group(3)) % P))
ks = sorted(sel)
pts = {}
for b in ks:
    (v1, a1), (v2, a2) = sel[b]
    for q in [((a1 + c3) % P, a2), ((a2 + c3) % P, a1),
              ((a1 + c3) % P, (-a2) % P), ((a2 + c3) % P, (-a1) % P)]:
        if (q[1] * q[1] - pow(q[0], 3, P) - BC) % P == 0:
            pts[b] = q
            break


def run(onset):
    def pol(u, roots, S=onset):
        if u in S:
            return 1 if 1 in roots else roots[0]
        x = F.witp[u]
        return x if x in roots else roots[0]
    return F.run(pol)


def formulas(x1, y1, x2, y2, x3, y3):
    A = ((x2 - x1)**2 * (x3 + x1 + x2 + K) - (y2 - y1)**2) % P
    B = ((y3 + y1) * (x2 - x1) - (y2 - y1) * (x1 - x3)) % P
    return A, B


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    okA = totA = 0
    lines = []
    for t in range(R):
        i, j = rng.sample(ks, 2)
        val, conf, _ = run({i, j})
        need = [12186, 16742, 14853, 24908, 22162, 30213]
        if any(val[v] is None for v in need):
            lines.append(f"  pair (X{i},X{j}): coordinates not all determined -- skipped")
            continue
        x1, y1, x2, y2, x3, y3 = (val[v] for v in need)
        totA += 1
        Pi, Pj = pts[i], pts[j]
        circ = {((x1 + c3) % P, y1), ((x2 + c3) % P, y2)}
        ptsok = circ == {Pi, Pj}
        A, B = formulas(x1, y1, x2, y2, x3, y3)
        fok = (A == val[35389] % P) and (B == val[6671] % P)
        T = ((x3 + c3) % P, y3)
        S = eadd(Pi, Pj)
        pred_zero = (S == T)
        act_zero = (A == 0 and B == 0)
        iff = (pred_zero == act_zero)
        good = ptsok and fok and iff
        okA += good
        lines.append(f"  pair (X{i},X{j}): coords==ladder pts:{ptsok}  "
                     f"formulas match circuit:{fok}  "
                     f"(A==B==0) <=> (Pi+Pj==T): {iff} [both {act_zero}]  "
                     f"conflicts={len(conf)}  {'OK' if good else 'MISMATCH'}")
    print(f"PART A: {okA}/{totA} random selector pairs behave exactly as "
          f"`P_i (+) P_j =? T` on the curve")
    for l in lines:
        print(l)

    # PART B: three selectors on -> (x1,y1) released; set it to the partial sum
    okB = totB = 0
    linesB = []
    for t in range(4):
        trio = rng.sample(ks, 3)
        val, conf, _ = run(set(trio))
        if val[12186] is not None:
            linesB.append(f"  trio {trio}: (x1,y1) not released -- skipped")
            continue
        # which two are absorbed into the accumulator?  try all 3 choices
        hit = None
        for drop in range(3):
            rest = [trio[k] for k in range(3) if k != drop]
            S = eadd(pts[rest[0]], pts[rest[1]])
            if S is None:
                continue
            pre = {12186: (S[0] - c3) % P, 16742: S[1]}

            def pol(u, roots, SS=set(trio)):
                if u in SS:
                    return 1 if 1 in roots else roots[0]
                x = F.witp[u]
                return x if x in roots else roots[0]
            v2, c2, _ = F.run(pol, preassign=pre)
            if len(c2) <= 3:
                hit = (rest, len(c2), c2[:6])
                break
        totB += 1
        if hit:
            okB += 1
            linesB.append(f"  trio {trio}: setting (x1,y1) := P_{hit[0][0]} (+) P_{hit[0][1]} "
                          f"leaves {hit[1]} violated atoms {hit[2]} "
                          f"(the single remaining rung)")
        else:
            linesB.append(f"  trio {trio}: no partial-sum assignment closed the new gadget")
    print(f"PART B: {okB}/{totB} triples close the newly created rung with the EC partial sum")
    for l in linesB:
        print(l)
    json.dump({'A_ok': okA, 'A_tot': totA, 'B_ok': okB, 'B_tot': totB},
              open(os.path.join(HERE, 'mech_results.json'), 'w'))


if __name__ == '__main__':
    main()
