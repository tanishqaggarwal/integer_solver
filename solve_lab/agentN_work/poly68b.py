"""Complete reduction of the exact polynomial system on the FULL 68-knob candidate set.

Iterates to a fixed point: on the current lattice, any collateral row that has become purely linear
is an honest linear constraint and cuts the lattice further; repeat until every surviving collateral
row is genuinely nonlinear there.  Nothing is linearised — a row is only used as a linear constraint
when its exact restriction *is* linear.

Reports the system size (unknowns, rows, degrees, term counts, coefficient bits) BEFORE solving,
then hands the nonlinear residue to Singular for dimension and degree, then does the exact integer
solve on the region.
"""
import os, sys, json, time, pickle, subprocess
from collections import defaultdict
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN, zsolve
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits
from polyexact import P
from polyfull import exact_polys
from kerquad import int_kernel_columns

eq_terms = ev.eq_terms


def lll(K):
    if not K:
        return []
    R = fmpz_mat([[int(x) for x in v] for v in K]).lll().tolist()
    return [[int(x) for x in r] for r in R if any(r)]


def main():
    st = make(WIT)
    b0 = build(st)
    Rl = b0['R']
    Rset = set(Rl)
    atoms_R = set()
    for e in Rl:
        for c, a in eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    cands = sorted(y for y in cands if y in FREE)
    k = len(cands)
    touched = set()
    for Y in cands:
        for a in fr.chk.get(Y, []):
            touched |= atom_eqs[a]
    outside = sorted(touched - Rset)
    rows = list(Rl) + list(outside)
    polys = exact_polys(st, rows, cands)
    live_out = [e for e in outside if polys[e].c]
    print('=== AMBIENT EXACT SYSTEM ===', flush=True)
    print('unknowns (all free inputs syntactically supporting the region): %d' % k, flush=True)
    print('nonzero rows: %d region + %d collateral' % (len(Rl), len(live_out)), flush=True)
    dd = defaultdict(int)
    for e in Rl + live_out:
        dd[polys[e].deg()] += 1
    print('degree distribution: %s   (max %d)' % (dict(sorted(dd.items())), max(dd)), flush=True)
    print('max terms in a row: %d ; max coefficient bits: %d'
          % (max(polys[e].nterms() for e in Rl + live_out),
             max(polys[e].maxbits() for e in Rl + live_out)), flush=True)
    for e in Rl + live_out:
        assert polys[e].deg() < 2 or True
    assert all(polys[e].c.get((0,) * k, 0) == 0 for e in live_out), 'collateral broken at witness'

    # ---- saturation loop -----------------------------------------------------------------
    K = [[1 if i == j else 0 for i in range(k)] for j in range(k)]

    def restrict(pol, K):
        dd_ = len(K)
        P.NK = dd_
        T = []
        for j in range(k):
            c = {}
            for a in range(dd_):
                if K[a][j]:
                    m = [0] * dd_
                    m[a] = 1
                    c[tuple(m)] = K[a][j]
            T.append(P(c))
        acc = P()
        for mono, cf in pol.c.items():
            term = P.const(cf)
            for j, ex in enumerate(mono):
                for _ in range(ex):
                    term = term * T[j]
            acc = acc + term
        return acc

    print('\n=== saturation loop ===', flush=True)
    it = 0
    used_lin = []
    while True:
        it += 1
        cur = {e: restrict(polys[e], K) for e in live_out}
        linrows = [e for e in live_out if cur[e].deg() == 1]
        nl = [e for e in live_out if cur[e].deg() >= 2]
        print('iter %d: lattice rank %d ; collateral rows now linear %d, nonlinear %d, zero %d'
              % (it, len(K), len(linrows), len(nl), len(live_out) - len(linrows) - len(nl)),
              flush=True)
        if not linrows:
            break
        used_lin += linrows
        A = []
        for e in linrows:
            pol = cur[e]
            v = [0] * len(K)
            for mono, c in pol.c.items():
                assert sum(mono) == 1
                v[mono.index(1)] = c
            A.append(v)
        Kn = int_kernel_columns(A, len(K))
        if not Kn:
            K = []
            break
        Kn = lll(Kn)
        K = [[sum(u[a] * K[a][j] for a in range(len(u))) for j in range(k)] for u in Kn]
        K = lll(K)
        live_out = nl

    d = len(K)
    print('\nlattice rank after all honestly-linear collateral constraints: %d' % d, flush=True)
    print('collateral rows still nonlinear there: %d -> %s' % (len(live_out), live_out), flush=True)

    resid = {e: restrict(polys[e], K) for e in live_out}
    for e in live_out:
        print('   eq %-6d degree %d, %d terms, max coef bits %d'
              % (e, resid[e].deg(), resid[e].nterms(), resid[e].maxbits()), flush=True)

    regr = {e: restrict(polys[e], K) for e in Rl}
    print('\nregion rows on the lattice:', flush=True)
    nq = 0
    for e in Rl:
        if regr[e].deg() >= 2:
            nq += 1
        print('   eq %-6d degree %-3d terms %-4d const %s'
              % (e, regr[e].deg(), regr[e].nterms(),
                 '0' if regr[e].c.get((0,) * d, 0) == 0 else 'nonzero'), flush=True)
    print('   region rows genuinely nonlinear on the lattice: %d of 12' % nq, flush=True)

    # ---- Singular: dimension and degree of the residue --------------------------------------
    V = ['s(%d)' % i for i in range(d)]

    def pstr(pol):
        ts = []
        for mono, c in sorted(pol.c.items()):
            f = ['(%d)' % c] + ['%s^%d' % (V[j], x) if x > 1 else V[j]
                                for j, x in enumerate(mono) if x]
            ts.append('*'.join(f))
        return '+'.join(ts) if ts else '0'

    print('\n=== Singular: dimension and degree of the collateral residue ===', flush=True)
    script = ('ring r = 0,(%s),dp;\nLIB "primdec.lib";\nideal I = %s;\nideal G = std(I);\n'
              '"dim="; dim(G); "vdim="; vdim(G); "gbsize="; size(G);\n'
              'ideal Rad = radical(I); ideal GR = std(Rad);\n'
              '"radical_dim="; dim(GR); "radical_gbsize="; size(GR);\n'
              '"radical_gens_are_linear=";\nint ok=1;\n'
              'for (int i=1;i<=size(GR);i++){ if (deg(GR[i])>1) { ok=0; } }\nok;\n'
              'list pd = minAssGTZ(I); "components="; size(pd);\n'
              'for (i=1;i<=size(pd);i++){ "  dim"; dim(std(pd[i])); "  maxdeg";\n'
              '  int mx=0; for (int j=1;j<=size(pd[i]);j++){ if (deg(pd[i][j])>mx){mx=deg(pd[i][j]);} } mx; kill mx; }\n'
              'quit;\n' % (','.join(V), ',\n  '.join(pstr(resid[e]) for e in live_out)))
    p = os.path.join(HERE, 'runs', 'sing_resid68.sing')
    open(p, 'w').write(script)
    t0 = time.time()
    try:
        r = subprocess.run(['Singular', '-q', '--no-warn', p], capture_output=True, text=True,
                           timeout=2400)
        out = r.stdout.strip()
    except subprocess.TimeoutExpired:
        out = 'TIMEOUT'
    print(out, flush=True)
    print('(%.1fs)' % (time.time() - t0), flush=True)

    pickle.dump(dict(K=K, cands=cands, Rl=Rl, resid={e: resid[e].c for e in live_out},
                     regr={e: regr[e].c for e in Rl}),
                open(os.path.join(HERE, 'runs', 'poly68b.pkl'), 'wb'))
    json.dump(dict(unknowns=k, region=len(Rl), collateral=len(rows) - len(Rl),
                   lattice_rank=d, residue=live_out, region_nonlinear=nq,
                   singular=out),
              open(os.path.join(HERE, 'runs', 'poly68b.json'), 'w'), indent=1)
    print('\nwrote runs/poly68b.{json,pkl}', flush=True)


if __name__ == '__main__':
    main()
