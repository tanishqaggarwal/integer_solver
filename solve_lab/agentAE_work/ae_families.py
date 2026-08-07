#!/usr/bin/env python3
"""agent AE -- the structured-key family sweep.

Each family is a set F of scalars; the run decides (probabilistically, see
STRUCTURED_KEYS.md) whether k0 in F.  Every family is realised as an interval
[lo, lo+2^R) mod N against a target point Q_fam, so one kangaroo run per family.

A HIT is never reported from here: candidates are re-verified in Python bignum
and written to HIT_*.json for a second, independent implementation to confirm.
"""
import sys, os, json, time, math
import ae_lib as L

N = L.N; T = L.T; G = L.G

def two_sided(c, R):
    """interval of width 2^R centred on c"""
    return (c - (1 << (R - 1))) % N

def fam_interval(name, Qfam, lo, R, note=''):
    return dict(kind='interval', name=name, Q=Qfam, lo=lo % N, R=R, note=note)

def build(tier, R_const=48, R_orbit=52, R_window=40):
    fams = []
    A = lambda *x: fams.append(fam_interval(*x))

    if tier == 'head':
        A('mag_small_R64', T, 0, 64, 'k0 < 2^64')

    if tier == 'orbit':
        R = R_orbit
        A('mag_top', T, (N - (1 << R)) % N, R, 'N - k0 < 2^R  (i.e. -k0 small)')
        P1 = L.phi(T); P2 = L.phi(P1)
        A('orb_lam',   P1, 0, R, 'lam*k0 < 2^R')
        A('orb_lam_neg', L.neg(P1), 0, R, '-lam*k0 < 2^R')
        A('orb_lam2',  P2, 0, R, 'lam^2*k0 < 2^R')
        A('orb_lam2_neg', L.neg(P2), 0, R, '-lam^2*k0 < 2^R')

    if tier == 'const':
        R = R_const
        ones = (1 << 256) - 1
        C = []
        for d in (2, 3, 4, 5, 6, 7, 8, 10, 16, 100):
            C.append(('N_over_%d' % d, N // d))
        C.append(('2N_over_3', (2 * N) // 3))
        C.append(('3N_over_4', (3 * N) // 4))
        for d in (3, 5, 7, 11, 13):
            C.append(('inv_%d' % d, pow(d, -1, N)))
        for e in (16, 32, 64, 96, 128, 160, 192, 224, 250, 254, 255):
            C.append(('pow2_%d' % e, (1 << e) % N))
        C.append(('two256_modN', (1 << 256) % N))
        C.append(('ones256_modN', ones % N))
        C.append(('p_modN', L.p % N))
        C.append(('lam', L.lam))
        C.append(('lam2', pow(L.lam, 2, N)))
        C.append(('neg_lam', (-L.lam) % N))
        C.append(('neg_lam2', (-pow(L.lam, 2, N)) % N))
        C.append(('beta_modN', L.beta % N))
        # rational multiples of the all-ones constant: 0x5555.., 0xaaaa.., 0x3333.., 0xcccc..
        C.append(('ones_over_3', (ones // 3) % N))
        C.append(('ones_2over3', (2 * (ones // 3)) % N))
        C.append(('ones_over_5', (ones // 5) % N))
        C.append(('ones_4over5', (4 * (ones // 5)) % N))
        C.append(('ones_over_15', (ones // 15) % N))
        C.append(('ones_over_17', (ones // 17) % N))
        C.append(('ones_over_257', (ones // 257) % N))
        # dedupe centres that fall in the same 2^R window
        seen = []
        for nm, c in C:
            dup = None
            for nm2, c2 in seen:
                dd = (c - c2) % N
                if dd < (1 << R) or (N - dd) < (1 << R): dup = nm2; break
            if dup: print('  (skip %s: same 2^%d window as %s)' % (nm, R, dup)); continue
            seen.append((nm, c))
            A('c_' + nm, T, two_sided(c, R), R, 'k0 within 2^%d of %s' % (R - 1, nm))

    if tier == 'window':
        R = R_window
        inv2 = pow(2, -1, N)
        P = T
        for s in range(256):
            # target 2^-s T has log 2^-s k0 ; interval centred on 0 covers +-2^(R-1)
            A('win_s%03d' % s, P, (N - (1 << (R - 1))) % N, R,
              'k0 = a*2^%d mod N with |a| < 2^%d' % (s, R - 1))
            P = L.mul(inv2, P)
    return fams

def run_tier(tier, threads=1, kpt=None, log2max=None, tablebits=21, seed=20250807,
             R_const=48, R_orbit=52, R_window=40):
    fams = build(tier, R_const, R_orbit, R_window)
    print('tier %s: %d families' % (tier, len(fams)))
    out = []
    hits = []
    for i, f in enumerate(fams):
        R = f['R']
        Q = L.sub(f['Q'], L.mulG(f['lo'])) if f['lo'] else f['Q']
        assert L.oncurve(Q)
        k = kpt if kpt else max(32, min(1024, 1 << max(5, (R // 2) - 12)))
        lm = log2max if log2max else int(R / 2.0 + 3)
        t0 = time.time()
        res = L.run_kangaroo('fam_' + f['name'], Q, R, threads=threads, kpt=k,
                             log2max=lm, seed=seed + i, tablebits=tablebits, quiet=True)
        good = L.check_cands(res, Q)
        rec = dict(name=f['name'], note=f['note'], R=R, lo=str(f['lo']),
                   jumps=res.get('jumps'), dps=res.get('dps'), dpbits=res['dpbits'],
                   log2max=lm, K=res['K'], rc=res['rc'], secs=round(time.time() - t0, 1),
                   dxzero=res.get('dxzero'), done=res.get('done'))
        if res.get('dps') and res.get('jumps'):
            rec['dp_over_closedform'] = round(res['dps'] / (res['jumps'] / 2.0 ** res['dpbits']), 4)
        rec['hit'] = bool(good)
        if good:
            ks = [(f['lo'] + g) % N for g in good]
            ks = [x for x in ks if L.mulG(x) == f['Q']]
            rec['candidate_k'] = [str(x) for x in ks]
            hits.append(rec)
            print('*** HIT in family %s: %r' % (f['name'], rec['candidate_k']))
            json.dump(rec, open('HIT_%s.json' % f['name'], 'w'), indent=1)
        out.append(rec)
        ok = 'OK' if res['rc'] == 0 else 'RC%d' % res['rc']
        print('  %-22s R=%2d jumps=%-12s dps=%-8s dp/cf=%s %5.1fs %s %s' %
              (f['name'], R, rec['jumps'], rec['dps'], rec.get('dp_over_closedform'),
               rec['secs'], ok, 'HIT' if good else 'miss'))
        sys.stdout.flush()
        try: os.remove(os.path.join(L.HERE, 'in_fam_%s.txt' % f['name']))
        except OSError: pass
    json.dump(out, open('res_%s.json' % tier, 'w'), indent=1)
    nbad = sum(1 for r in out if r['rc'] != 0)
    print('TIER %s COMPLETE: %d families, %d nonzero-rc, %d hits' % (tier, len(out), nbad, len(hits)))
    return out

if __name__ == '__main__':
    tier = sys.argv[1]
    kw = {}
    for arg in sys.argv[2:]:
        k, v = arg.split('='); kw[k] = int(v)
    run_tier(tier, **kw)
