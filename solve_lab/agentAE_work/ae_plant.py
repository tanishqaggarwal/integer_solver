#!/usr/bin/env python3
"""agent AE -- planted-key calibration of the kangaroo engine.

Every plant is a key k drawn uniformly from [0,2^R); a PASS requires the engine
to return a candidate that an independent Python bignum recomputation confirms
satisfies k*G == Q, and that equals the planted k.  The measured statistic is
jumps/sqrt(2^R), which is compared to the closed-form expectation ~2.
"""
import sys, random, statistics, json, os
import ae_lib as L

def calib(R, n, threads, kpt, seed0=1000, log2max=None):
    rng = random.Random(seed0 ^ R)
    ratios = []; fails = 0; wrong = 0; dpchk = []
    sq = 2.0 ** (R / 2.0)
    for t in range(n):
        k = rng.randrange(1, 1 << R)
        Q = L.mulG(k)
        res = L.run_kangaroo('plant%d_%d' % (R, t), Q, R, threads=threads, kpt=kpt,
                             seed=rng.randrange(1, 1 << 62),
                             log2max=log2max or int(R / 2.0 + 4), quiet=True)
        good = L.check_cands(res, Q)
        if not good:
            fails += 1
            print('  plant %d R=%d: NO HIT  jumps=%s' % (t, R, res.get('jumps')))
            continue
        if k not in good:
            wrong += 1
            print('  plant %d R=%d: hit but != planted k  (%r vs %d)' % (t, R, good, k))
            continue
        ratios.append(res['jumps'] / sq)
        if res.get('dps'):
            dpchk.append(res['dps'] / (res['jumps'] / float(1 << res['dpbits'])))
    return ratios, fails, wrong, dpchk

if __name__ == '__main__':
    R = int(sys.argv[1]); n = int(sys.argv[2])
    threads = int(sys.argv[3]); kpt = int(sys.argv[4])
    ratios, fails, wrong, dpchk = calib(R, n, threads, kpt)
    print('R=%d n=%d threads=%d kpt=%d  K=%d' % (R, n, threads, kpt, threads * kpt))
    print('  recovered exactly : %d / %d   (no-hit %d, wrong %d)' % (len(ratios), n, fails, wrong))
    if ratios:
        print('  jumps / sqrt(2^R): mean %.3f  median %.3f  min %.3f  max %.3f' %
              (statistics.mean(ratios), statistics.median(ratios), min(ratios), max(ratios)))
        rs = sorted(ratios)
        print('  quantiles 50/75/90/max: %.2f %.2f %.2f %.2f' %
              (rs[len(rs) // 2], rs[int(.75 * len(rs))], rs[min(len(rs) - 1, int(.9 * len(rs)))], rs[-1]))
    if dpchk:
        print('  DP count / closed form (jumps/2^dpbits): mean %.4f  (expect 1.0)' % statistics.mean(dpchk))
    json.dump(dict(R=R, n=n, threads=threads, kpt=kpt, ratios=ratios, fails=fails, wrong=wrong,
                   dpchk=dpchk), open('calib_R%d.json' % R, 'w'))
