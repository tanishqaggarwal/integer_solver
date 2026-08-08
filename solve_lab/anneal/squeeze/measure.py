#!/usr/bin/env python3
"""measure.py -- ONE 256-bit modular multiplication, every variant, measured.

Reports the triple the hardware actually cares about:

    (logical qubits, max clique, |J| dynamic range in bits)

with  physical ~ logical * max(1, clique/6)  (Pegasus calibration from embed.py).

Usage:  python3 measure.py <group> [group ...]
        groups: base reduce kara toom carry sq all
Results accumulate in results.json.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmqb import MMQB                                            # noqa: E402
from mm import build_modmul                                      # noqa: E402

P = 2 ** 256 - 2 ** 32 - 977
S = 256
HW = 4400
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results.json')


def phys(v, c):
    return v * max(1.0, c / 6.0)


def modmul(mult='schoolbook', leaf=32, red='naf', mode='wallace', square=False,
           chunk=16, dadda_height=2, s=S, p=P):
    """marginal cost of ONE modular multiplication: A, B, C already exist."""
    Q = MMQB(chunk=chunk, mode=mode, dadda_height=dadda_height)
    A = Q.mkword('A', s, lambda wv: wv['_a'])
    B = A if square else Q.mkword('B', s, lambda wv: wv['_b'])
    C = Q.mkword('C', s, lambda wv: wv['_c'])
    base = Q.n
    t0 = time.time()
    build_modmul(Q, p, A, B, C, mult=mult, leaf=leaf, red=red)
    Q.finalize()
    st = Q.stats()
    return dict(vars=Q.n - base, total_vars=Q.n, clique=st['max_clique'],
                jbits=st['dynamic_range_bits'], couplers=st['couplers'],
                ands=st['and_vars'], carries=st['carry_bits'],
                words=st['word_bits'] - (s if square else 2 * s) - s,
                and_lookups=Q.and_lookups, and_hits=Q.and_hits,
                phys=phys(Q.n - base, st['max_clique']), secs=round(time.time() - t0, 1))


def load():
    if os.path.exists(RES):
        return json.load(open(RES))
    return {}


def save(d):
    json.dump(d, open(RES, 'w'), indent=1, sort_keys=True)


HDR = f"{'variant':>44} {'logical':>10} {'K':>5} {'|J|':>6} {'physical':>12} {'vs4400':>8} {'AND':>8} {'carry':>9}"


def show(key, r):
    print(f"{key:>44} {r['vars']:10,d} {r['clique']:5d} 2^{r['jbits']:<4d} "
          f"{r['phys']:12,.0f} {r['phys']/HW:7,.0f}x {r['ands']:8,d} {r['carries']:9,d}")


def run(key, **kw):
    d = load()
    if key in d:
        show(key, d[key])
        return d[key]
    r = modmul(**kw)
    d = load()
    d[key] = r
    save(d)
    show(key, r)
    return r


def main(groups):
    print(HDR)
    print("-" * len(HDR))
    if 'base' in groups or 'all' in groups:
        # 'unary' is omitted at s=256 on purpose: its cliques are ~500 wide, so
        # the penalty expansion is ~c^2 per column and the build exhausts memory
        # before it can be measured.  See the 'small' group for what it costs.
        for mode in ('binary', 'wallace', 'dadda'):
            run(f"BASELINE school/quotient/{mode}", mult='schoolbook',
                red='quotient', mode=mode)
    if 'reduce' in groups or 'all' in groups:
        for mode in ('binary', 'wallace'):
            for red in ('quotient', 'naf', 'fold'):
                run(f"school/{red}/{mode}", mult='schoolbook', red=red, mode=mode)
    if 'kara' in groups or 'all' in groups:
        for mode in ('wallace', 'binary'):
            for leaf in (128, 64, 32, 16, 8, 4):
                run(f"karatsuba(leaf={leaf})/naf/{mode}", mult='karatsuba',
                    leaf=leaf, red='naf', mode=mode)
    if 'toom' in groups or 'all' in groups:
        for mode in ('wallace', 'binary'):
            for leaf in (86, 64, 32, 16, 8):
                run(f"toom3(leaf={leaf})/naf/{mode}", mult='toom3',
                    leaf=leaf, red='naf', mode=mode)
    if 'carry' in groups or 'all' in groups:
        for mode in ('wallace', 'dadda', 'unary'):
            run(f"school/naf/{mode}", mult='schoolbook', red='naf', mode=mode)
        for h in (2, 3, 4, 6):
            run(f"school/naf/dadda(h={h})", mult='schoolbook', red='naf',
                mode='dadda', dadda_height=h)
        for ch in (4, 8, 16, 64, 1024):
            run(f"school/naf/binary(chunk={ch})", mult='schoolbook', red='naf',
                mode='binary', chunk=ch)
    if 'hybrid' in groups or 'all' in groups:
        for leaf in (20, 24, 28, 40, 48):
            run(f"karatsuba(leaf={leaf})/naf/wallace", mult='karatsuba',
                leaf=leaf, red='naf', mode='wallace')
        for top, lf in ((86, 32), (86, 16), (86, 24), (64, 32), (64, 16), (64, 24)):
            run(f"toom3(>{top})+karatsuba({lf})/naf/wallace",
                mult=[('toom3', top), ('karatsuba', lf)], red='naf', mode='wallace')
        for top, lf in ((86, 32), (86, 16), (64, 24)):
            run(f"toom3(>{top})+karatsuba({lf})/naf/binary",
                mult=[('toom3', top), ('karatsuba', lf)], red='naf', mode='binary')
        for top, lf in ((128, 32), (128, 16)):
            run(f"karatsuba(>{top})+toom3({lf})/naf/wallace",
                mult=[('karatsuba', top), ('toom3', lf)], red='naf', mode='wallace')
    if 'small' in groups or 'all' in groups:
        from crossover import PM
        for s in (32, 64):
            for mode in ('binary', 'wallace', 'dadda', 'unary'):
                run(f"s={s} school/naf/{mode}", mult='schoolbook', red='naf',
                    mode=mode, s=s, p=PM[s])
    if 'sq' in groups or 'all' in groups:
        for mode in ('wallace', 'binary'):
            run(f"SQUARE school/naf/{mode}", mult='schoolbook', red='naf',
                mode=mode, square=True)
            run(f"SQUARE karatsuba(leaf=16)/naf/{mode}", mult='karatsuba',
                leaf=16, red='naf', mode=mode, square=True)


if __name__ == '__main__':
    main(sys.argv[1:] or ['all'])
