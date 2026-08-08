#!/usr/bin/env python3
"""tradeoff.py -- the runs-vs-V curve for every sound multi-execution scheme.

Reads the measured atoms in multirun/pieces256.json (built at the REAL prime
p = 2^256 - 2^32 - 977) and produces:

  * V_comb(mu, w)   qubits for a comb over mu free scalar bits, window width w
  * V_min(mu)       optimised over w
  * the D1 / S1 trade-off table  (they coincide -- see FINDINGS.md Prop. 4)
  * the Pollard-rho and BSGS crossovers
  * the V < 4400 verdict

Writes multirun/tradeoff.json.
"""
import json, math, os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(_HERE, 'pieces256.json')))

HW = 4400            # D-Wave Advantage2, Zephyr, 4400 qubits
HW_ADV = 5760        # D-Wave Advantage, Pegasus
WMAX = 12            # largest window width measured


def V_comb(mode, mu, w):
    """qubits for a comb over mu free scalar bits with window width w.

    M = ceil(mu/w) windows.  Each window is 2 one-hot table look-ups (x and y);
    consecutive windows are joined by one affine EC addition; the last
    accumulator is compared with T.  Every term is a measured fragment.
    """
    r = D[mode]
    if w > WMAX:
        return None
    M = -(-mu // w)
    return M * 2 * r[f'sel{w}']['vars'] + (M - 1) * r['add']['vars'] + r['final']['vars']


def V_min(mode, mu):
    best = None
    for w in range(1, WMAX + 1):
        v = V_comb(mode, mu, w)
        if v is not None and (best is None or v < best[0]):
            best = (v, w)
    return best


def couplers_comb(mode, mu, w):
    r = D[mode]
    M = -(-mu // w)
    return M * 2 * r[f'sel{w}']['couplers'] + (M - 1) * r['add']['couplers'] + r['final']['couplers']


# ------------------------------------------------------------------ tables
def main():
    out = {'note': 'all costs measured at the real p = 2^256-2^32-977', 'curve': {}}

    MUS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96,
           112, 120, 128, 129, 144, 160, 176, 192, 208, 224, 240, 256]

    for mode in ('binary', 'wallace'):
        print("=" * 108)
        print(f"THE RUNS-vs-SIZE CURVE   (mode={mode}; real p; V = qubits per QPU run)")
        print("  D1 exhaustive interval split  AND  S1 sampling/vOW collision search"
              " lie on this SAME curve")
        print("=" * 108)
        print(f"{'mu':>4} {'best w':>6} {'M':>4} {'V (qubits/run)':>16} {'couplers/run':>15} "
              f"{'runs R=2^(256-mu)':>19} {'R*V':>10} {'fits 4400?':>11}")
        rows = []
        for mu in MUS:
            v, w = V_min(mode, mu)
            M = -(-mu // w)
            c = couplers_comb(mode, mu, w)
            R = 256 - mu
            prod = R + math.log2(v)
            fits = 'YES' if v <= HW else 'no'
            print(f"{mu:4d} {w:6d} {M:4d} {v:16,d} {c:15,d} {'2^%d' % R:>19} "
                  f"{'2^%.1f' % prod:>10} {fits:>11}")
            rows.append(dict(mu=mu, w=w, M=M, V=v, couplers=c, log2_runs=R,
                             log2_RV=round(prod, 2), fits_4400=(v <= HW)))
        out['curve'][mode] = rows
        print()

    # ---------------------------------------------------- the 4400 verdict
    print("=" * 108)
    print("VERDICT ON A 4,400-QUBIT BUDGET")
    print("=" * 108)
    for mode in ('binary', 'wallace'):
        r = D[mode]
        print(f"\n-- {mode} --")
        print(f"  one affine EC addition (s=256)      : {r['add']['vars']:>12,d} qubits "
              f"= {r['add']['vars']/HW:7.1f} x the 4,400-qubit budget")
        print(f"  one 256x256 modular multiplication  : {r['mul']['vars']:>12,d} qubits "
              f"= {r['mul']['vars']/HW:7.1f} x")
        print(f"  one balanced column of a modmul     : {r['col']['vars']:>12,.0f} qubits "
              f"= {r['col']['vars']/HW:7.3f} x   (fits -- but is not self-contained)")
        # largest mu with V <= 4400, table-only regime (M = 1, no EC addition)
        best_mu = 0
        for mu in range(1, WMAX + 1):
            v = 2 * r[f'sel{mu}']['vars'] + r['final']['vars']
            if v <= HW:
                best_mu = mu
        v0 = (2 * r[f'sel{best_mu}']['vars'] + r['final']['vars']) if best_mu else None
        if best_mu:
            print(f"  largest mu with V <= 4,400          : mu = {best_mu}  (V = {v0:,d}; "
                  f"pure table look-up, no group law in the QUBO)")
            print(f"    => runs required                  : 2^{256-best_mu}")
        else:
            print(f"  largest mu with V <= 4,400          : NONE -- even mu = 1 costs "
                  f"{2*r['sel1']['vars']+r['final']['vars']:,d} qubits")
        out.setdefault('verdict', {})[mode] = dict(
            add=r['add']['vars'], mul=r['mul']['vars'], col=r['col']['vars'],
            max_mu_under_4400=best_mu, V_at_max_mu=v0, log2_runs=256 - best_mu)

    # idealised table gadget: 2^mu one-hot bits + one 256-bit comparison
    ideal = max(mu for mu in range(1, 40) if (1 << mu) + 2 * 256 <= HW)
    print(f"\n  Even an IDEALISED table gadget (2^mu one-hot bits + a 512-bit compare,")
    print(f"  i.e. throwing away every compiler overhead) tops out at mu = {ideal} "
          f"for V <= 4,400,")
    print(f"  hence >= 2^{256-ideal} runs.")
    out['idealised_table_mu'] = ideal

    # -------------------------------- do the tiny ones ACTUALLY fit the chip?
    print()
    print("=" * 108)
    print("DO THE QUBOs THAT FIT BY QUBIT COUNT ALSO FIT BY DEGREE AND PRECISION?")
    print("  Pegasus degree 15 / Zephyr degree 20; usable coupler precision ~4-5 bits")
    print("=" * 108)
    print(f"{'mode':>8} {'mu':>3} {'V':>8} {'couplers':>10} {'avg deg':>8} "
          f"{'|J| bits (max frag)':>21} {'embeds on Zephyr?':>18}")
    fit = []
    for mode in ('binary', 'wallace'):
        r = D[mode]
        for mu in range(1, 7):
            v = 2 * r[f'sel{mu}']['vars'] + r['final']['vars']
            c = 2 * r[f'sel{mu}']['couplers'] + r['final']['couplers']
            deg = 2.0 * c / v
            jb = max(r[f'sel{mu}']['jbits'], r['final']['jbits'])
            ok = 'yes (deg ok)' if deg <= 20 else f'no, deg {deg:.0f} > 20'
            print(f"{mode:>8} {mu:3d} {v:8,d} {c:10,d} {deg:8.1f} {'2^%d' % jb:>21} {ok:>18}")
            fit.append(dict(mode=mode, mu=mu, V=v, couplers=c, avg_degree=round(deg, 1),
                            jbits=jb))
    out['hardware_fit'] = fit

    # ---------------------------------------------------- classical baselines
    print()
    print("=" * 108)
    print("CROSSOVERS AGAINST THE CLASSICAL BASELINES")
    print("=" * 108)
    print("  Pollard rho : 2^128 group operations, O(1) memory")
    print("  BSGS        : 2^128 group operations, 2^128 stored points")
    print()
    for mode in ('binary', 'wallace'):
        v128, w128 = V_min(mode, 128)
        v129, w129 = V_min(mode, 129)
        print(f"  {mode:>7}: R <= 2^128 needs mu >= 128  ->  V = {v128:,d} qubits/run "
              f"(w={w128})   [{v128/HW:,.0f} x hardware]")
        print(f"  {mode:>7}: strictly beating rho (mu = 129) ->  V = {v129:,d} qubits/run "
              f"(w={w129})   [{v129/HW:,.0f} x hardware]")
        out.setdefault('crossover', {})[mode] = dict(V_mu128=v128, w=w128, V_mu129=v129)

    # wall-clock-corrected crossover
    ANNEAL_RATE = 1e4       # samples/s, generous sustained figure for a QPU
    GROUPOP_RATE = 1e7      # EC group operations/s/core, ordinary CPU
    shift = math.log2(GROUPOP_RATE / ANNEAL_RATE)
    mu_t = 128 + shift
    print(f"\n  Wall-clock-corrected: at {ANNEAL_RATE:.0e} anneals/s vs "
          f"{GROUPOP_RATE:.0e} group ops/s/core,")
    print(f"  one anneal costs 2^{shift:.1f} group ops, so beating rho on TIME needs "
          f"mu >= {mu_t:.0f}:")
    for mode in ('binary', 'wallace'):
        v, w = V_min(mode, int(math.ceil(mu_t)))
        print(f"    {mode:>7}: V = {v:,d} qubits/run  [{v/HW:,.0f} x hardware]")
        out['crossover'][mode]['V_wallclock'] = v
    out['crossover']['wallclock_mu'] = mu_t

    json.dump(out, open(os.path.join(_HERE, 'tradeoff.json'), 'w'), indent=1)
    print("\nwrote multirun/tradeoff.json")


if __name__ == '__main__':
    main()
