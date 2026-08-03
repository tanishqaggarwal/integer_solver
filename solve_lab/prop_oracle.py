#!/usr/bin/env python3
"""Use CORRECT mod-P propagation (dynamic orientation, ModPEngine) as the oracle,
instead of the fixed-orientation v5 forward-eval. For several bit settings report:
 - how many atoms float (a valid state should float few / only genuine checks)
 - x_9770, x_3183, x_18274, x_17728
Compare to what the v5 forward-eval said. Key question: does correct propagation
give a consistent state for nonzero bits, and does x_9770 depend only on 22 bits?"""
import json, time
from modp import ModPEngine, make_base, P, NVARS
from propagate import load_atoms, atom_vars

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
TWIST = {1817, 30378, 40782, 44271}
WATCH = [9770, 3183, 18274, 17728]

def main():
    t0 = time.time()
    atoms = load_atoms()
    control = json.load(open('control_bits.json'))
    bits233 = [b for b in control if b not in set(BITS22)]
    print(f"loaded {len(atoms)} atoms ({time.time()-t0:.0f}s)", flush=True)

    def eval_prop(ones):
        eng = ModPEngine(atoms)
        # set control bits: ones->1, rest->0
        oset = set(ones)
        for b in control:
            eng.assign(b, 1 if b in oset else 0)
        eng.propagate()
        # zero-fill undetermined
        undet = sum(1 for v in range(NVARS) if eng.val[v] is None)
        for v in range(NVARS):
            if eng.val[v] is None: eng.val[v] = 0
        val = eng.val
        tw = ex = 0; exa = []
        for a, poly in enumerate(atoms):
            s = 0
            for m, c in poly.items():
                t = c % P
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            if s % P:
                if a in TWIST: tw += 1
                else:
                    ex += 1
                    if len(exa) < 8: exa.append(a)
        return tw, ex, exa, {w: val[w] for w in WATCH}, undet

    base_w = None
    for label, ones in [('all-0', []), ('233-bit '+str(bits233[0]), [bits233[0]]),
                        ('233-bit '+str(bits233[50]), [bits233[50]]),
                        ('22-bit '+str(BITS22[0]), [BITS22[0]]),
                        ('22-bit '+str(BITS22[5]), [BITS22[5]])]:
        tw, ex, exa, wv, undet = eval_prop(ones)
        if base_w is None: base_w = wv
        chg = [w for w in WATCH if wv[w] != base_w[w]]
        print(f"[{label}] twist={tw} EXTRA={ex} undet={undet} changed_watch={chg} {exa}", flush=True)
        print(f"    x_9770={wv[9770]} x_18274={wv[18274]}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
