#!/usr/bin/env python3
"""Trace the definition cascade of x_18274 and x_17728 to reveal the number
representation used by the 233-bit residue-load side. Also inspect what the
current best_partial sets on control bits and watch vars."""
import json
from confluent_eval5 import build5
from propagate import NVARS, atom_vars

def main():
    A, kind, info, seq, bestval, ncyc = build5()
    control = set(json.load(open('control_bits.json')))
    BITS22 = {1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116}

    # what does best set on control bits?
    best = json.load(open('best/best_partial_39019.json'))
    bv = {int(k[2:]): v for k, v in best.items()}
    cvals = {}
    for b in control:
        cvals[b] = bv.get(b, 0)
    setbits = [b for b in control if cvals[b] != 0]
    print(f"best_partial sets {len(setbits)} of {len(control)} control bits to nonzero")
    print(f"  those values: {sorted(set(cvals[b] for b in control))[:10]}")
    print(f"  set control bits: {sorted(setbits)[:40]}")

    def show(v, depth, seen):
        if v in seen or depth > 14:
            return
        seen.add(v)
        k = kind.get(v)
        pad = '  ' * depth
        if k == 'load':
            bit, cbx, lt = info[v]
            huge = [(c, m) for c, m in lt if len(m) <= 1 and abs(c) > 10**30]
            print(f"{pad}x_{v} = LOAD bit=x_{bit} cbx={cbx} terms={len(lt)} huge={[c for c,m in huge][:2]}")
            for c, m in lt:
                for x in m:
                    if x != bit: show(x, depth+1, seen)
        elif k == 'div':
            c, u, rest = info[v]
            print(f"{pad}x_{v} = DIV c={c} u=x_{u} rest_terms={len(rest)}")
            show(u, depth+1, seen)
            # show a couple rest vars
            rv = set()
            for cc, m in rest:
                rv.update(m)
            for x in sorted(rv)[:4]:
                show(x, depth+1, seen)
        elif k == 'gate':
            coef, terms = info[v]
            tv = set()
            for cc, m in terms: tv.update(m)
            bigc = [cc for cc, m in terms if abs(cc) > 10**30]
            print(f"{pad}x_{v} = GATE coef={coef} nterms={len(terms)} nvars={len(tv)} big={bigc[:2]} incontrol={sorted(tv&control)[:6]}")
            for x in sorted(tv)[:5]:
                show(x, depth+1, seen)
        else:
            tag = 'CTRL-BIT' if v in control else ('const' if k=='const' else 'input')
            print(f"{pad}x_{v} [{tag}]")

    print("\n===== x_18274 cascade =====")
    show(18274, 0, set())
    print("\n===== x_17728 cascade =====")
    show(17728, 0, set())

if __name__ == '__main__':
    main()
