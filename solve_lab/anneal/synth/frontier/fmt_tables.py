#!/usr/bin/env python3
"""fmt_tables.py -- read the frontier JSON checkpoints and emit markdown tables."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else {}


def cell(v):
    if v is None:
        return "   -   "
    r = v['rate']
    if v['hits'] == 0:
        return f" 0/{v['n']:<3d} E{v['best']:.0f}".ljust(11)
    return f"{v['hits']}/{v['n']} {r*100:.1f}%".ljust(11)


def solver_s_effort(data, solvers, sizes, efforts, title):
    print(f"\n### {title}\n")
    for eff in efforts:
        print(f"effort = {eff}")
        hdr = "  s   n_vars " + " ".join(f"{s:>12}" for s in solvers)
        print(hdr)
        for s in sizes:
            nv = None
            cells = []
            for solver in solvers:
                v = data.get(f"{s}|{eff}|{solver}")
                if v: nv = v['n_vars']
                cells.append(cell(v))
            print(f"{s:3d} {str(nv):>6}  " + " ".join(f"{c:>12}" for c in cells))
        print()


if __name__ == '__main__':
    core = load('core_baseline.json')
    solver_s_effort(core, ['sa', 'pt', 'tabu', 'sb'], [4, 5, 6, 7, 8],
                    ['low', 'mid', 'high', 'extreme'], "CORE baseline (wallace)")
    combos = load('combos.json')
    if combos:
        solver_s_effort(combos, ['sb_tabu', 'pt_wide'], [6, 7, 8],
                        ['high', 'extreme'], "COMBOS")
    clamp = load('clamp_a.json')
    if clamp:
        solver_s_effort(clamp, ['tabu', 'pt'], [4, 5, 6, 7, 8, 10, 12],
                        ['mid', 'high'], "CLAMP operand a (one operand known)")
    for mode in ('wallace', 'dadda'):
        sq = load(f'squeeze_{mode}.json')
        if sq:
            solver_s_effort(sq, ['tabu', 'sb'], [4, 5, 6, 7, 8],
                            ['mid', 'high'], f"SQUEEZE {mode}")
    wand = load('wand.json')
    if wand:
        print("\n### W_and / mode tuning (tabu)\n")
        print("  key: mode|W_and|s|effort|solver -> hits/n rate best")
        for k in sorted(wand):
            v = wand[k]
            print(f"  {v['mode']:8s} W_and={str(v['W_and']):>5} s={v['s']} "
                  f"{v['effort']:5s}: {v['hits']}/{v['n']} "
                  f"rate={v['rate']*100:5.1f}% best={v['best']:.0f}")
