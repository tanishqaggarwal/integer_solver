#!/usr/bin/env python3
"""sweep.py -- submodular/supermodular coupler split vs operand width, and the
identification of the supermodular core by variable kind."""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common


def kindname(k):
    return {'word': 'word', 'and': 'AND', 'carry': 'carry', 'adder': 'adder',
            'chunk': 'chunk'}.get(k, k)


def run():
    sizes = [8, 16, 32, 64, 128, 256]
    rows = []
    print(f"{'s':>4} {'vars':>8} {'couplers':>9} {'sub':>8} {'sup':>8} "
          f"{'sub%':>6} {'sup%':>6}  build(s)")
    for s in sizes:
        t = time.time()
        p = common.real_p256() if s == 256 else common.prime_of_bits(s)
        m = common.build_mm(s, p=p)
        Q = m['Q']
        sp = common.coupler_split(Q)
        dt = time.time() - t
        subpct = 100 * sp['sub'] / sp['total']
        row = dict(s=s, vars=Q.n, couplers=sp['total'], sub=sp['sub'],
                   sup=sp['sup'], subpct=subpct, build_s=dt,
                   sup_kinds=sp['sup_kinds'], sub_kinds=sp['sub_kinds'])
        rows.append(row)
        print(f"{s:>4} {Q.n:>8} {sp['total']:>9} {sp['sub']:>8} {sp['sup']:>8} "
              f"{subpct:>5.1f}% {100-subpct:>5.1f}%  {dt:>7.1f}")

    print("\nSUPERMODULAR couplers by variable-kind pair (who is the AND core):")
    for row in rows:
        s = row['s']
        tot = row['sup']
        items = sorted(row['sup_kinds'].items(), key=lambda kv: -kv[1])
        pretty = ", ".join(f"{'x'.join(kindname(x) for x in k)}={v}"
                            f"({100*v/tot:.0f}%)" for k, v in items)
        print(f"  s={s:>4}: {pretty}")

    print("\nSUBMODULAR couplers by variable-kind pair (the min-cut part):")
    for row in rows:
        s = row['s']
        tot = row['sub']
        items = sorted(row['sub_kinds'].items(), key=lambda kv: -kv[1])
        pretty = ", ".join(f"{'x'.join(kindname(x) for x in k)}={v}"
                            f"({100*v/tot:.0f}%)" for k, v in items)
        print(f"  s={s:>4}: {pretty}")

    def strkeys(d):
        return {'x'.join(k): v for k, v in d.items()}
    dump = [dict(r, sup_kinds=strkeys(r['sup_kinds']),
                 sub_kinds=strkeys(r['sub_kinds'])) for r in rows]
    with open(os.path.join(os.path.dirname(__file__), 'sweep.json'), 'w') as f:
        json.dump(dump, f, indent=2)
    return rows


if __name__ == '__main__':
    run()
