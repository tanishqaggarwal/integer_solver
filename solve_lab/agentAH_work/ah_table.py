#!/usr/bin/env python3
"""Score every AH closure with checker.py's own loader/evaluator and emit the landscape."""
import sys, os, json, glob, time
sys.set_int_max_str_digits(50_000_000)
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import checker
AH = os.path.dirname(os.path.abspath(__file__))

FOOT15 = [4573,7123,7469,9648,11854,16622,17726,21382,25539,28653,29437,31061,32894,32916,34517]

def main():
    t0 = time.time()
    codes, _ = checker.load_equations()
    NEQ = len(codes)
    print('[ah_table] loaded %d equations in %.1fs' % (NEQ, time.time()-t0), flush=True)
    rows = []
    for mp in sorted(glob.glob(os.path.join(AH, 'meta_*.json'))):
        m = json.load(open(mp))
        jp = m['json']
        if not os.path.exists(jp):
            continue
        v = checker.load_assignment(jp)
        fails = checker.evaluate_all(codes, v)
        m['score'] = NEQ-len(fails); m['nfail'] = len(fails); m['fails'] = fails
        m['foot_is_15'] = (fails == FOOT15)
        rows.append(m)
    rows.sort(key=lambda r: (r['n'], r['seed']))
    json.dump(rows, open(os.path.join(AH, 'landscape.json'), 'w'), indent=1)
    print()
    print('%-5s %-6s %6s %6s %-22s %7s %5s %6s %6s %s' %
          ('|S|','seed','score','nfail','reason','wall','out','nzatm','smpMiss','footprint'))
    for r in rows:
        fp = 'THE-15' if r['foot_is_15'] else ('%d eqs: %s' % (r['nfail'], r['fails'][:8]))
        print('%-5d %-6d %6d %6d %-22s %7.0f %5d %6d %6d  %s' %
              (r['n'], r['seed'], r['score'], r['nfail'], r['reason'], r['wall'],
               r['outer_reached'], r['nz_atoms'], r.get('sampled_root_misses', -1), fp))

main()
