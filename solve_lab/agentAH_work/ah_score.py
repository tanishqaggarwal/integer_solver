#!/usr/bin/env python3
"""Batch scorer: solve_lab/checker.py's OWN loader and evaluator, loaded once,
applied to many assignment JSONs.  Digit cap raised exactly as agentE_work/verifyE.py does.
Nothing here reimplements the scoring; `checker.load_equations` and `checker.evaluate_all`
are called verbatim."""
import sys, os, json, time
sys.set_int_max_str_digits(50_000_000)
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import checker

def main():
    paths = sys.argv[1:]
    t0 = time.time()
    codes, _ = checker.load_equations()
    n = len(codes)
    print('[ah_score] loaded %d equations in %.1fs' % (n, time.time()-t0), flush=True)
    out = {}
    for pth in paths:
        try:
            v = checker.load_assignment(pth)
        except Exception as e:
            print('%-40s LOAD-ERROR %r' % (os.path.basename(pth), e), flush=True)
            continue
        t1 = time.time()
        fails = checker.evaluate_all(codes, v)
        rec = {'path': pth, 'score': n-len(fails), 'total': n, 'nfail': len(fails),
               'fails': fails, 'eval_s': round(time.time()-t1, 1)}
        out[os.path.basename(pth)] = rec
        print('%-34s %5d/%d  nfail=%-4d %s' % (os.path.basename(pth), rec['score'], n,
              len(fails), fails[:20]), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
              'scores_%d.json' % int(time.time())), 'w'), indent=1)

main()
