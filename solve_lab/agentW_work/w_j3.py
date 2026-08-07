"""TASK 1a: the 21 triples O never reached, j=3, b<=2, exhaustive per triple.
Same method as O's fb_j3.py (exact integer solve_sparse), same setup, complement of O's 14."""
import sys, os, itertools, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

LOG = open(os.path.join(S.HERE, 'w_j3.log'), 'w', buffering=1)
def say(*a): print(*a, file=LOG, flush=True)

CAP = float(os.environ.get('WCAP', '86400'))
# O's own T_COMPENSATION.md lists exactly these 14 triples as completed at b<=2.
# NB: that is 14 of the 15 triples containing eq12231 -- [12231,22044,29125] was NOT reached.
O14 = [[12231,12270,12350],[12231,12270,14584],[12231,12270,18673],[12231,12270,22044],
       [12231,12270,29125],[12231,12350,14584],[12231,12350,18673],[12231,12350,22044],
       [12231,12350,29125],[12231,14584,18673],[12231,14584,22044],[12231,14584,29125],
       [12231,18673,22044],[12231,18673,29125]]
TODO = [P for P in itertools.combinations(S.FAIL, 3) if list(P) not in O14]
assert len(TODO) == 21, len(TODO)
say('knobs %d rows %d SAT %d FAIL %d' % (len(S.KNOB), len(S.names), len(S.SAT), len(S.FAIL)))
say('triples to do: %d (the ones O did not reach)' % len(TODO))
for P in TODO: say('   %s' % list(P))
say('')

t0 = time.time(); nsolve = 0; win = None; done = []
for P in TODO:
    if win or time.time() - t0 > CAP: break
    tp = time.time()
    if S.solve(list(P)) is None:
        nsolve += 1
        say(' triple %s: infeasible on its own -- skipped' % list(P)); done.append((P,'infeasible-alone')); continue
    nsolve += 1
    s = S.solve(S.SAT + list(P)); nsolve += 1
    if s is not None:
        win = S.price(s, 'buy %s break 0' % list(P), log=LOG, tagfile='j3'); done.append((P,'b=0 feasible')); break
    hit = False
    for r in S.SAT:
        s = S.solve([e for e in S.SAT if e != r] + list(P)); nsolve += 1
        if s is not None:
            win = S.price(s, 'buy %s break eq%s' % (list(P), r), log=LOG, tagfile='j3'); hit = True; break
    if hit:
        done.append((P,'b=1 feasible')); break
    for r1, r2 in itertools.combinations(S.SAT, 2):
        if time.time() - t0 > CAP: break
        s = S.solve([e for e in S.SAT if e != r1 and e != r2] + list(P)); nsolve += 1
        if s is not None:
            win = S.price(s, 'buy %s break eq%s,eq%s' % (list(P), r1, r2), log=LOG, tagfile='j3'); hit = True; break
    done.append((P, 'b<=2 exhausted' if not hit else 'b=2 feasible'))
    say(' triple %s: b<=2 %s  (%d solves so far, %.0fs, total %.0fs)'
        % (list(P), 'FEASIBLE' if hit else 'exhausted, none', nsolve, time.time()-tp, time.time()-t0))
    if hit: break

say('')
say('triples completed at b<=2: %d of %d attempted' % (len([1 for _,s in done if 'exhaust' in s or 'infeas' in s]), len(TODO)))
for P, stt in done: say('   %s : %s' % (list(P), stt))
say('total solves %d, elapsed %.0fs, improvement: %s' % (nsolve, time.time()-t0, win is not None))
json.dump([[list(P), stt] for P, stt in done], open(os.path.join(S.HERE,'w_j3_done.json'),'w'), indent=1)
say('DONE')
