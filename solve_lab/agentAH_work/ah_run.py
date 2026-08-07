#!/usr/bin/env python3
"""AH: instrumented driver for the fleet's construction+closure routine.

The construction is NOT reimplemented.  `t_close2wj.close` is reproduced here verbatim
(same passes, same order, same guards) with four additions and no algorithmic change:

  1. a hard wall-clock DEADLINE (SIGALRM) so that a run that does not finish is reported
     as a STALL/TIMEOUT with a reason, never as "failed to close";
  2. a `safe` snapshot taken at the top of every outer iteration -- always a
     guard-accepted state -- which is what gets dumped if the deadline fires mid-pass;
  3. an explicit exit-reason string for every way out of the loop;
  4. telemetry counting how many of the two-wire root searches returned an empty
     result that was SAMPLED rather than EXHAUSTIVE (a sampled miss is a statement
     about the solver, not about the instance).

Usage: ah_run.py <tag> <n> <seed> <outer_max> <budget_sec>
"""
import os, sys, json, time, signal, collections, io
sys.set_int_max_str_digits(50_000_000)
AH = os.path.dirname(os.path.abspath(__file__))
T  = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)

tag, n, seed, outer_max, budget = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), \
                                  int(sys.argv[4]), int(sys.argv[5])

# ---------------------------------------------------------------- memory guard
# The fleet's `t_close2wj.joint_rootsets` expands a residue set with
#   rs = set(b for b in range(m) if b % ma in rs)
# where m is the atom's handle cofactor -- up to 16,595,977 in this instance -- and then
# `out.extend(...)` over it.  At |S|=24 seed 101 that reached 9.7 GB RSS and was killed by the
# kernel OOM killer (dmesg: "Killed process 30381 (python3) ... anon-rss:9727464kB").  A hard
# RLIMIT_AS turns that into a MemoryError this process reports as MEMORY_BLOWUP instead of an
# out-of-memory event that can take down another agent's job.
import resource
MEMCAP_GB = float(os.environ.get('AH_MEMCAP_GB', '3'))
_lim = int(MEMCAP_GB*1024**3)
resource.setrlimit(resource.RLIMIT_AS, (_lim, _lim))

t_import = time.time()
import t_close2wj as J
import t_close2w as C
t_import = time.time() - t_import

E = J.E; SL = J.SL; SHIFT = J.SHIFT; p = J.p; NV = J.NV; M = J.M
relift = J.relift; vars_of = J.vars_of; atomvalvars = J.atomvalvars
influences = J.influences; nzcount = J.nzcount; assignment = J.assignment
ORIENT = J.ORIENT; T1 = J.T1; T2 = J.T2; solve_group3 = J.solve_group3
from math import gcd
import itertools

# ---------------------------------------------------------------- cost instrumentation
# No algorithmic change: `nzcount` (the global guard, = relift + a full forward evaluation)
# is the routine's dominant cost, so it is wrapped in a counter.  The wrapper is installed
# in solve_group3's OWN globals dict so that guard evaluations made inside the fleet's
# single-wire pass are counted too.
COST = {'nz': 0, 'sg3': 0, 'sg3_t': 0.0, 'sg3_maxcombo': 0, 'capped': 0}
_nz_orig = nzcount

def nzcount(vv):
    COST['nz'] += 1
    return _nz_orig(vv)

_g3 = solve_group3.__globals__
_g3['nzcount'] = nzcount
J.nzcount = nzcount

CAP = int(os.environ.get('AH_CAP', '0'))   # 0 = fleet routine verbatim

# ---------------------------------------------------------------- exact fast root enumeration
# AH_FASTROOTS=1 swaps closeS3's brute-force `rootset_pp` for ah_roots, which returns the SAME
# SET via agent T's own t_poly.roots_pp.  This is an enumeration change, not an algorithm change;
# ah_roots_selftest.py checks set equality against the brute-force original on the instance's own
# prime powers.  Needed because one brute-force call can cost 1.7e7 big-int evaluations.
# ---------------------------------------------------------------- closure-ORDER knobs (controls)
# Neither changes the algorithm or the guards; both change only the ORDER in which the same
# candidate shifts are offered.  Used as the "different closure order" control: if a score
# ceiling is real it must survive both.
RNDSEED  = int(os.environ.get('AH_RND', '0'))      # 0 = the fleet's random.Random(20260807)
WIREORD  = os.environ.get('AH_WIREORD', 'desc')    # 'desc' = the fleet's -len(atoms) ordering
if RNDSEED:
    import random as _r
    J.rnd = _r.Random(RNDSEED)
    C.rnd = _r.Random(RNDSEED)
FASTROOTS = int(os.environ.get('AH_FASTROOTS', '0'))
AHR = None
if FASTROOTS:
    import ah_roots as AHR
    assert _g3 is _g3['roots_c'].__globals__, 'roots_c does not share solve_group3 globals'
    assert _g3 is _g3['rootset_pp'].__globals__, 'rootset_pp does not share those globals'
    _g3['rootset_pp'] = AHR.make(_g3['peval'])

def solve_group3_capped(vv, V, w, gen, base):
    """CONTROL VARIANT ONLY (enabled with AH_CAP>0).  Identical to the fleet's
    solve_group3 except that the CRT product over the violated atoms' root sets is
    truncated to the first CAP combinations.  Used to test whether a run that does not
    finish is limited by this enumeration (routine) or by the instance."""
    fitc = _g3['fitc']; roots_c = _g3['roots_c']; crt_list = _g3['crt_list']
    Rs = []
    for a in V:
        Cc = fitc(vv, a, w, gen)
        if Cc is None:
            return None
        rs = roots_c(Cc, abs(SL[a])//p)
        if not rs:
            return None
        Rs.append((rs, abs(SL[a])//p))
    tried = 0
    for combo in itertools.product(*[r for r, _ in Rs]):
        tried += 1
        if tried > CAP:
            COST['capped'] += 1
            return None
        t = crt_list([(r, c) for r, (_, c) in zip(combo, Rs)])
        if not t:
            continue
        old = vv[w]; vv[w] = old+p*t
        n = nzcount(vv)
        if n < base:
            return t
        vv[w] = old
    return None

_SG3 = solve_group3_capped if CAP else solve_group3

def SG3(vv, V, w, gen, base):
    COST['sg3'] += 1
    t0 = time.time()
    r = _SG3(vv, V, w, gen, base)
    COST['sg3_t'] += time.time()-t0
    return r

class Deadline(Exception):
    pass

def _alarm(sig, frame):
    raise Deadline()

signal.signal(signal.SIGALRM, _alarm)

LOG = open(os.path.join(AH, 'log_%s.txt' % tag), 'w')
SAMPLED_MISS = [0]      # empty two-wire root set found by SAMPLING  (solver gap)
EXHAUST_MISS = [0]      # empty two-wire root set found EXHAUSTIVELY (real obstruction on that pair)

def log(s):
    if 'NO JOINT ROOT' in s:
        if '(sampled)' in s:
            SAMPLED_MISS[0] += 1
        else:
            EXHAUST_MISS[0] += 1
    LOG.write(s + '\n'); LOG.flush()

# ------------------------------------------------------------------ verbatim close(), instrumented
def close(S, outer_max):
    """t_close2wj.close, unchanged except for the deadline/snapshot/reason instrumentation."""
    v, isl, valn = assignment(set(S), ORIENT); v[24468] = T1; v[18956] = T2
    vv = [0]*NV
    for k, x in v.items():
        vv[k] = x
    stats = {'relift_rounds': 0}
    for rd in range(60):
        bad = relift(vv)
        if not bad:
            break
        r = E.run(vv); fx = 0
        for a in bad:
            i = E.residx[a]; cur = r[i]; sm = abs(SL[a])
            if cur % p:
                continue
            imm = [q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old = vv[w]; vv[w] = old+p; d = E.run(vv)[i]-cur; vv[w] = old
                if d == 0:
                    continue
                gg = gcd(d, sm)
                if cur % gg:
                    continue
                mm = sm//gg
                t = (-(cur//gg))*pow((d//gg) % mm, -1, mm) % mm if mm > 1 else 0
                vv[w] = old+p*t; fx += 1; break
        stats['relift_rounds'] = rd+1
        log('   [pre] relift round %d: %d bad atoms, %d fixed, t=%.0fs'
            % (rd, len(bad), fx, time.time()-T0))
        if fx == 0:
            break
    stats['prelift_s'] = round(time.time()-T0, 1)
    log('   [pre] fixpoint done in %.0fs (%d rounds)'
        % (stats['prelift_s'], stats['relift_rounds']))
    gen = 0
    TGT = ('x24468', 'x18956')
    safe = vv[:]                      # last guard-consistent state
    reason = 'OUTER_MAX'
    outer_done = 0
    trace = []
    try:
        for outer in range(outer_max):
            safe = vv[:]
            outer_done = outer
            base = nzcount(vv); r = E.run(vv); gen += 1
            viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a]
                    and r[E.residx[a]] % abs(SL[a]) != 0]
            hl0 = [a for a in E.res if r[E.residx[a]] and a not in SL]
            trace.append({'outer': outer, 'nz': base, 'viol': len(viol), 'hl': len(hl0),
                          't': round(time.time()-T0, 1)})
            log('outer %d: global nonzero %d, violated c-conditions %d, nonzero handle-less %d'
                % (outer, base, len(viol), len(hl0)))
            if not viol and not hl0:
                reason = 'CLOSED_NO_VIOL'; break
            if not viol:
                if J.handleless_pass(vv, base, log):
                    continue
                if any(J.joint_pair(vv, a, base, log) for a in hl0):
                    continue
                if J.forced_exact_pass(vv, hl0, log):
                    continue
                log('   handle-less atoms remain and nothing moves them -> stop')
                reason = 'STALL_HANDLELESS'; break
            wires = collections.defaultdict(list)
            for a in viol:
                for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT) |
                          set(q for q in atomvalvars[a] if q in SHIFT)):
                    wires[w].append(a)
            prog = 0
            _wl = sorted(wires.items(),
                         key=(lambda kv: len(kv[1])) if WIREORD == 'asc'
                         else (lambda kv: -len(kv[1])))
            log('   single-wire pass over %d wires (max atoms/wire %d)'
                % (len(_wl), max([len(x[1]) for x in _wl], default=0)))
            _hb = time.time()
            for _wi, (w, ats) in enumerate(_wl):
                V = [a for a in ats if influences(vv, a, w)]
                if not V:
                    continue
                t = SG3(vv, V, w, gen, base)
                if t:
                    prog += 1; base = nzcount(vv); gen += 1
                if time.time()-_hb > 120:
                    _hb = time.time()
                    log('     [hb] wire %d/%d  nz-guard evals=%d  sg3 calls=%d  '
                        'sg3 time=%.0fs  accepted=%d  t=%.0fs'
                        % (_wi, len(_wl), COST['nz'], COST['sg3'], COST['sg3_t'],
                           prog, time.time()-T0))
            if prog:
                log('   single-wire pass: %d accepted' % prog); continue
            if J.handleless_pass(vv, nzcount(vv), log):
                continue
            log('   single-wire STALLED -> joint two-wire pass')
            r = E.run(vv)
            viol = [a for a in SL if r[E.residx[a]] != 0 and SL[a]
                    and r[E.residx[a]] % abs(SL[a]) != 0
                    and not any(t in a for t in TGT)]
            hl = [a for a in E.res if r[E.residx[a]] and a not in SL]
            if not viol and not hl:
                log('   only the two TARGET congruences remain -- CLOSED')
                reason = 'CLOSED_TARGETS_ONLY'; break
            if not any(J.joint_pair(vv, a, base, log) for a in viol+hl):
                if hl and J.forced_exact_pass(vv, hl, log):
                    gen += 1; continue
                log('   joint two-wire pass also stalled -> stop')
                reason = 'STALL_TWOWIRE'; break
            gen += 1
    except MemoryError:
        vv[:] = safe
        reason = 'MEMORY_BLOWUP'
        try:
            log('   *** MemoryError at the %.1f GB cap in outer %d -- rolled back to the '
                'last guard-consistent state; this is a BLOWUP OF THE ROUTINE, not a '
                'failure to close' % (MEMCAP_GB, outer_done))
        except Exception:
            pass
    except Deadline:
        vv[:] = safe
        reason = 'TIMEOUT'
        log('   *** DEADLINE %ds reached in outer %d -- rolled back to the last '
            'guard-consistent state; this is a STALL OF THE ROUTINE, not a failure to close'
            % (budget, outer_done))
    signal.alarm(0)
    relift(vv); r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    return vv, nz, reason, outer_done, trace, stats

if __name__ == '__main__':
    import random
    r7 = random.Random(seed)
    if n == 2 and seed == 7:
        S = [24601, 2081]                    # the historical |S|=2 control
    else:
        S = r7.sample(M['live'], n)
    assert len(set(S)) == n and set(S) <= set(M['live'])
    log('S(seed=%d,n=%d) = %s' % (seed, n, S))
    T0 = time.time()
    signal.alarm(budget)
    vv, nz, reason, outer_done, trace, stats = close(S, outer_max)
    wall = time.time()-T0
    out = os.path.join(AH, 'close_%s.json' % tag)
    json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]}, open(out, 'w'))
    meta = {'tag': tag, 'n': n, 'seed': seed, 'S': sorted(S), 'reason': reason,
            'outer_reached': outer_done, 'outer_max': outer_max, 'budget': budget,
            'wall': round(wall, 1), 'import_s': round(t_import, 1),
            'nz_atoms': len(nz), 'n_atoms': len(E.res), 'nz_list': nz[:20],
            'sampled_root_misses': SAMPLED_MISS[0], 'exhaustive_root_misses': EXHAUST_MISS[0],
            'trace': trace, 'relift_rounds': stats['relift_rounds'],
            'cost': COST, 'cap': CAP, 'fastroots': FASTROOTS,
            'rndseed': RNDSEED, 'wireord': WIREORD, 'memcap_gb': MEMCAP_GB,
            'peak_rss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'root_stats': (AHR.STATS if AHR else None),
            'json': out}
    json.dump(meta, open(os.path.join(AH, 'meta_%s.json' % tag), 'w'), indent=1)
    log('DONE tag=%s n=%d seed=%d reason=%s wall=%.1f nz_atoms=%d sampled_miss=%d exh_miss=%d'
        % (tag, n, seed, reason, wall, len(nz), SAMPLED_MISS[0], EXHAUST_MISS[0]))
    print('DONE tag=%s n=%d seed=%d reason=%s wall=%.1f nz_atoms=%d sampled_miss=%d exh_miss=%d'
          % (tag, n, seed, reason, wall, len(nz), SAMPLED_MISS[0], EXHAUST_MISS[0]), flush=True)
    LOG.close()
