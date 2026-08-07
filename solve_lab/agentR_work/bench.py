#!/usr/bin/env python3
"""Scaling benchmark: how far up the prime size can automated reasoning solve the
ladder-fold problem?  Records solver statistics, checkpointed to runs/bench.json."""
import sys, os, json, time, subprocess, resource
import sibling, z3enc, witness

OUT = 'runs/bench.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
def save(): json.dump(res, open(OUT, 'w'), indent=1)

def smt2(d, logic, pin=None, path=None):
    from z3 import Solver
    C, bits = z3enc.build(d, logic, pin=pin)
    s = Solver(); s.add(C)
    txt = s.to_smt2()
    if path: open(path, 'w').write(txt)
    return txt, len(C)

def run_z3(path, tmo):
    t = time.time()
    r = subprocess.run(['z3', '-T:%d' % tmo, '-st', path], capture_output=True, text=True)
    o = r.stdout
    st = {}
    for ln in o.splitlines():
        ln = ln.strip().strip('()')
        for kk in ('conflicts', 'decisions', 'propagations', 'memory', 'restarts',
                   'arith-', 'rlimit', 'max-memory', 'num-checks'):
            if ln.startswith(kk):
                p = ln.split()
                if len(p) >= 2: st[p[0]] = p[-1]
    ans = 'unknown'
    for a in ('unsat', 'sat', 'timeout', 'unknown'):
        if a in o.split('\n')[0] if o else False: ans = a; break
    first = o.strip().split('\n')[0] if o.strip() else 'ERR'
    return {'answer': first, 'time': round(time.time() - t, 2), 'stats': st}

def run_cvc5(path, tmo):
    try:
        import cvc5
    except ImportError:
        return {'answer': 'cvc5-missing'}
    t = time.time()
    code = ('import cvc5,sys\n'
            's=cvc5.Solver()\n'
            's.setOption("tlimit","%d")\n'
            's.setOption("stats","true")\n'
            'import cvc5.pythonic\n' % (tmo * 1000))
    # cvc5 python has no file parser in all builds; use the C++ binary if present
    return {'answer': 'cvc5-no-binary'}

def do(m, logic, pinned, tmo=600):
    key = 'm%d_%s_%s' % (m, logic, 'pin' if pinned else 'free')
    if key in res: return res[key]
    d = sibling.instance(m)
    witness.witness(d, d['k'])          # encoding is provably satisfiable
    path = 'encodings/%s.smt2' % key
    t0 = time.time()
    txt, nc = smt2(d, logic, pin=(d['k'] if pinned else None), path=path)
    enc_t = time.time() - t0
    r = run_z3(path, tmo)
    r.update(m=m, logic=logic, pinned=pinned, n=d['n'], p=d['p'], k=d['k'],
             constraints=nc, smt2_bytes=len(txt), encode_time=round(enc_t, 2), tmo=tmo)
    res[key] = r; save()
    print('%-18s %-8s %8.1fs  %s  clauses/constraints=%d bytes=%d' %
          (key, r['answer'], r['time'], r['stats'].get('conflicts', '-'), nc, len(txt)), flush=True)
    return r

if __name__ == '__main__':
    for m in (8, 10, 12, 16, 20, 24, 28, 32):
        for logic in ('bv', 'int'):
            for pinned in (True, False):
                do(m, logic, pinned, tmo=int(sys.argv[1]) if len(sys.argv) > 1 else 300)
