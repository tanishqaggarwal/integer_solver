"""Sensitivity scan: perturb each free input, ripple, record which CHECK atoms move."""
import pickle, collections, heapq, time, sys, json
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256-2**32-977
NV = 38748

checks = [a for a in range(len(polys)) if a not in atom_out]
print(f'{len(checks)} check atoms', file=sys.stderr)
freeinp = [x for x in range(NV) if x not in definer]
print(f'{len(freeinp)} free inputs', file=sys.stderr)

def check_vals(v):
    return {a: evalpoly(polys[a], v) for a in checks}

if __name__ == '__main__':
    v0 = H.load_assignment('../best/new_instance_partial_39022.json')
    base = check_vals(v0)
    nzbase = {a:c for a,c in base.items() if c}
    print('baseline nonzero checks:', nzbase.keys(), file=sys.stderr)
    t0=time.time()
    out = {}
    sample = freeinp if len(sys.argv)<2 else freeinp[:int(sys.argv[1])]
    for i, f in enumerate(sample):
        v = list(v0); 
        try:
            ripple(v, {f: v0[f] + 1})
        except Exception as e:
            out[f] = 'ERR'; continue
        cv = check_vals(v)
        moved = {a: cv[a]-base[a] for a in checks if cv[a] != base[a]}
        out[f] = moved
        if i % 500 == 0:
            print(f'{i}/{len(sample)} t={time.time()-t0:.0f}s', file=sys.stderr)
    pickle.dump(out, open('sens.pkl','wb'))
    print(f'done in {time.time()-t0:.0f}s', file=sys.stderr)
    cnt = collections.Counter()
    for f, m in out.items():
        if m == 'ERR': continue
        for a in m: cnt[a]+=1
    print('checks most-affected:', cnt.most_common(15))
