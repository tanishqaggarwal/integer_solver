"""Activate the dormant MUX select x_38144 to break the pin chain on x_12186.

x_2121 = x_38144*x_13636 and x_8211 = x_13636*(1 - x_38144).  With x_38144 = 0 the setter routes
x_29524 <- x_22152 (pinned to CONST).  With x_38144 = 1 that link dies and x_29524 <- x_11648,
a FREE input -- so x_22649, x_12186, x_14853 and x_7068 all become free and A = atom 22229 can be
driven to 0.  The price is atom 32288 (x_38144 - 0), a bare pin with no handle.
"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def stage(v, tag, show=0):
    nz = allnz(v); f = H.evaluate(CODES, v)
    print(f'[{tag}] atoms={len(nz)}  EQ {len(CODES)-len(f)}/{len(CODES)} ({len(f)} failing)')
    if show:
        for a in nz[:show]:
            print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:90]}')
    return f, nz


if __name__ == '__main__':
    print('equation cost of the bare pin atom 32288 (x_38144 - 0):', len(a2e.get(32288, [])))
    for a in (16203, 16205, 16207, 22770, 22775, 36691, 37754, 22767):
        print(f'   atom {a}: {len(a2e.get(a,[]))} eqs   {src[a][:80]}')
    v = H.load_assignment('../best/new_instance_partial_39026.json')
    stage(v, '0 baseline 39,026')
    ripple(v, {38144: 1})
    stage(v, '1 dormant select activated', show=25)
    # now x_29524 = x_8253 + x_11648 with x_11648 free; drive x_12186 -> K1 so that A = 0
    K1 = v[2099]
    need = K1 - v[12186]
    print(f'\nneed to shift x_12186 by {str(need)[:30]}...')
    # find the 1:1 driver of x_12186 in this configuration
    NV = 38748
    boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
    nb = [x for x in range(NV) if x not in definer and x not in boolv]
    drv = None
    for f in nb:
        w = list(v); ripple(w, {f: v[f] + 1})
        if w[12186] - v[12186] == 1:
            drv = f; break
    print('1:1 driver of x_12186 =', drv)
    if drv is not None:
        ripple(v, {drv: v[drv] + need, 14853: v[14853] + need, 7068: K1})
        print('  x_12186 == K1 ?', v[12186] == K1, '  x_7068 == K1 ?', v[7068] == K1,
              '  u =', v[29322])
        print('  atom 22229 (A) =', evalpoly(polys[22229], v))
        stage(v, '2 A driven to zero', show=25)
    H.save_assignment(v, 'dormsolve_out.json')
