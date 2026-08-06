"""S11 step 67: the FREE-INPUT SUPPORT of every atom, as bitsets.

Forward-mode AD gives one Jacobian column per free input.  To know whether a
closure is COLUMN-CLOSED -- whether the rows I put in a system depend on any free
input I left out -- I need the transpose question: which free inputs reach a given
check?  Propagating a bitset along the same topological order answers it for every
atom at once, in one pass.

Set-union over-approximates (two paths can cancel), so `supp(c) subset U` is a
sound proof that U is closed for c, while `supp(c) - U` is only a candidate list
that still has to be tested with a real Jacobian column.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P


def build(vm, definer=None, ORDER=None, FREE=None):
    """Return (idx, freelist, supp) where supp[atom] is a bitset over freelist."""
    if definer is None:
        definer = L.definer
    if ORDER is None:
        ORDER = ad.ORDER
    if FREE is None:
        FREE = set(t for t in range(L.NVARS) if t not in definer)
    freelist = sorted(FREE)
    idx = {u: i for i, u in enumerate(freelist)}
    vs = [0] * L.NVARS
    for u in freelist:
        vs[u] = 1 << idx[u]
    for t in ORDER:
        a = definer[t]
        if ad.dpart(a, t, vm) % P == 0:
            vs[t] = 0
            continue
        m = 0
        for w in L.avars[a]:
            if w == t:
                continue
            if ad.dpart(a, w, vm) % P:
                m |= vs[w]
        vs[t] = m
    return idx, freelist, vs


def atom_supp(a, vm, vs):
    m = 0
    for w in L.avars[a]:
        if ad.dpart(a, w, vm) % P:
            m |= vs[w]
    return m


if __name__ == '__main__':
    import time
    v = L.load(os.path.join(HERE, 'mod9118_0.json'))
    vm = [x % P for x in v]
    t0 = time.time()
    idx, freelist, vs = build(vm)
    print(f'{len(freelist)} free inputs; bitset pass {time.time()-t0:.1f}s')
    for a in (21617, 29539):
        m = atom_supp(a, vm, vs)
        print(f'  atom a{a}: free-input support {bin(m).count("1")}')
