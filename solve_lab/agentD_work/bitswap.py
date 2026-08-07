"""Swap the active gating bit in a bank (turn old off, new on AND set its two
pinned variables to their constants), re-solve advice, read the coordinates."""
import json, sys, time
import dlib as L
import engine2 as E
import adv3
P = L.P
COORD = {'x1': 12186, 'y1': 16742, 'x2': 14853, 'y2': 24908, 'x3': 22162, 'y3': 30213,
         's9192': 9192, 'A': 35389, 'B': 6671}
t = json.load(open('table.json'))
banks = json.load(open('banks.json'))
base = L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json')
st0 = E.St(base)


def entry(b):
    return [(x, C) for a, x, C in t[str(b)]]


def swap(old, new, sweeps=10):
    st = st0.clone()
    seeds = {old: 0, new: 1}
    for x, C in entry(new):
        seeds[x] = C
    st.apply(seeds)
    adv3.sweep(st, rounds=sweeps)
    return st


if __name__ == '__main__':
    bank = sys.argv[2] if len(sys.argv) > 2 else 'bank1'
    old = 24601 if bank == 'bank1' else 2081
    print('base', st0.score, {k: st0.v[u] % P for k, u in COORD.items()})
    for i in banks[bank][:8]:
        if i == old:
            continue
        st = swap(old, i)
        c = [C % P for x, C in entry(i)]
        d = {k: st.v[u] % P for k, u in COORD.items()}
        print(f'bit x_{i}: score={st.score} nz={len(st.nz())}')
        print('    consts', c)
        print('    coords', d)
